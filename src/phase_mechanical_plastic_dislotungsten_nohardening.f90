! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author David Cereceda, Lawrence Livermore National Laboratory
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Crystal plasticity model for bcc metals, especially tungsten, with fixed dislocation density
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) dislotungsten_nohardening

  type :: tParameters
    real(pREAL) :: &
      D = 1.0_pREAL                                                                                 !< grain size
    real(pREAL),               allocatable, dimension(:) :: &
      b_sl, &                                                                                       !< magnitude of Burgers vector [m]
      i_sl, &                                                                                       !< Adj. parameter for distance between 2 forest dislocations
      tau_Peierls, &                                                                                !< Peierls stress
      !* mobility law parameters
      Q_s, &                                                                                        !< activation energy for glide [J]
      p, &                                                                                          !< p-exponent in glide velocity
      q, &                                                                                          !< q-exponent in glide velocity
      B, &                                                                                          !< friction coefficient
      h, &                                                                                          !< height of the kink pair
      w, &                                                                                          !< width of the kink pair
      omega                                                                                         !< attempt frequency for kink pair nucleation
    real(pREAL),               allocatable, dimension(:,:) :: &
      h_sl_sl, &                                                                                    !< slip resistance from slip activity
      forestProjection
    real(pREAL),               allocatable, dimension(:,:,:) :: &
      P_sl, &
      P_nS_pos, &
      P_nS_neg
    integer :: &
      sum_N_sl                                                                                      !< total number of active slip system
    character(len=:),          allocatable               :: &
      isotropic_bound
    character(len=pSTRLEN), allocatable, dimension(:) :: &
      output
    character(len=:),          allocatable, dimension(:) :: &
      systems_sl
  end type tParameters                                                                              !< container type for internal constitutive parameters

  type :: tIndexDotState
    integer, dimension(2) :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tIndexDotState

  type :: tDislotungstenNohardeningState
    real(pREAL), dimension(:,:), pointer :: &
      rho_mob, &
      rho_dip, &
      gamma_sl
  end type tDislotungstenNohardeningState

  type :: tDislotungstenNohardeningDependentState
    real(pREAL), dimension(:,:), allocatable :: &
      Lambda_sl, &
      tau_pass
  end type tDislotungstenNohardeningDependentState

!--------------------------------------------------------------------------------------------------
! containers for parameters and state
  type(tParameters),                  allocatable, dimension(:) :: param
  type(tIndexDotState),               allocatable, dimension(:) :: indexDotState
  type(tDisloTungstenNohardeningState),          allocatable, dimension(:) :: state
  type(tDisloTungstenNohardeningDependentState), allocatable, dimension(:) :: dependentState

contains


!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_dislotungsten_nohardening_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph, i, &
    Nmembers, &
    sizeState, sizeDotState, &
    startIndex, endIndex
  integer,    dimension(:), allocatable :: &
    N_sl
  real(pREAL),dimension(:), allocatable :: &
    f_edge, &                                                                                       !< edge character fraction of total dislocation density
    rho_mob_0, &                                                                                    !< initial dislocation density
    rho_dip_0                                                                                       !< initial dipole density
    real(pREAL), dimension(:,:), allocatable :: &
    a_nS                                                                                            !< non-Schmid coefficients
  character(len=:), allocatable :: &
    refs, &
    extmsg, &
    nonSchmid_model
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech, &
    pl


  myPlasticity = plastic_active('dislotungsten_nohardening')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:dislotungsten_nohardening init  -+>>>'

  print'(/,1x,a)', 'D. Cereceda et al., International Journal of Plasticity 78:242–256, 2016'
  print'(  1x,a)', 'https://doi.org/10.1016/j.ijplas.2015.09.002'

  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(size(phases)))
  allocate(indexDotState(size(phases)))
  allocate(state(size(phases)))
  allocate(dependentState(size(phases)))
  extmsg = ''

  do ph = 1, size(phases)
    if (.not. myPlasticity(ph)) cycle

    associate(prm => param(ph), &
              stt => state(ph), dst => dependentState(ph), &
              idx_dot => indexDotState(ph))

    phase => phases%get_dict(ph)
    mech  => phase%get_dict('mechanical')
    pl  => mech%get_dict('plastic')

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)
    refs = config_listReferences(pl,indent=3)
    if (len(refs) > 0) print'(/,1x,a)', refs

#if defined (__GFORTRAN__)
    prm%output = output_as1dStr(pl)
#else
    prm%output = pl%get_as1dStr('output',defaultVal=emptyStrArray)
#endif

    prm%isotropic_bound = pl%get_asStr('isotropic_bound',defaultVal='isostrain')

!--------------------------------------------------------------------------------------------------
! slip related parameters
    N_sl = pl%get_as1dInt('N_sl',defaultVal=emptyIntArray)
    prm%sum_N_sl = sum(abs(N_sl))
    slipActive: if (prm%sum_N_sl > 0) then
      prm%P_sl = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph))
      prm%systems_sl = crystal_labels_slip(N_sl,phase_lattice(ph))

      a_nS = pl%get_as2dReal('a_non-Schmid',defaultVal=reshape(emptyRealArray,[0,0]))
      nonSchmid_model = pl%get_asStr('non-Schmid_model', defaultVal='G')
      prm%P_nS_pos = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph), &
                                               nonSchmidCoefficients=a_nS,sense=+1, &
                                               nonSchmid_model=nonSchmid_model)
      prm%P_nS_neg = crystal_SchmidMatrix_slip(N_sl,phase_lattice(ph),phase_cOverA(ph), &
                                               nonSchmidCoefficients=a_nS,sense=-1, &
                                               nonSchmid_model=nonSchmid_model)

      prm%h_sl_sl = crystal_interaction_SlipBySlip(N_sl,pl%get_as1dReal('h_sl-sl'), &
                                                   phase_lattice(ph))

      prm%D = pl%get_asReal('D')

      f_edge          = pl%get_as1dReal('f_edge',      requiredChunks=N_sl, &
                                                       defaultVal=[(0.5_pREAL, i=1,size(N_sl))])
      rho_mob_0       = pl%get_as1dReal('rho_mob_0',   requiredChunks=N_sl)
      rho_dip_0       = pl%get_as1dReal('rho_dip_0',   requiredChunks=N_sl)
      prm%b_sl        = pl%get_as1dReal('b_sl',        requiredChunks=N_sl)
      prm%Q_s         = pl%get_as1dReal('Q_s',         requiredChunks=N_sl)
      prm%i_sl        = pl%get_as1dReal('i_sl',        requiredChunks=N_sl)
      prm%tau_Peierls = pl%get_as1dReal('tau_Peierls', requiredChunks=N_sl)
      prm%p           = pl%get_as1dReal('p_sl',        requiredChunks=N_sl)
      prm%q           = pl%get_as1dReal('q_sl',        requiredChunks=N_sl)
      prm%h           = pl%get_as1dReal('h',           requiredChunks=N_sl)
      prm%w           = pl%get_as1dReal('w',           requiredChunks=N_sl)
      prm%omega       = pl%get_as1dReal('omega',       requiredChunks=N_sl)
      prm%B           = pl%get_as1dReal('B',           requiredChunks=N_sl)

      prm%forestProjection = spread(          f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_edge (N_sl,phase_lattice(ph),phase_cOverA(ph)) &
                           + spread(1.0_pREAL-f_edge,1,prm%sum_N_sl) &
                           * crystal_forestProjection_screw(N_sl,phase_lattice(ph),phase_cOverA(ph))

      ! sanity checks
      if (any(rho_mob_0        <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_mob_0'
      if (any(rho_dip_0        <  0.0_pREAL)) extmsg = trim(extmsg)//' rho_dip_0'
      if (any(prm%b_sl         <= 0.0_pREAL)) extmsg = trim(extmsg)//' b_sl'
      if (any(prm%Q_s          <= 0.0_pREAL)) extmsg = trim(extmsg)//' Q_s'
      if (any(prm%tau_Peierls  <  0.0_pREAL)) extmsg = trim(extmsg)//' tau_Peierls'
      if (any(prm%B            <  0.0_pREAL)) extmsg = trim(extmsg)//' B'

    else slipActive
      rho_mob_0 = emptyRealArray
      rho_dip_0 = emptyRealArray
      allocate(prm%b_sl, &
               prm%i_sl, &
               prm%tau_Peierls, &
               prm%Q_s, &
               prm%p, &
               prm%q, &
               prm%B, &
               prm%h, &
               prm%w, &
               prm%omega, &
               source = emptyRealArray)
      allocate(prm%forestProjection(0,0))
      allocate(prm%h_sl_sl         (0,0))
    end if slipActive

!--------------------------------------------------------------------------------------------------
! allocate state arrays
    Nmembers = count(material_ID_phase == ph)
    sizeDotState = size(['rho_mob ','rho_dip ','gamma_sl']) * prm%sum_N_sl
    sizeState = sizeDotState

    call phase_allocateState(plasticState(ph),Nmembers,sizeState,sizeDotState,0)
    deallocate(plasticState(ph)%dotState) ! ToDo: remove dotState completely

!--------------------------------------------------------------------------------------------------
! state aliases and initialization
    startIndex = 1
    endIndex   = prm%sum_N_sl
    idx_dot%rho_mob = [startIndex,endIndex]
    stt%rho_mob => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_mob = spread(rho_mob_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_rho'

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%rho_dip = [startIndex,endIndex]
    stt%rho_dip => plasticState(ph)%state(startIndex:endIndex,:)
    stt%rho_dip = spread(rho_dip_0,2,Nmembers)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_rho',defaultVal=1.0_pREAL)

    startIndex = endIndex + 1
    endIndex   = endIndex + prm%sum_N_sl
    idx_dot%gamma_sl = [startIndex,endIndex]
    stt%gamma_sl => plasticState(ph)%state(startIndex:endIndex,:)
    plasticState(ph)%atol(startIndex:endIndex) = pl%get_asReal('atol_gamma',defaultVal=1.0e-6_pREAL)
    if (any(plasticState(ph)%atol(startIndex:endIndex) < 0.0_pREAL)) extmsg = trim(extmsg)//' atol_gamma'

    allocate(dst%Lambda_sl(prm%sum_N_sl,Nmembers), source=0.0_pREAL)
    allocate(dst%tau_pass (prm%sum_N_sl,Nmembers), source=0.0_pREAL)

    end associate

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg))

  end do

end function plastic_dislotungsten_nohardening_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate plastic velocity gradient and its tangent.
!--------------------------------------------------------------------------------------------------
pure module subroutine dislotungsten_nohardening_LpAndItsTangent(Lp,dLp_dMp, &
                                                     Mp,ph,en)
  real(pREAL), dimension(3,3),     intent(out) :: &
    Lp                                                                                              !< plastic velocity gradient
  real(pREAL), dimension(3,3,3,3), intent(out) :: &
    dLp_dMp                                                                                         !< derivative of Lp with respect to the Mandel stress
  real(pREAL), dimension(3,3),      intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                          intent(in) :: &
    ph, &
    en

  integer :: &
    i,k,l,m,n
  real(pREAL) :: &
    T                                                                                               !< temperature
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    dot_gamma, ddot_gamma_dtau
  real(pREAL), dimension(3,3,param(ph)%sum_N_sl) :: &
    P_nS


  T = thermal_T(ph,en)
  Lp = 0.0_pREAL
  dLp_dMp = 0.0_pREAL

  associate(prm => param(ph))

    call kinetics(Mp,T,ph,en, dot_gamma,ddot_gamma_dtau)
    P_nS = merge(prm%P_nS_pos,prm%P_nS_neg, spread(spread(dot_gamma,1,3),2,3)>0.0_pREAL)            ! faster than 'merge' in loop
    do i = 1, prm%sum_N_sl
      Lp = Lp + dot_gamma(i)*prm%P_sl(1:3,1:3,i)
      forall (k=1:3,l=1:3,m=1:3,n=1:3) &
        dLp_dMp(k,l,m,n) = dLp_dMp(k,l,m,n) &
                         + ddot_gamma_dtau(i) * prm%P_sl(k,l,i) * P_nS(m,n,i)
    end do

  end associate

end subroutine dislotungsten_nohardening_LpAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief Keep the dislocation densities fixed and accumulate plastic shear.
!--------------------------------------------------------------------------------------------------
module function dislotungsten_nohardening_dotState(Mp,ph,en) result(dotState)

  real(pREAL), dimension(3,3),  intent(in) :: &
    Mp                                                                                              !< Mandel stress
  integer,                      intent(in) :: &
    ph, &
    en
  real(pREAL), dimension(plasticState(ph)%sizeDotState) :: &
    dotState

  associate(dot_rho_mob => dotState(indexDotState(ph)%rho_mob(1):indexDotState(ph)%rho_mob(2)), &
            dot_rho_dip => dotState(indexDotState(ph)%rho_dip(1):indexDotState(ph)%rho_dip(2)), &
            dot_gamma   => dotState(indexDotState(ph)%gamma_sl(1):indexDotState(ph)%gamma_sl(2)))

    dot_rho_mob = 0.0_pREAL
    dot_rho_dip = 0.0_pREAL
    call kinetics(Mp,thermal_T(ph,en),ph,en,dot_gamma)
    dot_gamma = abs(dot_gamma)

  end associate

end function dislotungsten_nohardening_dotState


!--------------------------------------------------------------------------------------------------
!> @brief Calculate derived quantities from state.
!--------------------------------------------------------------------------------------------------
module subroutine dislotungsten_nohardening_dependentState(ph,en)

  integer, intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    Lambda_sl_inv


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    dst%tau_pass(:,en) = elastic_mu(ph,en,prm%isotropic_bound)*prm%b_sl &
                       * sqrt(math_clip(matmul(prm%h_sl_sl, &
                              math_clip(stt%rho_mob(:,en),left=0.0_pREAL) + &
                              math_clip(stt%rho_dip(:,en),left=0.0_pREAL)), left=0.0_pREAL))

    Lambda_sl_inv = 1.0_pREAL/prm%D &
                  + sqrt(math_clip(matmul(prm%forestProjection, &
                         math_clip(stt%rho_mob(:,en),left=0.0_pREAL) + &
                         math_clip(stt%rho_dip(:,en),left=0.0_pREAL)), left=0.0_pREAL))/prm%i_sl
    dst%Lambda_sl(:,en) = math_clip(Lambda_sl_inv**(-1.0_pREAL), left = prm%w + prm%b_sl)

  end associate

end subroutine dislotungsten_nohardening_dependentState


!--------------------------------------------------------------------------------------------------
!> @brief Write results to HDF5 output file.
!--------------------------------------------------------------------------------------------------
module subroutine plastic_dislotungsten_nohardening_result(ph,group)

  integer,          intent(in) :: ph
  character(len=*), intent(in) :: group

  integer :: ou


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do ou = 1,size(prm%output)

      select case(trim(prm%output(ou)))

        case('rho_mob')
          call result_writeDataset(stt%rho_mob,group,trim(prm%output(ou)), &
                                   'mobile dislocation density','1/m²',prm%systems_sl)
        case('rho_dip')
          call result_writeDataset(stt%rho_dip,group,trim(prm%output(ou)), &
                                   'dislocation dipole density','1/m²',prm%systems_sl)
        case('gamma_sl')
          call result_writeDataset(stt%gamma_sl,group,trim(prm%output(ou)), &
                                   'plastic shear','1',prm%systems_sl)
        case('Lambda_sl')
          call result_writeDataset(dst%Lambda_sl,group,trim(prm%output(ou)), &
                                   'mean free path for slip','m',prm%systems_sl)
        case('tau_pass')
          call result_writeDataset(dst%tau_pass,group,trim(prm%output(ou)), &
                                   'threshold stress for slip','Pa',prm%systems_sl)
      end select

    end do

  end associate

end subroutine plastic_dislotungsten_nohardening_result


!--------------------------------------------------------------------------------------------------
!> @brief Calculate shear rates on slip systems, their derivatives with respect to resolved
!         stress, and the resolved stress.
!> @details Derivatives and resolved stress are calculated only optionally.
! NOTE: Contrary to common convention, here the result (i.e. intent(out)) variables have to be put
! at the end since some of them are optional.
!--------------------------------------------------------------------------------------------------
pure subroutine kinetics(Mp,T,ph,en, &
                         dot_gamma,ddot_gamma_dtau,tau)

  real(pREAL), dimension(3,3),                           intent(in) :: &
    Mp                                                                                              !< Mandel stress
  real(pREAL),                                           intent(in) :: &
    T                                                                                               !< temperature
  integer,                                               intent(in) :: &
    ph, &
    en

  real(pREAL), dimension(param(ph)%sum_N_sl),           intent(out) :: &
    dot_gamma
  real(pREAL), dimension(param(ph)%sum_N_sl), optional, intent(out) :: &
    ddot_gamma_dtau, &
    tau
  real(pREAL), dimension(param(ph)%sum_N_sl) :: &
    effectiveLength, &
    tau_pos, tau_neg, tau_eff
  real(pREAL) :: &
    StressRatio_i, StressRatio_p_i, StressRatio_pm1_i, &
    deltaH_factor_i, deltaH_pow_i, &
    BoltzmannRatio_i, b_rho_i, &
    t_n_i, t_k_i, dtn_i, dtk_i
  integer :: i


  associate(prm => param(ph), stt => state(ph), dst => dependentState(ph))

    do i = 1, prm%sum_N_sl
      tau_pos(i) = math_tensordot(Mp, prm%P_nS_pos(1:3,1:3,i))
      tau_neg(i) = math_tensordot(Mp, prm%P_nS_neg(1:3,1:3,i))
      tau_eff(i) = max(0.0_pREAL, max(tau_pos(i), tau_neg(i)) - dst%tau_pass(i,en))
    end do

    if (present(tau)) tau = tau_eff

    effectiveLength = math_clip(dst%Lambda_sl(:,en) - prm%w, left = prm%b_sl)

    dot_gamma = 0.0_pREAL
    if (present(ddot_gamma_dtau)) ddot_gamma_dtau = 0.0_pREAL

    do i = 1, prm%sum_N_sl
      if (tau_eff(i) > tol_math_check) then
        BoltzmannRatio_i = prm%Q_s(i)/(K_B*T)
        b_rho_i = max(0.0_pREAL, stt%rho_mob(i,en)) * prm%b_sl(i)

        StressRatio_i = tau_eff(i)/prm%tau_Peierls(i)
        StressRatio_p_i = StressRatio_i ** prm%p(i)
        deltaH_factor_i = max(0.0_pREAL, 1.0_pREAL - StressRatio_p_i)

        t_n_i = prm%b_sl(i)*exp(BoltzmannRatio_i * (deltaH_factor_i ** prm%q(i))) &
              / (prm%omega(i)*effectiveLength(i))
        t_k_i = effectiveLength(i) * prm%B(i) /(2.0_pREAL*prm%b_sl(i)*tau_eff(i))

        if (tau_pos(i) > tau_neg(i)) then
          dot_gamma(i) = b_rho_i * prm%h(i)/(t_n_i + t_k_i)
        else
          dot_gamma(i) = -b_rho_i * prm%h(i)/(t_n_i + t_k_i)
        end if

        if (present(ddot_gamma_dtau)) then
          if (StressRatio_p_i < 1.0_pREAL) then
            StressRatio_pm1_i = StressRatio_i**(prm%p(i)-1.0_pREAL)
            if (deltaH_factor_i > 0.0_pREAL .and. prm%q(i) > 1.0_pREAL) then
              deltaH_pow_i = deltaH_factor_i**(prm%q(i) - 1.0_pREAL)
            else if (deltaH_factor_i > 0.0_pREAL) then
              deltaH_pow_i = 1.0_pREAL
            else
              deltaH_pow_i = 0.0_pREAL
            end if
            dtn_i = -1.0_pREAL * t_n_i * BoltzmannRatio_i * prm%p(i) * prm%q(i) * deltaH_pow_i &
                  * StressRatio_pm1_i / prm%tau_Peierls(i)
          else
            dtn_i = 0.0_pREAL
          end if
          dtk_i = -1.0_pREAL * t_k_i / tau_eff(i)
          ddot_gamma_dtau(i) = -1.0_pREAL * dot_gamma(i) * (dtn_i + dtk_i) / (t_n_i + t_k_i)
        end if
      end if
    end do

  end associate

end subroutine kinetics

end submodule dislotungsten_nohardening
