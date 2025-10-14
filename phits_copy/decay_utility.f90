module decay_utility

  use udm_Parameter
  use udm_Utility
  implicit none
  private ! Functions and variables are set to private by default.
  public :: n_body_decay, decay_with_data
  logical         , save :: is_loaded=.false.
  integer         , allocatable, save ::  kf_in(:)
  integer         , allocatable, save :: n_mode(:)
  integer         , allocatable, save ::  n_out(:,:)
  double precision, allocatable, save ::     Br(:,:)
  integer         , allocatable, save :: kf_out(:,:,:)
  integer, save :: n_in
  contains

  !=======================================================================
  double precision function betabar(a,b)
    implicit none
    double precision a, b
    double precision tmp
    tmp=1d0-2d0*(a+b)+(a-b)**2
    if(tmp>=0d0) then
      betabar = sqrt(tmp)
    else
      betabar = -1d0
    endif
    return
  end function betabar

  !=======================================================================
  double precision function random(min,max)
    double precision min, max
    random = (max-min)*(1d0-get_random_0to1()) + min
    return
  end function random

  !=======================================================================
  function T(gamma, gamma_beta, COSTH, COSPHI, SINPHI)
    implicit none
    double precision :: T(4,4)
    double precision :: Bz(4,4), Ry(4,4), Rz(4,4)
    double precision gamma
    double precision gamma_beta
    double precision COSTH
    double precision SINTH
    double precision COSPHI
    double precision SINPHI

    SINTH = sqrt(dabs(1d0-COSTH**2))

    Bz(1,:) = (/ gamma     ,  0d0,  0d0,  gamma_beta /)
    Bz(2,:) = (/ 0d0       ,  1d0,  0d0,  0d0        /)
    Bz(3,:) = (/ 0d0       ,  0d0,  1d0,  0d0        /)
    Bz(4,:) = (/ gamma_beta,  0d0,  0d0,  gamma      /)

    Ry(1,:) = (/ 1d0,  0d0  ,  0d0,  0d0   /)
    Ry(2,:) = (/ 0d0,  COSTH,  0d0,  SINTH /)
    Ry(3,:) = (/ 0d0,  0d0  ,  1d0,  0d0   /)
    Ry(4,:) = (/ 0d0, -SINTH,  0d0,  COSTH /)

    Rz(1,:) = (/ 1d0,  0d0   ,  0d0   ,  0d0 /)
    Rz(2,:) = (/ 0d0,  COSPHI, -SINPHI,  0d0 /)
    Rz(3,:) = (/ 0d0,  SINPHI,  COSPHI,  0d0 /)
    Rz(4,:) = (/ 0d0,  0d0   ,  0d0   ,  1d0 /)

    T = matmul(Ry,Bz)
    T = matmul(Rz,T)
    return
  end function T


  !=======================================================================
  logical function n_body_decay(kf0,kin0,kfs)
    ! Isotropic decay in CM frame
    ! ---------------------------
    ! kf0  : kf-codes of the decaying particle
    ! kin0 : kinetic energy of the decaying particles at Lab-frame
    ! kfs  : kf-codes of the final states
  !=======================================================================
    implicit none
    integer kf0
    double precision kin0
    integer,intent(in) :: kfs(:)

    ! ------------------------------
    integer i,j,k
    integer n
    double precision m0 ! initial mass (M -> m0)
    double precision, allocatable :: m(:) ! final state masses
    double precision, allocatable :: r(:) ! normalized intermediate masses, Q^2/m0^2
    double precision, allocatable :: w(:) ! normalized final state masses, m^2/m0^2
    double precision, allocatable :: costh(:)
    double precision, allocatable :: sinth(:)
    double precision, allocatable :: cosphi(:)
    double precision, allocatable :: sinphi(:)
    double precision, allocatable :: final_states(:,:)
    double precision, allocatable :: rmin(:)
    double precision, allocatable :: rmax(:)
    double precision costhmin, costhmax
    double precision phimin, phimax
    double precision :: M_PI = 3.1415926535d0
    logical r_are_accepted, pass_all_condition
    double precision integrand,iMAX,integrand_max,temp,x,y,absp,sign_
    double precision :: p_before(4), p_after(4)
    double precision gamma,gamma_beta,COSTH_,COSPHI_,SINPHI_
    double precision :: Tn(4,4)
    double precision ener
    double precision mi_sum
    integer nTRY, nTRY_MAX, nPASS

    n_body_decay=.false.

    ! ------------------------------
    ! check
    if(size(kfs)<2) then
      print*,"[Caution] n_body_decay: size(kfs)<2"
      return
    endif

    ! ------------------------------
    ! initialize
    integrand_max=1d0
    n = size(kfs)      ! final states number
    m0 = get_mass(kf0) ! initial_mass (decay particle mass)
    ! -----
    allocate(m(n))
    do i=1,n
      m(i)=get_mass(kfs(i))
    enddo

    ! ------------------------------
    ! check
    mi_sum=0d0
    do i=1,n
      mi_sum=mi_sum+m(i)
    enddo
    if(mi_sum >= m0) then
      deallocate(m)
      return
    endif

    ! ------------------------------
    allocate(w(n))
    do i=1,n
      w(i)=m(i)**2/m0**2
    enddo
    ! -----
    allocate(r(n))
    allocate(rmin(n))
    allocate(rmax(n))
    allocate(costh(n))
    allocate(sinth(n))
    allocate(cosphi(n))
    allocate(sinphi(n))
    ! ------------------------------

    ! ------------------------------
    do i=2, n-1
      rmax(i)=1d0
      rmin(i)=0d0
      do j=i+1, n
        rmax(i)=rmax(i)-sqrt(w(j))
      enddo
      do j=1, i
        rmin(i)=rmin(i)+sqrt(w(j))
      enddo
      rmax(i)=rmax(i)**2
      rmin(i)=rmin(i)**2
    enddo
    costhmin = -1d0
    costhmax =  1d0
    phimin = 0d0
    phimax = 2d0*M_PI

    ! ------------------------------
    ! Generate r variables.
    if(n==2) then
      r_are_accepted=.true.
    else
      r_are_accepted=.false.
    endif

    r(1)=m(1)**2/m0**2
    r(n)=1d0 ! m0**2/m0**2

    ! ============================================================
    nTRY=0
    nTRY_MAX=50000000
    nPASS=0
    iMAX=0d0
    do while (.not. r_are_accepted)

      nTRY=nTRY+1
      if(nTRY>nTRY_MAX) then
        print*,"Failed: n_body_decay. Maximum number of attempts exceeded."
        ! print*,m0,n
        ! do i=1,n
        !   print*,"->",kfs(i),m(i)
        ! enddo
        return
      endif

      do i=2, n-1
        r(i)=random(rmin(i), rmax(i))
      enddo

      ! condition 1-1 ------------------
      pass_all_condition = .true.
      do i=3, n-1
        if(sqrt(r(i)) < sqrt(r(i-1)) + sqrt(w(i))) pass_all_condition = .false.
      enddo
      if(.not. pass_all_condition) cycle
      ! condition 1-2 ------------------
      do i=1,n
        if(i==1)then
          j   =i+1
        else
          j   =i
        endif
        if(betabar(r(j-1)/r(j), w(j)/r(j))<0) pass_all_condition = .false.
      enddo
      if(.not. pass_all_condition) cycle
      ! condition 2 ------------------
      integrand=1d0
      do i=2,n
        integrand = integrand*betabar(r(i-1)/r(i), w(i)/r(i))
      enddo
      if(integrand>iMAX) iMAX=integrand
      if(nPASS>1000 .or. nTRY>nTRY_MAX/5) integrand_max=iMAX*20d0
      if(integrand/integrand_max < get_random_0to1()) then
        nPASS=nPASS+1
        cycle
      endif

      r_are_accepted = .true.
    enddo
    ! ============================================================




    ! ------------------------------
    ! Generate angles in each CM-frame.
    costh(1)=1d0
    sinth(1)=0d0
    cosphi(1)=1d0
    sinphi(1)=0d0

    do i=2,n
      x = random(costhmin, costhmax)
      y = random(phimin  , phimax  )
      costh(i)=x
      sinth(i)=sqrt(1d0-x*x)
      cosphi(i)=cos(y)
      sinphi(i)=sin(y)
    enddo

    ! ------------------------------
    ! Define 4-momenta in each CM-frame.
    allocate(final_states(n,4))
    do i=1,n
      if(i==1)then
        j   =i+1
        sign_= 1d0
      else
        j   =i
        sign_=-1d0
      endif

      absp = 0.5d0 * m0 * sqrt(r(j)) * betabar(r(j-1)/r(j), w(j)/r(j))
      p_before(1) = 0.5d0*m0*sqrt(r(j)) * (1d0 + sign_ * (r(j-1)-w(j))/r(j) )
      p_before(2) = sign_ * absp * sinth(j) * cosphi(j)
      p_before(3) = sign_ * absp * sinth(j) * sinphi(j)
      p_before(4) = sign_ * absp * costh(j)

      ! -----
      ! Tn(j)
      Tn=0d0
      Tn(1,1)=1d0
      Tn(2,2)=1d0
      Tn(3,3)=1d0
      Tn(4,4)=1d0
      do k=j,n-1
        gamma      = 0.5d0 * sqrt(r(k+1)/r(k)) * (1d0 + r(k)/r(k+1) - w(k+1)/r(k+1))
        gamma_beta = 0.5d0 * sqrt(r(k+1)/r(k)) * betabar(r(k)/r(k+1), w(k+1)/r(k+1))
        COSTH_ = costh(k+1)
        COSPHI_ = cosphi(k+1)
        SINPHI_ = sinphi(k+1)
        Tn = matmul(T(gamma, gamma_beta, COSTH_, COSPHI_, SINPHI_), Tn)
      enddo
      p_after = matmul(Tn,p_before)

      ! ---------------------------------------
      ! Boost the 4-momenta into the Lab-frame.
      if(kin0>0d0) then
        ener=kin0+m0
        absp=sqrt(dabs(ener**2-m0**2))
        gamma      = ener/m0
        gamma_beta = absp/m0
        !                T(gamma, gamma_beta, COSTH, COSPHI, SINPHI)
        p_after = matmul(T(gamma, gamma_beta, 1d0  , 1d0   , 0d0   ), p_after)
      endif

      ! -----
      final_states(i,1)=p_after(1)
      final_states(i,2)=p_after(2)
      final_states(i,3)=p_after(3)
      final_states(i,4)=p_after(4)
    enddo

    ! ! ------------------------------
    ! ! check conservation law
    ! do i=1,4
    !   temp=0d0
    !   do j=1,n
    !     temp=temp+final_states(j,i)
    !   enddo
    !   print*,i,temp
    ! enddo

    ! ------------------------------
    ! Fill the boosted 4-momenta.
    call initialize_udm_event_info
    set_final_state_number=n
    do i=1,n
      set_kf(i)=kfs(i)
      set_Total_Energy_in_MeV(i)=final_states(i,1)
      set_Px_in_MeV(i)          =final_states(i,2)
      set_Py_in_MeV(i)          =final_states(i,3)
      set_Pz_in_MeV(i)          =final_states(i,4)
      if(.not. set_Total_Energy_in_MeV(i)-get_mass(set_kf(i))>0) then
        deallocate(m)
        deallocate(w)
        deallocate(r)
        deallocate(rmin)
        deallocate(rmax)
        deallocate(costh)
        deallocate(sinth)
        deallocate(cosphi)
        deallocate(sinphi)
        deallocate(final_states)
        return
      endif
    enddo

    ! ------------------------------
    ! everything OK
    call fill_final_state_for_decay
    n_body_decay=.true.
    deallocate(m)
    deallocate(w)
    deallocate(r)
    deallocate(rmin)
    deallocate(rmax)
    deallocate(costh)
    deallocate(sinth)
    deallocate(cosphi)
    deallocate(sinphi)
    deallocate(final_states)
    return

  end function n_body_decay

  !=======================================================================
  subroutine load_data
  !=======================================================================
    implicit none
    integer iii
    integer i_in, i_mode, i_out
    integer i1, i2
    character(len=400) :: decay_path, particle_name, line, ltmp
    logical exist_decay_file
    double precision BrSum,tmp

    ! -------------------
    if(is_loaded) return
    is_loaded = .true.

    ! -------------------
    allocate( kf_in(27))
    allocate(n_mode(27))
    allocate( n_out(27, 1100))
    allocate(    Br(27, 1100))
    allocate(kf_out(27, 1100, 15))

    ! -------------------
    ! load data
    i_in=0
    n_in=0
    do iii=1,2
      if(iii==1) decay_path=trim(udm_phits_path_fun())//"/data/decay/decay_data_1.dat"
      if(iii==2) decay_path=trim(udm_phits_path_fun())//"/data/decay/decay_data_2.dat"
      inquire(file=decay_path, exist=exist_decay_file)

      ! ----------
      if(exist_decay_file) then
        print*,"Loaded: "//trim(decay_path)
      else
        print*,trim(decay_path)//" does not exist. so the following kf's decays are ignored."
        if(iii==1) print*,"kf: 15,-15,113,213,-213,223,20213,-20213,313,-313,323,-323,333,411,-411,421,-421,431,-431"
        if(iii==2) print*,"kf: 511,-511,521,-521,531,-531,541,-541"
        exit
      endif

      ! ----------
      open(1,file=trim(decay_path),status="old")
      do while (.true.)
        read(1,*) particle_name, i1, i2
        if(trim(particle_name)=="END") exit
        ! ----------
        i_in=i_in+1
        n_in=n_in+1
        kf_in (i_in)=i1
        n_mode(i_in)=i2
        ! ----------
        BrSum=0d0
        do i_mode=1, n_mode(i_in)
          read(1,"(a)") line
          read(line,*) i1, Br(i_in, i_mode), tmp, n_out(i_in, i_mode)
          read(line,*) i1, Br(i_in, i_mode), tmp, n_out(i_in, i_mode), kf_out(i_in, i_mode, 1:n_out(i_in, i_mode))
          BrSum=BrSum+Br(i_in, i_mode)
        enddo
        ! ----------
        do i_mode=1, n_mode(i_in)
          if(iii==1) Br(i_in, i_mode)=Br(i_in, i_mode)/BrSum ! normalize only data-1.
        enddo
        ! ----------
      enddo
      close(1)
    enddo

    return

  end subroutine load_data

  !=======================================================================
  subroutine decay_with_data(kf0, kin0)
  !=======================================================================
    implicit none
    integer kf0
    double precision kin0
    logical :: data_exist = .false.
    integer i_in, i_out, i_mode
    double precision BrSum, random
    integer, allocatable :: kfs(:)
    logical is_success
    integer i

    ! ----------
    call load_data

    ! ----------
    ! get i_in
    do i_in=1,n_in
      if(kf0==kf_in(i_in)) then
        data_exist = .true.
        exit
      endif
    enddo
    if(.not. data_exist) return

    ! ----------
    ! choose decay mode
    BrSum=0d0
    random=get_random_0to1()
    do i_mode=1, n_mode(i_in)
      BrSum = BrSum+Br(i_in, i_mode)
      if(BrSum > random) then
        allocate(kfs(n_out(i_in, i_mode)))
        do i_out=1, n_out(i_in, i_mode)
          kfs(i_out)=kf_out(i_in, i_mode, i_out)
        enddo
        is_success=n_body_decay(kf0, kin0, kfs)
        deallocate(kfs)
        if(is_success) then
          return
        else
          ! print*,kf0, kin0, BrSum, random
          exit
        endif
      endif
    enddo

    ! ----------
    ! try sampling 5 times more
    do i=1,5
      if(abs(kf0)== 511 .or. abs(kf0)== 521 .or. abs(kf0)== 531 .or. abs(kf0)== 541) then
        i_mode=get_random_int(451,950)
      ! elseif(abs(kf0)== 411) then
      !   i_mode=get_random_int(1,10)
      else
        i_mode=1
      endif
      allocate(kfs(n_out(i_in, i_mode)))
      do i_out=1, n_out(i_in, i_mode)
        kfs(i_out)=kf_out(i_in, i_mode, i_out)
      enddo
      is_success=n_body_decay(kf0, kin0, kfs)
      deallocate(kfs)
      if(is_success) return
    enddo
    ! ----------

  end subroutine decay_with_data

!=======================================================================
end module decay_utility
