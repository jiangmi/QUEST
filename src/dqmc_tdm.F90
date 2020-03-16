module DQMC_TDM
#include "dqmc_include.h"

  ! 11/22/2017 (comments on DQMC_TDM_GetKFTold):
  ! DQMC_TDM_GetKFTold does calculates the orbital (in one unit cell) dependent
  ! FT of tdm quantities; but DQMC_TDM_GetKFT does not so be careful for new one

  ! Also, if compare SzSz(k, b1,b2,tau) and SxSx(k, b1,b2,tau) with 
  ! chizz and chixx computed in DQMC_TDM_Chi_q_orbital
  ! they are equal for b1=b2 but chizz = 2*SzSz(k, b1,b2,tau) for b1!=b2
  ! the reason is that DQMC_TDM_Chi_q_orbital uses F(class), which treats
  ! site pair (0,1) and (1,0) the same thing but DQMC_TDM_GetKFTold treats them
  ! as different orbital dependent quantities, namely 
  ! SzSz(k, b1,b2,tau) different from SzSz(k, b2,b1,tau) althought they are equal  

  ! 08/16/2015:
  ! Add features of conductivity calculation using the standard definition
  ! instead of Simone's tricks using link correlation
  ! Also add variable flagcond: specify in input file if to compute conductivity

  ! 01/04/2016:
  ! Add new FT routines, old ones consider intersite correlation for general unit cell
  ! which is not necessary. Now old FT routines only useful for computing SelfEnergy

  ! New FT routines have coincidence with currDs routines  
  ! So they can be separately computed via input file parameters

  use DQMC_UTIL
  use DQMC_STRUCT
 ! use LAPACK_MOD
 ! use BLAS_MOD
  use DQMC_GTAU

  implicit none 
  
  !
  ! This module is designed for the computation of time dependent 
  ! measurement (TDM), which requires the unequal time Green's 
  ! function Gtau (up and down).
  !
  ! Two measurements are taken in this module: Green's function (G) and Chi
  ! (chi) function.  Both are stored in a 3 dimensional array. 
  ! The first dimension is for space, the second is for time, and the 
  ! last dimension is for bins.
  !


  type tdmarray

     integer              :: n
     integer              :: nclass
     integer              :: nk
     integer              :: np

     integer, pointer     :: D(:,:)
     integer, pointer     :: F(:)

     integer, pointer     :: phase(:,:)
     complex(wp), pointer :: ftk(:,:) 
     complex(wp), pointer :: ftw(:,:) 

     real(wp),    pointer :: values(:,:,:)
     complex(wp), pointer :: valueskold(:,:,:) ! old FT, now only for selfE
     complex(wp), pointer :: valueskold_iw0(:,:)
     real(wp),    pointer :: valuesk(:,:,:,:)  ! (kx,ky,tau,bin)

     real(wp),    pointer :: tlink(:,:)

     real(wp),    pointer :: vecClass(:,:)
     character(label_len), pointer :: clabel(:)

  end type tdmarray

  ! Index of the array varaiables
  integer, parameter  :: NTDMARRAY = 11
  integer, parameter  :: IGFUN = 1
  integer, parameter  :: IGFUP = 2
  integer, parameter  :: IGFDN = 3
  integer, parameter  :: ISPXX = 4
  integer, parameter  :: ISPZZ = 5
  integer, parameter  :: IDENS = 6
  integer, parameter  :: IPAIRs = 7
  integer, parameter  :: IPAIRd = 8
  integer, parameter  :: ICOND  = 9
  integer, parameter  :: ICONDup = 10
  integer, parameter  :: ICONDdn = 11
  integer             :: ntdm = NTDMARRAY ! if needed, decrease it for not cond, d-wave sus etc.

  ! Index of the array varaiables
  character(len=9), parameter :: &
            pname(NTDMARRAY) = (/ &
                  "Gfun     ", &
                  "Gfun up  ", &
                  "Gfun dn  ", &
                  "SxSx     ", &
                  "SzSz     ", &
                  "Den-Den  ", &
                  "s-wave   ", &
                  "d-wave   ", &
                  "Cond     ", &
                  "Cond_up  ", &
                  "Cond_dn  " /)

  type TDM
     integer  :: L
     integer  :: nbin
     integer  :: avg
     integer  :: err
     integer  :: idx    !label bin

     integer  :: tmp
     integer  :: cnt

     logical  :: compute = .false.
 
     ! input file specifies if computing different correlation functions
     integer  :: flags(NTDMARRAY) 

     ! input file specifies if computing FT of correlation functions
     integer  :: flagsFT(NTDMARRAY)

     ! # of k points for FT (0:L/2,0:L/2) = 1+L/2
     integer  :: NkFT
     integer  :: norb      ! number of orbitals in unit cell

     ! Instead of using FT transformation, 
     ! manually compute chi_q(tau,2,bin) and chi_q_iw0(2,bin)
     ! for 2 special K=(0,0) and (pi,pi)
     real(wp), pointer :: chiqxx(:,:,:)    
     real(wp), pointer :: chiqzz(:,:,:)
     real(wp), pointer :: chiqxx_iw0(:,:)  
     real(wp), pointer :: chiqzz_iw0(:,:)
     ! In case of two sublattices, can also
     ! compute the chi_q=0 for two sublattices individually
     real(wp), pointer :: chiqxx_sub1(:,:)
     real(wp), pointer :: chiqzz_sub1(:,:)
     real(wp), pointer :: chiqxx_iw0_sub1(:)
     real(wp), pointer :: chiqzz_iw0_sub1(:)
     real(wp), pointer :: chiqxx_sub2(:,:)
     real(wp), pointer :: chiqzz_sub2(:,:)
     real(wp), pointer :: chiqxx_iw0_sub2(:)
     real(wp), pointer :: chiqzz_iw0_sub2(:)

     ! composite Simpson's iw=0 summation for spinxx and spinzz
     real(wp), pointer :: chixx_r_orb_iw0(:,:)  
     real(wp), pointer :: chizz_r_orb_iw0(:,:)

     real(wp) :: dtau
     real(wp), pointer :: sgn(:)
     type(tdmarray), pointer :: properties(:)

     ! record the average of all local sum_r G(r, tau) for average N(w)
     real(wp), pointer :: GtauAvg(:,:), GtupAvg(:,:), GtdnAvg(:,:)

     ! record the average of all local sum_r spin-xx(r, tau) for average spin-xx susceptibility
     real(wp), pointer :: spinxxAvg(:,:)
     real(wp), pointer :: spinzzAvg(:,:)

     ! record the average of all local sum_r pair(r, tau) for average pair susceptibility
     real(wp), pointer :: swaveAvg(:,:)

     ! correlated and uncorrelated/non-vertex d-wave susceptibility
     ! and Gamma_d = Pd^-1 - Pd0^-1
     ! final index denotes orbital; for hubbard, just 1;
     ! for PAM, Pd tensor can be abcd indices so that c,f orbitals gives 16 terms
     ! For now, only compute P_ffff components, which is expected to dominant
     integer :: NPd  ! No. of pairing susceptibility components
     real(wp), pointer :: Pd0tau(:,:,:)
     real(wp), pointer :: Pdtau(:,:,:)
     real(wp), pointer :: Pd0(:,:)
     real(wp), pointer :: Pd(:,:)
     real(wp), pointer :: Gammad(:,:)
     real(wp), pointer :: Gd_Pd0(:,:)

     ! Fourier transform matrix for bosonic and fermionic fields
     complex(wp), pointer :: ftwfer(:,:), ftwbos(:,:)

     ! curr-curr(qx=0,qy;iwn=0) is estimated by linear extrapolation of two smallest qy
     ! iwn=0 component is real,  (q,3,bin,3), first 3 denotes up, dn and total
     ! second 3: nospline, qwspline, wqspline for iwn=0
     real(wp), pointer :: Dsqx(:,:,:,:)   
     real(wp), pointer :: Dsqy(:,:,:,:)    

     ! pi* (<-Kx>-curr-curr(qx=0,qy->0;iwn=0))
     real(wp), pointer :: Ds(:,:,:)        ! (3,bin,3) 

     ! 08/15/2015
     ! used for conductivity, d-wave pairing sus etc.
     integer, ALLOCATABLE     :: rt(:), lf(:), up(:), dn(:)
     complex*16, ALLOCATABLE  :: hopup(:,:), hopdn(:,:)

     ! 07/17/2019
     ! obtain the cartesian coordinates of sites, used for in-plane d-wave pairing sus
     ! and chiq(pi,pi) manually
     real(wp), pointer :: cartpos(:,:)

     ! 11/19/2015
     ! used for self-energy, (natom, natom, L, nk, 3)
     ! natom for unit cell, nk for cluster(cell) K values, 3 for G, Gup, Gdn
     complex(wp), pointer :: GkwAvg(:,:,:,:,:), GkwErr(:,:,:,:,:)
     complex(wp), pointer :: SEkavg(:,:,:,:,:), SEkerr(:,:,:,:,:)

     ! 5/9/2019
     ! sigma(r,tau) and sigma(r,iw=0) with (natom, natom, L, 3)
     complex(wp), pointer :: GrwAvg(:,:,:,:), GrwErr(:,:,:,:)
     complex(wp), pointer :: SEravg(:,:,:,:), SErerr(:,:,:,:)

  end type TDM
 
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDM_Init(model, L, dtau, T1, nBin, S, Gwrap, flags, flagsFT)!, splinew0)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM. 
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout)   :: T1      ! time dependent measurement
    integer, intent(in)         :: L       ! No of time slice
    integer, intent(in)         :: nBin    ! No of Bins
    integer, intent(in)         :: model   ! Hubbard or PAM ...
    integer, intent(in)         :: flags(NTDMARRAY), flagsFT(NTDMARRAY)
    !character(len=2),intent(in) :: splinew0
    real(wp), intent(in)        :: dtau
    type(Struct), intent(in)    :: S
    type(GeomWrap), intent(in)  :: Gwrap

    ! ... local variables ...
    integer     :: i,j,nClass

    ! ... Executable ...

    T1%L      =  L
    T1%dtau   =  dtau
    T1%nbin   =  nBin
    nclass    =  S%nClass

    T1%tmp    =  nBin + 1
    T1%avg    =  nBin + 1
    T1%err    =  nBin + 2
    T1%idx    =  1   ! initial value, increase in DQMC_TDM_Avg

    T1%compute  = .true.
    T1%flags    = flags
    T1%flagsFT  = flagsFT

    allocate(T1%sgn(nBin+2))
    T1%sgn   = ZERO

    ! # of k points for FT, only for square lattice !!!
    T1%norb = Gwrap%lattice%natom     
    T1%NkFT = int(sqrt(real(S%nSite/T1%norb)))/2 

    allocate(T1%cartpos(3, 0:S%nSite-1))
    T1%cartpos = Gwrap%lattice%cartpos

    call  DQMC_TDM_InitFTw(T1)
    ntdm = sum(T1%flags)         ! how many quantities to compute
    allocate(T1%properties(NTDMARRAY))

    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         call DQMC_TDM_InitProp(T1, S, Gwrap, i)
       endif
       if (T1%flagsFT(i)==1) then
         call DQMC_TDM_InitPropFT(T1, Gwrap, i)
       endif
    enddo

    if (T1%flags(IGFUN) == 1) then
      allocate(T1%GtauAvg(0:T1%L-1,T1%err))
      T1%GtauAvg = 0.0_wp
    endif
    if (T1%flags(IGFUP) == 1) then
      allocate(T1%GtupAvg(0:T1%L-1,T1%err))
      T1%GtupAvg = 0.0_wp
    endif
    if (T1%flags(IGFDN) == 1) then
      allocate(T1%GtdnAvg(0:T1%L-1,T1%err))
      T1%GtdnAvg = 0.0_wp
    endif

    if (T1%flags(ISPXX) == 1) then
      allocate(T1%spinxxAvg(0:T1%L-1,T1%err))
      allocate(T1%chixx_r_orb_iw0(nClass, T1%err))
      T1%spinxxAvg = 0.0_wp
    endif
    if (T1%flagsFT(ISPXX) == 1) then
      allocate(T1%chiqxx(0:T1%L-1, 2, T1%err))
      allocate(T1%chiqxx_iw0(2, T1%err))

      ! for two sublattices individually
      allocate(T1%chiqxx_sub1(0:T1%L-1, T1%err))
      allocate(T1%chiqxx_iw0_sub1(T1%err))
      allocate(T1%chiqxx_sub2(0:T1%L-1, T1%err))
      allocate(T1%chiqxx_iw0_sub2(T1%err))
    endif

    if (T1%flags(ISPZZ) == 1) then
      allocate(T1%spinzzAvg(0:T1%L-1,T1%err))
      allocate(T1%chizz_r_orb_iw0(nClass, T1%err))
      T1%spinzzAvg = 0.0_wp
    endif
    if (T1%flagsFT(ISPZZ) == 1) then
      allocate(T1%chiqzz(0:T1%L-1, 2, T1%err))
      allocate(T1%chiqzz_iw0(2, T1%err))

      ! for two sublattices individually
      allocate(T1%chiqzz_sub1(0:T1%L-1, T1%err))
      allocate(T1%chiqzz_iw0_sub1(T1%err))
      allocate(T1%chiqzz_sub2(0:T1%L-1, T1%err))
      allocate(T1%chiqzz_iw0_sub2(T1%err))
    endif

    if (T1%flags(IPAIRs) == 1) then
      allocate(T1%swaveAvg(0:T1%L-1,T1%err))
      T1%swaveAvg = 0.0_wp
    endif

    if (T1%flags(IPAIRd) == 1) then
      ! Get No. of pairing susceptibility components
      select case (model)
        ! hubbard square
        case (0)
          T1%NPd = 1
        ! PAM square
        ! In principle, Pd tensor can be abcd indices so that c,f orbitals gives 16 terms
        ! For now, only compute P_ffff components, which is expected to dominant
        case (1)
          T1%NPd = 1
      end select

      allocate(T1%Pdtau(0:T1%L-1, 1:T1%err, T1%NPd))
      allocate(T1%Pd0tau(0:T1%L-1, 1:T1%err, T1%NPd))
      allocate(T1%Pd(1:T1%err, T1%NPd))
      allocate(T1%Pd0(1:T1%err, T1%NPd))
      allocate(T1%Gammad(1:T1%err, T1%NPd))
      allocate(T1%Gd_Pd0(1:T1%err, T1%NPd))

      T1%Pdtau  = 0.0_wp
      T1%Pd0tau = 0.0_wp
      T1%Pd  = 0.0_wp
      T1%Pd0 = 0.0_wp
      T1%Gammad = 0.0_wp
      T1%Gd_Pd0 = 0.0_wp
    endif

    ! used for conductivity, d-wave paring sus etc.
    ! IMPORTANT: the indices of Gtau(1:nsites) and hamilt(0:nsites-1) are different
    ! Note difference from rt etc. in dqmc_hamilt.F90
    ! So switch the values to from 1 to nSite, instead of 0 to nSite-1
    ! should be correct for any unit cell definition
    allocate(T1%rt(1:S%nSite))
    allocate(T1%lf(1:S%nSite))
    allocate(T1%up(1:S%nSite))
    allocate(T1%dn(1:S%nSite))
    allocate(T1%hopup(1:S%nSite,1:S%nSite))
    allocate(T1%hopdn(1:S%nSite,1:S%nSite))

    do i = 1, S%nSite
       T1%rt(i) = Gwrap%Hamilt%rt(i-1)+1
       T1%lf(i) = Gwrap%Hamilt%lf(i-1)+1
       T1%up(i) = Gwrap%Hamilt%up(i-1)+1
       T1%dn(i) = Gwrap%Hamilt%dn(i-1)+1
       do j = 1, S%nSite
          T1%hopup(i,j) = Gwrap%Hamilt%hopup(i-1,j-1)
          T1%hopdn(i,j) = Gwrap%Hamilt%hopdn(i-1,j-1)
       enddo
    enddo

  end subroutine DQMC_TDM_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_InitProp(T1, S, Gwrap, iprop)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    ! Initialize contents of tdmarray for iprop
    !
    ! Arguments
    ! =========
    type(TDM), intent(inout) :: T1
    type(Struct), intent(in)  :: S
    type(GeomWrap), intent(in):: Gwrap
    integer, intent(in)       :: iprop

    integer :: np, npp, nclass

    select case (iprop)

       case(IGFUN, IGFUP, IGFDN)

          nclass = S%nClass
          np     = Gwrap%lattice%natom

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  nclass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftw    => T1%ftwfer
          T1%properties(iprop)%phase  => S%gf_phase
          T1%properties(iprop)%clabel  => S%clabel
          T1%properties(iprop)%vecClass  => S%vecClass
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))

       case(ISPXX, ISPZZ, IDENS, IPAIRs, IPAIRd)

          nclass = S%nClass
          np     = Gwrap%lattice%natom

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  S%nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftw    => T1%ftwbos
          T1%properties(iprop)%phase  => S%chi_phase
          T1%properties(iprop)%clabel  => S%clabel
          T1%properties(iprop)%vecClass  => S%vecClass
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))

       case(ICOND, ICONDup, ICONDdn)

          ! originally nClass=1 for only q=0 components of conductivity
          ! modify nClass=1 to be same as other quantities
          nclass = S%nClass
          np     = Gwrap%lattice%natom

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%phase  => S%chi_phase
          T1%properties(iprop)%clabel  => S%clabel
          T1%properties(iprop)%vecClass  => S%vecClass
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))

        ! used in meas for average, for conductivity, renormalized by lattice size
        ! Below two lines are original for nClass=1 due to only q=0 components
        !  allocate(T1%properties(iprop)%F(nclass))
        !  T1%properties(iprop)%F(1) = S%nSite*2   ! *2 is for averaging x and y directions
     end select

     T1%properties(iprop)%values  = 0.0_wp

  end subroutine DQMC_TDM_InitProp

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_InitPropFT(T1, Gwrap, iprop)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    ! Initialize contents of tdmarray for iprop
    !
    ! Arguments
    ! =========
    type(TDM), intent(inout) :: T1
    type(GeomWrap), intent(in):: Gwrap
    integer, intent(in)       :: iprop

    integer :: nk, np, npp

    select case (iprop)

       case(IGFUN, IGFUP, IGFDN)

          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2
          nk = Gwrap%RecipLattice%nclass_k
          T1%properties(iprop)%nk     =  nk
          T1%properties(iprop)%ftk    => Gwrap%RecipLattice%FourierC
          allocate(T1%properties(iprop)%valueskold(nk*npp,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valueskold_iw0(nk*npp,T1%err))
          allocate(T1%properties(iprop)%valuesk(0:T1%NkFT,0:T1%NkFT,0:T1%L-1,T1%err))

       case(ISPXX, ISPZZ, IDENS, IPAIRs, IPAIRd)

          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2
          nk = Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%ftk    => Gwrap%GammaLattice%FourierC
          allocate(T1%properties(iprop)%valueskold(nk*npp,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valueskold_iw0(nk*npp,T1%err))
          allocate(T1%properties(iprop)%valuesk(0:T1%NkFT,0:T1%NkFT,0:T1%L-1,T1%err))

       case(ICOND, ICONDup, ICONDdn)

          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2
          nk = Gwrap%GammaLattice%nclass_k
          ! originally nClass=1 for only q=0 components of conductivity
          ! modify nClass=1 to be same as other quantities
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%ftk    => Gwrap%GammaLattice%FourierC
          allocate(T1%properties(iprop)%valueskold(nk*npp,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valueskold_iw0(nk*npp,T1%err))
          allocate(T1%properties(iprop)%valuesk(0:T1%NkFT,0:T1%NkFT,0:T1%L-1,T1%err))

        ! used in meas for average, for conductivity, renormalized by lattice
        ! size
        ! Below two lines are original for nClass=1 due to only q=0 components
        !  allocate(T1%properties(iprop)%F(nclass))
        !  T1%properties(iprop)%F(1) = S%nSite*2   ! *2 is for averaging x and y
        !  directions
     end select

     if(associated(T1%properties(iprop)%valueskold)) &
         T1%properties(iprop)%valueskold = 0.0_wp
         T1%properties(iprop)%valueskold_iw0 = 0.0_wp
     if(associated(T1%properties(iprop)%valuesk)) &
         T1%properties(iprop)%valuesk = 0.0_wp

  end subroutine DQMC_TDM_InitPropFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_InitFTw(T1)
    !
    ! Purpose
    ! =======
    ! Fill the matricess of fourier coefficients necessary to transform
    ! from imaginary time to imaginary frequency. Two matrices are
    ! needed: ftwbos for bosonic fields (spin, charge...) and 
    ! ftwfer for fermionic fields (particles)
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: T1      ! TDM to be freed

    integer :: iw, it, L
    real(wp) :: x, pi

    ! ... Executable ...

    if (.not.T1%compute) return

    L = T1%L
    pi = acos(-1.0_wp)

    allocate(T1%ftwbos(0:L-1,0:L-1))
    allocate(T1%ftwfer(0:L-1,0:L-1))

    do iw = 0, L-1
       do it = 0, L-1
          x = 2*it*iw*pi/L
          T1%ftwbos(iw,it) = T1%dtau * cmplx(cos(x),sin(x))
       enddo
    enddo

    do iw = 0, L-1
       do it = 0, L-1
          x = 2*it*(iw+0.5_wp)*pi/L
          T1%ftwfer(iw,it) = T1%dtau * cmplx(cos(x),sin(x))
       enddo
    enddo

  end subroutine DQMC_TDM_InitFTw

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Free(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM.
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: T1

    integer :: i

    ! ... Executable ...

    if (.not.T1%compute) return
    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         deallocate(T1%properties(i)%values)
         nullify(T1%properties(i)%D)
         nullify(T1%properties(i)%F)
         nullify(T1%properties(i)%ftw)
       endif
       if (T1%flagsFT(i)==1) then
         nullify(T1%properties(i)%ftk)
         deallocate(T1%properties(i)%valueskold)
         deallocate(T1%properties(i)%valueskold_iw0)
         deallocate(T1%properties(i)%valuesk)
       endif
    enddo

    deallocate(T1%ftwbos)
    deallocate(T1%ftwfer)
    deallocate(T1%properties)

    if (T1%flags(IGFUN) == 1) then
      deallocate(T1%GtauAvg)
    endif
    if (T1%flags(IGFUP) == 1) then
      deallocate(T1%GtupAvg)
    endif
    if (T1%flags(IGFDN) == 1) then
      deallocate(T1%GtdnAvg)
    endif
    if (T1%flags(ISPXX) == 1) then
      deallocate(T1%spinxxAvg)
      deallocate(T1%chixx_r_orb_iw0)
    endif
    if (T1%flagsFT(ISPXX) == 1) then
      deallocate(T1%chiqxx)
      deallocate(T1%chiqxx_iw0)
      ! for two sublattices individually
      deallocate(T1%chiqxx_sub1)
      deallocate(T1%chiqxx_iw0_sub1)
      deallocate(T1%chiqxx_sub2)
      deallocate(T1%chiqxx_iw0_sub2)
    endif
    if (T1%flags(ISPZZ) == 1) then
      deallocate(T1%spinzzAvg)
      deallocate(T1%chizz_r_orb_iw0)
    endif
    if (T1%flagsFT(ISPZZ) == 1) then
      deallocate(T1%chiqzz)
      deallocate(T1%chiqzz_iw0)
      ! for two sublattices individually
      deallocate(T1%chiqzz_sub1)
      deallocate(T1%chiqzz_iw0_sub1)
      deallocate(T1%chiqzz_sub2)
      deallocate(T1%chiqzz_iw0_sub2)
    endif
    if (T1%flags(IPAIRs) == 1) then
      deallocate(T1%swaveAvg)
    endif
    if (T1%flags(IPAIRd) == 1) then
      deallocate(T1%Pdtau)
      deallocate(T1%Pd0tau)
      deallocate(T1%Pd)
      deallocate(T1%Pd0)
      deallocate(T1%Gammad)
      deallocate(T1%Gd_Pd0)
    endif

    deallocate(T1%rt)
    deallocate(T1%lf)
    deallocate(T1%up)
    deallocate(T1%dn)
    deallocate(T1%hopup)
    deallocate(T1%hopdn)

  end subroutine DQMC_TDM_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Meas(T1, tau, model)
    !
    ! Purpose
    ! =======
    !    This subroutine fills the bin. It assumes that tau%A_up and,
    !    when necessary, tau%A_dn have been filled.
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout)    :: t1
    type(Gtau), intent(inout)   :: tau
    
    ! ... Local var ...
    integer  :: i, k, m, L, cnt, dt, i0, it, j0, jt, dtau, iprop
    integer, intent(in)  :: model
    real(wp) :: sgn, factor, chiv
    real(wp),pointer :: up0t(:,:), upt0(:,:), dn0t(:,:), dnt0(:,:)
    real(wp),pointer :: up00(:,:), uptt(:,:), dn00(:,:), dntt(:,:)
    real(wp), pointer :: values(:,:,:)

    if (.not.T1%compute) return
    ! ... executable ...
    cnt = 0
    L   = tau%L
    k   = mod(tau%north,2)
    m   = (tau%north-k) / 2

    upt0 => tau%upt0
    up0t => tau%up0t
    dnt0 => tau%dnt0
    dn0t => tau%dn0t
    up00 => tau%up00
    uptt => tau%uptt
    dn00 => tau%dn00
    dntt => tau%dntt

    blocks: do i0 = 1, tau%nb
       do dtau = 0, tau%nb-1
          it = mod(i0+dtau-1,tau%nb) + 1

          ! Stored value
          call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
          if (tau%comp_dn .or. .not.tau%neg_u) &
             call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          
          jt = tau%it_up; j0 = tau%i0_up

          call DQMC_TDM_Compute(T1, model, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM_Compute(T1, model, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
          enddo

          if (m .gt. 0) then
             call DQMC_Gtau_DumpA(tau, TAU_UP, it, i0)
             if (tau%comp_dn .or. .not.tau%neg_u) &
                call DQMC_Gtau_DumpA(tau, TAU_DN, it, i0)
          endif
          ! Increment index tau%it
          do dt = 1, m
             call DQMC_change_gtau_time(tau, TMINUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TMINUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif
             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM_Compute(T1, model, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
                
          enddo

       enddo

       cnt = cnt + 1
    enddo blocks
 
    if (i0 .ne. tau%nb+1) then
       write(*,*) "Up and down time slices are mismatched. Stop"
       stop
    endif

    sgn = tau%sgnup * tau%sgndn

    do iprop = 1, NTDMARRAY
       if (T1%flags(iprop)==1) then       
          values => T1%properties(iprop)%values
          do i = 1, T1%properties(iprop)%nClass
            factor = sgn/(T1%properties(iprop)%F(i)*cnt)
            values(i,0:L-1,T1%idx)   = values(i,0:L-1,T1%idx)   + factor*values(i,0:L-1,T1%tmp)

            !compute chixx_r0_iw0 with composite Simpson's rule:
            if (iprop==ISPXX) then
              chiv = 0.0
              call convert_to_iw0_real(values(i,0:L-1,T1%tmp), chiv, L, T1%dtau)
              T1%chixx_r_orb_iw0(i,T1%idx) = T1%chixx_r_orb_iw0(i,T1%idx) + factor*chiv  
            endif

            !compute chizz_r0_iw0 with composite Simpson's rule:
            if (iprop==ISPZZ) then
              chiv = 0.0
              call convert_to_iw0_real(values(i,0:L-1,T1%tmp), chiv, L, T1%dtau)
              T1%chizz_r_orb_iw0(i,T1%idx) = T1%chizz_r_orb_iw0(i,T1%idx) + factor*chiv  
            endif
          end do
          values(:,:,T1%tmp) = ZERO
       endif
    enddo

    ! compute chixx_q of f-electron at q=(0,0) and (pi,pi)
    if (T1%flagsFT(ISPXX) == 1) then
      factor = sgn/cnt
      T1%chiqxx(0:L-1,1:2,T1%idx) = T1%chiqxx(0:L-1,1:2,T1%idx) + factor*T1%chiqxx(0:L-1,1:2,T1%tmp)
      T1%chiqxx_sub1(0:L-1,T1%idx) = T1%chiqxx_sub1(0:L-1,T1%idx) + factor*T1%chiqxx_sub1(0:L-1,T1%tmp)
      T1%chiqxx_sub2(0:L-1,T1%idx) = T1%chiqxx_sub2(0:L-1,T1%idx) + factor*T1%chiqxx_sub2(0:L-1,T1%tmp)

      ! compute chi(iw=0) static chi using Composite Simpson's rule
      do i = 1,2
        chiv = 0.0
        call convert_to_iw0_real(T1%chiqxx(0:L-1,i,T1%tmp), chiv, L, T1%dtau)
        T1%chiqxx_iw0(i,T1%idx) = T1%chiqxx_iw0(i,T1%idx) + factor*chiv
      enddo

      chiv = 0.0
      call convert_to_iw0_real(T1%chiqxx_sub1(0:L-1,T1%tmp), chiv, L, T1%dtau)
      T1%chiqxx_iw0_sub1(T1%idx) = T1%chiqxx_iw0_sub1(T1%idx) + factor*chiv
      chiv = 0.0
      call convert_to_iw0_real(T1%chiqxx_sub2(0:L-1,T1%tmp), chiv, L, T1%dtau)
      T1%chiqxx_iw0_sub2(T1%idx) = T1%chiqxx_iw0_sub2(T1%idx) + factor*chiv

      T1%chiqxx(:,:,T1%tmp) = ZERO
      T1%chiqxx_sub1(:,T1%tmp) = ZERO
      T1%chiqxx_sub2(:,T1%tmp) = ZERO
      T1%chiqxx_iw0(:,T1%tmp) = ZERO
      T1%chiqxx_iw0_sub1(T1%tmp) = ZERO
      T1%chiqxx_iw0_sub2(T1%tmp) = ZERO
    endif

    ! compute chizz_q of f-electron at q=(0,0) and (pi,pi)
    if (T1%flagsFT(ISPZZ) == 1) then
      factor = sgn/cnt
      T1%chiqzz(0:L-1,1:2,T1%idx) = T1%chiqzz(0:L-1,1:2,T1%idx) + factor*T1%chiqzz(0:L-1,1:2,T1%tmp)
      T1%chiqzz_sub1(0:L-1,T1%idx) = T1%chiqzz_sub1(0:L-1,T1%idx) + factor*T1%chiqzz_sub1(0:L-1,T1%tmp)
      T1%chiqzz_sub2(0:L-1,T1%idx) = T1%chiqzz_sub2(0:L-1,T1%idx) + factor*T1%chiqzz_sub2(0:L-1,T1%tmp)

      ! compute chi(iw=0) static chi using Composite Simpson's rule
      do i = 1,2
        chiv = 0.0
        call convert_to_iw0_real(T1%chiqzz(0:L-1,i,T1%tmp), chiv, L, T1%dtau)
        T1%chiqzz_iw0(i,T1%idx) = T1%chiqzz_iw0(i,T1%idx) + factor*chiv
      enddo

      chiv = 0.0
      call convert_to_iw0_real(T1%chiqzz_sub1(0:L-1,T1%tmp), chiv, L, T1%dtau)
      T1%chiqzz_iw0_sub1(T1%idx) = T1%chiqzz_iw0_sub1(T1%idx) + factor*chiv
      chiv = 0.0
      call convert_to_iw0_real(T1%chiqzz_sub2(0:L-1,T1%tmp), chiv, L, T1%dtau)
      T1%chiqzz_iw0_sub2(T1%idx) = T1%chiqzz_iw0_sub2(T1%idx) + factor*chiv

      T1%chiqzz(:,:,T1%tmp) = ZERO
      T1%chiqzz_sub1(:,T1%tmp) = ZERO
      T1%chiqzz_sub2(:,T1%tmp) = ZERO
      T1%chiqzz_iw0(:,T1%tmp) = ZERO
      T1%chiqzz_iw0_sub1(T1%tmp) = ZERO
      T1%chiqzz_iw0_sub2(T1%tmp) = ZERO
    endif

    if (T1%flags(IPAIRd) == 1) then
      factor = sgn/(T1%properties(IPAIRd)%n*cnt)
      T1%Pdtau(0:T1%L-1, T1%idx, 1:T1%NPd) = T1%Pdtau(0:T1%L-1, T1%idx, 1:T1%NPd) &
                                    + T1%Pdtau(0:T1%L-1, T1%tmp, 1:T1%NPd) * factor

      do i = 1,T1%NPd
        chiv = 0.0
        call convert_to_iw0_real(T1%Pdtau(0:T1%L-1, T1%tmp, i), chiv, T1%L, T1%dtau)
        T1%Pd(T1%idx, i) = T1%Pd(T1%idx, i) + factor*chiv
      enddo

      T1%Pdtau(:, T1%tmp, :) = ZERO
      T1%Pd(T1%tmp, :) = ZERO
    endif

    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

  end subroutine DQMC_TDM_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Compute(T1, model, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, it, i0)
    !
    ! Purpose
    ! =======
    !    This subroutine assembles the time dependent properties
    !    starting from the 1-body Green's function
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout)    :: T1
    real(wp), intent(in)        :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)        :: dnt0(:,:), dn0t(:,:)
    real(wp), intent(in)        :: up00(:,:), uptt(:,:)
    real(wp), intent(in)        :: dn00(:,:), dntt(:,:)
    integer, intent(in)         :: it, i0, model
 
    ! ... Local scalar ...

    character(label_len) :: label
    integer  :: i, j, k, dt, dt1, dt2, b1, b2
    real*8   :: a,b,c,d,x,y,z
    real*8   :: x1,y1,z1,x2,y2,z2
    real(wp), pointer :: value1(:), value2(:)
    real(wp) :: factor
    real(wp), dimension(1:5) :: vec
    real(wp) :: val1, val2

    ! ... Executable ...
    if (.not.T1%compute) return

    ! Below dt1 and dt2's length switch for two cases
    ! 2*factor=0.5 for double counting of site loops
    ! factor=0.25 in Gfun for average over up and dn,
    ! which only applies for it != i0, 
    dt = it - i0
    if (dt .gt. 0) then
       ! it > i0
       dt1  =  dt
       dt2  =  T1%L - dt
       factor  = 0.25d0
    elseif (dt .lt. 0) then
       ! it < i0
       dt1  =  dt + T1%L
       dt2  = -dt
       factor = -0.25d0
    else
       dt1 = 0
       dt2 = 0
       factor = 0.5d0
    endif

    if (dt .ne. 0) then
       ! value2:
       ! for dt>0: G(dt2,0) = G(beta-dt1,0) = -G(-dt1,0) = -G(0,dt1)
       !     dt<0: G(dt2,0) = G(-dt,0) = G(beta-dt1,0) = -G(-dt1,0) = -G(0,dt1)
       ! So rules for value1 --> value2:
       ! upt0 <--> -up0t, dnt0 <--> -dn0t
       ! uptt <--> up00, dntt <--> dn00
       ! If G's appear as a pair product, do not need since the product 
       ! remains the same even switching it. 
       ! In this product case, switch the spatial sites instead and average by 2
       ! below *0.5_wp is for averaging over switching sites i,j

     if (T1%flags(IGFUN) == 1) then
       value1  => T1%properties(IGFUN)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IGFUN)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(IGFUN)%n
          do j = 1, T1%properties(IGFUN)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IGFUN)%D(i,j)
             value1(k)  = value1(k) + factor*(upt0(i,j) + dnt0(i,j))
             value2(k)  = value2(k) - factor*(up0t(i,j) + dn0t(i,j))
          end do
       end do
     endif

     if (T1%flags(IGFUP) == 1) then
       value1  => T1%properties(IGFUP)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IGFUP)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(IGFUP)%n
          do j = 1, T1%properties(IGFUP)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IGFUP)%D(i,j)
             value1(k)  = value1(k) + 2*factor*upt0(i,j)
             value2(k)  = value2(k) - 2*factor*up0t(i,j)
          end do
       end do
      endif

     if (T1%flags(IGFDN) == 1) then
       value1  => T1%properties(IGFDN)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IGFDN)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(IGFDN)%n
          do j = 1, T1%properties(IGFDN)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IGFDN)%D(i,j)
             value1(k)  = value1(k) + 2*factor*dnt0(i,j)
             value2(k)  = value2(k) - 2*factor*dn0t(i,j)
          end do
       end do
     endif

     if (T1%flags(ISPXX) == 1) then
       value1  => T1%properties(ISPXX)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ISPXX)%values(:, dt2, T1%tmp)
       do i = 1,  T1%properties(ISPXX)%n
          do j = 1,  T1%properties(ISPXX)%n
             ! k is the distance index of site i and site j
             k = T1%properties(ISPXX)%D(i,j)
             ! SxSx = (Si+ * Sj + Sj+ * Si)
             ! Si+ = c^+_i,up * c_i,dn

             val1 = up0t(j,i)*dnt0(i,j)
             val2 = up0t(i,j)*dnt0(j,i)
             value1(k)  = value1(k) - val1 
             value2(k)  = value2(k) - val2 
           !  value1(k)  = value1(k) - (up0t(j,i)*dnt0(i,j) &
           !       + up0t(i,j)*dnt0(j,i))/2
           !  value2(k)  = value2(k) - (up0t(i,j)*dnt0(j,i) &
           !       + up0t(j,i)*dnt0(i,j))/2

             ! compute chi_q of f-electron at q=(0,0) and (pi,pi)
             if (T1%flagsFT(ISPXX) == 1) then
               ! get the cartesian coordinates of site i
               ! Note i-1 accounts for the different convention of labelling sites  
               x1 = T1%cartpos(1,i-1)
               y1 = T1%cartpos(2,i-1)
               z1 = T1%cartpos(3,i-1)
               x2 = T1%cartpos(1,j-1)
               y2 = T1%cartpos(2,j-1)
               z2 = T1%cartpos(3,j-1)

               ! chi_ff calculation
               if (abs(z1-1.)<1.e-4 .and. abs(z2-1.)<1.e-4) then
                 ! 1) sum both inter- and intra-site even for multi-site unit cell
                 ! q=(0,0)
                 T1%chiqxx(dt1,1,T1%tmp) = T1%chiqxx(dt1,1,T1%tmp) - val1
                 T1%chiqxx(dt2,1,T1%tmp) = T1%chiqxx(dt2,1,T1%tmp) - val2

                 ! q=(pi,pi), note the sign accouting for exp(separation*pi)
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   T1%chiqxx(dt1,2,T1%tmp) = T1%chiqxx(dt1,2,T1%tmp) - val1
                   T1%chiqxx(dt2,2,T1%tmp) = T1%chiqxx(dt2,2,T1%tmp) - val2
                 else
                   T1%chiqxx(dt1,2,T1%tmp) = T1%chiqxx(dt1,2,T1%tmp) + val1
                   T1%chiqxx(dt2,2,T1%tmp) = T1%chiqxx(dt2,2,T1%tmp) + val2
                 endif

                 ! 2) bipartite lattice can compute q=0 individually for each sublat
                 ! separation limitation:
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   ! one sublattice
                   if (mod(int(x1+y1),2)==0) then
                     T1%chiqxx_sub1(dt1,T1%tmp) = T1%chiqxx_sub1(dt1,T1%tmp) - val1
                     T1%chiqxx_sub1(dt2,T1%tmp) = T1%chiqxx_sub1(dt2,T1%tmp) - val2

                   ! the other sublattice
                   else
                     T1%chiqxx_sub2(dt1,T1%tmp) = T1%chiqxx_sub2(dt1,T1%tmp) - val1
                     T1%chiqxx_sub2(dt2,T1%tmp) = T1%chiqxx_sub2(dt2,T1%tmp) - val2
                   endif
                 endif
               endif
             endif
          end do
       end do
     endif

     if (T1%flags(ISPZZ) == 1) then
       value1  => T1%properties(ISPZZ)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ISPZZ)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(ISPZZ)%n
          do j = 1, T1%properties(ISPZZ)%n
            ! k is the distance index of site i and site j
             k = T1%properties(ISPZZ)%D(i,j)

             val1 = (up0t(j,i)*upt0(i,j) + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )*0.5_wp
             val2 = (up0t(i,j)*upt0(j,i) + dn0t(i,j)*dnt0(j,i) - (uptt(j,j)-dntt(j,j))*(up00(i,i)-dn00(i,i)) )*0.5_wp
             value1(k)  = value1(k) - val1
             value2(k)  = value2(k) - val2

             ! compute chi_q of f-electron at q=(0,0) and (pi,pi)
             if (T1%flagsFT(ISPZZ) == 1) then
               ! get the cartesian coordinates of site i
               ! Note i-1 accounts for the different convention of labelling sites  
               x1 = T1%cartpos(1,i-1)
               y1 = T1%cartpos(2,i-1)
               z1 = T1%cartpos(3,i-1)
               x2 = T1%cartpos(1,j-1)
               y2 = T1%cartpos(2,j-1)
               z2 = T1%cartpos(3,j-1)

               ! chi_ff calculation
               if (abs(z1-1.)<1.e-4 .and. abs(z2-1.)<1.e-4) then
                 ! 1) sum both inter- and intra-site even for multi-site unit
                 ! cell
                 ! q=(0,0)
                 T1%chiqzz(dt1,1,T1%tmp) = T1%chiqzz(dt1,1,T1%tmp) - val1
                 T1%chiqzz(dt2,1,T1%tmp) = T1%chiqzz(dt2,1,T1%tmp) - val2

                 ! q=(pi,pi), note the sign accouting for exp(separation*pi)
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   T1%chiqzz(dt1,2,T1%tmp) = T1%chiqzz(dt1,2,T1%tmp) - val1
                   T1%chiqzz(dt2,2,T1%tmp) = T1%chiqzz(dt2,2,T1%tmp) - val2
                 else
                   T1%chiqzz(dt1,2,T1%tmp) = T1%chiqzz(dt1,2,T1%tmp) + val1
                   T1%chiqzz(dt2,2,T1%tmp) = T1%chiqzz(dt2,2,T1%tmp) + val2
                 endif

                 ! 2) bipartite lattice can compute q=0 individually for each
                 ! sublat
                 ! separation limitation:
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   ! one sublattice
                   if (mod(int(x1+y1),2)==0) then
                     T1%chiqzz_sub1(dt1,T1%tmp) = T1%chiqzz_sub1(dt1,T1%tmp) - val1
                     T1%chiqzz_sub1(dt2,T1%tmp) = T1%chiqzz_sub1(dt2,T1%tmp) - val2

                   ! the other sublattice
                   else
                     T1%chiqzz_sub2(dt1,T1%tmp) = T1%chiqzz_sub2(dt1,T1%tmp) - val1
                     T1%chiqzz_sub2(dt2,T1%tmp) = T1%chiqzz_sub2(dt2,T1%tmp) - val2
                   endif
                 endif
               endif
             endif
          end do
       end do
     endif

     if (T1%flags(IDENS) == 1) then
       value1  => T1%properties(IDENS)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IDENS)%values(:, dt2, T1%tmp)
       do i = 1, T1%properties(IDENS)%n
          do j = 1, T1%properties(IDENS)%n
            ! k is the distance index of site i and site j
             k = T1%properties(IDENS)%D(i,j)
             value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
               + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)+dntt(i,i))*(up00(j,j)+dn00(j,j)) )*0.5_wp
             value2(k)  = value2(k) - (up0t(i,j)*upt0(j,i) &
               + dn0t(i,j)*dnt0(j,i) - (uptt(j,j)+dntt(j,j))*(up00(i,i)+dn00(i,i)) )*0.5_wp
          end do
       end do
     endif

     if (T1%flags(IPAIRs) == 1) then
       value1  => T1%properties(IPAIRs)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IPAIRs)%values(:, dt2, T1%tmp)
       do i = 1,  T1%properties(IPAIRs)%n
          do j = 1,  T1%properties(IPAIRs)%n
             ! Delta^+_i = c^+_i,up * c^+_i,dn 
             ! So rules for value1 --> value2:
             ! upt0 <--> -up0t, dnt0 <--> -dn0t
             k = T1%properties(IPAIRs)%D(i,j)
             value1(k)  = value1(k) + upt0(i,j)*dnt0(i,j) *0.5_wp 
             value2(k)  = value2(k) + up0t(i,j)*dn0t(i,j) *0.5_wp
          end do
       end do
     endif

     if (T1%flags(IPAIRd) == 1) then
       value1  => T1%properties(IPAIRd)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IPAIRd)%values(:, dt2, T1%tmp)

       ! Following the general def of P_abcd in WeiWu's PRX (2015)
       ! See also RTS's paper in 1989 PRB
       ! D_i*D_j = sum_{dd'} Gup(i,j)*Gdn(i+d,j+d')
       do i = 1,  T1%properties(IPAIRd)%n
          do j = 1,  T1%properties(IPAIRd)%n
             ! Delta^+_i = c^+_i,up * c^+_i,dn 
             ! So rules for value1 --> value2:
             ! upt0 <--> -up0t, dnt0 <--> -dn0t
             k = T1%properties(IPAIRd)%D(i,j)
             vec = T1%properties(IPAIRd)%vecClass(k,:)

             ! get the cartesian coordinates of site i
             ! Note i-1 accounts for the different convention of labelling sites  
             z1 = T1%cartpos(3,i-1)
             z2 = T1%cartpos(3,j-1)

             ! 4 neighbors of each site so that totally 16 terms
             ! with different phase factors for d-wave pairing
             ! record all 16 possible Gdn(i+d,j+d') as a below
             a = dnt0(T1%rt(i), T1%rt(j)) - dnt0(T1%rt(i), T1%up(j))  &
                +dnt0(T1%rt(i), T1%lf(j)) - dnt0(T1%rt(i), T1%dn(j))  &
                -dnt0(T1%up(i), T1%rt(j)) + dnt0(T1%up(i), T1%up(j))  &
                -dnt0(T1%up(i), T1%lf(j)) + dnt0(T1%up(i), T1%dn(j))  &
                +dnt0(T1%lf(i), T1%rt(j)) - dnt0(T1%lf(i), T1%up(j))  &
                +dnt0(T1%lf(i), T1%lf(j)) - dnt0(T1%lf(i), T1%dn(j))  &
                -dnt0(T1%dn(i), T1%rt(j)) + dnt0(T1%dn(i), T1%up(j))  &
                -dnt0(T1%dn(i), T1%lf(j)) + dnt0(T1%dn(i), T1%dn(j))  

             b = dn0t(T1%rt(i), T1%rt(j)) - dn0t(T1%rt(i), T1%up(j))  &
                +dn0t(T1%rt(i), T1%lf(j)) - dn0t(T1%rt(i), T1%dn(j))  &
                -dn0t(T1%up(i), T1%rt(j)) + dn0t(T1%up(i), T1%up(j))  &
                -dn0t(T1%up(i), T1%lf(j)) + dn0t(T1%up(i), T1%dn(j))  &
                +dn0t(T1%lf(i), T1%rt(j)) - dn0t(T1%lf(i), T1%up(j))  &
                +dn0t(T1%lf(i), T1%lf(j)) - dn0t(T1%lf(i), T1%dn(j))  &
                -dn0t(T1%dn(i), T1%rt(j)) + dn0t(T1%dn(i), T1%up(j))  &
                -dn0t(T1%dn(i), T1%lf(j)) + dn0t(T1%dn(i), T1%dn(j))

             ! *0.25 or /4 accounts for the convention for definition
             ! See 1989 PRB paper: Numerical study of 2D Hubbard model
             a = a*0.5_wp*0.25
             b = b*0.5_wp*0.25
             value1(k)  = value1(k) + upt0(i,j)*a
             value2(k)  = value2(k) + up0t(i,j)*b

             ! only compute in-plane d-wave pairing susceptibility (square lattice)
             select case (model)
               ! Hubbard
               case (0)
                 T1%Pdtau(dt1, T1%tmp, 1) = T1%Pdtau(dt1, T1%tmp, 1) + upt0(i,j)*a
                 T1%Pdtau(dt2, T1%tmp, 1) = T1%Pdtau(dt2, T1%tmp, 1) + up0t(i,j)*b

               ! PAM (P_ffff only temporarily)
               case (1)
                 if (abs(z1-1.d0)<1.d-6 .and. abs(z2-1.d0)<1.d-6) then
                   T1%Pdtau(dt1, T1%tmp, 1) = T1%Pdtau(dt1, T1%tmp, 1) + upt0(i,j)*a
                   T1%Pdtau(dt2, T1%tmp, 1) = T1%Pdtau(dt2, T1%tmp, 1) + up0t(i,j)*b
                 endif
             end select
          end do
       end do
     endif

     if (T1%flags(ICOND) == 1) then
       ! J-J correlation is not following Simone's trick of link correlation here
       ! here use standard definition <sum_ij jx(i,tau)*jx(j,0)>
       ! Only for q=0 component
       ! IMPORTANT: the indices of Gtau(1:nsites) and T1ilt(0:nsites-1) are different 
       ! e.g. T1%hopup(i-1,T1%rt(i-1))

       ! Dec.9, 2015
       ! It is interesting to compute spin-dependent conductivity
       ! See RTS email Oct.20, 2013 for further projects

       ! First compute cond from spin up
       value1  => T1%properties(ICONDup)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ICONDup)%values(:, dt2, T1%tmp)

       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n

            k = T1%properties(ICOND)%D(i,j)

            a = T1%hopup(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            c = T1%hopup(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))

            ! up*up terms (note two ways for contraction!)           
            value1(k) = value1(k) - a* &
                        (-upt0(i,T1%rt(j))*up0t(j,T1%rt(i)) - upt0(T1%rt(i),j)*up0t(T1%rt(j),i) &
                         +upt0(i,j)*up0t(T1%rt(j),T1%rt(i)) + upt0(T1%rt(i),T1%rt(j))*up0t(j,i))*0.5_wp
            value1(k) = value1(k) - a* &
                        (uptt(i,T1%rt(i)) - uptt(T1%rt(i),i))*(up00(j,T1%rt(j)) - up00(T1%rt(j),j))*0.5_wp
            ! up*dn terms
            value1(k) = value1(k) - c* &
                        (uptt(i,T1%rt(i)) - uptt(T1%rt(i),i))*(dn00(j,T1%rt(j)) - dn00(T1%rt(j),j))*0.5_wp

            ! Below assuming isotropic system
            ! for averaging over x and y directions
            ! for anisotropic systems, need to modify this part
!            a = T1%hopup(i,T1%up(i))*T1%hopup(j,T1%up(j))
!            c = T1%hopup(i,T1%up(i))*T1%hopdn(j,T1%up(j))

            ! up*up terms (note two ways for contraction!)
!            value1(k) = value1(k) - a* &
!                        (-upt0(i,T1%up(j))*up0t(j,T1%up(i)) - upt0(T1%up(i),j)*up0t(T1%up(j),i) &
!                         +upt0(i,j)*up0t(T1%up(j),T1%up(i)) + upt0(T1%up(i),T1%up(j))*up0t(j,i))*0.5_wp
!            value1(k) = value1(k) - a* &
!                        ( uptt(i,T1%up(i))*up00(j,T1%up(j)) - uptt(T1%up(i),i)*up00(T1%up(j),j) &
!                         -uptt(i,T1%up(i))*up00(T1%up(j),j) + uptt(T1%up(i),i)*up00(j,T1%up(j)))*0.5_wp
            ! up*dn terms
!            value1(k) = value1(k) - c* &
!                        ( uptt(i,T1%up(i))*dn00(j,T1%up(j)) + uptt(T1%up(i),i)*dn00(T1%up(j),j) &
!                         -uptt(i,T1%up(i))*dn00(T1%up(j),j) - uptt(T1%up(i),i)*dn00(j,T1%up(j)))*0.5_wp
  
            ! value2: see the rules at the beginning of routine
            a = T1%hopup(j,T1%rt(j))*T1%hopup(i,T1%rt(i))
            c = T1%hopup(j,T1%rt(j))*T1%hopdn(i,T1%rt(i))
            ! up*up terms
            value2(k) = value2(k) - a* &
                        (-upt0(j,T1%rt(i))*up0t(i,T1%rt(j)) - upt0(T1%rt(j),i)*up0t(T1%rt(i),j) &
                         +upt0(j,i)*up0t(T1%rt(i),T1%rt(j)) + upt0(T1%rt(j),T1%rt(i))*up0t(i,j))*0.5_wp
            value2(k) = value2(k) - a* &
                        (uptt(j,T1%rt(j)) - uptt(T1%rt(j),j))*(up00(i,T1%rt(i)) - up00(T1%rt(i),i))*0.5_wp
            ! up*dn terms
            value2(k) = value2(k) - c* &
                        (uptt(j,T1%rt(j)) - uptt(T1%rt(j),j))*(dn00(i,T1%rt(i)) - dn00(T1%rt(i),i))*0.5_wp

!            a = T1%hopup(j,T1%up(j))*T1%hopup(i,T1%up(i))
!            c = T1%hopup(j,T1%up(j))*T1%hopdn(i,T1%up(i))
            ! up*up terms
!            value2(k) = value2(k) - a* &
!                        (-upt0(j,T1%up(i))*up0t(i,T1%up(j)) - upt0(T1%up(j),i)*up0t(T1%up(i),j) &
!                         +upt0(j,i)*up0t(T1%up(i),T1%up(j)) + upt0(T1%up(j),T1%up(i))*up0t(i,j))*0.5_wp
!            value2(k) = value2(k) - a* &
!                        ( uptt(j,T1%up(j))*up00(i,T1%up(i)) - uptt(T1%up(j),j)*up00(T1%up(i),i) &
!                         -uptt(j,T1%up(j))*up00(T1%up(i),i) + uptt(T1%up(j),j)*up00(i,T1%up(i)))*0.5_wp
            ! up*dn terms
!            value2(k) = value2(k) - c* &
!                        ( uptt(j,T1%up(j))*dn00(i,T1%up(i)) + uptt(T1%up(j),j)*dn00(T1%up(i),i) &
!                         -uptt(j,T1%up(j))*dn00(T1%up(i),i) - uptt(T1%up(j),j)*dn00(i,T1%up(i)))*0.5_wp
          end do
       end do

      ! Compute cond for down spin
       value1  => T1%properties(ICONDdn)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ICONDdn)%values(:, dt2, T1%tmp)

       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n

            k = T1%properties(ICOND)%D(i,j)

            b = T1%hopdn(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            d = T1%hopdn(i,T1%rt(i))*T1%hopup(j,T1%rt(j))

            ! dn*dn terms (note two ways for contraction!)
            value1(k) = value1(k) - b* &
                        (-dnt0(i,T1%rt(j))*dn0t(j,T1%rt(i)) - dnt0(T1%rt(i),j)*dn0t(T1%rt(j),i) &
                         +dnt0(i,j)*dn0t(T1%rt(j),T1%rt(i)) + dnt0(T1%rt(i),T1%rt(j))*dn0t(j,i))*0.5_wp
            value1(k) = value1(k) - b* &
                        (dntt(i,T1%rt(i)) - dntt(T1%rt(i),i))*(dn00(j,T1%rt(j)) - dn00(T1%rt(j),j))*0.5_wp
            ! dn*up terms
            value1(k) = value1(k) - d* &
                        (dntt(i,T1%rt(i)) - dntt(T1%rt(i),i))*(up00(j,T1%rt(j)) - up00(T1%rt(j),j))*0.5_wp

            ! Below assuming isotropic system
            ! for averaging over x and y directions
            ! for anisotropic systems, need to modify this part
!            b = T1%hopdn(i,T1%up(i))*T1%hopdn(j,T1%up(j))
!            d = T1%hopdn(i,T1%up(i))*T1%hopup(j,T1%up(j))

            ! dn*dn terms (note two ways for contraction!)
!            value1(k) = value1(k) - b* &
!                        (-dnt0(i,T1%up(j))*dn0t(j,T1%up(i)) - dnt0(T1%up(i),j)*dn0t(T1%up(j),i) &
!                         +dnt0(i,j)*dn0t(T1%up(j),T1%up(i)) + dnt0(T1%up(i),T1%up(j))*dn0t(j,i))*0.5_wp
!            value1(k) = value1(k) - b* &
!                        ( dntt(i,T1%up(i))*dn00(j,T1%up(j)) - dntt(T1%up(i),i)*dn00(T1%up(j),j) &
!                         -dntt(i,T1%up(i))*dn00(T1%up(j),j) + dntt(T1%up(i),i)*dn00(j,T1%up(j)))*0.5_wp
            ! dn*up terms
!            value1(k) = value1(k) - d* &
!                        ( dntt(i,T1%up(i))*up00(j,T1%up(j)) + dntt(T1%up(i),i)*up00(T1%up(j),j) &
!                         -dntt(i,T1%up(i))*up00(T1%up(j),j) - dntt(T1%up(i),i)*up00(j,T1%up(j)))*0.5_wp

            ! value2: see the rules at the beginning of routine
            b = T1%hopdn(j,T1%rt(j))*T1%hopdn(i,T1%rt(i))
            d = T1%hopdn(j,T1%rt(j))*T1%hopup(i,T1%rt(i))

            ! dn*dn terms (note two ways for contraction!)
            value2(k) = value2(k) - b* &
                        (-dnt0(j,T1%rt(i))*dn0t(i,T1%rt(j)) - dnt0(T1%rt(j),i)*dn0t(T1%rt(i),j) &
                         +dnt0(j,i)*dn0t(T1%rt(i),T1%rt(j)) + dnt0(T1%rt(j),T1%rt(i))*dn0t(i,j))*0.5_wp
            value2(k) = value2(k) - b* &
                        (dntt(j,T1%rt(j)) - dntt(T1%rt(j),j))*(dn00(i,T1%rt(i)) - dn00(T1%rt(i),i))*0.5_wp
            ! dn*up terms
            value2(k) = value2(k) - d* &
                        (dntt(j,T1%rt(j)) - dntt(T1%rt(j),j))*(up00(i,T1%rt(i)) - up00(T1%rt(i),i))*0.5_wp

            ! for averaging over x and y directions
!            b = T1%hopdn(j,T1%up(j))*T1%hopdn(i,T1%up(i))
!            d = T1%hopdn(j,T1%up(j))*T1%hopup(i,T1%up(i))

            ! dn*dn terms (note two ways for contraction!)
!            value2(k) = value2(k) - b* &
!                        (-dnt0(j,T1%up(i))*dn0t(i,T1%up(j)) - dnt0(T1%up(j),i)*dn0t(T1%up(i),j) &
!                         +dnt0(j,i)*dn0t(T1%up(i),T1%up(j)) + dnt0(T1%up(j),T1%up(i))*dn0t(i,j))*0.5_wp
!            value2(k) = value2(k) - b* &
!                        ( dntt(j,T1%up(j))*dn00(i,T1%up(i)) - dntt(T1%up(j),j)*dn00(T1%up(i),i) &
!                         -dntt(j,T1%up(j))*dn00(T1%up(i),i) + dntt(T1%up(j),j)*dn00(i,T1%up(i)))*0.5_wp
            ! dn*up terms
!            value2(k) = value2(k) - d* &
!                        ( dntt(j,T1%up(j))*up00(i,T1%up(i)) + dntt(T1%up(j),j)*up00(T1%up(i),i) &
!                         -dntt(j,T1%up(j))*up00(T1%up(i),i) - dntt(T1%up(j),j)*up00(i,T1%up(i)))*0.5_wp
          end do
       end do

       ! sum for conventional total conductivity
       T1%properties(ICOND)%values(:,dt1,T1%tmp) = &
           T1%properties(ICONDup)%values(:,dt1,T1%tmp) + T1%properties(ICONDdn)%values(:,dt1,T1%tmp)
       T1%properties(ICOND)%values(:,dt2,T1%tmp) = &
           T1%properties(ICONDup)%values(:,dt2,T1%tmp) + T1%properties(ICONDdn)%values(:,dt2,T1%tmp)

     endif

    else

    if (dt1/=0) then
      write(*,*) "dt should be 0, error in tdm.F90"
      return
    endif

     if (T1%flags(IGFUN) == 1) then
       value1  => T1%properties(IGFUN)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(IGFUN)%n
          do j = 1, T1%properties(IGFUN)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IGFUN)%D(i,j)
             value1(k)  = value1(k) + factor*(upt0(i,j) + dnt0(i,j))
          end do
       end do
     endif

     if (T1%flags(IGFUP) == 1) then
       value1  => T1%properties(IGFUP)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(IGFUP)%n
          do j = 1, T1%properties(IGFUP)%n
            ! k is the distance index of site i and site j
             k = T1%properties(IGFUP)%D(i,j)
             value1(k)  = value1(k) + 2*factor*upt0(i,j)
          end do
       end do
     endif

     if (T1%flags(IGFDN) == 1) then
       value1  => T1%properties(IGFDN)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(IGFDN)%n
          do j = 1, T1%properties(IGFDN)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IGFDN)%D(i,j)
             value1(k)  = value1(k) + 2*factor*dnt0(i,j)
          end do
       end do
     endif

     if (T1%flags(ISPXX) == 1) then
       value1  => T1%properties(ISPXX)%values(:, dt1, T1%tmp)
       do i = 1,  T1%properties(ISPXX)%n
          do j = 1,  T1%properties(ISPXX)%n
             ! k is the distance index of site i and site j
             k = T1%properties(ISPXX)%D(i,j)
 
             val1 = (up0t(j,i)*dnt0(i,j) + up0t(i,j)*dnt0(j,i))
             value1(k)  = value1(k) - val1

             ! compute chi_q of f-electron at q=(0,0) and (pi,pi)
             if (T1%flagsFT(ISPXX) == 1) then
               ! get the cartesian coordinates of site i
               ! Note i-1 accounts for the different convention of labelling
               ! sites  
               x1 = T1%cartpos(1,i-1)
               y1 = T1%cartpos(2,i-1)
               z1 = T1%cartpos(3,i-1)
               x2 = T1%cartpos(1,j-1)
               y2 = T1%cartpos(2,j-1)
               z2 = T1%cartpos(3,j-1)

               ! chi_ff calculation
               if (abs(z1-1.)<1.e-4 .and. abs(z2-1.)<1.e-4) then
                 ! 1) sum both inter- and intra-site even for multi-site unit cell
                 ! q=(0,0)
                 T1%chiqxx(dt1,1,T1%tmp) = T1%chiqxx(dt1,1,T1%tmp) - val1

                 ! q=(pi,pi), note the sign accouting for exp(separation*pi)
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   T1%chiqxx(dt1,2,T1%tmp) = T1%chiqxx(dt1,2,T1%tmp) - val1
                 else
                   T1%chiqxx(dt1,2,T1%tmp) = T1%chiqxx(dt1,2,T1%tmp) + val1
                 endif

                 ! 2) bipartite lattice can compute q=0 individually for each
                 ! sublat
                 ! separation limitation:
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   ! one sublattice
                   if (mod(int(x1+y1),2)==0) then
                     T1%chiqxx_sub1(dt1,T1%tmp) = T1%chiqxx_sub1(dt1,T1%tmp) - val1
                   ! the other sublattice
                   else
                     T1%chiqxx_sub2(dt1,T1%tmp) = T1%chiqxx_sub2(dt1,T1%tmp) - val1
                   endif
                 endif
               endif
             endif
          end do
       end do
     endif

     if (T1%flags(ISPZZ) == 1) then
       value1  => T1%properties(ISPZZ)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(ISPZZ)%n
          do j = 1, T1%properties(ISPZZ)%n
             ! k is the distance index of site i and site j
             k = T1%properties(ISPZZ)%D(i,j)

             val1 = (up0t(j,i)*upt0(i,j) + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )
             value1(k)  = value1(k) - val1

             ! compute chi_q of f-electron at q=(0,0) and (pi,pi)
             if (T1%flagsFT(ISPZZ) == 1) then
               ! get the cartesian coordinates of site i
               ! Note i-1 accounts for the different convention of labelling
               ! sites  
               x1 = T1%cartpos(1,i-1)
               y1 = T1%cartpos(2,i-1)
               z1 = T1%cartpos(3,i-1)
               x2 = T1%cartpos(1,j-1)
               y2 = T1%cartpos(2,j-1)
               z2 = T1%cartpos(3,j-1)

               ! chi_ff calculation
               if (abs(z1-1.)<1.e-4 .and. abs(z2-1.)<1.e-4) then
                 ! 1) sum both inter- and intra-site even for multi-site unit
                 ! cell
                 ! q=(0,0)
                 T1%chiqzz(dt1,1,T1%tmp) = T1%chiqzz(dt1,1,T1%tmp) - val1

                 ! q=(pi,pi), note the sign accouting for exp(separation*pi)
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   T1%chiqzz(dt1,2,T1%tmp) = T1%chiqzz(dt1,2,T1%tmp) - val1
                 else
                   T1%chiqzz(dt1,2,T1%tmp) = T1%chiqzz(dt1,2,T1%tmp) + val1
                 endif

                 ! 2) bipartite lattice can compute q=0 individually for each
                 ! sublat
                 ! separation limitation:
                 if (mod(int(x1-x2+y1-y2),2)==0) then
                   ! one sublattice
                   if (mod(int(x1+y1),2)==0) then
                     T1%chiqzz_sub1(dt1,T1%tmp) = T1%chiqzz_sub1(dt1,T1%tmp) - val1
                   ! the other sublattice
                   else
                     T1%chiqzz_sub2(dt1,T1%tmp) = T1%chiqzz_sub2(dt1,T1%tmp) - val1
                   endif
                 endif
               endif
            endif
          end do
       end do
     endif

     if (T1%flags(IDENS) == 1) then
       value1  => T1%properties(IDENS)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(IDENS)%n
          do j = 1, T1%properties(IDENS)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IDENS)%D(i,j)
             value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
               + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)+dntt(i,i))*(up00(j,j)+dn00(j,j)) )
          end do
       end do
     endif

     if (T1%flags(IPAIRs) == 1) then
       value1  => T1%properties(IPAIRs)%values(:, dt1, T1%tmp)
       do i = 1,  T1%properties(IPAIRs)%n
          do j = 1,  T1%properties(IPAIRs)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IPAIRs)%D(i,j)
             value1(k)  = value1(k) + upt0(i,j)*dnt0(i,j) 
          end do
       end do
     endif

     if (T1%flags(IPAIRd) == 1) then
       value1  => T1%properties(IPAIRd)%values(:, dt1, T1%tmp)

       ! D_i*D_j = sum_{dd'} Gup(i,j)*Gdn(i+d,j+d')
       do i = 1,  T1%properties(IPAIRd)%n
          do j = 1,  T1%properties(IPAIRd)%n
             k = T1%properties(IPAIRd)%D(i,j)

             ! get the cartesian coordinates of site i
             ! Note i-1 accounts for the different convention of labelling sites  
             z = T1%cartpos(3,i-1)
             vec = T1%properties(IPAIRd)%vecClass(k,:)

             ! 4 neighbors of each site so that totally 16 terms                                      
             ! with different phase factors for d-wave pairing
             ! record all 16 possible Gdn(i+d,j+d') as a below                                        
             a = dnt0(T1%rt(i), T1%rt(j)) - dnt0(T1%rt(i), T1%up(j))  &                               
                +dnt0(T1%rt(i), T1%lf(j)) - dnt0(T1%rt(i), T1%dn(j))  &                               
                -dnt0(T1%up(i), T1%rt(j)) + dnt0(T1%up(i), T1%up(j))  &                               
                -dnt0(T1%up(i), T1%lf(j)) + dnt0(T1%up(i), T1%dn(j))  &                               
                +dnt0(T1%lf(i), T1%rt(j)) - dnt0(T1%lf(i), T1%up(j))  &                               
                +dnt0(T1%lf(i), T1%lf(j)) - dnt0(T1%lf(i), T1%dn(j))  &                               
                -dnt0(T1%dn(i), T1%rt(j)) + dnt0(T1%dn(i), T1%up(j))  & 
                -dnt0(T1%dn(i), T1%lf(j)) + dnt0(T1%dn(i), T1%dn(j))    
               
             ! *0.25 or /4 accounts for the convention for definition                                 
             ! See 1989 PRB paper: Numerical study of 2D Hubbard model                                
             a = a*0.25                                                                        
             value1(k) = value1(k) + upt0(i,j)*a
                                            
             ! only compute in-plane d-wave pairing susceptibility (square lattice)
             select case (model)
               ! Hubbard
               case (0)
                 T1%Pdtau(dt1, T1%tmp, 1) = T1%Pdtau(dt1, T1%tmp, 1) + upt0(i,j)*a

               ! PAM (P_ffff only temporarily)
               case (1)
                 if (abs(z1-1.d0)<1.d-6 .and. abs(z2-1.d0)<1.d-6) then
                   T1%Pdtau(dt1, T1%tmp, 1) = T1%Pdtau(dt1, T1%tmp, 1) + upt0(i,j)*a
                 endif
             end select
          end do
       end do
     endif

     if (T1%flags(ICOND) == 1) then
       
       ! First compute cond from spin up
       value1  => T1%properties(ICONDup)%values(:, dt1, T1%tmp)

       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n

            k = T1%properties(ICOND)%D(i,j)

            a = T1%hopup(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            c = T1%hopup(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            ! up*up terms (note two ways for contraction!)
            value1(k) = value1(k) - a* &
                        (-upt0(i,T1%rt(j))*up0t(j,T1%rt(i)) - upt0(T1%rt(i),j)*up0t(T1%rt(j),i) &
                         +upt0(i,j)*up0t(T1%rt(j),T1%rt(i)) + upt0(T1%rt(i),T1%rt(j))*up0t(j,i))
            value1(k) = value1(k) - a* &
                        (uptt(i,T1%rt(i)) - uptt(T1%rt(i),i))*(up00(j,T1%rt(j)) - up00(T1%rt(j),j))*0.5_wp
            ! up*dn terms
            value1(k) = value1(k) - c* &
                        (uptt(i,T1%rt(i)) - uptt(T1%rt(i),i))*(dn00(j,T1%rt(j)) - dn00(T1%rt(j),j))*0.5_wp

!            a = T1%hopup(i,T1%up(i))*T1%hopup(j,T1%up(j))
!            c = T1%hopup(i,T1%up(i))*T1%hopdn(j,T1%up(j))
            ! up*up terms (note two ways for contraction!)
!            value1(k) = value1(k) - a* &
!                        (-upt0(i,T1%up(j))*up0t(j,T1%up(i)) - upt0(T1%up(i),j)*up0t(T1%up(j),i) &
!                         +upt0(i,j)*up0t(T1%up(j),T1%up(i)) + upt0(T1%up(i),T1%up(j))*up0t(j,i))
!            value1(k) = value1(k) - a* &
!                        ( uptt(i,T1%up(i))*up00(j,T1%up(j)) - uptt(T1%up(i),i)*up00(T1%up(j),j) &
!                         -uptt(i,T1%up(i))*up00(T1%up(j),j) + uptt(T1%up(i),i)*up00(j,T1%up(j)))
            ! up*dn terms
!            value1(k) = value1(k) - c* &
!                        ( uptt(i,T1%up(i))*dn00(j,T1%up(j)) + uptt(T1%up(i),i)*dn00(T1%up(j),j) &
!                         -uptt(i,T1%up(i))*dn00(T1%up(j),j) - uptt(T1%up(i),i)*dn00(j,T1%up(j)))
          end do
       end do

       ! Compute cond for down spin
       value1  => T1%properties(ICONDdn)%values(:, dt1, T1%tmp)
       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n

            k = T1%properties(ICOND)%D(i,j)

            b = T1%hopdn(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            d = T1%hopdn(i,T1%rt(i))*T1%hopup(j,T1%rt(j))

            ! dn*dn terms (note two ways for contraction!)
            value1(k) = value1(k) - b* &
                        (-dnt0(i,T1%rt(j))*dn0t(j,T1%rt(i)) - dnt0(T1%rt(i),j)*dn0t(T1%rt(j),i) &
                         +dnt0(i,j)*dn0t(T1%rt(j),T1%rt(i)) + dnt0(T1%rt(i),T1%rt(j))*dn0t(j,i))
            value1(k) = value1(k) - b* &
                        (dntt(i,T1%rt(i)) - dntt(T1%rt(i),i))*(dn00(j,T1%rt(j)) - dn00(T1%rt(j),j))*0.5_wp
            ! dn*up terms
            value1(k) = value1(k) - d* &
                        (dntt(i,T1%rt(i)) - dntt(T1%rt(i),i))*(up00(j,T1%rt(j)) - up00(T1%rt(j),j))*0.5_wp

!            b = T1%hopdn(i,T1%up(i))*T1%hopdn(j,T1%up(j))
!            d = T1%hopdn(i,T1%up(i))*T1%hopup(j,T1%up(j))

            ! dn*dn terms (note two ways for contraction!)
!            value1(k) = value1(k) - b* &
!                        (-dnt0(i,T1%up(j))*dn0t(j,T1%up(i)) - dnt0(T1%up(i),j)*dn0t(T1%up(j),i) &
!                         +dnt0(i,j)*dn0t(T1%up(j),T1%up(i)) + dnt0(T1%up(i),T1%up(j))*dn0t(j,i))
!            value1(k) = value1(k) - b* &
!                        ( dntt(i,T1%up(i))*dn00(j,T1%up(j)) - dntt(T1%up(i),i)*dn00(T1%up(j),j) &
!                         -dntt(i,T1%up(i))*dn00(T1%up(j),j) + dntt(T1%up(i),i)*dn00(j,T1%up(j)))
            ! dn*up terms
!            value1(k) = value1(k) - d* &
!                        ( dntt(i,T1%up(i))*up00(j,T1%up(j)) + dntt(T1%up(i),i)*up00(T1%up(j),j) &
!                         -dntt(i,T1%up(i))*up00(T1%up(j),j) - dntt(T1%up(i),i)*up00(j,T1%up(j)))
          end do
       end do

       ! sum for conventional total conductivity
       T1%properties(ICOND)%values(:,dt1,T1%tmp) = &
            T1%properties(ICONDup)%values(:,dt1,T1%tmp) + T1%properties(ICONDdn)%values(:,dt1,T1%tmp)

     endif

    endif

  end subroutine DQMC_TDM_Compute

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Avg(T1, model)
    !
    ! Purpose
    ! =======
    !    This subroutine average properties in a bin and
    !    increment the bin count (idx).
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer  :: nl, idx, i, j, k, it, z1, z2
    integer, intent(in)  :: model
    real(wp) :: factor, fac, a
    real(wp), allocatable :: value1(:), value2(:)

    ! ... Executable ...
    if (.not.T1%compute) return
    idx    = T1%idx
    factor = ONE/T1%cnt
    
    ! /2.0 accounts for two orbitals per site
    ! Note that this is only for PAM
    fac = factor/(T1%properties(ISPXX)%n/2.d0)  

    allocate(value1(1:T1%properties(IGFUN)%nClass))
    allocate(value2(1:T1%properties(IGFUN)%nClass))

    ! Compute average on Green's function
    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         nl = T1%properties(i)%nClass
         do j = 0, T1%L-1
           call dscal(nl, factor, T1%properties(i)%values(:,j,idx), 1)
         enddo
       endif
    enddo

    if (T1%flags(ISPXX) == 1) then
      do i = 1, T1%properties(ISPXX)%nClass 
        T1%chixx_r_orb_iw0(i,idx) = T1%chixx_r_orb_iw0(i,idx) * fac
      enddo
    endif

    if (T1%flags(ISPZZ) == 1) then
      do i = 1, T1%properties(ISPZZ)%nClass 
        T1%chizz_r_orb_iw0(i,idx) = T1%chizz_r_orb_iw0(i,idx) * fac
      enddo
    endif

    if (T1%flagsFT(ISPXX) == 1) then
      T1%chiqxx(0:T1%L-1,1:2,idx) = T1%chiqxx(0:T1%L-1,1:2,idx) * fac
      T1%chiqxx_sub1(0:T1%L-1,idx) = T1%chiqxx_sub1(0:T1%L-1,idx) * fac
      T1%chiqxx_sub2(0:T1%L-1,idx) = T1%chiqxx_sub2(0:T1%L-1,idx) * fac
      T1%chiqxx_iw0(1:2,idx) = T1%chiqxx_iw0(1:2,idx) * fac
      T1%chiqxx_iw0_sub1(idx) = T1%chiqxx_iw0_sub1(idx) * fac
      T1%chiqxx_iw0_sub2(idx) = T1%chiqxx_iw0_sub2(idx) * fac 
    endif

    if (T1%flagsFT(ISPZZ) == 1) then
      T1%chiqzz(0:T1%L-1,1:2,idx) = T1%chiqzz(0:T1%L-1,1:2,idx) * fac
      T1%chiqzz_sub1(0:T1%L-1,idx) = T1%chiqzz_sub1(0:T1%L-1,idx) * fac
      T1%chiqzz_sub2(0:T1%L-1,idx) = T1%chiqzz_sub2(0:T1%L-1,idx) * fac
      T1%chiqzz_iw0(1:2,idx) = T1%chiqzz_iw0(1:2,idx) * fac
      T1%chiqzz_iw0_sub1(idx) = T1%chiqzz_iw0_sub1(idx) * fac
      T1%chiqzz_iw0_sub2(idx) = T1%chiqzz_iw0_sub2(idx) * fac
    endif

    ! record the average of all local sum_r G(r, tau) for average N(w)
    if (T1%flags(IGFUN) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IGFUN)%n    ! avg over n sites
            k = T1%properties(IGFUN)%D(i,i)  ! local quantity
            T1%GtauAvg(j, T1%idx) = T1%GtauAvg(j, T1%idx) &
                   + T1%properties(IGFUN)%values(k,j,idx) / T1%properties(IGFUN)%n
         end do
      enddo
    endif
    if (T1%flags(IGFUP) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IGFUP)%n
            k = T1%properties(IGFUP)%D(i,i)
            T1%GtupAvg(j, T1%idx) = T1%GtupAvg(j, T1%idx) &
                   + T1%properties(IGFUP)%values(k,j,idx) / T1%properties(IGFUP)%n
         end do
      enddo
    endif
    if (T1%flags(IGFDN) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IGFDN)%n
            k = T1%properties(IGFDN)%D(i,i)
            T1%GtdnAvg(j, T1%idx) = T1%GtdnAvg(j, T1%idx) &
                   + T1%properties(IGFDN)%values(k,j,idx) / T1%properties(IGFDN)%n
         end do
      end do
    endif

    ! record the average of all local sum_r spin-xx(r, tau) for average spin-xx susceptibility
    if (T1%flags(ISPXX) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(ISPXX)%n    ! avg over n sites
            k = T1%properties(ISPXX)%D(i,i)  ! local quantity
            T1%spinxxAvg(j, T1%idx) = T1%spinxxAvg(j, T1%idx) &
                   + T1%properties(ISPXX)%values(k,j,idx) / T1%properties(ISPXX)%n
         end do
      enddo
    endif

    ! record the average of all local sum_r spin-zz(r, tau) for average spin-zz susceptibility
    if (T1%flags(ISPZZ) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(ISPZZ)%n    ! avg over n sites
            k = T1%properties(ISPZZ)%D(i,i)  ! local quantity
            T1%spinzzAvg(j, T1%idx) = T1%spinzzAvg(j, T1%idx) &
                   + T1%properties(ISPZZ)%values(k,j,idx) / T1%properties(ISPZZ)%n
         end do
      enddo
    endif

    ! record the average of all local sum_r pair(r, tau) for average pair
    ! susceptibility
    if (T1%flags(IPAIRs) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IPAIRs)%n    ! avg over n sites
            k = T1%properties(IPAIRs)%D(i,i)  ! local quantity
            T1%swaveAvg(j, T1%idx) = T1%swaveAvg(j, T1%idx) &
                   + T1%properties(IPAIRs)%values(k,j,idx) / T1%properties(IPAIRs)%n
         end do
      enddo
    endif

    ! Correlated d-wave susceptibility Pd
    if (T1%flags(IPAIRd) == 1) then
       T1%Pdtau(0:T1%L-1, T1%idx, :) = T1%Pdtau(0:T1%L-1, T1%idx, :) * factor
       T1%Pd(T1%idx, :) = T1%Pd(T1%idx, :) * factor
    endif

    ! compute uncorrelated or non-vertex d-wave susceptibility Pd^bar
    ! D_i*D_j = sum_{dd'} <Gup(i,j)>*<Gdn(i+d,j+d')>
    ! Pd0 = int^beta_0 sum_ij D_i*D_j
    if (T1%flags(IPAIRd) == 1) then
       do it = 0, T1%L-1
         ! Note that precisely value1 (value2) should use
         ! Gup and Gdn respectively. Because sometimes the code
         ! input file does not specify computing them
         ! Here assume that G = Gup = Gdn
         value1 = T1%properties(IGFUN)%values(:, it, idx)
         value2 = T1%properties(IGFUN)%values(:, it, idx)

         a = 0.0_wp
         do i = 1,  T1%properties(IPAIRd)%n
            do j = 1,  T1%properties(IPAIRd)%n
               ! 4 neighbors of each site so that totally 16 terms
               ! with different phase factors for d-wave pairing
               ! record all 16 possible Gdn(i+d,j+d') as a below
               ! *0.25 or /4 is convention, see the computation of Pd

               k = T1%properties(IPAIRd)%D(T1%rt(i), T1%rt(j))
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%rt(i), T1%up(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%rt(i), T1%lf(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%rt(i), T1%dn(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%up(i), T1%rt(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%up(i), T1%up(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%up(i), T1%lf(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%up(i), T1%dn(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%lf(i), T1%rt(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%lf(i), T1%up(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%lf(i), T1%lf(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%lf(i), T1%dn(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%dn(i), T1%rt(j))                                        
               a = a + value2(k)
               k = T1%properties(IPAIRd)%D(T1%dn(i), T1%up(j))                                        
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%dn(i), T1%lf(j))
               a = a - value2(k)
               k = T1%properties(IPAIRd)%D(T1%dn(i), T1%dn(j))
               a = a + value2(k)

               k = T1%properties(IPAIRd)%D(i,j)

               ! get the cartesian coordinates of site i
               ! Note i-1 accounts for the different convention of labelling sites
               z1 = T1%cartpos(3,i-1)
               z2 = T1%cartpos(3,j-1)

               ! 6/18/2019:
               ! *2.0 is only aimed to agree with RTScode and Pd at U=0
               ! The reason may be that computing Pd in TDM_Compute assume
               ! G(i,j) = G(j,i) because (i,j) and (j,i) have same class k

               ! only compute in-plane d-wave pairing susceptibility (square lattice)
               select case (model)
                 ! Hubbard
                 case (0)
                   T1%Pd0tau(it, idx, 1) = T1%Pd0tau(it, idx, 1) + &
                                           value1(k)*a/(4.0*T1%properties(IPAIRd)%n)*2.0

                 ! PAM (P_ffff only temporarily)
                 case (1)
                   if (abs(z1-1.d0)<1.d-6 .and. abs(z2-1.d0)<1.d-6) then
                     T1%Pd0tau(it, idx, 1) = T1%Pd0tau(it, idx, 1) + &
                                             value1(k)*a/(4.0*T1%properties(IPAIRd)%n)*2.0
                   endif
               end select
            end do
         end do
         !write(*,*) T1%Pd0tau(it)
       enddo

       do i = 1,T1%NPd
         call convert_to_iw0_real(T1%Pd0tau(0:T1%L-1, idx, i), T1%Pd0(idx, i), T1%L, T1%dtau)

         ! calculate the d-wave pairing vertex Gammad
         T1%Gammad(T1%idx,i) = 1.0/T1%Pd(T1%idx,i) - 1.0/T1%Pd0(idx,i)

         ! product of d-wave pairing vertex Gammad and Pd0
         ! superconducting instability is signified by Gd_Pd0 -> -1
         ! See PRB 86, 184506 (2012)
         T1%Gd_Pd0(T1%idx,i) = T1%Gammad(T1%idx,i)*T1%Pd0(idx,i)
       enddo
    endif

    T1%sgn(idx) = T1%sgn(idx)*factor
    T1%cnt = 0
    T1%idx = T1%idx + 1

  end subroutine DQMC_TDM_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_GetErr(T1)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine compute the error in tdm using the jackknife
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer   :: i, j, k, iprop
    integer   :: nproc, n, avg, err, mpi_err
    real(wp)  :: sum_sgn, sgn(T1%nBin), y(T1%nBin), data(T1%nBin)
    real(wp)  :: average, error

#   ifdef _QMC_MPI
      real(wp), pointer :: binptr(:,:), aveptr(:,:), errptr(:,:)
      real(wp), pointer :: bins(:), aves(:), errs(:)
#   endif

    if (.not.T1%compute) return
    ! ... Executable ...
    nproc  = qmc_sim%size
    n      = T1%nbin
    avg    = T1%avg
    err    = T1%err

    if (nproc .eq. 1) then

       data = T1%sgn(1:n)
       call DQMC_JackKnife(n, T1%sgn(avg), T1%sgn(err), data , &
            y, sgn, sum_sgn)

       do iprop = 1, NTDMARRAY
         if (T1%flags(iprop)==1) then
           do i = 1, T1%properties(iprop)%nClass
             do j = 0, T1%L-1
                data =  T1%properties(iprop)%values(i, j, 1:n)
                call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
                T1%properties(iprop)%values(i, j, avg) = average
                T1%properties(iprop)%values(i, j, err) = error
             enddo
           end do
         endif
       enddo

       if (T1%flags(ISPXX) == 1) then
         do i = 1, T1%properties(ISPXX)%nClass
           data =  T1%chixx_r_orb_iw0(i,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chixx_r_orb_iw0(i,avg) = average
           T1%chixx_r_orb_iw0(i,err) = error
         enddo
       endif

       if (T1%flagsFT(ISPXX) == 1) then
         do i = 1, 2
           do k = 0, T1%L-1
             data =  T1%chiqxx(k,i,1:n)
             call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
             T1%chiqxx(k,i,avg) = average
             T1%chiqxx(k,i,err) = error
           enddo

           data =  T1%chiqxx_iw0(i,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqxx_iw0(i,avg) = average
           T1%chiqxx_iw0(i,err) = error
         enddo

         ! for two sublattices
         do k = 0, T1%L-1
           data =  T1%chiqxx_sub1(k,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqxx_sub1(k,avg) = average
           T1%chiqxx_sub1(k,err) = error
         enddo
         do k = 0, T1%L-1
           data =  T1%chiqxx_sub2(k,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqxx_sub2(k,avg) = average
           T1%chiqxx_sub2(k,err) = error
         enddo

         data =  T1%chiqxx_iw0_sub1(1:n)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%chiqxx_iw0_sub1(avg) = average
         T1%chiqxx_iw0_sub1(err) = error

         data =  T1%chiqxx_iw0_sub2(1:n)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%chiqxx_iw0_sub2(avg) = average
         T1%chiqxx_iw0_sub2(err) = error
       endif

       if (T1%flags(ISPZZ) == 1) then
         do i = 1, T1%properties(ISPZZ)%nClass
           data =  T1%chizz_r_orb_iw0(i,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chizz_r_orb_iw0(i,avg) = average
           T1%chizz_r_orb_iw0(i,err) = error
         enddo
       endif

       if (T1%flagsFT(ISPZZ) == 1) then
         do i = 1, 2
           do k = 0, T1%L-1
             data =  T1%chiqzz(k,i,1:n)
             call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
             T1%chiqzz(k,i,avg) = average
             T1%chiqzz(k,i,err) = error
           enddo

           data =  T1%chiqzz_iw0(i,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqzz_iw0(i,avg) = average
           T1%chiqzz_iw0(i,err) = error
         enddo

         ! for two sublattices
         do k = 0, T1%L-1
           data =  T1%chiqzz_sub1(k,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqzz_sub1(k,avg) = average
           T1%chiqzz_sub1(k,err) = error
         enddo
         do k = 0, T1%L-1
           data =  T1%chiqzz_sub2(k,1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%chiqzz_sub2(k,avg) = average
           T1%chiqzz_sub2(k,err) = error
         enddo

         data =  T1%chiqzz_iw0_sub1(1:n)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%chiqzz_iw0_sub1(avg) = average
         T1%chiqzz_iw0_sub1(err) = error

         data =  T1%chiqzz_iw0_sub2(1:n)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%chiqzz_iw0_sub2(avg) = average
         T1%chiqzz_iw0_sub2(err) = error
       endif

       if (T1%flags(IGFUN) == 1) then
         do j = 0, T1%L-1
           data =  T1%GtauAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%GtauAvg(j, avg) = average
           T1%GtauAvg(j, err) = error
         enddo
       endif
       if (T1%flags(IGFUP) == 1) then
         do j = 0, T1%L-1
           data =  T1%GtupAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%GtupAvg(j, avg) = average
           T1%GtupAvg(j, err) = error
         enddo
       endif
       if (T1%flags(IGFDN) == 1) then
         do j = 0, T1%L-1
           data =  T1%GtdnAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%GtdnAvg(j, avg) = average
           T1%GtdnAvg(j, err) = error
         enddo
       endif

       if (T1%flags(ISPXX) == 1) then
         do j = 0, T1%L-1
           data =  T1%spinxxAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%spinxxAvg(j, avg) = average
           T1%spinxxAvg(j, err) = error
         enddo
       endif

       if (T1%flags(ISPZZ) == 1) then
         do j = 0, T1%L-1
           data =  T1%spinzzAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%spinzzAvg(j, avg) = average
           T1%spinzzAvg(j, err) = error
         enddo
       endif

       if (T1%flags(IPAIRs) == 1) then
         do j = 0, T1%L-1
           data =  T1%swaveAvg(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%swaveAvg(j, avg) = average
           T1%swaveAvg(j, err) = error
         enddo
       endif

       if (T1%flags(IPAIRd) == 1) then
         do j = 0, T1%L-1
           data =  T1%Pdtau(j, 1:n, 1)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%Pdtau(j, avg, 1) = average
           T1%Pdtau(j, err, 1) = error
         enddo
         do j = 0, T1%L-1
           data =  T1%Pd0tau(j, 1:n, 1)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%Pd0tau(j, avg, 1) = average
           T1%Pd0tau(j, err, 1) = error
         enddo

         data =  T1%Pd(1:n, 1)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%Pd(avg, 1) = average
         T1%Pd(err, 1) = error

         data =  T1%Pd0(1:n, 1)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%Pd0(avg, 1) = average
         T1%Pd0(err, 1) = error

         data =  T1%Gammad(1:n, 1)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%Gammad(avg, 1) = average
         T1%Gammad(err, 1) = error

         data =  T1%Gd_Pd0(1:n, 1)
         call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
         T1%Gd_Pd0(avg, 1) = average
         T1%Gd_Pd0(err, 1) = error
       endif

    else

       mpi_err = 0

#      ifdef _QMC_MPI
          
          ! comments below similar to dqmc_phy0.F90

          ! note for MPI, each processor has only ONE bin          
          ! The process below for computing err is JackKnife similar to non-MPI case
          ! See also DQMC_SignJackKnife_Real in dqmc_util.F90
          ! Note DQMC_SignJackKnife_Real differs DQMC_JackKnife_Real, so computing y_i does not have /(n-1)

          !    y_i = (sum(x)-x_i)/sgn_i
          !    The JackKnife variance of X with sign is defined as 
          !     
          !          n-1  
          !    sqrt(----- *sum(y_i-avg_y)^2))
          !           n
          ! 
          !    where avg_y = sum(y)/n

          ! compute sum(x)          
          call mpi_allreduce(T1%sgn(1), T1%sgn(avg), 1, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                n = T1%properties(iprop)%nClass * T1%L
                call mpi_allreduce(binptr, aveptr, n, mpi_double, &
                   mpi_sum, mpi_comm_world, mpi_err)
            endif
          enddo

          if (T1%flags(ISPXX) == 1) then
            n = T1%properties(ISPXX)%nClass
            bins => T1%chixx_r_orb_iw0(:,1)
            aves => T1%chixx_r_orb_iw0(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
          endif

          if (T1%flagsFT(ISPXX) == 1) then
            n = T1%L
            do k = 1,2
              bins => T1%chiqxx(:,k,1)
              aves => T1%chiqxx(:,k,avg)             
              call mpi_allreduce(bins, aves, n, mpi_double, &
                 mpi_sum, mpi_comm_world, mpi_err)
            enddo

            ! for two sublattices
            bins => T1%chiqxx_sub1(:,1)
            aves => T1%chiqxx_sub1(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            bins => T1%chiqxx_sub2(:,1)
            aves => T1%chiqxx_sub2(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            bins => T1%chiqxx_iw0(1:2,1)
            aves => T1%chiqxx_iw0(1:2,avg)
            call mpi_allreduce(bins, aves, 2, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            call mpi_allreduce(T1%chiqxx_iw0_sub1(1), T1%chiqxx_iw0_sub1(avg), 1, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
            call mpi_allreduce(T1%chiqxx_iw0_sub2(1), T1%chiqxx_iw0_sub2(avg), 1, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
          endif

          if (T1%flags(ISPZZ) == 1) then
            n = T1%properties(ISPZZ)%nClass
            bins => T1%chizz_r_orb_iw0(:,1)
            aves => T1%chizz_r_orb_iw0(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
          endif

          if (T1%flagsFT(ISPZZ) == 1) then
            n = T1%L
            do k = 1,2
              bins => T1%chiqzz(:,k,1)
              aves => T1%chiqzz(:,k,avg)
              call mpi_allreduce(bins, aves, n, mpi_double, &
                 mpi_sum, mpi_comm_world, mpi_err)
            enddo

            ! for two sublattices
            bins => T1%chiqzz_sub1(:,1)
            aves => T1%chiqzz_sub1(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            bins => T1%chiqzz_sub2(:,1)
            aves => T1%chiqzz_sub2(:,avg)
            call mpi_allreduce(bins, aves, n, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            bins => T1%chiqzz_iw0(1:2,1)
            aves => T1%chiqzz_iw0(1:2,avg)
            call mpi_allreduce(bins, aves, 2, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)

            call mpi_allreduce(T1%chiqzz_iw0_sub1(1), T1%chiqzz_iw0_sub1(avg), 1, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
            call mpi_allreduce(T1%chiqzz_iw0_sub2(1), T1%chiqzz_iw0_sub2(avg), 1, mpi_double, &
               mpi_sum, mpi_comm_world, mpi_err)
          endif

          if (T1%flags(IGFUN) == 1) then
             bins => T1%GtauAvg(:, 1)
             aves => T1%GtauAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(IGFUP) == 1) then
             bins => T1%GtupAvg(:, 1)
             aves => T1%GtupAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(IGFDN) == 1) then
             bins => T1%GtdnAvg(:, 1)
             aves => T1%GtdnAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(ISPXX) == 1) then
             bins => T1%spinxxAvg(:, 1)
             aves => T1%spinxxAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(ISPZZ) == 1) then
             bins => T1%spinzzAvg(:, 1)
             aves => T1%spinzzAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(IPAIRs) == 1) then
             bins => T1%swaveAvg(:, 1)
             aves => T1%swaveAvg(:, avg)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif
          if (T1%flags(IPAIRd) == 1) then
             bins => T1%Pdtau(:, 1, 1)
             aves => T1%Pdtau(:, avg, 1)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)

             bins => T1%Pd0tau(:, 1, 1)
             aves => T1%Pd0tau(:, avg, 1)
             call mpi_allreduce(bins, aves, T1%L, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)

             call mpi_allreduce(T1%Pd(1), T1%Pd(avg), 1, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
             call mpi_allreduce(T1%Pd0(1), T1%Pd0(avg), 1, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
             call mpi_allreduce(T1%Gammad(1), T1%Gammad(avg), 1, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
             call mpi_allreduce(T1%Gd_Pd0(1), T1%Gd_Pd0(avg), 1, mpi_double, &
                     mpi_sum, mpi_comm_world, mpi_err)
          endif

          ! Compute y_i, note original binned values will be updated
          ! non-MPI version of code would not update binned values
          ! as DQMC_SignJackKnife does not change data
          ! Note that ave is now storing sum(x)
          T1%sgn(1)   = (T1%sgn(avg) - T1%sgn(1)) / dble(nproc - 1)
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                binptr = (aveptr - binptr) / dble(nproc - 1)
                binptr =  binptr / T1%sgn(1)
            endif
          enddo

          if (T1%flags(ISPXX) == 1) then
            bins => T1%chixx_r_orb_iw0(:,1)
            aves => T1%chixx_r_orb_iw0(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)
          endif

          if (T1%flagsFT(ISPXX) == 1) then
            binptr => T1%chiqxx(:,:,1)
            aveptr => T1%chiqxx(:,:,avg) 
            binptr = (aveptr - binptr) / dble(nproc - 1)
            binptr =  binptr / T1%sgn(1)

            bins => T1%chiqxx_iw0(:,1)
            aves => T1%chiqxx_iw0(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)

            bins => T1%chiqxx_sub1(:,1)
            aves => T1%chiqxx_sub1(:,avg)
            bins = (aves - bins) / dble(nproc - 1)                                                    
            bins =  bins / T1%sgn(1)

            bins => T1%chiqxx_sub2(:,1)
            aves => T1%chiqxx_sub2(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)

            T1%chiqxx_iw0_sub1(1) =  (T1%chiqxx_iw0_sub1(avg) &
                   - T1%chiqxx_iw0_sub1(1))/dble(nproc - 1)/T1%sgn(1)
            T1%chiqxx_iw0_sub2(1) =  (T1%chiqxx_iw0_sub2(avg) &
                   - T1%chiqxx_iw0_sub2(1))/dble(nproc - 1)/T1%sgn(1)
          endif

          if (T1%flags(ISPZZ) == 1) then
            bins => T1%chizz_r_orb_iw0(:,1)
            aves => T1%chizz_r_orb_iw0(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)
          endif

          if (T1%flagsFT(ISPZZ) == 1) then
            binptr => T1%chiqzz(:,:,1)
            aveptr => T1%chiqzz(:,:,avg)
            binptr = (aveptr - binptr) / dble(nproc - 1)
            binptr =  binptr / T1%sgn(1)

            bins => T1%chiqzz_iw0(:,1)
            aves => T1%chiqzz_iw0(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)

            bins => T1%chiqzz_sub1(:,1)
            aves => T1%chiqzz_sub1(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)

            bins => T1%chiqzz_sub2(:,1)
            aves => T1%chiqzz_sub2(:,avg)
            bins = (aves - bins) / dble(nproc - 1)
            bins =  bins / T1%sgn(1)

            T1%chiqzz_iw0_sub1(1) =  (T1%chiqzz_iw0_sub1(avg) &
                   - T1%chiqzz_iw0_sub1(1))/dble(nproc - 1)/T1%sgn(1)
            T1%chiqzz_iw0_sub2(1) =  (T1%chiqzz_iw0_sub2(avg) &
                   - T1%chiqzz_iw0_sub2(1))/dble(nproc - 1)/T1%sgn(1)
          endif

          if (T1%flags(IGFUN) == 1) then
             bins => T1%GtauAvg(:, 1)
             aves => T1%GtauAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(IGFUP) == 1) then
             bins => T1%GtupAvg(:, 1)
             aves => T1%GtupAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(IGFDN) == 1) then
             bins => T1%GtdnAvg(:, 1)
             aves => T1%GtdnAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(ISPXX) == 1) then
             bins => T1%spinxxAvg(:, 1)
             aves => T1%spinxxAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(ISPZZ) == 1) then
             bins => T1%spinzzAvg(:, 1)
             aves => T1%spinzzAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(IPAIRs) == 1) then
             bins => T1%swaveAvg(:, 1)
             aves => T1%swaveAvg(:, avg)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)
          endif
          if (T1%flags(IPAIRd) == 1) then
             bins => T1%Pdtau(:, 1, 1)
             aves => T1%Pdtau(:, avg, 1)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)

             bins => T1%Pd0tau(:, 1, 1)
             aves => T1%Pd0tau(:, avg, 1)
             bins = (aves - bins) / dble(nproc - 1)
             bins =  bins / T1%sgn(1)

             T1%Pd(1, 1) = (T1%Pd(avg, 1) - T1%Pd(1, 1)) / dble(nproc - 1)
             T1%Pd(1, 1) =  T1%Pd(1, 1) / T1%sgn(1)
             T1%Pd0(1, 1) = (T1%Pd0(avg, 1) - T1%Pd0(1, 1)) / dble(nproc - 1)
             T1%Pd0(1, 1) =  T1%Pd0(1, 1) / T1%sgn(1)
             T1%Gammad(1, 1) = (T1%Gammad(avg, 1) - T1%Gammad(1, 1)) / dble(nproc - 1)
             T1%Gammad(1, 1) =  T1%Gammad(1, 1) / T1%sgn(1)
             T1%Gd_Pd0(1, 1) = (T1%Gd_Pd0(avg, 1) - T1%Gd_Pd0(1, 1)) / dble(nproc - 1)
             T1%Gd_Pd0(1, 1) =  T1%Gd_Pd0(1, 1) / T1%sgn(1)
          endif

          ! Compute avg_y = avg_x = sum_x/sum_sgn
          ! because in previous compute sum(x) step
          ! Recall that aveptr is still for storing sum(x)
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                aveptr => T1%properties(iprop)%values(:,:,avg)
                aveptr =  aveptr / T1%sgn(avg) 
            endif
          enddo

          if (T1%flags(ISPXX) == 1) then
            aves => T1%chixx_r_orb_iw0(:,avg)
            aves =  aves / T1%sgn(avg)
          endif

          if (T1%flagsFT(ISPXX) == 1) then
            aveptr => T1%chiqxx(:,:,avg)
            aveptr =  aveptr / T1%sgn(avg)

            aves => T1%chiqxx_iw0(:,avg)
            aves =  aves / T1%sgn(avg)

            aves => T1%chiqxx_sub1(:,avg)
            aves =  aves / T1%sgn(avg)
            aves => T1%chiqxx_sub2(:,avg)
            aves =  aves / T1%sgn(avg)

            T1%chiqxx_iw0_sub1(avg) = T1%chiqxx_iw0_sub1(avg) / T1%sgn(avg)
            T1%chiqxx_iw0_sub2(avg) = T1%chiqxx_iw0_sub2(avg) / T1%sgn(avg)
          endif

          if (T1%flags(ISPZZ) == 1) then
            aves => T1%chizz_r_orb_iw0(:,avg)
            aves =  aves / T1%sgn(avg)
          endif

          if (T1%flagsFT(ISPZZ) == 1) then
            aveptr => T1%chiqzz(:,:,avg)
            aveptr =  aveptr / T1%sgn(avg)

            aves => T1%chiqzz_iw0(:,avg)
            aves =  aves / T1%sgn(avg)

            aves => T1%chiqzz_sub1(:,avg)
            aves =  aves / T1%sgn(avg)
            aves => T1%chiqzz_sub2(:,avg)
            aves =  aves / T1%sgn(avg)

            T1%chiqzz_iw0_sub1(avg) = T1%chiqzz_iw0_sub1(avg) / T1%sgn(avg)
            T1%chiqzz_iw0_sub2(avg) = T1%chiqzz_iw0_sub2(avg) / T1%sgn(avg)
          endif

          if (T1%flags(IGFUN) == 1) then
             aves => T1%GtauAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(IGFUP) == 1) then
             aves => T1%GtupAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(IGFDN) == 1) then
             aves => T1%GtdnAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(ISPXX) == 1) then
             aves => T1%spinxxAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(ISPZZ) == 1) then
             aves => T1%spinzzAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(IPAIRs) == 1) then
             aves => T1%swaveAvg(:, avg)
             aves =  aves / T1%sgn(avg)
          endif
          if (T1%flags(IPAIRd) == 1) then
             aves => T1%Pdtau(:, avg, 1)
             aves =  aves / T1%sgn(avg)

             aves => T1%Pd0tau(:, avg, 1)
             aves =  aves / T1%sgn(avg)

             T1%Pd (avg, 1) = T1%Pd (avg, 1) / T1%sgn(avg)
             T1%Pd0(avg, 1) = T1%Pd0(avg, 1) / T1%sgn(avg)
             T1%Gammad(avg, 1) = T1%Gammad(avg, 1) / T1%sgn(avg)
             T1%Gd_Pd0(avg, 1) = T1%Gd_Pd0(avg, 1) / T1%sgn(avg)
          endif

          T1%sgn(avg)  = T1%sgn(avg) / dble(nproc)

          ! Compute error: sum(y_i-avg_y)^2
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                !note here binned values are y_i
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                errptr => T1%properties(iprop)%values(:,:,err)
                n = T1%properties(iprop)%nClass * T1%L
                call mpi_allreduce((binptr-aveptr)**2, errptr, n, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
                errptr = sqrt(errptr * dble(nproc-1)/dble(nproc))
            endif
          enddo

          if (T1%flags(ISPXX) == 1) then
            n = T1%properties(ISPXX)%nClass
            bins => T1%chixx_r_orb_iw0(:,1)
            aves => T1%chixx_r_orb_iw0(:,avg)
            errs => T1%chixx_r_orb_iw0(:,err)
            call mpi_allreduce((bins-aves)**2, errs, n, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif

          if (T1%flagsFT(ISPXX) == 1) then
            n = 2*T1%L
            binptr => T1%chiqxx(:,:,1)
            aveptr => T1%chiqxx(:,:,avg)
            errptr => T1%chiqxx(:,:,err)
            call mpi_allreduce((binptr-aveptr)**2, errptr, n, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errptr = sqrt(errptr * dble(nproc-1)/dble(nproc))

            bins => T1%chiqxx_iw0(:,1)
            aves => T1%chiqxx_iw0(:,avg)
            errs => T1%chiqxx_iw0(:,err)
            call mpi_allreduce((bins-aves)**2, errs, 2, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            bins => T1%chiqxx_sub1(:,1)
            aves => T1%chiqxx_sub1(:,avg)
            errs => T1%chiqxx_sub1(:,err)
            call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            bins => T1%chiqxx_sub2(:,1)
            aves => T1%chiqxx_sub2(:,avg)
            errs => T1%chiqxx_sub2(:,err)
            call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            call mpi_allreduce((T1%chiqxx_iw0_sub1(1)-T1%chiqxx_iw0_sub1(avg))**2, &
              T1%chiqxx_iw0_sub1(err), 1, mpi_double, mpi_sum, mpi_comm_world, mpi_err)
            T1%chiqxx_iw0_sub1(err) = sqrt(T1%chiqxx_iw0_sub1(err) &
                                      * dble(nproc-1)/dble(nproc))

            call mpi_allreduce((T1%chiqxx_iw0_sub2(1)-T1%chiqxx_iw0_sub2(avg))**2, &
              T1%chiqxx_iw0_sub2(err), 1, mpi_double, mpi_sum, mpi_comm_world, mpi_err)
            T1%chiqxx_iw0_sub2(err) = sqrt(T1%chiqxx_iw0_sub2(err) &
                                      * dble(nproc-1)/dble(nproc))
          endif

          if (T1%flags(ISPZZ) == 1) then
            n = T1%properties(ISPZZ)%nClass
            bins => T1%chizz_r_orb_iw0(:,1)
            aves => T1%chizz_r_orb_iw0(:,avg)
            errs => T1%chizz_r_orb_iw0(:,err)
            call mpi_allreduce((bins-aves)**2, errs, n, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif

          if (T1%flagsFT(ISPZZ) == 1) then
            n = 2*T1%L
            binptr => T1%chiqzz(:,:,1)
            aveptr => T1%chiqzz(:,:,avg)
            errptr => T1%chiqzz(:,:,err)
            call mpi_allreduce((binptr-aveptr)**2, errptr, n, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errptr = sqrt(errptr * dble(nproc-1)/dble(nproc))

            bins => T1%chiqzz_iw0(:,1)
            aves => T1%chiqzz_iw0(:,avg)
            errs => T1%chiqzz_iw0(:,err)
            call mpi_allreduce((bins-aves)**2, errs, 2, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            bins => T1%chiqzz_sub1(:,1)
            aves => T1%chiqzz_sub1(:,avg)
            errs => T1%chiqzz_sub1(:,err)
            call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            bins => T1%chiqzz_sub2(:,1)
            aves => T1%chiqzz_sub2(:,avg)
            errs => T1%chiqzz_sub2(:,err)
            call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                mpi_sum, mpi_comm_world, mpi_err)
            errs = sqrt(errs * dble(nproc-1)/dble(nproc))

            call mpi_allreduce((T1%chiqzz_iw0_sub1(1)-T1%chiqzz_iw0_sub1(avg))**2, &
              T1%chiqzz_iw0_sub1(err), 1, mpi_double, mpi_sum, mpi_comm_world, mpi_err)
            T1%chiqzz_iw0_sub1(err) = sqrt(T1%chiqzz_iw0_sub1(err) &
                                      * dble(nproc-1)/dble(nproc))

            call mpi_allreduce((T1%chiqzz_iw0_sub2(1)-T1%chiqzz_iw0_sub2(avg))**2, &
              T1%chiqzz_iw0_sub2(err), 1, mpi_double, mpi_sum, mpi_comm_world, mpi_err)
            T1%chiqzz_iw0_sub2(err) = sqrt(T1%chiqzz_iw0_sub2(err) &
                                      * dble(nproc-1)/dble(nproc))
          endif

          if (T1%flags(IGFUN) == 1) then
             bins => T1%GtauAvg(:, 1)
             aves => T1%GtauAvg(:, avg)
             errs => T1%GtauAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(IGFUP) == 1) then
             bins => T1%GtupAvg(:, 1)
             aves => T1%GtupAvg(:, avg)
             errs => T1%GtupAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(IGFDN) == 1) then
             bins => T1%GtdnAvg(:, 1)
             aves => T1%GtdnAvg(:, avg)
             errs => T1%GtdnAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(ISPXX) == 1) then
             bins => T1%spinxxAvg(:, 1)
             aves => T1%spinxxAvg(:, avg)
             errs => T1%spinxxAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(ISPZZ) == 1) then
             bins => T1%spinzzAvg(:, 1)
             aves => T1%spinzzAvg(:, avg)
             errs => T1%spinzzAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(IPAIRs) == 1) then
             bins => T1%swaveAvg(:, 1)
             aves => T1%swaveAvg(:, avg)
             errs => T1%swaveAvg(:, err)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))
          endif
          if (T1%flags(IPAIRd) == 1) then
             bins => T1%Pdtau(:, 1, 1)
             aves => T1%Pdtau(:, avg, 1)
             errs => T1%Pdtau(:, err, 1)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))

             bins => T1%Pd0tau(:, 1, 1)
             aves => T1%Pd0tau(:, avg, 1)
             errs => T1%Pd0tau(:, err, 1)
             call mpi_allreduce((bins-aves)**2, errs, T1%L, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             errs = sqrt(errs * dble(nproc-1)/dble(nproc))

             call mpi_allreduce((T1%Pd(1, 1)-T1%Pd(avg, 1))**2, T1%Pd(err, 1), 1, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             T1%Pd(err, 1) = sqrt(T1%Pd(err, 1) * dble(nproc-1)/dble(nproc))

             call mpi_allreduce((T1%Pd0(1, 1)-T1%Pd0(avg, 1))**2, T1%Pd0(err, 1), 1, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             T1%Pd0(err, 1) = sqrt(T1%Pd0(err, 1) * dble(nproc-1)/dble(nproc))

             call mpi_allreduce((T1%Gammad(1, 1)-T1%Gammad(avg, 1))**2, T1%Gammad(err, 1), 1, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             T1%Gammad(err, 1) = sqrt(T1%Gammad(err, 1) * dble(nproc-1)/dble(nproc))

             call mpi_allreduce((T1%Gd_Pd0(1, 1)-T1%Gd_Pd0(avg, 1))**2, T1%Gd_Pd0(err, 1), 1, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
             T1%Gd_Pd0(err, 1) = sqrt(T1%Gd_Pd0(err, 1) * dble(nproc-1)/dble(nproc))
          endif

#      endif

    endif

  end subroutine DQMC_TDM_GetErr

  !--------------------------------------------------------------------!
  ! print all tdm quantities

  subroutine DQMC_TDM_Print(T1, ofile, OPT, OPT1, OPT2)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(in)   :: T1                 ! T1
    integer, intent(in)     :: OPT

    integer             :: i, j, t, iprop, OPT1, OPT2, b1, b2
    real(wp)            :: tmp(T1%L, 2)
    real(wp)            :: x, y
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title
    character(len=80)   :: ofile
    character(label_len) :: lab
    real(wp), dimension(1:5) :: vec

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
    enddo

    do iprop = 1, NTDMARRAY
      if (T1%flags(iprop)==1) then
        do i = 1, T1%properties(iprop)%nclass
          do j = 0, T1%L-1
             tmp(j+1, 1:2) = T1%properties(iprop)%values(i, j, T1%avg:T1%err)
          enddo
          title=pname(iprop)//" "//trim(adjustl(T1%properties(iprop)%clabel(i)))
          call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
          write(OPT,'(1x)')
        enddo
      endif
    enddo

    if (T1%flags(IPAIRd) == 1) then
       ! Pd(tau) and Pd
       do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%Pdtau(j, T1%avg:T1%err, 1)
       enddo
       title="Pd(tau)"
       call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
       write(OPT,'(1x)')

       write(OPT,"(a20,e16.8,e16.8)") 'Pd = ', T1%Pd(T1%avg, 1), T1%Pd(T1%err, 1)
       write(OPT,FMT_DBLINE)

       ! Pd0(tau) and Pd0
       do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%Pd0tau(j, T1%avg:T1%err, 1)
       enddo
       title="Pd0 (nonvertex) (tau)"
       call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
       write(OPT,'(1x)')

       write(OPT,"(a20,e16.8,e16.8)") 'Pd0 (nonvertex) = ', T1%Pd0(T1%avg, 1), T1%Pd0(T1%err, 1)
       write(OPT,FMT_DBLINE)

       write(OPT,"(a20,e16.8,e16.8)") 'accumulated Gammad = ', T1%Gammad(T1%avg, 1), T1%Gammad(T1%err, 1)
       x = 1.0/T1%Pd(T1%avg, 1) - 1.0/T1%Pd0(T1%avg, 1)
       write(OPT,"(a20,e16.8)") 'calculated Gammad = ', x
       write(OPT,FMT_DBLINE) 

       write(OPT,"(a20,e16.8,e16.8)") 'accumulated Gd*Pd0 = ', T1%Gd_Pd0(T1%avg, 1), T1%Gd_Pd0(T1%err, 1)
       y = x*T1%Pd0(T1%avg, 1)
       write(OPT,"(a20,e16.8)") 'calculated Gd*Pd0 = ', y
       write(OPT,FMT_DBLINE)
   endif

    ! print out chi_q_iw0 and chi_q(tau):
    if (T1%flagsFT(ISPXX)==1 .and. T1%flagsFT(ISPZZ)==1) then
      call DQMC_open_file('chi_q_iw0_'//adjustl(trim(ofile)),'replace', OPT1)
      call DQMC_open_file('chi_q_tau_'//adjustl(trim(ofile)),'replace', OPT2)

      write(OPT1,"(a)") &
             "chi_xx(q=0,iw=0)       err    chi_zz(q=0,iw=0)       err  &
              chi_xx(q=pi,iw=0)      err    chi_zz(q=pi,iw=0)      err  &
              chi_xx_sub1(q=0,iw=0)  err    chi_zz_sub1(q=0,iw=0)  err  &
              chi_xx_sub2(q=0,iw=0)  err    chi_zz_sub2(q=0,iw=0)  err  "

      write(OPT1,'(16(e16.8))') &
             T1%chiqxx_iw0(1,T1%avg), T1%chiqxx_iw0(1,T1%err), &
             T1%chiqzz_iw0(1,T1%avg), T1%chiqzz_iw0(1,T1%err), &
             T1%chiqxx_iw0(2,T1%avg), T1%chiqxx_iw0(2,T1%err), &
             T1%chiqzz_iw0(2,T1%avg), T1%chiqzz_iw0(2,T1%err), &
             T1%chiqxx_iw0_sub1(T1%avg), T1%chiqxx_iw0_sub1(T1%err), &
             T1%chiqzz_iw0_sub1(T1%avg), T1%chiqzz_iw0_sub1(T1%err), &
             T1%chiqxx_iw0_sub2(T1%avg), T1%chiqxx_iw0_sub2(T1%err), &
             T1%chiqzz_iw0_sub2(T1%avg), T1%chiqzz_iw0_sub2(T1%err)

      !write(OPT1,"(a)") "==================================================================================="
      !write(OPT1,"(a)") "   b   b   tau         chi_xx(q=0)           err       chi_zz(q=0)             err"
      do t = 0, T1%L-1
        write(OPT2,'(f10.5,8(f16.8))') t*T1%dtau, &
             T1%chiqxx(t,1,T1%avg), T1%chiqxx(t,1,T1%err), &
             T1%chiqzz(t,1,T1%avg), T1%chiqzz(t,1,T1%err), &
             T1%chiqxx(t,2,T1%avg), T1%chiqxx(t,2,T1%err), &
             T1%chiqzz(t,2,T1%avg), T1%chiqzz(t,2,T1%err), &
             T1%chiqxx_sub1(t,T1%avg), T1%chiqxx_sub1(t,T1%err), &
             T1%chiqzz_sub1(t,T1%avg), T1%chiqzz_sub1(t,T1%err), &
             T1%chiqxx_sub2(t,T1%avg), T1%chiqxx_sub2(t,T1%err), &
             T1%chiqzz_sub2(t,T1%avg), T1%chiqzz_sub2(t,T1%err)
      enddo
    endif

  end subroutine DQMC_TDM_Print

  !--------------------------------------------------------------------!
  ! Below print out tdm quantities separately (if needed)
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Print_local(T1, ofile, OPT1, OPT2, OPT3, OPT4, OPT5)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(in)   :: T1                 ! T1
    integer                 :: OPT1, OPT2, OPT3, OPT4, OPT5

    integer             :: i, j, b1, b2
    real(wp)            :: tmp(T1%L, 2)
    real(wp)            :: x, y
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title
    character(len=80)   :: ofile
    character(len=label_len)  :: lab
    real(wp), dimension(1:5) :: vec

    ! ... Executable ...
    if (.not.T1%compute .or. T1%flags(IGFUN) == 0) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
    enddo

    if (T1%flags(IGFUN) == 1) then
      call DQMC_open_file('Gr0_'//adjustl(trim(ofile)),'replace', OPT1)

      ! Print local G(tau)'s
      do i = 1, T1%properties(IGFUN)%nclass
        do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%properties(IGFUN)%values(i, j, T1%avg:T1%err)
        enddo
        title=pname(IGFUN)//" "//trim(adjustl(T1%properties(IGFUN)%clabel(i)))
        if (index(title, " 0.0000000   0.0000000   0.0000000") > 0) then
          call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
        endif
        ! write(OPT1,'(1x)')
      enddo

      ! Print average of local G(tau) for average N(w)
      title="average local G(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%GtauAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
    endif

    if (T1%flags(IGFUP) == 1) then
      title="average local Gup(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%GtupAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
    endif

    if (T1%flags(IGFDN) == 1) then
      title="average local Gdn(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%GtdnAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
    endif

    !############################################################################
    if (T1%flags(ISPXX) == 1) then
      call DQMC_open_file('spinxx_r0_'//adjustl(trim(ofile)),'replace', OPT2)

      ! Print local spin-xx correlation
      do i = 1, T1%properties(ISPXX)%nclass
        do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%properties(ISPXX)%values(i, j, T1%avg:T1%err)
        enddo
        title=pname(ISPXX)//" "//trim(adjustl(T1%properties(ISPXX)%clabel(i)))
        if (index(title, " 0.0000000   0.0000000   0.0000000") > 0) then
          call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT2)
        endif
        ! write(OPT1,'(1x)')
      enddo

      ! Print average of sum_r spin-xx(r,tau) for average spin-xx susceptibility
      title="average local spin-xx(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%spinxxAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT2)
    endif

    !############################################################################
    if (T1%flags(ISPZZ) == 1) then
      call DQMC_open_file('spinzz_r0_'//adjustl(trim(ofile)),'replace', OPT3)

      ! Print local spin-zz correlation
      do i = 1, T1%properties(ISPZZ)%nclass
        do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%properties(ISPZZ)%values(i, j, T1%avg:T1%err)
        enddo
        title=pname(ISPZZ)//" "//trim(adjustl(T1%properties(ISPZZ)%clabel(i)))
        if (index(title, " 0.0000000   0.0000000   0.0000000") > 0) then
          call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT3)
        endif
        ! write(OPT1,'(1x)')
      enddo

      ! Print average of sum_r spin-zz(r,tau) for average spin-zz susceptibility
      title="average local spin-zz(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%spinzzAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT3)
    endif

    !############################################################################
    if (T1%flags(IPAIRs) == 1) then
      call DQMC_open_file('swave_r0_'//adjustl(trim(ofile)),'replace', OPT4)

      ! Print local s-wave
      do i = 1, T1%properties(IPAIRs)%nclass
        do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%properties(IPAIRs)%values(i, j, T1%avg:T1%err)
        enddo
        title=pname(IPAIRs)//" "//trim(adjustl(T1%properties(IPAIRs)%clabel(i)))
        if (index(title, " 0.0000000   0.0000000   0.0000000") > 0) then
          call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT4)
        endif
        ! write(OPT1,'(1x)')
      enddo

      ! Print average of sum_r pair(r,tau) for average pair susceptibility
      title="average local swave(tau)"
      do j = 0, T1%L-1
         tmp(j+1, 1:2) = T1%swaveAvg(j, T1%avg:T1%err)
      enddo
      call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT4)
    endif

    !############################################################################
    ! print out chi_r0_iw0 and chi_r0(tau):
    if (T1%flags(ISPXX)==1 .and. T1%flags(ISPZZ)==1) then
      call DQMC_open_file('chi_r0_iw0_orb_'//adjustl(trim(ofile)),'replace', OPT5)
      write(OPT5,"(a)") "   b   b     chi_xx(r=0,iw=0)       err    chi_zz(r=0,iw=0)       err"

      do i = 1, T1%properties(ISPXX)%nclass
        vec = T1%properties(ISPXX)%vecClass(i,:)
        b1 = int(vec(1))
        b2 = int(vec(2))
        x  = vec(3)
        y  = vec(4)

        if (abs(x)<1.e-5 .and. abs(y)<1.e-5) then
          write(OPT5,'(2(i4),4(e16.8))') b1,b2, &
                 T1%chixx_r_orb_iw0(i,T1%avg), T1%chixx_r_orb_iw0(i,T1%err), &
                 T1%chizz_r_orb_iw0(i,T1%avg), T1%chizz_r_orb_iw0(i,T1%err)
        endif
      enddo
    endif

  end subroutine DQMC_TDM_Print_local

  !--------------------------------------------------------------------!
  ! Below three routines are old FT of tdm quantities
  ! 1/4/2016:
  ! old FTk considering intersite correlation within unit cells
  ! Now old FTk routines are only useful for SelfEnergy
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_GetKFTold(T1)

    type(tdm), intent(inout) :: T1

    integer              :: ip, it, n, nclass, np, npp, nk, ibin, j
    integer,     pointer :: class(:,:)
    complex(wp), pointer :: wgtftk(:,:)
    integer,     pointer :: phase(:,:)

    real(wp)   , pointer :: value(:)
    complex(wp), pointer :: valuek(:)
    complex(wp), pointer :: valuek_iw0

    if (.not.T1%compute) return
 
    !Loop over properties to Fourier transform
    do ip = 1, NTDMARRAY-1  ! -1 for exclusion of conductivity for FT
       if (T1%flags(ip)==1) then
         if (.not.associated(T1%properties(ip)%valueskold)) cycle

         ! Aliases
         n        =  T1%properties(ip)%n
         nclass   =  T1%properties(ip)%nclass
         np       =  T1%properties(ip)%np
         nk       =  T1%properties(ip)%nk
         class    => T1%properties(ip)%D
         wgtftk   => T1%properties(ip)%ftk
         phase    => T1%properties(ip)%phase
         npp      =  (np*(np+1))/2

         !Fourier transform each bin and average
         !For MPI run, nbin=1 so T1%avg=T1%nBin+1 = 2
         do ibin = T1%avg, 1, -1

            ! More aliases
            do it = 0, T1%L-1
               ! 11/20/2015: note here for MPI run
               ! binned values() (only 1 bin) are not for MC, which has been
               ! updated in JackKnife among procs in DQMC_TDM_GetErr
               ! avg value is already averaged over proc
               value  =>  T1%properties(ip)%values(:,it,ibin)
               valuek =>  T1%properties(ip)%valueskold(:,it,ibin)

               ! In util.F90:
               !Note valueskold has dimension (nk,npp), where npp=np(np+1)/2
               !npp will obtained in dqmc_GetFTk
               call dqmc_GetFTk(value, n, nclass, class, np, nk, wgtftk, phase, valuek)
            enddo

            ! FT convert to iw=0 quantity
            do j = 1,nk*npp
               valuek     =>  T1%properties(ip)%valueskold(j,:,ibin)
               valuek_iw0 =>  T1%properties(ip)%valueskold_iw0(j,ibin)
               call convert_to_iw0(valuek, valuek_iw0, T1%L, T1%dtau)
            enddo
         enddo ! Loop over bins
      endif
    enddo ! Loop over properties

  end subroutine DQMC_TDM_GetKFTold

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_GetErrKFTold(T1)

    use DQMC_MPI

    type(tdm), intent(inout) :: T1

    integer :: ip, it, nproc, n, i

    complex(wp), pointer  :: average(:), binval(:), error(:), temp(:)

    !Loop over properties to Fourier transform
    nproc = qmc_sim%size

    if (.not.T1%compute) return

    if (nproc .eq. 1) then

      do ip = 1, NTDMARRAY-1 ! -1 for exclusion of conductivity for FT
        if (T1%flags(ip)==1) then
          if (.not.associated(T1%properties(ip)%valueskold)) cycle

          do it = 0, T1%L-1

             !Note that valueskold(avg) is known from DQMC_TDM_GetKFT
             !in which FT from values(avg), here do not use JackKnife for simplicity
             !But valueskold(err) is unknown
             average  => T1%properties(ip)%valueskold(:,it,T1%avg)
             error    => T1%properties(ip)%valueskold(:,it,T1%err)

             do i = 1, T1%nbin
                binval => T1%properties(ip)%valueskold(:,it,i)
                error  = error + cmplx((real(average-binval))**2,(aimag(average-binval))**2)
             enddo 

             error = error* dble(T1%nbin-1)/dble(T1%nbin)
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)))
          enddo

          ! Then iw=0 component:
          average  => T1%properties(ip)%valueskold_iw0(:,T1%avg)
          error    => T1%properties(ip)%valueskold_iw0(:,T1%err)

          do i = 1, T1%nbin
             binval => T1%properties(ip)%valueskold_iw0(:,i)
             error  = error + cmplx((real(average-binval))**2,(aimag(average-binval))**2)
          enddo

          error = error* dble(T1%nbin-1)/dble(T1%nbin)
          error = cmplx(sqrt(real(error)),sqrt(aimag(error)))
        endif
      enddo 

    else

#  ifdef _QMC_MPI
               ! 11/20/2015: note here for MPI run
               ! binned values() (only 1 bin) are not for MC, which has been
               ! updated in JackKnife among procs in DQMC_TDM_GetErr
               ! avg value is already averaged over proc

      do ip = 1, NTDMARRAY-1
        if (T1%flags(ip)==1) then
          if (.not.associated(T1%properties(ip)%valueskold)) cycle
    
          n = T1%properties(ip)%nk * T1%properties(ip)%np*(T1%properties(ip)%np+1)/2          
          allocate(temp(n))

          do it = 0, T1%L-1
             !Note that avg value is already averaged over proc in DQMC_TDM_GetKFT
             !and binval is JackKnifed among proc
             average  => T1%properties(ip)%valueskold(:,it,T1%avg)
             binval   => T1%properties(ip)%valueskold(:,it,1)

             !Compute error: sum(y_i-avg_y)^2
             error  => T1%properties(ip)%valueskold(:,it,T1%err)
             temp   =  cmplx((real(average-binval))**2,(aimag(average-binval))**2)
             call mpi_allreduce(temp, error, n, mpi_double_complex, mpi_sum, mpi_comm_world, i)

             error  = error*dble(nproc-1)/dble(nproc)
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)))
          enddo

          ! Then iw=0 component:
          average  => T1%properties(ip)%valueskold_iw0(:,T1%avg)
          binval   => T1%properties(ip)%valueskold_iw0(:,1)

          !Compute error: sum(y_i-avg_y)^2
          error  => T1%properties(ip)%valueskold_iw0(:,T1%err)
          temp   = cmplx((real(average-binval))**2,(aimag(average-binval))**2)
          call mpi_allreduce(temp, error, n, mpi_double_complex, mpi_sum, mpi_comm_world, i)

          error  = error*dble(nproc-1)/dble(nproc)
          error = cmplx(sqrt(real(error)),sqrt(aimag(error)))

          deallocate(temp)
        endif
      enddo ! Loop over properties

#  endif

    endif


  end subroutine DQMC_TDM_GetErrKFTold

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_PrintKFTold(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT

    integer             :: i, j, k, ip, jp, iprop
    integer             :: np, npp
    complex(wp)         :: tmp(T1%L, 2)
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
       label(j) = adjustl(label(j))
    enddo

    do iprop = 1, NTDMARRAY-1 ! -1 for exclusion of conductivity for FT
       if (T1%flags(iprop)==1) then
         if (.not.associated(T1%properties(iprop)%valueskold)) cycle
         np = T1%properties(iprop)%np
         npp = (np*(np+1))/2
         do k = 1, T1%properties(iprop)%nk
            i = (k-1)*npp
            do ip = 1, np
               do jp = ip, np
                  i = i + 1
                  do j = 0, T1%L-1
                     tmp(j+1, 1:2) = T1%properties(iprop)%valueskold(i, j, T1%avg:T1%err)
                  enddo
                  write(title,'(A,i3,A,i3,A,i3,A)') 'k=',k,'   cell_site_pair=',ip-1,',',jp-1
                  title=pname(iprop)//" "//trim(adjustl(title))
                  !In util.F90, call DQMC_Print_Array for complex array
                  call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                  write(OPT,'(1x)')
               enddo
            enddo
         enddo
      endif
    enddo

    ! Then print out iw=0 component
    do iprop = 1, NTDMARRAY-1 ! -1 for exclusion of conductivity for FT
       if (T1%flags(iprop)==1) then
         if (.not.associated(T1%properties(iprop)%valueskold_iw0)) cycle
         np = T1%properties(iprop)%np
         npp = (np*(np+1))/2
         do k = 1, T1%properties(iprop)%nk
            i = (k-1)*npp                  
            write(OPT,'(A,i3,A)') pname(iprop)//' cell_site_pair   k=',k,', iw=0'

            do ip = 1, np
               do jp = ip, np
                  i = i + 1
                  tmp(1, 1:2) = T1%properties(iprop)%valueskold_iw0(i, T1%avg:T1%err)
                  write(OPT,"(I4,I4,e16.8,e16.8)") ip-1,jp-1, dble(tmp(1, 1:1)), dble(tmp(1, 2:2))
               enddo
            enddo
            write(OPT,FMT_DBLINE)
            write(OPT,'(1x)')
         enddo
      endif
    enddo

  end subroutine DQMC_TDM_PrintKFTold

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_SelfEnergy_r0(T1, tau)
    ! May.14, 2019: obtain sigma(r=0,iwm=0), here r=0 means within the
    ! same unit cell instead of intercell

    use DQMC_MPI

    type(TDM), intent(inout) :: T1
    type(gtau), intent(inout) :: tau

    real(wp),    allocatable  :: g0tau(:,:), tmp(:,:)
    complex(wp), allocatable  :: tdmg0(:,:,:), tdmg(:,:,:)
    complex(wp), allocatable  :: tdmg0w(:,:,:), tdmgw(:,:,:)
    complex(wp), pointer      :: avgGw(:,:,:), errGw(:,:,:)
    complex(wp), pointer      :: avgSE(:,:,:), errSE(:,:,:)
    complex(wp), allocatable  :: binSE(:,:,:),binGw(:,:,:)

    integer :: i, j, k, h, m, nproc, ia, ja, it, jt
    integer :: L, n, nclass, np, nt
    real(wp):: factor
    integer,     pointer :: class(:,:), ph(:,:)
    complex(wp), pointer :: ftw(:,:)

    integer, parameter  :: gflist(3) = (/IGFUN, IGFUP, IGFDN/)
    integer, parameter  :: splist(3) = (/  0, TAU_UP, TAU_DN/)

    if (.not.T1%compute) return

    nproc = qmc_sim%size

    L      =  T1%L
    n      =  T1%properties(IGFUN)%n
    nclass =  T1%properties(IGFUN)%nclass
    np     =  T1%properties(IGFUN)%np
    class  => T1%properties(IGFUN)%D
    ftw    => T1%properties(IGFUN)%ftw
    ph     => T1%properties(IGFUN)%phase

    nt = n / np
    !npp = np*(np+1)/2

    !non-interacting green's function
    allocate(tmp(nclass,0:L-1))
    allocate(g0tau(n,n))
    allocate(tdmg0(np,np,0:L-1))
    allocate(tdmg (np,np,0:L-1))

    !frequency dependent green's functions
    allocate(tdmg0w(np,np,0:L-1))
    allocate(tdmgw (np,np,0:L-1))

    !self energy: (np, np, L, 3)
    ! natom np for unit cell, 3 for G, Gup, Gdn
    allocate(T1%SEravg(np,np,0:L-1,3))
    allocate(T1%SErerr(np,np,0:L-1,3))
    allocate(binSE(np,np,0:L-1))

    allocate(T1%GrwAvg(np,np,0:L-1,3))
    allocate(T1%GrwErr(np,np,0:L-1,3))
    allocate(binGw(np,np,0:L-1))

    do h = 1, 1

      if (T1%flags(gflist(h))==1) then
       ! Get G0 and G_avg in real space for the non-interacting system
       tmp = 0.0_wp
       do m = 0, L-1
          call dqmc_Gtau_GetG0(n, tau, splist(h), m, g0tau)                                          
          do i = 1, n
             do j = 1, n
                k = class(i,j)
                tmp(k,m) = tmp(k,m) + g0tau(i,j)                                                  
             enddo                                                                                   
          enddo                                                                                      
       enddo
       do k = 1, T1%properties(h)%nClass
          tmp(k,:) = tmp(k,:) / T1%properties(IGFUN)%F(k)     
       enddo                                                                                          

        do m = 0,L-1
          do ia = 1, np
            do ja = 1, np
              !sum over translations
              do it = 1, nt
                 do jt = 1, nt
                   !Find atom which is the translation of "ja" by "it"
                   !Note that jt=it now because we only want within the same
                   !unit cell, which is different from sigma_k
                   i = (it-1) * np + ia
                   j = (jt-1) * np + ja
                   k = class(i,j)
                   tdmg0(ja,ia,m) = cmplx(tmp(k,m)) 
                   tdmg(ja,ia,m)  = cmplx(T1%properties(gflist(h))%values(k,m,T1%avg))
                 enddo
               enddo
            enddo
          enddo
        enddo

       ! Transform G0 from tau to iw, note for MPI, only ONE bin
       call convert_to_iwn(tdmg0, tdmg0w)
       call invertG(tdmg0w)

       ! Get self-energy and errorbar
       avgSE => T1%SEravg(:,:,0:L-1,h)
       errSE => T1%SErerr(:,:,0:L-1,h)

       avgGw => T1%GrwAvg(:,:,0:L-1,h)
       errGw => T1%GrwErr(:,:,0:L-1,h)

       ! Compute average self-energy; For MPI, avgSE is already averaged over proc
       call convert_to_iwn(tdmg, tdmgw)
       avgGw = tdmgw
       call invertG(tdmgw)
       avgSE = tdmg0w - tdmgw
 
       if (nproc .eq. 1) then

         do m = 1, T1%nbin
           ! Transform G from tau to iwn for bin "m"
           do ia = 1, np
             do ja = 1, np
               !sum over translations
               do it = 1, nt
                 do jt = 1, nt
                   !Find atom which is the translation of "ja" by "it"
                   i = (it-1) * np + ia
                   j = (jt-1) * np + ja
                   k = class(i,j)
                   tdmg(ja,ia,0:L-1) = cmplx(T1%properties(gflist(h))%values(k,0:L-1,m))
                 enddo
               enddo
             enddo
           enddo

           call convert_to_iwn(tdmg, tdmgw)

           ! collect G(r,w) for each bin
           binGw = tdmgw
           errGw = ZERO
           errGw = errGw + cmplx((real(binGw-avgGw))**2,(aimag(binGw-avgGw))**2)

           call invertG(tdmgw)

           ! Compute self-energy for bin
           binSE = tdmg0w - tdmgw
           errSE = ZERO
           errSE = errSE + cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
         enddo 

         errGw  = cmplx(sqrt(real(errGw)),sqrt(aimag(errGw))) * sqrt(dble(T1%nbin-1)/dble(T1%nbin))
         errSE  = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE))) * sqrt(dble(T1%nbin-1)/dble(T1%nbin))

       else

#  ifdef _QMC_MPI
         ! 11/20/2015: note here for MPI run
         ! binned values() (only 1 bin) are not for MC, which has been
         ! updated in JackKnife among procs in DQMC_TDM_GetErr
         ! avg value is already averaged over proc
         ! see DQMC_TDM_GetKFT

         m = np*np*L

         ! Compute self-energy for bin for each proc
         do ia = 1, np
           do ja = 1, np
             !sum over translations
             do it = 1, nt
               do jt = 1, nt
                 !Find atom which is the translation of "ja" by "it"
                 i = (it-1) * np + ia
                 j = (jt-1) * np + ja
                 k = class(i,j)
                 tdmg(ja,ia,0:L-1) = cmplx(T1%properties(gflist(h))%values(k,0:L-1,1))
               enddo
             enddo
           enddo
         enddo

         call convert_to_iwn(tdmg, tdmgw)
         binGw = tdmgw
         call invertG(tdmgw)
         binSE = tdmg0w - tdmgw

         !Compute error of G(k,w): sum(y_i-avg_y)^2
         tdmgw =  cmplx((real(binGw-avgGw))**2,(aimag(binGw-avgGw))**2)
         call mpi_allreduce(tdmgw, errGw, m, mpi_double_complex, mpi_sum, mpi_comm_world, i)

         errGw = cmplx(sqrt(real(errGw)),sqrt(aimag(errGw))) * sqrt(dble(nproc-1)/dble(nproc))

         !Compute error of self-energy: sum(y_i-avg_y)^2
         tdmgw =  cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
         call mpi_allreduce(tdmgw, errSE, m, mpi_double_complex, mpi_sum, mpi_comm_world, i)

         errSE = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE))) * sqrt(dble(nproc-1)/dble(nproc))
# endif
       endif
     endif   ! end if flags=1
    enddo    ! end loop h

    contains

       subroutine convert_to_iwn(tdmgtau, tdmgw)
          complex(wp), intent(in)  :: tdmgtau(np, np, 0:L-1)
          complex(wp), intent(out) :: tdmgw  (np, np, 0:L-1)
          ! Local variables
          complex(wp) :: valuetl(0:L-1), valuewl(0:L-1)
          integer     :: ipl, jpl
          complex(wp), parameter :: unum=(1.0_wp,0.0_wp), nil=(0.0_wp,0.0_wp)

          do ipl = 1, np
             do jpl = ipl, np
                valuetl = tdmgtau(jpl, ipl, 0:L-1)
                if (ipl .eq. jpl) valuetl(0) = valuetl(0) - 0.5_wp
                call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                tdmgw(ipl, jpl, 0:L-1) = valuewl
                if (ipl .ne. jpl) then
                   valuetl = conjg(tdmgtau(jpl, ipl, 0:L-1))
                   call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                   tdmgw(jpl, ipl, 0:L-1) = valuewl
                endif
             enddo
          enddo
       end subroutine convert_to_iwn

       subroutine invertG(tdmgw)
          complex(wp), target, intent(inout) :: tdmgw(np, np, 0:L-1)
          ! Local variables
          complex(wp) :: work(np)
          complex(wp), pointer :: gw(:,:)
          integer     :: ipiv(np)
          integer     :: iwl, info
          do iwl = 0, L-1
             ! Fill matrix. Note that G is complex symmetric. Not hermitian.
             gw => tdmgw(1:np, 1:np, iwl)
             call zgetrf(np, np, gw, np, ipiv, info)
             call zgetri(np, gw, np, ipiv, work, np, info)
          enddo
       end subroutine invertG

  end subroutine DQMC_TDM_SelfEnergy_r0

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_SelfEnergy_k(T1, tau)

    use DQMC_MPI

    type(TDM), intent(inout) :: T1
    type(gtau), intent(inout) :: tau

    real(wp),    allocatable  :: g0tau(:,:), tdmg0(:,:)
    complex(wp), allocatable  :: tdmg0k(:,:), tdmgk(:,:)
    complex(wp), allocatable  :: tdmg0kw(:,:,:), tdmgkw(:,:,:)
    complex(wp), pointer      :: avgGkw(:,:,:), errGkw(:,:,:)
    complex(wp), pointer      :: avgSE(:,:,:), errSE(:,:,:)
    complex(wp), allocatable  :: binSE(:,:,:),binGkw(:,:,:)

    integer :: i, j, k, h, m, nproc
    integer :: L, n, nclass, np, nk, npp, ia, ja, it, nt
    integer,     pointer :: class(:,:), ph(:,:)
    complex(wp), pointer :: ftk(:,:), ftw(:,:)

    integer, parameter  :: gflist(3) = (/IGFUN, IGFUP, IGFDN/)
    integer, parameter  :: splist(3) = (/  0, TAU_UP, TAU_DN/)

    if (.not.T1%compute) return

    nproc = qmc_sim%size

    L      =  T1%L
    n      =  T1%properties(IGFUN)%n
    nclass =  T1%properties(IGFUN)%nclass
    np     =  T1%properties(IGFUN)%np
    nk     =  T1%properties(IGFUN)%nk
    class  => T1%properties(IGFUN)%D
    ftk    => T1%properties(IGFUN)%ftk
    ftw    => T1%properties(IGFUN)%ftw
    ph     => T1%properties(IGFUN)%phase

    npp = np*(np+1)/2

    !non-interacting green's function
    allocate(g0tau(n,n))
    allocate(tdmg0(nclass,0:L-1))
    allocate(tdmg0k(nk*npp,0:L-1))

    !temporary storage
    allocate(tdmgk(npp,0:L-1))

    !frequency dependent green's functions
    allocate(tdmg0kw(np,np,0:L-1))
    allocate(tdmgkw(np,np,0:L-1))

    !self energy: (natom, natom, L, nk, 3)
    ! natom for unit cell, nk for cluster(cell) K values, 3 for G, Gup, Gdn
    allocate(T1%SEkavg(np,np,0:L-1,nk,3))
    allocate(T1%SEkerr(np,np,0:L-1,nk,3))
    allocate(binSE(np,np,0:L-1))

    allocate(T1%GkwAvg(np,np,0:L-1,nk,3))
    allocate(T1%GkwErr(np,np,0:L-1,nk,3))
    allocate(binGkw(np,np,0:L-1))

    do h = 1, 3

      if (T1%flags(gflist(h))==1) then
       ! Get G0 in real space for the non-interacting system
       tdmg0 = 0.0_wp
       do m = 0, L-1
          call dqmc_Gtau_GetG0(n, tau, splist(h), m, g0tau)
          do i = 1, n
             do j = 1, n
                k = class(i,j)
                tdmg0(k,m) = tdmg0(k,m) + g0tau(i,j) 
             enddo
          enddo
       enddo
       do k = 1, T1%properties(h)%nClass
          tdmg0(k,:) = tdmg0(k,:) / T1%properties(IGFUN)%F(k)
       enddo

       ! Get G in the k-space
       do m = 0, L-1
          call dqmc_getFTk(tdmg0(:,m), n, nclass, class, np, nk, ftk, ph, tdmg0k(:,m))
       enddo

       do k = 1, nk 

          avgSE => T1%SEkavg(:,:,0:L-1,k,h)
          errSE => T1%SEkerr(:,:,0:L-1,k,h)

          avgGkw => T1%GkwAvg(:,:,0:L-1,k,h)
          errGkw => T1%GkwErr(:,:,0:L-1,k,h)

          i = (k-1) * npp + 1
          j = k * npp

          ! Transform G0 from tau to iwn
          tdmgk = tdmg0k(i:j,0:L-1)
          call convert_to_iwn(tdmgk, tdmg0kw)
          call invertG(tdmg0kw)

          ! Transform G from tau to iwn for average, note for MPI, only ONE bin

          tdmgk = T1%properties(gflist(h))%valueskold(i:j,0:L-1,T1%avg)
          call convert_to_iwn(tdmgk, tdmgkw)

          !Compute average self-energy
          !For MPI, avgSE is already averaged over proc
          avgGkw = tdmgkw
          call invertG(tdmgkw)
          avgSE = tdmg0kw - tdmgkw
 
          if (nproc .eq. 1) then

            do m = 1, T1%nbin
               ! Transform G from tau to iwn for bin "m"
               tdmgk = T1%properties(gflist(h))%valueskold(i:j,0:L-1,m)
               call convert_to_iwn(tdmgk, tdmgkw)

               ! collect G(k,w) for each bin
               binGkw = tdmgkw
               errGkw = ZERO
               errGkw = errGkw + cmplx((real(binGkw-avgGkw))**2,(aimag(binGkw-avgGkw))**2)

               call invertG(tdmgkw)

               ! Compute self-energy for bin
               binSE = tdmg0kw - tdmgkw
               errSE = ZERO
               errSE = errSE + cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
            enddo 

            errGkw = cmplx(sqrt(real(errGkw)),sqrt(aimag(errGkw))) * sqrt(dble(T1%nbin-1)/dble(T1%nbin))
            errSE  = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE))) * sqrt(dble(T1%nbin-1)/dble(T1%nbin))

          else

#  ifdef _QMC_MPI
            ! 11/20/2015: note here for MPI run
            ! binned values() (only 1 bin) are not for MC, which has been
            ! updated in JackKnife among procs in DQMC_TDM_GetErr
            ! avg value is already averaged over proc
            ! see DQMC_TDM_GetKFT

            m = np*np*L

            ! Compute self-energy for bin for each proc
            tdmgk = T1%properties(gflist(h))%valueskold(i:j,0:L-1,1)
            call convert_to_iwn(tdmgk, tdmgkw)
            binGkw = tdmgkw
            call invertG(tdmgkw)
            binSE = tdmg0kw - tdmgkw

            !Compute error of G(k,w): sum(y_i-avg_y)^2
            tdmgkw =  cmplx((real(binGkw-avgGkw))**2,(aimag(binGkw-avgGkw))**2)
            call mpi_allreduce(tdmgkw, errGkw, m, mpi_double_complex, mpi_sum, mpi_comm_world, i)

            errGkw = cmplx(sqrt(real(errGkw)),sqrt(aimag(errGkw))) * sqrt(dble(nproc-1)/dble(nproc))

            !Compute error of self-energy: sum(y_i-avg_y)^2
            tdmgkw =  cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
            call mpi_allreduce(tdmgkw, errSE, m, mpi_double_complex, mpi_sum, mpi_comm_world, i)

            errSE = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE))) * sqrt(dble(nproc-1)/dble(nproc))
# endif
          endif
       enddo ! end loop k
     endif   ! end if flags=1
    enddo    ! end loop h

    contains

       subroutine convert_to_iwn(tdmgtau, tdmgw)
          complex(wp), intent(in)  :: tdmgtau(npp, 0:L-1)
          complex(wp), intent(out) :: tdmgw(np, np, 0:L-1)
          ! Local variables
          complex(wp) :: valuetl(0:L-1), valuewl(0:L-1)
          integer     :: ipl, jpl, ijpl
          complex(wp), parameter :: unum=(1.0_wp,0.0_wp), nil=(0.0_wp,0.0_wp)
          ijpl = 0
          do ipl = 1, np
             do jpl = ipl, np
                ijpl = ijpl + 1
                valuetl = tdmgtau(ijpl, 0:L-1)

                ! why is there this correction for local equal-time G(k,tau)?
                if (ipl .eq. jpl) valuetl(0) = valuetl(0) - 0.5_wp

                call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                tdmgw(ipl, jpl, 0:L-1) = valuewl
                if (ipl .ne. jpl) then
                   valuetl = conjg(tdmgtau(ijpl, 0:L-1))
                   call zgemv('N', L, L, unum, ftw, L, valuetl, 1, nil, valuewl, 1)
                   tdmgw(jpl, ipl, 0:L-1) = valuewl
                endif
             enddo
          enddo
       end subroutine convert_to_iwn

     
       subroutine invertG(tdmgw)
          complex(wp), target, intent(inout) :: tdmgw(np, np, 0:L-1)
          ! Local variables
          complex(wp) :: work(np)
          complex(wp), pointer :: gw(:,:)
          integer     :: ipiv(np)
          integer     :: iwl, info
          do iwl = 0, L-1
             ! Fill matrix. Note that G is complex symmetric. Not hermitian.
             gw => tdmgw(1:np, 1:np, iwl)
             call zgetrf(np, np, gw, np, ipiv, info)
             call zgetri(np, gw, np, ipiv, work, np, info)
          enddo
       end subroutine invertG

  end subroutine DQMC_TDM_SelfEnergy_k

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_Print_SelfEnergy(T1, OPT)

    use DQMC_MPI

    type(TDM), intent(in)    :: T1
    integer, intent(in)      :: OPT

    character(len=10)    :: label(T1%L)
    character(len=50)    :: title

    complex(wp)          :: tmp(T1%L,2)
    integer, parameter  :: gflist(3) = (/IGFUN, IGFUP, IGFDN/)
    real(wp), parameter :: pi = 3.1415926535898

    integer             :: i, j, h, k

    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    !Fermionic Matsubara frequency w_n = (2n+1)*pi/beta
    do j = 1, T1%L
       write(label(j),'(f10.5)') (2*j-1)*pi/(T1%dtau*T1%L)
       label(j) = adjustl(label(j))
    enddo

    ! print self-energy(k,iwn)
    do h = 1, 3
      if (T1%flags(gflist(h))==1) then
        do k = 1, T1%properties(IGFUN)%nk
          do i = 1, T1%properties(IGFUN)%np
             do j = 1, T1%properties(IGFUN)%np
                tmp(1:T1%L,1) = T1%SEkavg(i,j,0:T1%L,k,h)
                tmp(1:T1%L,2) = T1%SEkerr(i,j,0:T1%L,k,h)
                write(title,'(A,i3)') trim(pname(gflist(h)))//" SelfEn k=", k
                write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   cell_site_pair=',i-1,',',j-1  ! site index starts from 0
                call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                write(OPT,'(1x)')
             enddo
          enddo
        enddo
      endif
    enddo

    ! print G(k,iwn)
    do h = 1, 3
      if (T1%flags(gflist(h))==1) then
        do k = 1, T1%properties(IGFUN)%nk
          do i = 1, T1%properties(IGFUN)%np
             do j = 1, T1%properties(IGFUN)%np
                tmp(1:T1%L,1) = T1%GkwAvg(i,j,0:T1%L,k,h)
                tmp(1:T1%L,2) = T1%GkwErr(i,j,0:T1%L,k,h)
                write(title,'(A,i3)') trim(pname(gflist(h)))//" k=", k
                write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   cell_site_pair=',i-1,',',j-1  ! site index starts from 0
                call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                write(OPT,'(1x)')
             enddo
          enddo
        enddo
      endif
    enddo

    ! print self-energy(r,iwn)
    do h = 1, 3
      if (T1%flags(gflist(h))==1) then
          do i = 1, T1%properties(IGFUN)%np
             do j = 1, T1%properties(IGFUN)%np
                tmp(1:T1%L,1) = T1%SEravg(i,j,0:T1%L,h)
                tmp(1:T1%L,2) = T1%SErerr(i,j,0:T1%L,h)
                write(title,'(A,i3)') trim(pname(gflist(h)))//" SelfEn r"
                write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   cell_site_pair=',i-1,',',j-1  ! site index starts from 0
                call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                write(OPT,'(1x)')
             enddo
          enddo
      endif
    enddo

    ! print G(r,iwn)
    do h = 1, 3
      if (T1%flags(gflist(h))==1) then
          do i = 1, T1%properties(IGFUN)%np
             do j = 1, T1%properties(IGFUN)%np
                tmp(1:T1%L,1) = T1%GrwAvg(i,j,0:T1%L,h)
                tmp(1:T1%L,2) = T1%GrwErr(i,j,0:T1%L,h)
                write(title,'(A,i3)') trim(pname(gflist(h)))//" r"
                write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   cell_site_pair=',i-1,',',j-1  ! site index starts from 0
                call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                write(OPT,'(1x)')
             enddo
          enddo
      endif
    enddo

  end subroutine DQMC_TDM_Print_SelfEnergy

  !--------------------------------------------------------------------!
  ! Below three routines for new FT of tdm quantities
  ! 1/4/2016:
  ! original FTk considering intersite correlation within unit cells
  ! Here for any unit cell, compute FT for (kx,ky) = (0:n/2, 0:n/2) 
  ! Now old FTk routines are only useful for SelfEnergy
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_GetKFT(T1, Hub)
    use DQMC_Hubbard
    ! For MPI run
    ! binned values() (only 1 bin) are not from MC, which has been
    ! updated in JackKnife among procs in DQMC_TDM_GetErr
    ! avg value is already averaged over proc, similar to DQMC_TDM_GetKFT

    type(tdm), intent(inout) :: T1
    type(Hubbard), intent(in) :: Hub

    integer              :: i, ip, kx, ky, n, nclass, ibin
    integer,     pointer :: F(:)
    real(wp),    pointer :: vec(:,:)
    real(wp),    pointer :: value(:), valuek(:)
    real(wp)             :: twopi, L, qx, qy, factor
    real(wp),    pointer :: ql(:,:,:)

    if (.not.T1%compute) return

    ! Aliases
    n        =  Hub%S%nsite  
    nclass   =  Hub%S%nclass
    F        => Hub%S%F
    vec      => Hub%S%vecClass
    twopi    =  2.d0*acos(-1.0_wp)
    L        =  sqrt(real(n/T1%norb))  ! for square lattice only!!!

    allocate(ql(0:T1%NkFT, 0:T1%NkFT, nclass))

    ! for averaging weighted summation of sum_ly exp(i*qy*ly) * curr-curr
    factor = 1.d0/n

    !Prepare FT coefficients
    do kx = 0, T1%NkFT
       do ky = 0, T1%NkFT
          do i = 1, nclass
             qx = kx*twopi/L
             qy = ky*twopi/L

             ql(kx,ky,i) = F(i) * dcos(qx*vec(i,1)+qy*vec(i,2))
          enddo
       enddo
    enddo
 
    !Loop over properties to Fourier transform
    do ip = 1, NTDMARRAY
       if (T1%flagsFT(ip)==1) then
         if (.not.associated(T1%properties(ip)%valuesk)) cycle

         do ibin = T1%avg, 1, -1
            do kx = 0, T1%NkFT
               do ky = 0, T1%NkFT
                  do i = 1, nclass
                     value  =>  T1%properties(ip)%values(i,:,ibin)
                     valuek =>  T1%properties(ip)%valuesk(kx,ky,:,ibin)
                     valuek = valuek + ql(kx,ky,i)*value
                  enddo
               enddo
            enddo
         enddo

         T1%properties(ip)%valuesk = T1%properties(ip)%valuesk * factor         
       endif
    enddo

  end subroutine DQMC_TDM_GetKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_GetErrKFT(T1)

    use DQMC_MPI

    type(tdm), intent(inout) :: T1

    integer :: ip, it, nproc, n, i

    real(wp), pointer  :: average(:,:), binval(:,:), error(:,:), temp(:,:)
 
    !Loop over properties to Fourier transform
    nproc = qmc_sim%size

    if (.not.T1%compute) return

    if (nproc .eq. 1) then

      do ip = 1, NTDMARRAY
        if (T1%flagsFT(ip)==1) then
          if (.not.associated(T1%properties(ip)%valuesk)) cycle

          do it = 0, T1%L-1

             !Note that valuesk(avg) is known from DQMC_TDM_GetKFT
             !in which FT from values(avg), here do not use JackKnife for simplicity
             !But valuesk(err) is unknown
             average  => T1%properties(ip)%valuesk(:,:,it,T1%avg)
             error    => T1%properties(ip)%valuesk(:,:,it,T1%err)

             do i = 1, T1%nbin
                binval => T1%properties(ip)%valuesk(:,:,it,i)
                error  = error + (average-binval)**2
             enddo 

             error  = error* dble(T1%nbin-1)/dble(T1%nbin)
          enddo
        endif
      enddo 

    else

#  ifdef _QMC_MPI
      ! 11/20/2015: note here for MPI run
      ! binned values() (only 1 bin) are not for MC, which has been
      ! updated in JackKnife among procs in DQMC_TDM_GetErr
      ! avg value is already averaged over proc
      ! see DQMC_TDM_GetKFT

      do ip = 1, NTDMARRAY
        if (T1%flagsFT(ip)==1) then
          if (.not.associated(T1%properties(ip)%valuesk)) cycle
    
          n = (T1%NkFT+1)**2
          allocate(temp(0:T1%NkFT,0:T1%NkFT))
          
          do it = 0, T1%L-1
             !Note that avg value is already averaged over proc in DQMC_TDM_GetKFT
             !and binval is JackKnifed among proc
             average  => T1%properties(ip)%valuesk(:,:,it,T1%avg)
             binval   => T1%properties(ip)%valuesk(:,:,it,1)

             !Compute error: sum(y_i-avg_y)^2
             error  => T1%properties(ip)%valuesk(:,:,it,T1%err)
             temp   =  (average-binval)**2

             call mpi_allreduce(temp, error, n, mpi_double, mpi_sum, mpi_comm_world, i)
             error  = error*dble(nproc-1)/dble(nproc)
          enddo

          deallocate(temp)
        endif
      enddo ! Loop over properties

#  endif

    endif


  end subroutine DQMC_TDM_GetErrKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_PrintKFT(T1, OPT, ofile, OPT1, OPT2)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM), intent(in)   :: T1                 ! T1
    integer                  :: OPT, OPT1, OPT2

    integer             :: j, kx, ky, iprop
    real(wp)            :: tmp(T1%L, 2)
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title
    character(len=slen)   :: ofile

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
       label(j) = adjustl(label(j))
    enddo

    do iprop = 1, NTDMARRAY
       if (T1%flagsFT(iprop)==1) then
         if (.not.associated(T1%properties(iprop)%valuesk)) cycle
         do kx = 0, T1%NkFT
            do ky = 0, T1%NkFT
               do j = 0, T1%L-1
                  tmp(j+1, 1:2) = T1%properties(iprop)%valuesk(kx, ky, j, T1%avg:T1%err)
               enddo
               write(title,'(a5,i3,a5,i3)') 'kx=',kx,'ky=',ky
               title=pname(iprop)//" "//trim(adjustl(title))
               call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
               write(OPT,'(1x)')

               !print out k=0 Gtau and s-wave pairing correlation for maxent
               if (iprop==IGFUN .and. kx==0 .and. ky==0) then
                  call DQMC_open_file('Gk0_'//adjustl(trim(ofile)),'replace', OPT1)
                  call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
               endif
               if (iprop==IPAIRs .and. kx==0 .and. ky==0) then
                  call DQMC_open_file('swave_k0_'//adjustl(trim(ofile)),'replace', OPT2)
                  call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT2)
               endif

            enddo
         enddo
      endif
    enddo

  end subroutine DQMC_TDM_PrintKFT

  !--------------------------------------------------------------------!                                                       
                                                                                                                               
  subroutine DQMC_TDM_PrintKFT_allBins(T1, OPT, ofile, OPT1, OPT2)
    use dqmc_mpi
    !
    ! Purpose                                                                                                                  
    ! =======
    !    This subroutine prints properties to file                                                                             
    !    to print out all binned Gtau(k) and s-wave pairing correlation(k) 
    !    for Sandvik's new analytical continuation method
    !
    ! Arguments
    ! =========                                                                                                                
    !
    type(TDM), intent(in)   :: T1                 ! T1                                                                         
    integer                 :: OPT, OPT1, OPT2                                                                                
                                                                                                                               
    integer             :: j, kx, ky, iprop, ibin
    real(wp)            :: tmp(T1%L, 2)                                                                                        
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title                                                                                               
    character(len=slen)   :: ofile                                                                                               
  
    ! ... Executable ...                                                                                                       
    if (.not.T1%compute) return
  
    if (qmc_sim%rank .ne. 0) return                                                                                            
    
    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
       label(j) = adjustl(label(j))                                                                                            
    enddo                                                                                                                      

    do iprop = 1, NTDMARRAY
       if (T1%flagsFT(iprop)==1) then
         if (.not.associated(T1%properties(iprop)%valuesk)) cycle                                                              
         do kx = 0, T1%NkFT
            do ky = 0, T1%NkFT                                                                                               
               do j = 0, T1%L-1                                                                                                
                  tmp(j+1, 1:2) = T1%properties(iprop)%valuesk(kx, ky, j, T1%avg:T1%err)                                       
               enddo
               write(title,'(a5,i3,a5,i3)') 'kx=',kx,'ky=',ky                                                                  
               title=pname(iprop)//" "//trim(adjustl(title))
               call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)                                    
               write(OPT,'(1x)')
            enddo
         enddo
      endif
    enddo

    !print out all binned Gtau(k) and s-wave pairing correlation(k) for Sandvik's new analytical continuation method                                              
    if (T1%flagsFT(IGFUN)==1) then                                                                  
       call DQMC_open_file('Gk_bin_'//adjustl(trim(ofile)),'replace', OPT1)                                            

       ! bin-by-bin, for each bin, print kx,ky first and then G(tau)
       do ibin = 1, T1%nBin
         do kx = 0, T1%NkFT
           do ky = 0, T1%NkFT                                                                                               
             write(OPT1,'(i3,i3)') kx, ky
             do j = 0, T1%L-1                                                                                                
                write(OPT1,'(e16.8)')  T1%properties(IGFUN)%valuesk(kx, ky, j, ibin)                                       
             enddo
             write(OPT1,'(1x)')
             enddo
          enddo
       enddo
    endif

    if (T1%flagsFT(IPAIRs)==1) then           
       call DQMC_open_file('swavek_bin_'//adjustl(trim(ofile)),'replace', OPT2)        

       ! bin-by-bin, for each bin, print kx,ky first and then s-wave(tau)                                                           
       do ibin = 1, T1%nBin 
         do kx = 0, T1%NkFT                                                                                                  
           do ky = 0, T1%NkFT                                                                                               
             write(OPT2,'(i3,i3)') kx, ky                                                                                      
             do j = 0, T1%L-1                                                                                                  
                write(OPT2,'(e16.8)')  T1%properties(IPAIRs)%valuesk(kx, ky, j, ibin)                                           
             enddo                                                                                                             
             write(OPT2,'(1x)')                                                                                                 
             enddo
          enddo
       enddo
    endif
  
  end subroutine DQMC_TDM_PrintKFT_allBins  

  !--------------------------------------------------------------------!
  ! Below three routines for curr-curr(qx=0,qy;iwn) is estimated by 
  ! linear extrapolation of two smallest qy
  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_currDs(T1,Hub)
    use DQMC_Hubbard
    ! curr-curr(qx=0,qy;iwn) is estimated by linear extrapolation of two smallest qy
    ! same as RTScode dishbpar2.f, only need cos real part for FT coefficients
    ! Whether to call this routine is specified in ggeom.F90

    ! For MPI run
    ! binned values() (only 1 bin) are not from MC, which has been
    ! updated in JackKnife among procs in DQMC_TDM_GetErr
    ! avg value is already averaged over proc, similar to DQMC_TDM_GetKFT

    type(tdm), intent(inout) :: T1
    type(Hubbard), intent(in) :: Hub

    integer              :: i, k, n, nclass, ibin, Nq
    integer,     pointer :: F(:)
    real(wp),    pointer :: vec(:,:)
    real(wp),    pointer :: value(:)
    real(wp),    pointer :: xtemp1(:), ytemp1(:), xtemp2(:), ytemp2(:)
    complex(wp), pointer :: valuew(:)
    real(wp)             :: twopi, L, q, factor
    real(wp)             :: ql

    if (.not.T1%compute .or. T1%flags(ICOND)==0) return

    ! Aliases
    n        =  Hub%S%nsite
    nclass   =  Hub%S%nclass
    F        => Hub%S%F
    vec      => Hub%S%vecClass
    L        =  sqrt(real(n))     ! for square lattice only!!!
    twopi    =  2.d0*acos(-1.0_wp)

    ! for averaging weighted summation of sum_ly exp(i*qy*ly) * curr-curr
    factor = 1.d0/n

    ! 3 smallest qx, qy. For MPI run, nbin=1 so T1%avg=T1%nBin+1 = 2
   ! final 3 denotes three cases: nospline, qwspline, wqspline for iwn=0
    Nq = 3
    allocate(T1%Dsqx(0:Nq,3,T1%err,3))
    allocate(T1%Dsqy(0:Nq,3,T1%err,3))
    allocate(T1%Ds  (3,T1%err,3))  ! <-Kx>-curr

    ! for usage in qwspline case
    allocate(xtemp1(0:T1%L-1))
    allocate(xtemp2(0:T1%L-1))
    allocate(ytemp1(0:T1%L-1))
    allocate(ytemp2(0:T1%L-1))
    ! totally L wn values, but actually only need iwn=0
    allocate(valuew(0:T1%L-1))  

    T1%Dsqx = 0.d0
    T1%Dsqy = 0.d0
    T1%Ds   = 0.d0

    !Get Dsqy for each bin and average
    do ibin = T1%avg, 1, -1

       do k = 0, Nq  ! including qy=0
          q = k*twopi/L

          !for usage in qwspline case
          xtemp1  = 0.d0
          xtemp2  = 0.d0
          ytemp1  = 0.d0
          ytemp2  = 0.d0
          do i = 1, nclass

!write(*,'(a6,i2,a3,f4.1,a3,f4.1,a3,i4)') "class",i," x=",vec(i,1)," y=",vec(i,2), " F=",F(i)

             !First do FT to qy: vec(i,2)=ly for (qx=0,qy)
             ql = F(i)*dcos(q*vec(i,2))  

             !only need iwn=0 component
             value => T1%properties(ICONDup)%values(i,0:T1%L-1,ibin)
             call DQMC_getFTw(value,valuew,T1%L,Hub%dtau,1,T1%L/2)
             T1%Dsqy(k,1,ibin,1) = T1%Dsqy(k,1,ibin,1) + ql*sum(value)*Hub%dtau
             T1%Dsqy(k,1,ibin,3) = T1%Dsqy(k,1,ibin,3) + ql*valuew(0)
             ytemp1 = ytemp1 + ql*value

             value => T1%properties(ICONDdn)%values(i,0:T1%L-1,ibin)
             call DQMC_getFTw(value,valuew,T1%L,Hub%dtau,1,T1%L/2)
             T1%Dsqy(k,2,ibin,1) = T1%Dsqy(k,2,ibin,1) + ql*sum(value)*Hub%dtau
             T1%Dsqy(k,2,ibin,3) = T1%Dsqy(k,2,ibin,3) + ql*valuew(0)
             ytemp2 = ytemp2 + ql*value
 
             !------------------------------------------------------------
             !FT to qx: vec(i,1)=lx for (qx,qy=0)
             ql = F(i)*dcos(q*vec(i,1)) 

             !only need iwn=0 component
             value => T1%properties(ICONDup)%values(i,0:T1%L-1,ibin)
             call DQMC_getFTw(value,valuew,T1%L,Hub%dtau,1,T1%L/2)
             T1%Dsqx(k,1,ibin,1) = T1%Dsqx(k,1,ibin,1) + ql*sum(value)*Hub%dtau
             T1%Dsqx(k,1,ibin,3) = T1%Dsqx(k,1,ibin,3) + ql*valuew(0)
             xtemp1 = xtemp1 + ql*value

             value => T1%properties(ICONDdn)%values(i,0:T1%L-1,ibin)
             call DQMC_getFTw(value,valuew,T1%L,Hub%dtau,1,T1%L/2)
             T1%Dsqx(k,2,ibin,1) = T1%Dsqx(k,2,ibin,1) + ql*sum(value)*Hub%dtau
             T1%Dsqx(k,2,ibin,3) = T1%Dsqx(k,2,ibin,3) + ql*valuew(0)
             xtemp2 = xtemp2 + ql*value
          enddo  ! loop over nclass
  
          !Get iwn=0 for the case qwspline
          call DQMC_getFTw(ytemp1,valuew,T1%L,Hub%dtau,1,T1%L/2)
          T1%Dsqy(k,1,ibin,2) = valuew(0)
          call DQMC_getFTw(ytemp2,valuew,T1%L,Hub%dtau,1,T1%L/2)
          T1%Dsqy(k,2,ibin,2) = valuew(0)
          call DQMC_getFTw(xtemp1,valuew,T1%L,Hub%dtau,1,T1%L/2)
          T1%Dsqx(k,1,ibin,2) = valuew(0)
          call DQMC_getFTw(xtemp2,valuew,T1%L,Hub%dtau,1,T1%L/2)
          T1%Dsqx(k,2,ibin,2) = valuew(0)

          !total=up+dn contributions
          T1%Dsqy(k,3,ibin,:) = T1%Dsqy(k,1,ibin,:) + T1%Dsqy(k,2,ibin,:) 
          T1%Dsqx(k,3,ibin,:) = T1%Dsqx(k,1,ibin,:) + T1%Dsqx(k,2,ibin,:)
  
       enddo  ! Loop over q

       ! averaging weighted summation over # of lattice sites n
       T1%Dsqy(:,1:3,ibin,:) = T1%Dsqy(:,1:3,ibin,:) * factor
       T1%Dsqx(:,1:3,ibin,:) = T1%Dsqx(:,1:3,ibin,:) * factor

       !curr-curr(qx=0,qy;iwn) is estimated by linear extrapolation of two smallest qy
       !22,23,24 denotes Kx, Kx_up and Kx_dn, from dqmc_phy0.F90
       T1%Ds(1,ibin,:) = Hub%P0%meas(23,ibin) - (2.d0*T1%Dsqy(1,1,ibin,:)-T1%Dsqy(2,1,ibin,:))
       T1%Ds(2,ibin,:) = Hub%P0%meas(24,ibin) - (2.d0*T1%Dsqy(1,2,ibin,:)-T1%Dsqy(2,2,ibin,:))
       T1%Ds(3,ibin,:) = T1%Ds(1,ibin,:) + T1%Ds(2,ibin,:)

    enddo     ! Loop over bins

  end subroutine DQMC_TDM_currDs

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_currDs_Err(T1)

    use DQMC_MPI

    type(tdm), intent(inout) :: T1

    integer :: nproc, i, j, s, Nq

    real(wp), pointer  :: average(:), binval(:), error(:), temp(:)

    if (.not.T1%compute .or. T1%flags(ICOND)==0) return

    nproc = qmc_sim%size

    if (nproc .eq. 1) then

       !Note that values(avg) is known from DQMC_TDM_GetErr, values(err) is unknown
       !in which FT from values(avg), here do not use JackKnife for simplicity

       ! loop over nospline, qwspline, wqspline
       do j = 1,3
   
          ! loop over spin up, dn, total
          do s = 1,3
             average  => T1%Dsqy(:,s,T1%avg,j)
             error    => T1%Dsqy(:,s,T1%err,j)
             do i = 1, T1%nbin
                binval => T1%Dsqy(:,s,i,j)
                error  =  error + (average-binval)**2
             enddo
             error  = error* dble(T1%nbin-1)/dble(T1%nbin)
             !--------------------------------------------------
             average  => T1%Dsqx(:,s,T1%avg,j)
             error    => T1%Dsqx(:,s,T1%err,j)
             do i = 1, T1%nbin
                binval => T1%Dsqx(:,s,i,j)
                error  =  error + (average-binval)**2
             enddo
             error  = error* dble(T1%nbin-1)/dble(T1%nbin)
             !--------------------------------------------------
             do i = 1, T1%nbin
                T1%Ds(s,T1%err,j) = T1%Ds(s,T1%err,j) + (T1%Ds(s,T1%avg,j) - T1%Ds(s,i,j))**2
             enddo
             T1%Ds(s,T1%err,j) = T1%Ds(s,T1%err,j)* dble(T1%nbin-1)/dble(T1%nbin)
          enddo
       enddo

    else

#   ifdef _QMC_MPI
      ! 12/27/2015: note here for MPI run
      ! binned values() (only 1 bin) are not for MC, which has been
      ! updated in JackKnife among procs in DQMC_TDM_GetErr
      ! avg value is already averaged over proc, see DQMC_TDM_GetKFT

      Nq = 3
      allocate(temp(0:Nq))  ! two lowest qy 

      ! loop over nospline, qwspline, wqspline
      do j = 1,3
         ! loop over spin up, dn, total
         do s = 1,3
            average  => T1%Dsqy(0:Nq,s,T1%avg,j)
            binval   => T1%Dsqy(0:Nq,s,1,j)
            ! Compute error: sum(y_i-avg_y)^2
            error  => T1%Dsqy(0:Nq,s,T1%err,j)
            temp   =  (average-binval)**2
            call mpi_allreduce(temp, error, Nq+1, mpi_double, mpi_sum, mpi_comm_world, i)

            error  = error*dble(nproc-1)/dble(nproc)
            !-------------------------------------------------------
            average  => T1%Dsqx(0:Nq,s,T1%avg,j)
            binval   => T1%Dsqx(0:Nq,s,1,j)
            ! Compute error: sum(y_i-avg_y)^2
            error  => T1%Dsqx(0:Nq,s,T1%err,j)
            temp   =  (average-binval)**2
            call mpi_allreduce(temp, error, Nq+1, mpi_double, mpi_sum, mpi_comm_world, i)

            error  = error*dble(nproc-1)/dble(nproc)
            !-------------------------------------------------------
            call mpi_allreduce((T1%Ds(s,T1%avg,j)-T1%Ds(s,1,j))**2, &
                          T1%Ds(s,T1%err,j), 1, mpi_double, mpi_sum, mpi_comm_world, i)
            T1%Ds(s,T1%err,j) = T1%Ds(s,T1%err,j)*dble(nproc-1)/dble(nproc)
         enddo
      enddo
      deallocate(temp)

#    endif

    endif

  end subroutine DQMC_TDM_currDs_Err

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM_currDs_Print(T1, ofile, OPT, Dsqy)
    use dqmc_mpi

    type(TDM), intent(in)   :: T1                 ! T1
    integer                  :: OPT
    integer,    intent(in)   :: Dsqy

    integer  :: i, j, Nq
    character(len=18) :: spline(3)
    character(len=slen) :: ofile

    if (.not.T1%compute .or. T1%flags(ICOND)==0) return
    if (qmc_sim%rank .ne. 0) return

    if (Dsqy>0) then
      call DQMC_open_file('current_'//adjustl(trim(ofile)),'replace', OPT)
    endif

    Nq = 3
    spline(1) = "nospline for iwn=0"
    spline(2) = "qwspline for iwn=0"
    spline(3) = "wqspline for iwn=0"

!--------------------------------------------------------------------------------------
  ! loop over nospline, qwspline, wqspline
  do j = 1,3
    write(OPT,'(a45)') "Ds/pi = <-Kx> - curr-curr(qx=0, qy->0; iwn=0)"
    write(OPT,'(a50)')"=================================================================="
    write(OPT,'(a18)') spline(j)
    write(OPT,'(a50)')"=================================================================="
    write(OPT,'(a3,e16.8,a2,e16.8)') "Ds=", T1%Ds(3,T1%avg,j), " ", T1%Ds(3,T1%err,j)
    write(OPT,'(a3,e16.8,a2,e16.8)') "UP ", T1%Ds(1,T1%avg,j), " ", T1%Ds(1,T1%err,j)
    write(OPT,'(a3,e16.8,a2,e16.8)') "DN ", T1%Ds(2,T1%avg,j), " ", T1%Ds(2,T1%err,j)
  enddo
  write(OPT,'(a30)') " "
  
  ! loop over nospline, qwspline, wqspline
  do j = 1,3
    write(OPT,'(a50)')"=================================================================="
    write(OPT,'(a43)') "qy*Ly/2*pi       curr-curr(qx=0, qy; iwn=0)"
    write(OPT,'(a50)')"=================================================================="
    write(OPT,'(a18)') spline(j)
    write(OPT,'(a50)')"=================================================================="

    write(OPT,'(a6)') "Dsqy"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqy(i,3,T1%avg,j), " ", T1%Dsqy(i,3,T1%err,j)
    enddo
    write(OPT,'(a50)')"------------------------------------------------------------------"

    write(OPT,'(a6)') "Dsqyup"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqy(i,1,T1%avg,j), " ", T1%Dsqy(i,1,T1%err,j)
    enddo
    write(OPT,'(a50)')"------------------------------------------------------------------"

    write(OPT,'(a6)') "Dsqydn"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqy(i,2,T1%avg,j), " ", T1%Dsqy(i,2,T1%err,j)
    enddo
    write(OPT,'(a50)')"=================================================================="

!--------------------------------------------------------------------------------------
    write(OPT,'(a43)') "qx*Lx/2*pi       curr-curr(qx, qy=0; iwn=0)"
    write(OPT,'(a50)')"=================================================================="

    write(OPT,'(a6)') "Dsqx"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqx(i,3,T1%avg,j), " ", T1%Dsqx(i,3,T1%err,j)
    enddo
    write(OPT,'(a50)')"------------------------------------------------------------------"

    write(OPT,'(a6)') "Dsqxup"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqx(i,1,T1%avg,j), " ", T1%Dsqx(i,1,T1%err,j)
    enddo
    write(OPT,'(a50)')"------------------------------------------------------------------"

    write(OPT,'(a6)') "Dsqxdn"
    do i = 0,Nq
       write(OPT,'(i2,a2,e16.8,a2,e16.8)') i, " ", T1%Dsqx(i,2,T1%avg,j), " ", T1%Dsqx(i,2,T1%err,j)
    enddo
    write(OPT,'(a50)') " "
  enddo
  
  end subroutine DQMC_TDM_currDs_Print

end module DQMC_TDM
