module DQMC_TDM1
#include "dqmc_include.h"

  ! 08/16/2015:
  ! Add features of conductivity calculation using the standard definition
  ! instead of Simone's tricks using link correlation
  ! Also add variable flagcond: specify in input file if to compute conductivity

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
     complex(wp), pointer :: valuesk(:,:,:)

     real(wp),    pointer :: tlink(:,:)

     character(label_len), pointer :: clabel(:)

  end type tdmarray

  ! Index of the array varaiables
  integer, parameter  :: NTDMARRAY = 8 
  integer, parameter  :: IGFUN = 1
  integer, parameter  :: IGFUP = 2
  integer, parameter  :: IGFDN = 3
  integer, parameter  :: ISPXX = 4
  integer, parameter  :: ISPZZ = 5
  integer, parameter  :: IDENS = 6
  integer, parameter  :: IPAIR = 7
  integer, parameter  :: ICOND = 8
  integer             :: ntdm = NTDMARRAY ! if needed, decrease it for not cond, d-wave sus etc.

  ! Index of the array varaiables
  character(len=12), parameter :: &
            pname(NTDMARRAY) = (/ &
                  "Gfun        ", &
                  "Gfun up     ", &
                  "Gfun dn     ", &
                  "SxSx        ", &
                  "SzSz        ", &
                  "Den-Den     ", &
                  "S-wave      ", &
                  "Conductivity" /)

  type TDM1
     integer  :: L
     integer  :: nbin
     integer  :: avg
     integer  :: err
     integer  :: idx

     integer  :: tmp
     integer  :: cnt

     logical  :: compute = .false.
 
     ! input file specifies if computing different correlation functions
     integer  :: flags(NTDMARRAY) 

     real(wp) :: dtau
     real(wp), pointer :: sgn(:)
     type(tdmarray), pointer :: properties(:)

     ! record the average of all local G(tau) for total N(w)
     real(wp), pointer :: GtupAve(:,:), GtdnAve(:,:)

     ! record Chi_FM and Chi_AF
     ! Here directly compute them within the code, instead of Fourier transform
     ! since sometimes for the general geometry, the computed Fourier components
     ! are somewhat incomplete
     real(wp), pointer :: ChiXX(:), ChiZZ(:)   ! : stores bins only

     ! Fourier transform matrix for bosonic and fermionic fields
     complex(wp), pointer :: ftwfer(:,:), ftwbos(:,:)

     ! 08/15/2015
     ! used for conductivity, d-wave paring sus etc.
     integer, ALLOCATABLE     :: rt(:), lf(:), up(:), dn(:)
     complex*16, ALLOCATABLE  :: hopup(:,:), hopdn(:,:)

  end type TDM1
  
contains

 !--------------------------------------------------------------------!
  
  subroutine DQMC_TDM1_Init(L, dtau, T1, nBin, S, Gwrap, flags)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    !    This subroutine initializes TDM1. 
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1      ! time dependent measurement
    integer, intent(in)       :: L       ! No of time slice
    integer, intent(in)       :: nBin    ! No of Bins
    integer, intent(in)       :: flags(NTDMARRAY)
    real(wp), intent(in)      :: dtau
    type(Struct), intent(in)  :: S
    type(GeomWrap), intent(in):: Gwrap

    ! ... local variables ...
    integer     :: i,j

    ! ... Executable ...

    T1%L      =  L
    T1%dtau   =  dtau
    T1%nbin   =  nBin

    T1%tmp    =  nBin + 1
    T1%avg    =  nBin + 1
    T1%err    =  nBin + 2
    T1%idx    =  1

    T1%compute  = .true.
    T1%flags    = flags

    ! Allocate storages
    allocate(T1%sgn(nBin+2))
    ! initialize values
    T1%sgn   = ZERO

    call  DQMC_TDM1_InitFTw(T1)
    ntdm = sum(T1%flags)         ! how many quantities to compute
!    write(*,*) "# of tdm quantities:", ntdm
    allocate(T1%properties(NTDMARRAY))

    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         call DQMC_TDM1_InitProp(T1, S, Gwrap, i)
       endif
    enddo

    if (T1%flags(IGFUP) == 1) then
      allocate(T1%GtupAve(0:T1%L-1,T1%err))
      T1%GtupAve = 0.0_wp
    endif
    if (T1%flags(IGFDN) == 1) then
      allocate(T1%GtdnAve(0:T1%L-1,T1%err))
      T1%GtdnAve = 0.0_wp
    endif

    allocate(T1%ChiXX(T1%err))
    allocate(T1%ChiZZ(T1%err))

    T1%ChiXX = 0.0_wp
    T1%ChiZZ = 0.0_wp
   
    ! IMPORTANT: the indices of Gtau(1:nsites) and hamilt(0:nsites-1) are different
    ! Note difference from rt etc. in dqmc_hamilt.F90
    allocate(T1%rt(1:S%nSite))
    allocate(T1%lf(1:S%nSite))
    allocate(T1%up(1:S%nSite))
    allocate(T1%dn(1:S%nSite))
    allocate(T1%hopup(1:S%nSite,1:S%nSite))
    allocate(T1%hopdn(1:S%nSite,1:S%nSite))
    ! switch the values to from 1 to nSite, instead of 0 to nSite-1
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

  end subroutine DQMC_TDM1_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_InitFTw(T1)
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
    type(TDM1), intent(inout) :: T1      ! TDM to be freed

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

  end subroutine DQMC_TDM1_InitFTw

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_InitProp(T1, S, Gwrap, iprop)
    use DQMC_Geom_Wrap
    !
    ! Purpose
    ! =======
    ! Initialize contents of tdmarray for iprop
    !
    ! Arguments
    ! =========
    type(TDM1), intent(inout) :: T1
    type(Struct), intent(in)  :: S
    type(GeomWrap), intent(in):: Gwrap
    integer, intent(in)       :: iprop

    integer :: nk, np, npp, nclass

    select case (iprop)

       case(IGFUN, IGFUP, IGFDN)

          nk     = Gwrap%RecipLattice%nclass_k
          nclass = S%nClass
          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  nclass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%nk     =  nk
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftk    => Gwrap%RecipLattice%FourierC
          T1%properties(iprop)%ftw    => T1%ftwfer
          T1%properties(iprop)%phase  => S%gf_phase
          T1%properties(iprop)%clabel  => S%clabel
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valuesk(nk*npp,0:T1%L-1,T1%err))

       case(ISPXX, ISPZZ, IDENS, IPAIR)

          nk = Gwrap%GammaLattice%nclass_k
          nclass = S%nClass
          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  S%nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%F      => S%F
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftk    => Gwrap%GammaLattice%FourierC
          T1%properties(iprop)%ftw    => T1%ftwbos
          T1%properties(iprop)%phase  => S%chi_phase
          T1%properties(iprop)%clabel  => S%clabel
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valuesk(nk*npp,0:T1%L-1,T1%err))

       case(ICOND)

          nk = Gwrap%GammaLattice%nclass_k
          nclass = 1
          np     = Gwrap%lattice%natom
          npp    = (np*(np+1))/2

          nullify(T1%properties(iprop)%tlink)
          T1%properties(iprop)%n      =  S%nSite
          T1%properties(iprop)%nclass =  nClass
          T1%properties(iprop)%D      => S%D
          T1%properties(iprop)%nk     =  Gwrap%GammaLattice%nclass_k
          T1%properties(iprop)%np     =  np
          T1%properties(iprop)%ftk    => Gwrap%GammaLattice%FourierC
          T1%properties(iprop)%ftw    => T1%ftwbos
          T1%properties(iprop)%phase  => S%chi_phase
          allocate(T1%properties(iprop)%values(nclass,0:T1%L-1,T1%err))
          allocate(T1%properties(iprop)%valuesk(nk*npp,0:T1%L-1,T1%err))

          allocate(T1%properties(iprop)%clabel(nclass))
          T1%properties(iprop)%clabel(1)=" "
        ! used in meas for average, for conductivity, renormalized by lattice size
          allocate(T1%properties(iprop)%F(nclass))
          T1%properties(iprop)%F(1) = S%nSite*2   ! *2 is for averaging x and y directions
     end select

     T1%properties(iprop)%values  = 0.0_wp
     if(associated(T1%properties(iprop)%valuesk)) &
         T1%properties(iprop)%valuesk = 0.0_wp

  end subroutine DQMC_TDM1_InitProp

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Free(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine frees TDM1.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1

    integer :: i

    ! ... Executable ...

    if (.not.T1%compute) return
    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         deallocate(T1%properties(i)%values)
         nullify(T1%properties(i)%D)
         nullify(T1%properties(i)%F)
         nullify(T1%properties(i)%ftk)
         nullify(T1%properties(i)%ftw)
         deallocate(T1%properties(i)%valuesk)
       endif
    enddo

    deallocate(T1%ftwbos)
    deallocate(T1%ftwfer)
    deallocate(T1%properties)

    if (T1%flags(IGFUP) == 1) then
      deallocate(T1%GtupAve)
    endif
    if (T1%flags(IGFDN) == 1) then
      deallocate(T1%GtdnAve)
    endif

    deallocate(T1%rt)
    deallocate(T1%lf)
    deallocate(T1%up)
    deallocate(T1%dn)
    deallocate(T1%hopup)
    deallocate(T1%hopdn)

  end subroutine DQMC_TDM1_Free

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Meas(T1, tau)
    !
    ! Purpose
    ! =======
    !    This subroutine fills the bin. It assumes that tau%A_up and,
    !    when necessary, tau%A_dn have been filled.
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)   :: t1
    type(Gtau), intent(inout)   :: tau
    
    ! ... Local var ...
    integer  :: i, k, m, L, cnt, dt, i0, it, j0, jt, dtau, iprop
    real(wp) :: sgn, factor
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

          call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
          ! Decrement index tau%it. If north is even do only north/2-1 decrements.
          do dt = 1, m-1+k
             call DQMC_change_gtau_time(tau, TPLUS, TAU_UP)
             if (tau%comp_dn) then
                call DQMC_change_gtau_time(tau, TPLUS, TAU_DN)
             elseif (.not.tau%neg_u) then
                call DQMC_Gtau_CopyUp(tau)
             endif

             jt = tau%it_up; j0 = tau%i0_up
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
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
             call DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, jt, j0)
                
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
          do it = 0, L-1
            do i = 1, T1%properties(iprop)%nClass
              factor = sgn/(T1%properties(iprop)%F(i)*cnt)
              values(i,it,T1%idx)   = values(i,it,T1%idx)   + factor*values(i,it,T1%tmp)
            end do
          end do
          values(:,:,T1%tmp)   = ZERO
       endif
    enddo

    T1%sgn(T1%idx) =  T1%sgn(T1%idx) + sgn
    T1%cnt = T1%cnt + 1

  end subroutine DQMC_TDM1_Meas

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Compute(T1, upt0, up0t, dnt0, dn0t, up00, uptt, dn00, dntt, it, i0)
    !
    ! Purpose
    ! =======
    !    This subroutine assembles the time dependent properties
    !    starting from the 1-body Green's function
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout)    :: T1
    real(wp), intent(in)         :: up0t(:,:), upt0(:,:)
    real(wp), intent(in)         :: dnt0(:,:), dn0t(:,:)
    real(wp), intent(in)         :: up00(:,:), uptt(:,:)
    real(wp), intent(in)         :: dn00(:,:), dntt(:,:)
    integer, intent(in)          :: it, i0
 
    ! ... Local scalar ...

    integer  :: i, j, k, dt, dt1, dt2
    real*8   :: a,b,c,d  
    real(wp), pointer :: value1(:), value2(:)
    real(wp) :: factor


    ! ... Executable ...
    if (.not.T1%compute) return

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
       ! value2 is for dt>0: G(beta-t,0) = -G(-t,0) = -G(0,t)
       !               dt<0: G(-t,0) = -G(beta-t,0)
       ! Some rules for value1 --> value2:
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
             ! SxSx = (Si+ * Sj+Si * Sj+)
             value1(k)  = value1(k) - (up0t(j,i)*dnt0(i,j) &
                  + up0t(i,j)*dnt0(j,i))/2
             value2(k)  = value2(k) - (up0t(i,j)*dnt0(j,i) &
                  + up0t(j,i)*dnt0(i,j))/2
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
             value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
               + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )*0.5_wp
             value2(k)  = value2(k) - (up0t(i,j)*upt0(j,i) &
               + dn0t(i,j)*dnt0(j,i) - (uptt(j,j)-dntt(j,j))*(up00(i,i)-dn00(i,i)) )*0.5_wp
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

     if (T1%flags(IPAIR) == 1) then
       value1  => T1%properties(IPAIR)%values(:, dt1, T1%tmp)
       value2  => T1%properties(IPAIR)%values(:, dt2, T1%tmp)
       do i = 1,  T1%properties(IPAIR)%n
          do j = 1,  T1%properties(IPAIR)%n
             ! someone decided there were two equivalent terms and 
             ! is averaging them for smaller error bars-  Eg if you wanted
             ! <A(t)A(0)> and you knew there was time reversal symmetry you
             ! might use 0.5* { <A(t)A(0) + A(0)A(t) }
             k = T1%properties(IPAIR)%D(i,j)
             value1(k)  = value1(k) + upt0(i,j)*dnt0(i,j)*0.5_wp 
             value2(k)  = value2(k) + upt0(j,i)*dnt0(j,i)*0.5_wp
          end do
       end do
     endif

     if (T1%flags(ICOND) == 1) then
       value1  => T1%properties(ICOND)%values(:, dt1, T1%tmp)
       value2  => T1%properties(ICOND)%values(:, dt2, T1%tmp)

       ! J-J correlation is not following Simone's trick of link correlation here
       ! here use standard definition <sum_ij jx(i,tau)*jx(j,0)>
       ! Only for q=0 component
       ! IMPORTANT: the indices of Gtau(1:nsites) and T1ilt(0:nsites-1) are different 
       ! e.g. T1%hopup(i-1,T1%rt(i-1))

       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n
            a = T1%hopup(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            b = T1%hopdn(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            c = T1%hopup(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            d = T1%hopdn(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            ! up*up terms (note two ways for contraction!)           
            value1(1) = value1(1) - a* &
                        (-upt0(i,T1%rt(j))*up0t(j,T1%rt(i)) - upt0(T1%rt(i),j)*up0t(T1%rt(j),i) &
                         +upt0(i,j)*up0t(T1%rt(j),T1%rt(i)) + upt0(T1%rt(i),T1%rt(j))*up0t(j,i))*0.5_wp
            value1(1) = value1(1) - a* &
                        ( uptt(i,T1%rt(i))*up00(j,T1%rt(j)) - uptt(T1%rt(i),i)*up00(T1%rt(j),j) &
                         -uptt(i,T1%rt(i))*up00(T1%rt(j),j) + uptt(T1%rt(i),i)*up00(j,T1%rt(j)))*0.5_wp
            ! dn*dn terms (note two ways for contraction!)
            value1(1) = value1(1) - b* &
                        (-dnt0(i,T1%rt(j))*dn0t(j,T1%rt(i)) - dnt0(T1%rt(i),j)*dn0t(T1%rt(j),i) &
                         +dnt0(i,j)*dn0t(T1%rt(j),T1%rt(i)) + dnt0(T1%rt(i),T1%rt(j))*dn0t(j,i))*0.5_wp
            value1(1) = value1(1) - b* &
                        ( dntt(i,T1%rt(i))*dn00(j,T1%rt(j)) - dntt(T1%rt(i),i)*dn00(T1%rt(j),j) &
                         -dntt(i,T1%rt(i))*dn00(T1%rt(j),j) + dntt(T1%rt(i),i)*dn00(j,T1%rt(j)))*0.5_wp
            ! up*dn terms
            value1(1) = value1(1) - c* &
                        ( uptt(i,T1%rt(i))*dn00(j,T1%rt(j)) + uptt(T1%rt(i),i)*dn00(T1%rt(j),j) &
                         -uptt(i,T1%rt(i))*dn00(T1%rt(j),j) - uptt(T1%rt(i),i)*dn00(j,T1%rt(j)))*0.5_wp
            ! dn*up terms
            value1(1) = value1(1) - d* &
                        ( dntt(i,T1%rt(i))*up00(j,T1%rt(j)) + dntt(T1%rt(i),i)*up00(T1%rt(j),j) &
                         -dntt(i,T1%rt(i))*up00(T1%rt(j),j) - dntt(T1%rt(i),i)*up00(j,T1%rt(j)))*0.5_wp

            a = T1%hopup(i,T1%up(i))*T1%hopup(j,T1%up(j))
            b = T1%hopdn(i,T1%up(i))*T1%hopdn(j,T1%up(j))
            c = T1%hopup(i,T1%up(i))*T1%hopdn(j,T1%up(j))
            d = T1%hopdn(i,T1%up(i))*T1%hopup(j,T1%up(j))
            ! up*up terms (note two ways for contraction!)
            value1(1) = value1(1) - a* &
                        (-upt0(i,T1%up(j))*up0t(j,T1%up(i)) - upt0(T1%up(i),j)*up0t(T1%up(j),i) &
                         +upt0(i,j)*up0t(T1%up(j),T1%up(i)) + upt0(T1%up(i),T1%up(j))*up0t(j,i))*0.5_wp
            value1(1) = value1(1) - a* &
                        ( uptt(i,T1%up(i))*up00(j,T1%up(j)) - uptt(T1%up(i),i)*up00(T1%up(j),j) &
                         -uptt(i,T1%up(i))*up00(T1%up(j),j) + uptt(T1%up(i),i)*up00(j,T1%up(j)))*0.5_wp
            ! dn*dn terms (note two ways for contraction!)
            value1(1) = value1(1) - b* &
                        (-dnt0(i,T1%up(j))*dn0t(j,T1%up(i)) - dnt0(T1%up(i),j)*dn0t(T1%up(j),i) &
                         +dnt0(i,j)*dn0t(T1%up(j),T1%up(i)) + dnt0(T1%up(i),T1%up(j))*dn0t(j,i))*0.5_wp
            value1(1) = value1(1) - b* &
                        ( dntt(i,T1%up(i))*dn00(j,T1%up(j)) - dntt(T1%up(i),i)*dn00(T1%up(j),j) &
                         -dntt(i,T1%up(i))*dn00(T1%up(j),j) + dntt(T1%up(i),i)*dn00(j,T1%up(j)))*0.5_wp
            ! up*dn terms
            value1(1) = value1(1) - c* &
                        ( uptt(i,T1%up(i))*dn00(j,T1%up(j)) + uptt(T1%up(i),i)*dn00(T1%up(j),j) &
                         -uptt(i,T1%up(i))*dn00(T1%up(j),j) - uptt(T1%up(i),i)*dn00(j,T1%up(j)))*0.5_wp
            ! dn*up terms
            value1(1) = value1(1) - d* &
                        ( dntt(i,T1%up(i))*up00(j,T1%up(j)) + dntt(T1%up(i),i)*up00(T1%up(j),j) &
                         -dntt(i,T1%up(i))*up00(T1%up(j),j) - dntt(T1%up(i),i)*up00(j,T1%up(j)))*0.5_wp
  
            ! value2: see the rules at the beginning of routine
            a = T1%hopup(j,T1%rt(j))*T1%hopup(i,T1%rt(i))
            b = T1%hopdn(j,T1%rt(j))*T1%hopdn(i,T1%rt(i))
            c = T1%hopup(j,T1%rt(j))*T1%hopdn(i,T1%rt(i))
            d = T1%hopdn(j,T1%rt(j))*T1%hopup(i,T1%rt(i))
            ! up*up terms
            value2(1) = value2(1) - a* &
                        (-upt0(j,T1%rt(i))*up0t(i,T1%rt(j)) - upt0(T1%rt(j),i)*up0t(T1%rt(i),j) &
                         +upt0(j,i)*up0t(T1%rt(i),T1%rt(j)) + upt0(T1%rt(j),T1%rt(i))*up0t(i,j))*0.5_wp
            value2(1) = value2(1) - a* &
                        ( uptt(j,T1%rt(j))*up00(i,T1%rt(i)) - uptt(T1%rt(j),j)*up00(T1%rt(i),i) &
                         -uptt(j,T1%rt(j))*up00(T1%rt(i),i) + uptt(T1%rt(j),j)*up00(i,T1%rt(i)))*0.5_wp
            ! dn*dn terms (note two ways for contraction!)
            value2(1) = value2(1) - b* &
                        (-dnt0(j,T1%rt(i))*dn0t(i,T1%rt(j)) - dnt0(T1%rt(j),i)*dn0t(T1%rt(i),j) &
                         +dnt0(j,i)*dn0t(T1%rt(i),T1%rt(j)) + dnt0(T1%rt(j),T1%rt(i))*dn0t(i,j))*0.5_wp
            value2(1) = value2(1) - b* &
                        ( dntt(j,T1%rt(j))*dn00(i,T1%rt(i)) - dntt(T1%rt(j),j)*dn00(T1%rt(i),i) &
                         -dntt(j,T1%rt(j))*dn00(T1%rt(i),i) + dntt(T1%rt(j),j)*dn00(i,T1%rt(i)))*0.5_wp
            ! up*dn terms
            value2(1) = value2(1) - c* &
                        ( uptt(j,T1%rt(j))*dn00(i,T1%rt(i)) + uptt(T1%rt(j),j)*dn00(T1%rt(i),i) &
                         -uptt(j,T1%rt(j))*dn00(T1%rt(i),i) - uptt(T1%rt(j),j)*dn00(i,T1%rt(i)))*0.5_wp
            ! dn*up terms
            value2(1) = value2(1) - d* &
                        ( dntt(j,T1%rt(j))*up00(i,T1%rt(i)) + dntt(T1%rt(j),j)*up00(T1%rt(i),i) &
                         -dntt(j,T1%rt(j))*up00(T1%rt(i),i) - dntt(T1%rt(j),j)*up00(i,T1%rt(i)))*0.5_wp

            a = T1%hopup(j,T1%up(j))*T1%hopup(i,T1%up(i))
            b = T1%hopdn(j,T1%up(j))*T1%hopdn(i,T1%up(i))
            c = T1%hopup(j,T1%up(j))*T1%hopdn(i,T1%up(i))
            d = T1%hopdn(j,T1%up(j))*T1%hopup(i,T1%up(i))
            ! up*up terms
            value2(1) = value2(1) - a* &
                        (-upt0(j,T1%up(i))*up0t(i,T1%up(j)) - upt0(T1%up(j),i)*up0t(T1%up(i),j) &
                         +upt0(j,i)*up0t(T1%up(i),T1%up(j)) + upt0(T1%up(j),T1%up(i))*up0t(i,j))*0.5_wp
            value2(1) = value2(1) - a* &
                        ( uptt(j,T1%up(j))*up00(i,T1%up(i)) - uptt(T1%up(j),j)*up00(T1%up(i),i) &
                         -uptt(j,T1%up(j))*up00(T1%up(i),i) + uptt(T1%up(j),j)*up00(i,T1%up(i)))*0.5_wp
            ! dn*dn terms (note two ways for contraction!)
            value2(1) = value2(1) - b* &
                        (-dnt0(j,T1%up(i))*dn0t(i,T1%up(j)) - dnt0(T1%up(j),i)*dn0t(T1%up(i),j) &
                         +dnt0(j,i)*dn0t(T1%up(i),T1%up(j)) + dnt0(T1%up(j),T1%up(i))*dn0t(i,j))*0.5_wp
            value2(1) = value2(1) - b* &
                        ( dntt(j,T1%up(j))*dn00(i,T1%up(i)) - dntt(T1%up(j),j)*dn00(T1%up(i),i) &
                         -dntt(j,T1%up(j))*dn00(T1%up(i),i) + dntt(T1%up(j),j)*dn00(i,T1%up(i)))*0.5_wp
            ! up*dn terms
            value2(1) = value2(1) - c* &
                        ( uptt(j,T1%up(j))*dn00(i,T1%up(i)) + uptt(T1%up(j),j)*dn00(T1%up(i),i) &
                         -uptt(j,T1%up(j))*dn00(T1%up(i),i) - uptt(T1%up(j),j)*dn00(i,T1%up(i)))*0.5_wp
            ! dn*up terms
            value2(1) = value2(1) - d* &
                        ( dntt(j,T1%up(j))*up00(i,T1%up(i)) + dntt(T1%up(j),j)*up00(T1%up(i),i) &
                         -dntt(j,T1%up(j))*up00(T1%up(i),i) - dntt(T1%up(j),j)*up00(i,T1%up(i)))*0.5_wp
          end do
       end do
     endif

    else

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
             value1(k)  = value1(k) - (up0t(j,i)*dnt0(i,j) &
                  + up0t(i,j)*dnt0(j,i))
          end do
       end do
     endif

     if (T1%flags(ISPZZ) == 1) then
       value1  => T1%properties(ISPZZ)%values(:, dt1, T1%tmp)
       do i = 1, T1%properties(ISPZZ)%n
          do j = 1, T1%properties(ISPZZ)%n
             ! k is the distance index of site i and site j
             k = T1%properties(ISPZZ)%D(i,j)
             value1(k)  = value1(k) - (up0t(j,i)*upt0(i,j) &
               + dn0t(j,i)*dnt0(i,j) - (uptt(i,i)-dntt(i,i))*(up00(j,j)-dn00(j,j)) )
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

     if (T1%flags(IPAIR) == 1) then
       value1  => T1%properties(IPAIR)%values(:, dt1, T1%tmp)
       do i = 1,  T1%properties(IPAIR)%n
          do j = 1,  T1%properties(IPAIR)%n
             ! k is the distance index of site i and site j
             k = T1%properties(IPAIR)%D(i,j)
             value1(k)  = value1(k) + upt0(i,j)*dnt0(i,j) 
          end do
       end do
     endif

     if (T1%flags(ICOND) == 1) then
       value1  => T1%properties(ICOND)%values(:, dt1, T1%tmp)
       do i = 1,  T1%properties(ICOND)%n
          do j = 1,  T1%properties(ICOND)%n
            a = T1%hopup(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            b = T1%hopdn(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            c = T1%hopup(i,T1%rt(i))*T1%hopdn(j,T1%rt(j))
            d = T1%hopdn(i,T1%rt(i))*T1%hopup(j,T1%rt(j))
            ! up*up terms (note two ways for contraction!)
            value1(1) = value1(1) - a* &
                        (-upt0(i,T1%rt(j))*up0t(j,T1%rt(i)) - upt0(T1%rt(i),j)*up0t(T1%rt(j),i) &
                         +upt0(i,j)*up0t(T1%rt(j),T1%rt(i)) + upt0(T1%rt(i),T1%rt(j))*up0t(j,i))
            value1(1) = value1(1) - a* &
                        ( uptt(i,T1%rt(i))*up00(j,T1%rt(j)) - uptt(T1%rt(i),i)*up00(T1%rt(j),j) &
                         -uptt(i,T1%rt(i))*up00(T1%rt(j),j) + uptt(T1%rt(i),i)*up00(j,T1%rt(j)))
            ! dn*dn terms (note two ways for contraction!)
            value1(1) = value1(1) - b* &
                        (-dnt0(i,T1%rt(j))*dn0t(j,T1%rt(i)) - dnt0(T1%rt(i),j)*dn0t(T1%rt(j),i) &
                         +dnt0(i,j)*dn0t(T1%rt(j),T1%rt(i)) + dnt0(T1%rt(i),T1%rt(j))*dn0t(j,i))
            value1(1) = value1(1) - b* &
                        ( dntt(i,T1%rt(i))*dn00(j,T1%rt(j)) - dntt(T1%rt(i),i)*dn00(T1%rt(j),j) &
                         -dntt(i,T1%rt(i))*dn00(T1%rt(j),j) + dntt(T1%rt(i),i)*dn00(j,T1%rt(j)))
            ! up*dn terms
            value1(1) = value1(1) - c* &
                        ( uptt(i,T1%rt(i))*dn00(j,T1%rt(j)) + uptt(T1%rt(i),i)*dn00(T1%rt(j),j) &
                         -uptt(i,T1%rt(i))*dn00(T1%rt(j),j) - uptt(T1%rt(i),i)*dn00(j,T1%rt(j)))
            ! dn*up terms
            value1(1) = value1(1) - d* &
                        ( dntt(i,T1%rt(i))*up00(j,T1%rt(j)) + dntt(T1%rt(i),i)*up00(T1%rt(j),j) &
                         -dntt(i,T1%rt(i))*up00(T1%rt(j),j) - dntt(T1%rt(i),i)*up00(j,T1%rt(j)))

            a = T1%hopup(i,T1%up(i))*T1%hopup(j,T1%up(j))
            b = T1%hopdn(i,T1%up(i))*T1%hopdn(j,T1%up(j))
            c = T1%hopup(i,T1%up(i))*T1%hopdn(j,T1%up(j))
            d = T1%hopdn(i,T1%up(i))*T1%hopup(j,T1%up(j))
            ! up*up terms (note two ways for contraction!)
            value1(1) = value1(1) - a* &
                        (-upt0(i,T1%up(j))*up0t(j,T1%up(i)) - upt0(T1%up(i),j)*up0t(T1%up(j),i) &
                         +upt0(i,j)*up0t(T1%up(j),T1%up(i)) + upt0(T1%up(i),T1%up(j))*up0t(j,i))
            value1(1) = value1(1) - a* &
                        ( uptt(i,T1%up(i))*up00(j,T1%up(j)) - uptt(T1%up(i),i)*up00(T1%up(j),j) &
                         -uptt(i,T1%up(i))*up00(T1%up(j),j) + uptt(T1%up(i),i)*up00(j,T1%up(j)))
            ! dn*dn terms (note two ways for contraction!)
            value1(1) = value1(1) - b* &
                        (-dnt0(i,T1%up(j))*dn0t(j,T1%up(i)) - dnt0(T1%up(i),j)*dn0t(T1%up(j),i) &
                         +dnt0(i,j)*dn0t(T1%up(j),T1%up(i)) + dnt0(T1%up(i),T1%up(j))*dn0t(j,i))
            value1(1) = value1(1) - b* &
                        ( dntt(i,T1%up(i))*dn00(j,T1%up(j)) - dntt(T1%up(i),i)*dn00(T1%up(j),j) &
                         -dntt(i,T1%up(i))*dn00(T1%up(j),j) + dntt(T1%up(i),i)*dn00(j,T1%up(j)))
            ! up*dn terms
            value1(1) = value1(1) - c* &
                        ( uptt(i,T1%up(i))*dn00(j,T1%up(j)) + uptt(T1%up(i),i)*dn00(T1%up(j),j) &
                         -uptt(i,T1%up(i))*dn00(T1%up(j),j) - uptt(T1%up(i),i)*dn00(j,T1%up(j)))
            ! dn*up terms
            value1(1) = value1(1) - d* &
                        ( dntt(i,T1%up(i))*up00(j,T1%up(j)) + dntt(T1%up(i),i)*up00(T1%up(j),j) &
                         -dntt(i,T1%up(i))*up00(T1%up(j),j) - dntt(T1%up(i),i)*up00(j,T1%up(j)))
          end do
       end do
     endif

    endif

  end subroutine DQMC_TDM1_Compute

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Avg(T1)
    !
    ! Purpose
    ! =======
    !    This subroutine average properties in a bin and
    !    increment the bin count (idx).
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer  :: nl, idx, i, j, k
    real(wp) :: factor

    ! ... Executable ...
    if (.not.T1%compute) return
    idx    = T1%idx
    factor = ONE/T1%cnt

    ! Compute average on Green's function
    do i = 1, NTDMARRAY
       if (T1%flags(i)==1) then
         nl = T1%properties(i)%nClass * T1%L
         call dscal(nl, factor, T1%properties(i)%values(:,0,idx), 1)
       endif
    enddo

    ! get the sum of all local G(tau) for total N(w)
    if (T1%flags(IGFUP) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IGFUP)%n
            k = T1%properties(IGFUP)%D(i,i)
               T1%GtupAve(j, T1%idx) = T1%GtupAve(j, T1%idx) &
                   + T1%properties(IGFUP)%values(k,j,idx) / T1%properties(IGFUP)%n
         end do
      enddo
    endif
    if (T1%flags(IGFDN) == 1) then
      do j = 0, T1%L-1
         do i = 1, T1%properties(IGFDN)%n
            k = T1%properties(IGFDN)%D(i,i)
               T1%GtdnAve(j, T1%idx) = T1%GtdnAve(j, T1%idx) &
                   + T1%properties(IGFDN)%values(k,j,idx) / T1%properties(IGFDN)%n
         end do
      end do
    endif

    T1%sgn(idx) = T1%sgn(idx)*factor
    T1%cnt = 0
    T1%idx = T1%idx + 1

  end subroutine DQMC_TDM1_Avg

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetErr(T1)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine compute the error in tdm using the jackknife
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(inout) :: T1                 ! T1

    ! ... local scalar ...
    integer   :: i, j, iprop
    integer   :: nproc, n, avg, err, mpi_err
    real(wp)  :: sum_sgn, sgn(T1%nBin), y(T1%nBin), data(T1%nBin)
    real(wp)  :: average, error

#   ifdef _QMC_MPI
      real(wp), pointer :: binptr(:,:), aveptr(:,:), errptr(:,:)
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

       if (T1%flags(IGFUP) == 1) then
         do j = 0, T1%L-1
           data =  T1%GtupAve(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%GtupAve(j, avg) = average
           T1%GtupAve(j, err) = error
         enddo
       endif
       if (T1%flags(IGFDN) == 1) then
         do j = 0, T1%L-1
           data =  T1%GtdnAve(j, 1:n)
           call DQMC_SignJackKnife(n, average, error, data, y, sgn, sum_sgn)
           T1%GtdnAve(j, avg) = average
           T1%GtdnAve(j, err) = error
         enddo
       endif
    else

       mpi_err = 0

#      ifdef _QMC_MPI
          
          ! comments below similar to dqmc_phy0.F90

          ! note for MPI, each processor has only ONE bin          
          ! The process below for computing err is JackKnife similar to non-MPI case
          ! See also DQMC_SignJackKnife_Real in dqmc_util.F90

          !    y_i = (sum(x)-x_i)/sgn_i
          !    The JackKnife variance of X with sign is defined as 
          !     
          !          n-1  
          !    sqrt(----- *sum(y_i-avg_y)^2))
          !           n
          ! 
          !    where avg_y = sum(y)/n
          
          call mpi_allreduce(T1%sgn(1), T1%sgn(avg), 1, mpi_double, &
             mpi_sum, mpi_comm_world, mpi_err)

          !Average properties
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                n = T1%properties(iprop)%nClass * T1%L
                call mpi_allreduce(binptr, aveptr, n, mpi_double, &
                   mpi_sum, mpi_comm_world, mpi_err)
            endif
          enddo

          ! Compute y_i
          T1%sgn(1)   = (T1%sgn(avg) - T1%sgn(1)) / dble(nproc - 1)
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                binptr = (aveptr - binptr) / dble(nproc - 1)
                binptr =  binptr / T1%sgn(1)
            endif
          enddo

          ! Compute avg = sum_x/sum_sgn
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                aveptr => T1%properties(iprop)%values(:,:,avg)
                aveptr =  aveptr / T1%sgn(avg) 
            endif
          enddo
          T1%sgn(avg)  = T1%sgn(avg) / dble(nproc)

          ! Compute error: sum(y_i-avg_y)^2
          do iprop = 1, NTDMARRAY
            if (T1%flags(iprop)==1) then
                binptr => T1%properties(iprop)%values(:,:,1)
                aveptr => T1%properties(iprop)%values(:,:,avg)
                errptr => T1%properties(iprop)%values(:,:,err)
                n = T1%properties(iprop)%nClass * T1%L
                call mpi_allreduce((binptr-aveptr)**2, errptr, n, mpi_double, &
                    mpi_sum, mpi_comm_world, mpi_err)
                errptr = sqrt(errptr * dble(nproc-1)/dble(nproc))
            endif
          enddo

#      endif

    endif

  end subroutine DQMC_TDM1_GetErr

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Print(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT

    integer             :: i, j, iprop
    real(wp)            :: tmp(T1%L, 2)
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title

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

  end subroutine DQMC_TDM1_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Print1(T1, OPT1)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT1!, OPT2

    integer             :: i, j
    real(wp)            :: tmp(T1%L, 2)
    character(len=10)   :: label(T1%L)
    character(len=slen) :: title

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    do j = 1, T1%L
       write(label(j),'(f10.5)') (j-1)*T1%dtau
    enddo

    ! Print average of local G(tau) for total N(w)
!    do j = 0, T1%L-1
!       tmp(j+1, 1:2) = T1%GtupAve(j, T1%avg:T1%err)
!    enddo
!    call DQMC_Print_Array(0, T1%L, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)

!    do j = 0, T1%L-1
!       tmp(j+1, 1:2) = T1%GtdnAve(j, T1%avg:T1%err)
!    enddo
!    call DQMC_Print_Array(0, T1%L, label, tmp(:, 1:1), tmp(:, 2:2), OPT2)

    ! Print local G(tau)'s
    if (T1%flags(IGFUN) == 1) then
      do i = 1, T1%properties(IGFUN)%nclass
        do j = 0, T1%L-1
          tmp(j+1, 1:2) = T1%properties(IGFUN)%values(i, j, T1%avg:T1%err)
        enddo
        title=pname(IGFUN)//" "//trim(adjustl(T1%properties(IGFUN)%clabel(i)))
        if (index(title, " 0.000   0.000   0.000") > 0) then
          call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT1)
        endif
     ! write(OPT1,'(1x)')
      enddo
    endif

!    do i = 1, T1%properties(IGFDN)%nclass
!      do j = 0, T1%L-1
!        tmp(j+1, 1:2) = T1%properties(IGFDN)%values(i, j, T1%avg:T1%err)
!      enddo
!      title=pname(IGFDN)//" "//trim(adjustl(T1%properties(IGFDN)%clabel(i)))
!      if (index(title, " 0.0000  0.0000") > 0) then
!        call DQMC_Print_Array(0, T1%L, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT2)
!      endif
     ! write(OPT2,'(1x)')
!    enddo

  end subroutine DQMC_TDM1_Print1

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_Chi_Print(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints temporal sum of correlation for computing Chi
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
    integer, intent(in)      :: OPT
    integer      :: nn

    integer             :: i, j, iprop
    real(wp)            :: tmp(T1%properties(ISPXX)%nclass, 2)
    character(len=30)   :: label(T1%properties(ISPXX)%nclass)
    character(len=slen) :: title

    nn = T1%properties(ISPXX)%nclass

    ! ... Executable ...
    if (.not.T1%compute) return

    if (qmc_sim%rank .ne. 0) return

    if (T1%flags(ISPXX) == 1) then
      do j = 1, nn
         write(label(j),*) trim(adjustl(T1%properties(ISPXX)%clabel(j)))
      enddo
      do i = 1, nn
         ! sum over all tau's components 
         tmp(i, 1) = sum(T1%properties(ISPXX)%values(i, 0:T1%L-1, T1%avg))
         ! Note that tmp(i,2) is in fact wrong (it is sum of error)
         tmp(i, 2) = sum(T1%properties(ISPXX)%values(i, 0:T1%L-1, T1%err))
      enddo
      title = pname(ISPXX)
      call DQMC_Print_RealArray(0, nn, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
      write(OPT,'(1x)')
    endif

    if (T1%flags(ISPZZ) == 1) then
      do j = 1, nn
         write(label(j),*) trim(adjustl(T1%properties(ISPZZ)%clabel(j)))
      enddo
      do i = 1, nn
         ! sum over all tau's components 
         tmp(i, 1) = sum(T1%properties(ISPZZ)%values(i, 0:T1%L-1, T1%avg))
         ! Note that tmp(i,2) is in fact wrong (it is sum of error)
         tmp(i, 2) = sum(T1%properties(ISPZZ)%values(i, 0:T1%L-1, T1%err))
      enddo
      title = pname(ISPZZ)
      call DQMC_Print_RealArray(0, nn, title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
      write(OPT,'(1x)')
    endif

  end subroutine DQMC_TDM1_Chi_Print

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetKFT(T1)

    type(tdm1), intent(inout) :: T1

    integer :: ip, it, n, nclass, np, nk, ibin
    integer, pointer :: class(:,:)
    complex(wp), pointer :: wgtftk(:,:)
    integer, pointer :: phase(:,:)

    real(wp), pointer :: value(:)
    complex(wp), pointer :: valuet(:)

    if (.not.T1%compute) return
 
    !Loop over properties to Fourier transform
    do ip = 1, NTDMARRAY
       if (T1%flags(ip)==1) then
         if (.not.associated(T1%properties(ip)%valuesk)) cycle

         ! Aliases
         n        =  T1%properties(ip)%n
         nclass   =  T1%properties(ip)%nclass
         np       =  T1%properties(ip)%np
         nk       =  T1%properties(ip)%nk
         class    => T1%properties(ip)%D
         wgtftk   => T1%properties(ip)%ftk
         phase    => T1%properties(ip)%phase

         !Fourier transform each bin and average
         do ibin = T1%avg, 1, -1

            ! More aliases
            do it = 0, T1%L-1
               value  =>  T1%properties(ip)%values(:,it,ibin)
               valuet =>  T1%properties(ip)%valuesk(:,it,ibin)
               call dqmc_getFTk(value, n, nclass, class, np, nk, wgtftk, phase, valuet)
            enddo

         enddo ! Loop over bins
      endif
    enddo ! Loop over properties

  end subroutine DQMC_TDM1_GetKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_GetErrKFT(T1)

    use DQMC_MPI

    type(tdm1), intent(inout) :: T1

    integer :: ip, it, n, nproc, i

    complex(wp), pointer  :: average(:), binval(:), error(:), temp(:)
 
    !Loop over properties to Fourier transform
    nproc = qmc_sim%size

    if (.not.T1%compute) return

    if (nproc .eq. 1) then

      do ip = 1, NTDMARRAY
        if (T1%flags(ip)==1) then
          if (.not.associated(T1%properties(ip)%valuesk)) cycle

          do it = 0, T1%L-1

             average  => T1%properties(ip)%valuesk(:,it,T1%avg)
             error    => T1%properties(ip)%valuesk(:,it,T1%err)

             !Fourier transform each bin and average
             do i = 1, T1%nbin

                binval => T1%properties(ip)%valuesk(:,it,i)
                error  = error  +  cmplx((real(average-binval))**2,(aimag(average-binval))**2)

             enddo ! Loop over bins

             error  = (T1%nbin-1)*error/T1%nbin
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)))
  
          enddo
        endif
      enddo ! Loop over properties

    else

      do ip = 1, NTDMARRAY
        if (T1%flags(ip)==1) then
          if (.not.associated(T1%properties(ip)%valuesk)) cycle
    
          n = T1%properties(ip)%nk * T1%properties(ip)%np
          allocate(temp(n))
          
          do it = 0, T1%L-1
             average  => T1%properties(ip)%valuesk(:,it,T1%avg)
             error    => T1%properties(ip)%valuesk(:,it,T1%err)
             binval   => T1%properties(ip)%valuesk(:,it,1)
             temp     =  cmplx((real(average-binval))**2,(aimag(average-binval))**2)
#            ifdef _QMC_MPI
             call mpi_allreduce(temp, error, n, mpi_double, mpi_sum, mpi_comm_world, i)
#            endif
             error  = (nproc-1)*error/nproc
             error = cmplx(sqrt(real(error)),sqrt(aimag(error)))
          enddo

          deallocate(temp)
        endif
      enddo ! Loop over properties

    endif


  end subroutine DQMC_TDM1_GetErrKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_PrintKFT(T1, OPT)
    use dqmc_mpi
    !
    ! Purpose
    ! =======
    !    This subroutine prints properties to file
    !
    ! Arguments
    ! =========
    !
    type(TDM1), intent(in)   :: T1                 ! T1
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

    do iprop = 1, NTDMARRAY
       if (T1%flags(iprop)==1) then
         if (.not.associated(T1%properties(iprop)%valuesk)) cycle
         np = T1%properties(iprop)%np
         npp = (np*(np+1))/2
         do k = 1, T1%properties(iprop)%nk
            i = (k-1)*npp
            do ip = 1, np
               do jp = ip, np
                  i = i + 1
                  do j = 0, T1%L-1
                     tmp(j+1, 1:2) = T1%properties(iprop)%valuesk(i, j, T1%avg:T1%err)
                  enddo
                  write(title,'(A,i3,A,i3,A,i3,A)') 'k=',k,'   pair=',ip,',',jp
                  title=pname(iprop)//" "//trim(adjustl(title))
                  call DQMC_Print_Array(0, T1%L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                  write(OPT,'(1x)')
               enddo
            enddo
         enddo
      endif
    enddo

  end subroutine DQMC_TDM1_PrintKFT

  !--------------------------------------------------------------------!

  subroutine DQMC_TDM1_SelfEnergy(T1, tau, OPT)

    use DQMC_MPI

    type(TDM1), intent(in)    :: T1
    type(gtau), intent(inout) :: tau
    integer, intent(in)       :: OPT

    real(wp),    allocatable  :: g0tau(:,:), tdmg0(:,:)
    complex(wp), allocatable  :: tdmg0k(:,:), tdmgk(:,:)
    complex(wp), allocatable  :: tdmg0kw(:,:,:), tdmgkw(:,:,:)
    complex(wp), allocatable  :: avgSE(:,:,:), errSE(:,:,:), binSE(:,:,:)

    integer :: i, j, k, h, m
    integer :: L, n, nclass, np, nk, npp
    integer,     pointer :: class(:,:), ph(:,:)
    complex(wp), pointer :: ftk(:,:), ftw(:,:)

    complex(wp)          :: tmp(T1%L,2)
    character(len=10)    :: label(T1%L)
    character(len=50)    :: title

    integer, parameter  :: gflist(3) = (/IGFUN, IGFUP, IGFDN/)
    integer, parameter  :: splist(3) = (/  0, TAU_UP, TAU_DN/)
    real(wp), parameter :: pi = 3.1415926535898

    if (.not.T1%compute) return

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

    do j = 1, L
       write(label(j),'(f10.5)') (2*j-1)*pi/(T1%dtau*L)
       label(j) = adjustl(label(j))
    enddo

    !non-interacting green's function
    allocate(g0tau(n,n))
    allocate(tdmg0(nclass,0:L-1))
    allocate(tdmg0k(nk*npp,0:L-1))

    !temporary storage
    allocate(tdmgk(npp,0:L-1))

    !frequency dependent green's functions
    allocate(tdmg0kw(np,np,0:L-1))
    allocate(tdmgkw(np,np,0:L-1))

    !self energy
    allocate(avgSE(np,np,0:L-1))
    allocate(errSE(np,np,0:L-1))
    allocate(binSE(np,np,0:L-1))

    do h = 1, 3

       ! Get G for the non-interacting system
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
 
          i = (k-1) * npp + 1
          j = k * npp

          ! Transform G0 from tau to iwn
          tdmgk = tdmg0k(i:j,0:L-1)
          call convert_to_iwn(tdmgk, tdmg0kw)
          call invertG(tdmg0kw)

          ! Transform G from tau to iwn for average
          tdmgk = T1%properties(gflist(h))%valuesk(i:j,0:L-1,T1%avg)
          call convert_to_iwn(tdmgk, tdmgkw)
          call invertG(tdmgkw)

          !Compute average self-energy
          avgSE = tdmg0kw - tdmgkw

          errSE = ZERO
          do m = 0, T1%nbin-1

             ! Transform G from tau to iwn for bin "m"
             tdmgk = T1%properties(gflist(h))%valuesk(i:j,0:L-1,m+1)
             call convert_to_iwn(tdmgk, tdmgkw)
             call invertG(tdmgkw)

             ! Compute self-energy for bin
             binSE = tdmg0kw - tdmgkw
             if (qmc_sim%size .eq. 1) &
                errSE = errSE + cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
  
          enddo 

          if (qmc_sim%size .gt. 1) then

             m = qmc_sim%size
             ! Reuse tdmgkw for temporary storage
             tdmgkw =  cmplx((real(binSE-avgSE))**2,(aimag(binSE-avgSE))**2)
#            ifdef _QMC_MPI
             call mpi_allreduce(tdmgkw, errSE, n, mpi_double_complex, mpi_sum, mpi_comm_world, i)
#            endif

          endif

          errSE = cmplx(sqrt(real(errSE)),sqrt(aimag(errSE))) * sqrt(dble(m-1)/m)

          ! Take care of printing
          if (qmc_sim%rank .eq. 0) then
             do i = 1, np
                do j = 1, np
                   tmp(1:L,1) = avgSE(i,j,0:L-1)
                   tmp(1:L,2) = errSE(i,j,0:L-1)
                   write(title,'(A,i3)') trim(pname(gflist(h)))//" SelfEn k=", k
                   write(title,'(A,i3,A,i3)') trim(adjustl(title))//'   pair=',i,',',j
                   call DQMC_Print_Array(0, L , title, label, tmp(:, 1:1), tmp(:, 2:2), OPT)
                   write(OPT,'(1x)')
                enddo
             enddo
          endif

       enddo

    enddo

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


  end subroutine DQMC_TDM1_SelfEnergy

  !--------------------------------------------------------------------!

end module DQMC_TDM1
