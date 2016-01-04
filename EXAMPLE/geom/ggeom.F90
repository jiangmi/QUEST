program dqmc_ggeom

  ! 10/25/2015:
  ! add flags for specifying in input file if to compute tdm quantities 

  use dqmc_cfg
  use dqmc_geom_wrap
  use dqmc_hubbard
  use dqmc_mpi
  use dqmc_tdm1

  implicit none

  integer :: values(8)
  character(3) :: mons(12)

  real                :: t1, t2, t3, t4, t5
  type(config)        :: cfg
  type(Hubbard)       :: Hub
  type(GeomWrap)      :: Gwrap
  type(tdm1)          :: tm
  type(Gtau)          :: tau
  character(len=slen) :: gfile
  logical             :: tformat
  integer             :: na, nt, nkt, nkg, i, j, k, slice, nhist, comp_tdm
  integer             :: nBin, nIter, ierr  
  character(len=50)   :: ofile  
  integer             :: OPT,OPT1,OPT2  !,OPT3,OPT4
  !integer             :: HSF_output_file_unit
  integer             :: symmetries_output_file_unit
  integer             :: FLD_UNIT, TDM_UNIT
  real(wp)            :: randn(1)
  integer, pointer    :: flags(:)      ! 10 possible tdm1 quantities
  integer             :: nflag
  integer             :: ntry2, FTphy0, FTtdm, SelfE, Dsqy

  call cpu_time(t1)  

  !Count the number of processors
  call DQMC_MPI_Init(qmc_sim, PLEVEL_1)
 
  !Read input
  call DQMC_Read_Config(cfg)

  !Get output file name header
  call CFG_Get(cfg, "ofile", ofile)

  !Get general geometry input
  call CFG_Get(cfg, "gfile", gfile)

  !Save whether to use refinement for G used in measurements.
  call CFG_Get(cfg, "nhist", nhist)

  !Get some input variables
  call CFG_Get(cfg, "ntry2" , ntry2)
  call CFG_Get(cfg, "tdm"   , comp_tdm)
  call CFG_Get(cfg, "FTphy0", FTphy0)
  call CFG_Get(cfg, "FTtdm" , FTtdm)
  call CFG_Get(cfg, "SelfE" , SelfE)
  call CFG_Get(cfg, "Dsqy"  , Dsqy)

  !if (nhist > 0) then
  !   call DQMC_open_file(adjustl(trim(ofile))//'.HSF.stream','unknown', HSF_output_file_unit)
  !endif

  call DQMC_open_file(adjustl(trim(ofile))//'.geometry','unknown', symmetries_output_file_unit)
  !Determines type of geometry file
  call DQMC_Geom_Read_Def(Hub%S, gfile, tformat)   ! dqmc_struct.F90, not really used since tableformat=false
  if (.not.tformat) then
     ! In dqmc_geom_wrap.F90:
     !If free format fill gwrap
     call DQMC_Geom_Fill(Gwrap, gfile, cfg, symmetries_output_file_unit)
     !Transfer info in Hub%S
     call DQMC_Geom_Init(Gwrap,Hub%S,cfg)
  endif
  call DQMC_Geom_Print(Hub%S, symmetries_output_file_unit) ! dqmc_struct.F90

  ! Initialize the rest data (including seed for MPI processes), dqmc_hubbard.F90
  call DQMC_Hub_Config(Hub, cfg, Gwrap)

  ! Initialize time dependent properties if comp_tdm > 0
  if (comp_tdm > 0) then
     ! If to compute 8 possible tdm1 quantities
     call CFG_Get(cfg, "flags", nflag, flags)
     call DQMC_Gtau_Init(Hub, tau)
     call DQMC_TDM1_Init(Hub%L, Hub%dtau, tm, Hub%P0%nbin, Hub%S, Gwrap, flags) 
  endif

  call cpu_time(t2)

  ! print date/time info into output
  if (qmc_sim%rank == qmc_sim%aggr_root) then
    mons = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    call date_and_time(VALUES=values)

    write(*,*) "Initialization time:",  t2-t1, "(second)"
    write(*,*) "Initialization time:",  (t2-t1)/60, "(minutes)"
    write(*,*) "Initialization time:",  (t2-t1)/3600, "(hours)"
    write(*,*) "Initialization time:",  (t2-t1)/3600/24, "(days)"
  endif

  ! Warmup sweep
  do i = 1, Hub%nWarm
     if (qmc_sim%rank == qmc_sim%aggr_root) then
       if (mod(i, 10)==0) write(*,'(A,i6,1x,i3)')' Warmup Sweep, nwrap  : ', i, Hub%G_up%nwrap
     endif
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)   ! NO_MEAS0 = -1, parameter defined in dqmc_hubbard.F90
     call DQMC_Hub_Sweep2(Hub, ntry2)     ! sweep2 is for global move update
  end do

  call cpu_time(t3)
  if (qmc_sim%rank == qmc_sim%aggr_root) then
    write(*,*) "Warmup time:",  t3-t2, "(second)"
    write(*,*) "Warmup time:",  (t3-t2)/60, "(minutes)"
    write(*,*) "Warmup time:",  (t3-t2)/3600, "(hours)"
    write(*,*) "Warmup time:",  (t3-t2)/3600/24, "(days)"
  endif

  ! We divide all the measurement into nBin, each having nPass/nBin pass
  nBin   = Hub%P0%nBin
  nIter  = Hub%nPass / Hub%tausk / nBin
  if (nIter > 0) then
     do i = 1, nBin
        do j = 1, nIter
           ! only measure once every tausk warmup steps
           do k = 1, Hub%tausk
              call DQMC_Hub_Sweep(Hub, NO_MEAS0)
              call DQMC_Hub_Sweep2(Hub, Hub%nTry)
           enddo

           ! Fetch a random slice for measurement 
           call ran0(1, randn, Hub%seed)
           slice = ceiling(randn(1)*Hub%L)
           write(*,'(a,3i6)') ' Measurement Sweep, bin, iter, slice : ', i, j, slice

           if (comp_tdm > 0) then
              ! Compute full Green's function 
              call DQMC_Gtau_LoadA(tau, TAU_UP, slice, Hub%G_up%sgn)
              call DQMC_Gtau_LoadA(tau, TAU_DN, slice, Hub%G_dn%sgn)
              ! Measure equal-time properties
              call DQMC_Hub_FullMeas(Hub, tau%nnb, tau%A_up, tau%A_dn, tau%sgnup, tau%sgndn)
              ! Measure time-dependent properties
              call DQMC_TDM1_Meas(tm, tau)
           else if (comp_tdm == 0) then
              call DQMC_Hub_Meas(Hub, slice)
           endif

           if (qmc_sim%rank == qmc_sim%aggr_root) then
             write(*,*)"Meas #",(i-1)*nIter + j
           endif
           !Write fields 
           !if (nhist > 0) call DQMC_Hub_Output_HSF(Hub, .false., slice, HSF_output_file_unit)
        end do

        ! Accumulate results for each bin
        call DQMC_Phy0_Avg(Hub%P0)
        call DQMC_tdm1_Avg(tm)

        if (Hub%meas2) then
           if(Hub%P2%diagonalize)then
             call DQMC_Phy2_Avg(Hub%P2, Hub%S)
           else
             call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
           endif
        end if

     end do
  endif  ! for if nIter  = Hub%nPass/Hub%tausk/nBin > 0

  !Read HSF configurations from file if no sweep was perfomed
  if (Hub%nWarm + Hub%nPass == 0) then
     Hub%nMeas = -1
     call DQMC_count_records(Hub%npass, FLD_UNIT)
     nIter = Hub%npass / nbin
     do i = 1, nBin
        do j = 1, nIter / qmc_sim%aggr_size
           call DQMC_Hub_Input_HSF(Hub, .false., slice, FLD_UNIT)
           call DQMC_Hub_Init_Vmat(Hub)
           if (comp_tdm > 0) then
              ! Compute full Green's function - if fullg is on -
              call DQMC_Gtau_LoadA(tau, TAU_UP, slice, Hub%G_up%sgn)
              call DQMC_Gtau_LoadA(tau, TAU_DN, slice, Hub%G_dn%sgn)
              ! Measure equal-time properties. Pass gtau in case fullg was computed.
              call DQMC_Hub_FullMeas(Hub, tau%nb, &
                 tau%A_up, tau%A_dn, tau%sgnup, tau%sgndn)
              ! Measure time-dependent properties. Reuses fullg when possible. 
              call DQMC_TDM1_Meas(tm, tau)
           else if (comp_tdm == 0) then
              call DQMC_Hub_Meas(Hub, slice)
           endif
        enddo

        call DQMC_Phy0_Avg(Hub%P0)
        call DQMC_TDM1_Avg(tm)

        if (Hub%meas2) then
           if(Hub%P2%diagonalize)then
             call DQMC_Phy2_Avg(Hub%P2, Hub%S)
           else
             call DQMC_Phy2_Avg(Hub%P2, Hub%S%W)
           endif
        end if

     enddo
  endif ! for if (Hub%nWarm + Hub%nPass == 0)

  !Compute average and error
  call DQMC_Phy0_GetErr(Hub%P0)
  call DQMC_TDM1_GetErr(tm)
  if (Hub%meas2) then
     call DQMC_Phy2_GetErr(Hub%P2)
  end if

  call cpu_time(t4)
  if (qmc_sim%rank == qmc_sim%aggr_root) then
    write(*,*) "Meas time:",  t4-t3, "(second)"
    write(*,*) "Meas time:",  (t4-t3)/60, "(minutes)"
    write(*,*) "Meas time:",  (t4-t3)/3600, "(hours)"
    write(*,*) "Meas time:",  (t4-t3)/3600/24, "(days)"
  endif

  ! ============================================================================================
  ! Print results at root
  if (qmc_sim%rank == qmc_sim%aggr_root) then

    call DQMC_open_file(adjustl(trim(ofile))//'.out', 'unknown', OPT)
    if (comp_tdm > 0) then
      call DQMC_open_file(adjustl(trim(ofile))//'.tdm.out','unknown', TDM_UNIT)
      call DQMC_open_file('G_'//adjustl(trim(ofile)),'replace', OPT1)
      if (Dsqy > 0) then
        call DQMC_open_file('current_'//adjustl(trim(ofile)),'replace', OPT2)
      endif
    endif

    write(OPT,'(a14,a3,a1,i2,a2,i2,a1,i2,a1,i2,a2,i4)') &
               "STARTING JOB: ", mons(values(2)), " ", values(3), "  ", &
               values(5), ":", values(6), ":", values(7), "  ", values(1)
    write(OPT,*) "============================================================================"
    write(OPT,*) "Initialization time:",  t2-t1, "(second)"
    write(OPT,*) "Initialization time:",  (t2-t1)/60, "(minutes)"
    write(OPT,*) "Initialization time:",  (t2-t1)/3600, "(hours)"
    write(OPT,*) "Initialization time:",  (t2-t1)/3600/24, "(days)"
    write(OPT,*) "-------------------------------------------------"
    write(OPT,*) "Warmup time:",  t3-t2, "(second)"
    write(OPT,*) "Warmup time:",  (t3-t2)/60, "(minutes)"
    write(OPT,*) "Warmup time:",  (t3-t2)/3600, "(hours)"
    write(OPT,*) "Warmup time:",  (t3-t2)/3600/24, "(days)"
    write(OPT,*) "-------------------------------------------------"
    write(OPT,*) "Meas time:",  t4-t3, "(second)"
    write(OPT,*) "Meas time:",  (t4-t3)/60, "(minutes)"
    write(OPT,*) "Meas time:",  (t4-t3)/3600, "(hours)"
    write(OPT,*) "Meas time:",  (t4-t3)/3600/24, "(days)"
    write(OPT,*) "============================================================================"

  endif ! end if (qmc_sim%rank == qmc_sim%aggr_root)

  !Print P0 results: determine if qmc_sim%rank==0 in subroutines
  call DQMC_Hub_OutputParam(Hub, OPT)
  call DQMC_Phy0_Print(Hub%P0, Hub%S, OPT)

  !Print tdm and G(tau) local: determine if qmc_sim%rank==0 in subroutines
  if (comp_tdm > 0) then
    call DQMC_TDM1_Print(tm, TDM_UNIT)
    call DQMC_TDM1_Print_localGtau(tm, OPT1)
  endif

! ==== curr-curr(qx=0,qy;iwn=0) is estimated by linear extrapolation of two smallest qy ======
  !Direct access to binned data; no need to be in the loop above
  if (Dsqy > 0) then
    call DQMC_TDM1_currDs(tm,Hub)  ! use Hub%S and Hub%dtau
    call DQMC_TDM1_currDs_Err(tm)
    call DQMC_TDM1_currDs_Print(tm,OPT2)
  endif

! ==============  Fourier transform ============================================
  !Compute Fourier transform
  !Direct access to binned data; no need to be in the loop above
  if (FTphy0 > 0) then
    !Aliases for Fourier transform
    na  =  Gwrap%lattice%natom
    nt  =  Gwrap%lattice%ncell
    nkt =  Gwrap%RecipLattice%nclass_k
    nkg =  Gwrap%GammaLattice%nclass_k

    call DQMC_phy0_GetFT(Hub%P0, Hub%S%D, Hub%S%gf_phase, Gwrap%RecipLattice%FourierC, &
                         Gwrap%GammaLattice%FourierC, nkt, nkg, na, nt)
    call DQMC_Phy0_GetErrFT(Hub%P0)
  endif

  !In DQMC_TDM1_GetKFT determines if comp_tdm>0
  if (FTtdm > 0)  then
    call DQMC_TDM1_GetKFT(tm)
    call DQMC_TDM1_GetErrKFT(tm)
  endif
  if (SelfE > 0)  then
    call DQMC_TDM1_SelfEnergy(tm, tau)
  endif

  !Printing: determine if qmc_sim%rank==0 in subroutines
  if (FTphy0 > 0) then
     !Print info on k-points and construct clabel, in dqmc_geom_wrap.F90
     call DQMC_Print_HeaderFT(Gwrap, OPT, .true.)   ! k grid for Green function
     call DQMC_Print_HeaderFT(Gwrap, OPT, .false.)  ! k grid for spin/charge correlation
     call DQMC_Phy0_PrintFT(Hub%P0, na, nkt, nkg, OPT)
  endif
  if (FTtdm > 0) then
     call DQMC_TDM1_PrintKFT(tm, TDM_UNIT)
  endif
  if (SelfE > 0)  then
     call DQMC_TDM1_Print_SelfEnergy(tm, TDM_UNIT)
  endif

! =========== End Fourier transform ============================================

  ! Clean up the used storage
  call DQMC_TDM1_Free(tm)
  call DQMC_Hub_Free(Hub)
  call DQMC_Config_Free(cfg)
  
  call cpu_time(t5)

  if (qmc_sim%rank == qmc_sim%aggr_root) then
    write(*,*) "Total time:",  t5-t1, "(second)"
    write(*,*) "Total time:",  (t5-t1)/60, "(minutes)"
    write(*,*) "Total time:",  (t5-t1)/3600, "(hours)"
    write(*,*) "Total time:",  (t5-t1)/3600/24, "(days)"

    write(OPT,*) "Total time:",  t5-t1, "(second)"
    write(OPT,*) "Total time:",  (t5-t1)/60, "(minutes)"
    write(OPT,*) "Total time:",  (t5-t1)/3600, "(hours)"
    write(OPT,*) "Total time:",  (t5-t1)/3600/24, "(days)"
  endif

  call DQMC_MPI_Final(qmc_sim)

  close(symmetries_output_file_unit)

end program dqmc_ggeom

