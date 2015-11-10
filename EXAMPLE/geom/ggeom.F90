program dqmc_ggeom

  ! 10/25/2015:
  ! add flags for specifying in input file if to compute tdm quantities 

  use dqmc_cfg
  use dqmc_geom_wrap
  use dqmc_hubbard
  use dqmc_mpi
  use dqmc_tdm1

  implicit none

  real                :: t1, t2, t3, t4, t5
  type(config)        :: cfg
  type(Hubbard)       :: Hub
  type(GeomWrap)      :: Gwrap
  type(tdm1)          :: tm
  type(Gtau)          :: tau
  character(len=slen) :: gfile
  logical             :: tformat
  integer             :: na, nt, nkt, nkg, i, j, k, slice, nhist, comp_tdm
  integer             :: nBin, nIter  
  character(len=50)   :: ofile  
  integer             :: OPT,OPT1 !,OPT2,OPT3,OPT4
  !integer             :: HSF_output_file_unit
  integer             :: symmetries_output_file_unit
  integer             :: FLD_UNIT, TDM_UNIT
  real(wp)            :: randn(1)
  integer, pointer    :: flags(:)      ! 8 possible tdm1 quantities
  integer             :: nflag
  character(len=1)    :: ranklabel

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
  call DQMC_Hub_Config(Hub, cfg)

  ! Initialize time dependent properties if comp_tdm > 0
  call CFG_Get(cfg, "tdm", comp_tdm)
  if (comp_tdm > 0) then
     ! If to compute 8 possible tdm1 quantities
     call CFG_Get(cfg, "flags", nflag, flags)
     call DQMC_Gtau_Init(Hub, tau)
     call DQMC_TDM1_Init(Hub%L, Hub%dtau, tm, Hub%P0%nbin, Hub%S, Gwrap, flags) 
  endif

  call cpu_time(t2)

  if (qmc_sim%rank == qmc_sim%aggr_root) then
    call DQMC_open_file(adjustl(trim(ofile))//'.out', 'unknown', OPT)

    write(*,*) "Initialization time:",  t2-t1, "(second)"
    write(*,*) "Initialization time:",  (t2-t1)/60, "(minutes)"
    write(*,*) "Initialization time:",  (t2-t1)/3600, "(hours)"
    write(*,*) "Initialization time:",  (t2-t1)/3600/24, "(days)"
    write(OPT,*) "Initialization time:",  t2-t1, "(second)"
    write(OPT,*) "Initialization time:",  (t2-t1)/60, "(minutes)"
    write(OPT,*) "Initialization time:",  (t2-t1)/3600, "(hours)"
    write(OPT,*) "Initialization time:",  (t2-t1)/3600/24, "(days)"
  endif

  ! Warmup sweep
  do i = 1, Hub%nWarm
     if (qmc_sim%rank == qmc_sim%aggr_root) then
       if (mod(i, 10)==0) write(*,'(A,i6,1x,i3)')' Warmup Sweep, nwrap  : ', i, Hub%G_up%nwrap
     endif
     call DQMC_Hub_Sweep(Hub, NO_MEAS0)   ! NO_MEAS0 = -1, parameter defined in dqmc_hubbard.F90
     call DQMC_Hub_Sweep2(Hub, Hub%nTry)  ! sweep2 is for global move update
  end do

  call cpu_time(t3)
  if (qmc_sim%rank == qmc_sim%aggr_root) then
    write(*,*) "Warmup time:",  t3-t2, "(second)"
    write(*,*) "Warmup time:",  (t3-t2)/60, "(minutes)"
    write(*,*) "Warmup time:",  (t3-t2)/3600, "(hours)"
    write(*,*) "Warmup time:",  (t3-t2)/3600/24, "(days)"

    write(OPT,*) "Warmup time:",  t3-t2, "(second)"
    write(OPT,*) "Warmup time:",  (t3-t2)/60, "(minutes)"
    write(OPT,*) "Warmup time:",  (t3-t2)/3600, "(hours)"
    write(OPT,*) "Warmup time:",  (t3-t2)/3600/24, "(days)"
  endif

  ! We divide all the measurement into nBin,
  ! each having nPass/nBin pass.
  nBin   = Hub%P0%nBin
  nIter  = Hub%nPass / Hub%tausk / nBin
  if (nIter > 0) then
     do i = 1, nBin
        do j = 1, nIter
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
  endif

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
  endif

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

    write(OPT,*) "Meas time:",  t4-t3, "(second)"
    write(OPT,*) "Meas time:",  (t4-t3)/60, "(minutes)"
    write(OPT,*) "Meas time:",  (t4-t3)/3600, "(hours)"
    write(OPT,*) "Meas time:",  (t4-t3)/3600/24, "(days)"
    write(OPT,*) " "
  endif

  ! Prepare output file
 if (qmc_sim%rank == qmc_sim%aggr_root) then
   if (comp_tdm > 0) then
      call DQMC_open_file(adjustl(trim(ofile))//'.tdm.out','unknown', TDM_UNIT)
   endif

   ! Print computed results
   call DQMC_Hub_OutputParam(Hub, OPT)
   call DQMC_Phy0_Print(Hub%P0, Hub%S, OPT)
   call DQMC_TDM1_Print(tm, TDM_UNIT)

   ! G(tau) local
   if (comp_tdm > 0) then
     call DQMC_open_file('G_'//adjustl(trim(ofile)),'replace', OPT1)
     call DQMC_TDM1_Print1(tm, OPT1)
   endif
 endif

  !Aliases for Fourier transform
  na  =  Gwrap%lattice%natom
  nt  =  Gwrap%lattice%ncell
  nkt =  Gwrap%RecipLattice%nclass_k
  nkg =  Gwrap%GammaLattice%nclass_k
 
  !Print info on k-points and construct clabel
!  call DQMC_Print_HeaderFT(Gwrap, OPT, .true.)
!  call DQMC_Print_HeaderFT(Gwrap, OPT, .false.)

  !Compute Fourier transform
!  call DQMC_phy0_GetFT(Hub%P0, Hub%S%D, Hub%S%gf_phase, Gwrap%RecipLattice%FourierC, &
!       Gwrap%GammaLattice%FourierC, nkt, nkg, na, nt)
!  call DQMC_Phy0_GetErrFt(Hub%P0)
!  call DQMC_Phy0_PrintFT(Hub%P0, na, nkt, nkg, OPT)

  !Compute Fourier transform and error for TDM's
!  call DQMC_TDM1_GetKFT(tm)
!  call DQMC_TDM1_GetErrKFT(tm)
!  call DQMC_TDM1_PrintKFT(tm, TDM_UNIT)

  !Compute and print the self-energy
!  call DQMC_TDM1_SelfEnergy(tm, tau, TDM_UNIT)

!  if(Hub%P2%compute)then
!     if(Hub%P2%diagonalize)then
        !Obtain waves from diagonalization
!        call DQMC_Phy2_GetIrrep(Hub%P2, Hub%S)
        !Get error for waves
!        call DQMC_Phy2_GetErrIrrep(Hub%P2, Hub%P0%G_fun, Hub%S)
        !Analyze symmetry of pairing modes
!        call DQMC_Phy2_WaveSymm(Hub%S,Hub%P2,Gwrap%SymmOp)
        !Print Pairing info
!        call dqmc_phy2_PrintSymm(Hub%S, Hub%P2, OPT)
!     else
!        call dqmc_phy2_print(Hub%P2, Hub%S%wlabel, OPT)
!     endif
!  endif

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

