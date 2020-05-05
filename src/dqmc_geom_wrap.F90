module DQMC_GEOM_WRAP
  use DQMC_GEOM_PARAM
  use DQMC_HAMILT
  use DQMC_SYMM
  use DQMC_LATT
  use DQMC_RECLATT
  use DQMC_BONDS
  use DQMC_Cfg
  use DQMC_STRUCT

  implicit none


  type GeomWrap
     type(lattice_t)       :: Lattice
     type(recip_lattice_t) :: RecipLattice
     type(recip_lattice_t) :: GammaLattice
     type(hamiltonian_t)   :: Hamilt
     type(bonds_t)         :: Bonds
     type(symm_operations) :: SymmOp
     type(pairing)         :: Pairs
  end type

  contains

  subroutine DQMC_Geom_Fill(gwrap, gfile, cfg, optsym, SOP)

    type(GeomWrap),    intent(inout)  :: gwrap
    type(config), intent(inout)       :: cfg
    character(len=slen), intent(in)   :: gfile
    integer, intent(in)               :: optsym, SOP
    logical                           :: found, connected

    inquire(file=gfile,exist=found)
    if(found)then
     !Check whether file is already connected
     inquire(file=gfile,opened=connected)
     if(connected)then
      !If connected retrieve the unit
      inquire(file=gfile,number=inpunit)
      rewind(inpunit)
     else
      !If not connected, find a free unit
      do inpunit=7,99
        inquire(unit=inpunit,opened=connected)
        if(.not.connected)exit
      enddo
      !open file
      open(file=gfile,unit=inpunit)
     endif
    else 
     call DQMC_Error("cannot open geom def file "//gfile, 0)
    endif

    !Scan geometry file to see which fields are specified, see dqmc_geom_param.F90
    call analyze_input
 
    !Initialize basic info about real space cluster, see dqmc_latt.F90
    call init_lattice(gwrap%Lattice)

    !Construct full lattice in dqmc_latt.F90
    call construct_lattice(gwrap%Lattice)
 
    !Initialize basic info about reciprocal lattice
    call init_recip_latt(gwrap%Lattice,gwrap%RecipLattice,.true.,cfg)
    call init_recip_latt(gwrap%Lattice,gwrap%GammaLattice,.false.,cfg)

    !construct reciprocal lattice shifted by k-point
    call construct_recip_lattice(gwrap%RecipLattice)

    !construct reciprocal lattice that includes Gamma
    call construct_recip_lattice(gwrap%GammaLattice)

    !Construct Hamiltonian in dqmc_hamilt.F90
    call construct_hamilt(gwrap%Hamilt,gwrap%Lattice,gwrap%RecipLattice,cfg)
    
    !Read point-symmetry (optional)
    call read_symm(gwrap%SymmOp)

    !pair-to-pair map of action of each symmetry in real space (also translations)
    call map_symm_lattice(gwrap%SymmOp,gwrap%Lattice, gwrap%Hamilt, optsym, SOP)

    !point-to-point map of action of each symmetry in reciprocal space
    call map_symm_recip_lattice(gwrap%SymmOp,gwrap%RecipLattice,.true.)
    call map_symm_recip_lattice(gwrap%SymmOp,gwrap%GammaLattice,.false.)

    !Group pairs of lattice points into classes in dqmc_symm.F90
    call construct_lattice_classes(gwrap%SymmOp,gwrap%Lattice)

    !Group k-points into classes
    call construct_recip_lattice_classes(gwrap%SymmOp,gwrap%RecipLattice,.true.)
    call construct_recip_lattice_classes(gwrap%SymmOp,gwrap%GammaLattice,.false.)

    ! Fill the matrix of fourier weights
    call DQMC_Fill_FourierC(Gwrap%RecipLattice, Gwrap%Lattice)
    call DQMC_Fill_FourierC(Gwrap%GammaLattice, Gwrap%Lattice)

    !Group hopping part of Hamiltonian in classes in dqmc_hamilt.F90
    call count_hop_class(gwrap%Lattice,gwrap%Hamilt)

    !Group separately local classes for U and mu in dqmc_hamilt.F90
    call count_local_classes(gwrap%Lattice,gwrap%Hamilt)

    !Assign phase for Green's function
    call assign_gf_phase(gwrap%Lattice,gwrap%RecipLattice%ktwist)

    !Read Bonds (optional input)   
    call read_bonds(gwrap%Bonds)

    if(gwrap%Bonds%initialized)then
       !bond-to-bond map of action of symmetry on bonds
       call map_symm_bonds(gwrap%Bonds,gwrap%SymmOp,gwrap%Lattice)

       !Group pairs of bonds into classes
       call construct_bond_classes(gwrap%Bonds,gwrap%SymmOp)

       !Map bonds throughout the entire lattice
       call construct_pairs(gwrap%Bonds,gwrap%Pairs,gwrap%Lattice)

       !pair-to-pair map of action of symmetry on bonds
       call map_symm_pairs(gwrap%Pairs, gwrap%SymmOp)

       !Group pairs of pairs into classes
       call construct_pair_classes(gwrap%Bonds,gwrap%Pairs)
    endif

    !Assign phase to each atom (optional input)
    call assign_phase(gwrap%Lattice)

    !Write some info
    !call write_files(gwrap%Lattice,gwrap%RecipLattice,gwrap%Hamilt)

  end subroutine DQMC_Geom_Fill

 !---------------------------------------------------------------------!

  subroutine DQMC_Geom_Init(gwrap, S, cfg)    

    type(GeomWrap), intent(in)    :: gwrap   
    type(Struct),   intent(inout) :: S       ! Struct
    type(config),   intent(inout) :: cfg

    ! ... local scalar ...
    integer  :: n                ! Order of matrix T and D 
    integer  :: i, j             ! Loop iterator
    integer  :: ic, ib, model
    real(wp), pointer  :: clab(:,:) => null()
    real(wp), pointer  :: tupvalue(:) => null()
    real(wp), pointer  :: tdnvalue(:) => null()
    integer, pointer   :: tmp(:,:) => null()

    character(label_len) :: label 

    ! ... Executable ...
    S%checklist=.false.

    call CFG_Get(cfg, "model" , model)  
    n        = gwrap%Lattice%nsites
    S%nSite  = n
    S%nCell  = gwrap%Lattice%ncell
    S%nGroup = gwrap%Hamilt%nlocclass
    S%Name   = 'General Geometry - Free Format'
    allocate(tmp(n,n))

    !Fill T
    call group_hopping(gwrap%Hamilt, n, S%n_t, tmp, tupvalue, tdnvalue)
    call DQMC_CCS_Compress(n,-1, tmp, S%T)
    S%checklist(STRUCT_ADJ)=.true.

    !Fill ckb
    tmp = 0
    do ib = 0, size(gwrap%hamilt%tckb,2)-1
       i = gwrap%hamilt%tckb(1,ib) + 1
       j = gwrap%hamilt%tckb(2,ib) + 1
       tmp(i,j) = gwrap%hamilt%tckb(3,ib)
       tmp(j,i) = gwrap%hamilt%tckb(3,ib)
    enddo
    call DQMC_CCS_Compress(n,-1, tmp, S%ckb)

    !Symmetry Classes, see construct_lattice_classes in dqmc_symm.F90
    S%nClass = gwrap%Lattice%nclass

    ! 12/27/2015:
    ! get site/orbital indices and its vector of nClass for various usage
    ! 1. can be used for phase for computing S_AF in plane 
    ! for bilayer cases, only include z=0 component
    ! 2. for general FT:
    ! For unit cell including more than 1 atom, e.g. with staggered pot,
    ! QUEST code only compute FT between cells while still real space for intracell
    ! here use relative vector between sites for general FT
    allocate(S%vecClass(S%nClass,5))  

    allocate(S%D(n,n))
    allocate(S%F(S%nClass))
    allocate(S%clabel(S%nClass))
    allocate(S%AFphase(S%nClass))
    allocate(S%pi0phase(S%nClass))

    !dqmc_latt.F90
    S%D(1:n,1:n) =  gwrap%Lattice%myclass(0: n - 1, 0: n - 1)
    S%F(:)       =  gwrap%Lattice%class_size(:)
    clab         => gwrap%Lattice%class_label
    do ic = 1, S%nClass       
      write(S%clabel(ic),'(2(i4),i5,3(f12.7))') (int(clab(ic,j)),j=4,5),S%F(ic),(clab(ic,j),j=1,3)    

      ! 12/27/2015:
      ! Get site/orbital indices and relative vector between them for nClass
      label = S%clabel(ic)
      read(label(1:4),*) S%vecClass(ic,1)
      read(label(5:8),*) S%vecClass(ic,2)
      read(label(14:25),*) S%vecClass(ic,3)
      read(label(26:37),*) S%vecClass(ic,4)
      read(label(38:49),*) S%vecClass(ic,5)
      !write(*,*) 'label=', label(14:25), label(26:37), label(38:49)

      ! decide (-1)^(dx+dy) for computing S_AF and S_CDW in plane
      ! for Hubbard model
      if (model==0) then
        if ( mod(int(abs(S%vecClass(ic,3)))+int(abs(S%vecClass(ic,4))),2) == 0) then  ! (-1)**(x+y)=1
          S%AFphase(ic) = 1.0
        else
          S%AFphase(ic) = -1.0
        endif
        write(*,*) 'AFphase = '
        write(*,'(a6,2(f4.1),a5,3(f5.1),a7,f5.1)') 'sites', S%vecClass(ic,1),S%vecClass(ic,2), &
              'r=',S%vecClass(ic,3),S%vecClass(ic,4),S%vecClass(ic,5),'phase=',S%AFphase(ic)
      endif

      ! mod(int(S%vecClass(ic,1)),2)==1 .and. mod(int(S%vecClass(ic,2)),2)==1 if
      ! for PAM model, which only need S_AF for f-electrons
      if (model==1 .or. model==2) then
        if (abs(S%vecClass(ic,5))<0.0001 .and. mod(int(S%vecClass(ic,1)),2)==1 .and. mod(int(S%vecClass(ic,2)),2)==1) then                 
          if ( mod(int(abs(S%vecClass(ic,3)))+int(abs(S%vecClass(ic,4))),2) == 0) then  ! (-1)**(x+y)=1
            S%AFphase(ic) = 1.0
          else
            S%AFphase(ic) = -1.0
          endif
          write(*,*) 'AFphase = '
          write(*,'(a6,2(f4.1),a5,3(f5.1),a7,f5.1)') 'sites', S%vecClass(ic,1),S%vecClass(ic,2), &
                ' r=',S%vecClass(ic,3),S%vecClass(ic,4),S%vecClass(ic,5),'phase=',S%AFphase(ic)
        else
          S%AFphase(ic) = 0.0
        endif
      endif

      ! model3: stack PAM with additional f-orbital
      ! need S_AF for two layers of f-electrons
      if (model==3) then
        if (abs(S%vecClass(ic,5))<0.0001 .and. ((mod(int(S%vecClass(ic,1)),3)==1 .and. mod(int(S%vecClass(ic,2)),3)==1) .or. &
                                                (mod(int(S%vecClass(ic,1)),3)==2 .and. mod(int(S%vecClass(ic,2)),3)==2))) then
          if ( mod(int(abs(S%vecClass(ic,3)))+int(abs(S%vecClass(ic,4))),2) == 0) then  ! (-1)**(x+y)=1
            S%AFphase(ic) = 1.0
          else
            S%AFphase(ic) = -1.0
          endif
          write(*,*) 'AFphase = '
          write(*,'(a6,2(f4.1),a5,3(f5.1),a7,f5.1)') 'sites', S%vecClass(ic,1),S%vecClass(ic,2), &
                ' r=',S%vecClass(ic,3),S%vecClass(ic,4),S%vecClass(ic,5),'phase=',S%AFphase(ic)
        else
          S%AFphase(ic) = 0.0
        endif
      endif

      ! for stacked two PAMs coupled with additional V12
      ! need S_AF for two layers of f-electrons
      if (model==4) then
        if (abs(S%vecClass(ic,5))<0.0001 .and. ((mod(int(S%vecClass(ic,1)),4)==1 .and. mod(int(S%vecClass(ic,2)),4)==1) .or. &
                                                (mod(int(S%vecClass(ic,1)),4)==2 .and. mod(int(S%vecClass(ic,2)),4)==2))) then
          if ( mod(int(abs(S%vecClass(ic,3)))+int(abs(S%vecClass(ic,4))),2) == 0) then  ! (-1)**(x+y)=1
            S%AFphase(ic) = 1.0
          else
            S%AFphase(ic) = -1.0
          endif
          write(*,*) 'AFphase = '
          write(*,'(a6,2(f4.1),a5,3(f5.1),a7,f5.1)') 'sites', S%vecClass(ic,1),S%vecClass(ic,2), &
                ' r=',S%vecClass(ic,3),S%vecClass(ic,4),S%vecClass(ic,5),'phase=',S%AFphase(ic)
        else
          S%AFphase(ic) = 0.0
        endif
      endif

      ! decide (-1)^x for computing S(pi,0) in plane
      ! mod(int(S%vecClass(ic,1)),2)==1 .and. mod(int(S%vecClass(ic,2)),2)==1 if
      ! for PAM model, which only need S_AF for f-electrons
      if (abs(S%vecClass(ic,5))<0.0001 .and. mod(int(S%vecClass(ic,1)),2)==1 .and. mod(int(S%vecClass(ic,2)),2)==1) then
        if ( mod(int(abs(S%vecClass(ic,3))),2) == 0) then  ! (-1)**dx=1
          S%pi0phase(ic) = 1.0
        else
          S%pi0phase(ic) = -1.0
        endif
      else
        S%pi0phase(ic) = 0.0
      endif

    enddo

    !store GF phase
    allocate(S%gf_phase(n,n))
    do i = 1, n
       do j = 1, n
          !if(abs(int(gwrap%lattice%gf_phase(i-1,j-1)))/=1)stop 'Problem with gf_phase'
          !S%gf_phase(i,j)=int(gwrap%Lattice%gf_phase(i-1,j-1))
          S%gf_phase(i, j) = gwrap%Lattice%gf_phase(i - 1, j - 1)
       enddo
    enddo
    allocate(S%chi_phase(n, n))
    S%chi_phase = 1
    S%checklist(STRUCT_CLASS)=.true.    

    allocate(S%map(n))
    S%map(1:n)=gwrap%Hamilt%mylocclass(0:n-1)
    
    !Fill B
    if (Found_Field(BONDS_F)) then
       tmp = 0
       S%n_b = gwrap%Pairs%nbond
       do ic = 0, size(gwrap%Pairs%nbondv)-1
          do ib = 1, gwrap%Pairs%nbondv(ic)
             i = gwrap%Pairs%bond_origin(ib,ic)
             j = gwrap%Pairs%bond_end(ib,ic)
             tmp(i+1,j+1) = gwrap%Pairs%bond_number(ib,ic)
          enddo
       enddo
       call DQMC_CCS_Compress(n,-1, tmp, S%B)

       !Store symmetry in S     
       allocate(S%class_b(S%n_b,S%n_b),S%size_b(gwrap%Pairs%nclass_p))
       S%class_b=gwrap%Pairs%myclass_p
       S%size_b=gwrap%Pairs%class_size_p  
       S%nClass_b=gwrap%Pairs%nclass_p
       S%checklist(STRUCT_BOND)=.true.    

       !Waves
       S%nWave=gwrap%Pairs%nWave
       allocate(S%wlabel(S%nWave))
       S%wlabel(:)=gwrap%Pairs%wave_label(:)
       if(Found_Field(PAIRS_F))then
          allocate(S%W(S%n_b,S%nWave))
          do i=1,S%nWave
             S%W(1:S%n_b,i)=gwrap%Pairs%bond_wgt(i,1:S%n_b)
          enddo
          S%checklist(STRUCT_WAVE)=.true.
       endif
    endif

    deallocate(tmp)
    
    if(Found_field(PHASE_F))then 
       allocate(S%P(n))
       S%P(1:n)=gwrap%Lattice%phase(0:n-1)
       S%checklist(STRUCT_PHASE)=.true.
    endif

    !Setting variables

    !Set variables that are otherwise read from main input
    !Note that mu_up becomes array for all orbitals, similar for others
    call CFG_Set(cfg,"n",n)
    call CFG_Set(cfg,"t_up",S%n_t,tupvalue)
    call CFG_Set(cfg,"t_dn",S%n_t,tdnvalue)
    call CFG_Set(cfg,"U",S%nGroup,gwrap%Hamilt%Uvalue)
    call CFG_Set(cfg,"mu_up",S%nGroup,gwrap%Hamilt%muupvalue)
    call CFG_Set(cfg,"mu_dn",S%nGroup,gwrap%Hamilt%mudnvalue)

    deallocate(tupvalue)
    deallocate(tdnvalue)

    S%checklist(STRUCT_INIT)=.true.

  end subroutine DQMC_Geom_Init

  !--------------------------------------------------------------------!

  subroutine DQMC_Print_HeaderFT(Gwrap, OPT, applytwist)
    use dqmc_mpi

    type(GeomWrap), intent(in)  :: Gwrap
    integer, intent(in)         :: OPT
    logical, intent(in) :: applytwist

    ! ... Local scalar ...
    integer :: na, nk, ii, i, ik, j, jj, nkpts
    integer, pointer :: myclass_k(:) 
    real*8, pointer  :: klist(:,:) 

    ! ... Executable ...

    if (qmc_sim%rank .ne. 0) return

    na   = Gwrap%Lattice%natom
    if(applytwist)then 
       nk         = Gwrap%RecipLattice%nclass_k
       nkpts      = Gwrap%RecipLattice%nkpts
       myclass_k => Gwrap%RecipLattice%myclass_k
       klist     => Gwrap%RecipLattice%klist
    else
       nk         = Gwrap%GammaLattice%nclass_k
       nkpts      = Gwrap%GammaLattice%nkpts
       myclass_k => Gwrap%GammaLattice%myclass_k
       klist     => Gwrap%GammaLattice%klist
    endif

    !Print general info about k-space
    if(applytwist)then
      write(OPT,'(A)')' Grid for Green''s function'
    else
      write(OPT,'(A)')' Grid for spin/charge correlations'
    endif

    write(OPT,'(A)')'  K-points'
    ii=0
    write(OPT,'(A)')'  Class'
    do i=1,nk
      jj=1
      do ik=1,nkpts
        if(myclass_k(ik)==i)then
          if(jj==1)then
            write(OPT,'(2x,i3,6x,3(f10.5))')i,(klist(ik,j),j=1,Gwrap%Lattice%ndim)
            jj=-1
          else
            write(OPT,'(11x,3(f10.5))')(klist(ik,j),j=1,Gwrap%Lattice%ndim)
          endif
        endif
      enddo
      write(OPT,*)
    enddo
    write(OPT,FMT_DBLINE)

  end subroutine DQMC_Print_HeaderFT

  !--------------------------------------------------------------------!

end module DQMC_GEOM_WRAP
