!*******************************************************************************
! Data Distribute, part of XNet 6, 1/31/12
!
! Needed for MPI execution.
! These routines broadcast the nuclear and network data between MPI tasks.
!
!*******************************************************************************

Module file_data
!=======================================================================
! This module contains data about the input and output files which must 
! be shared.
!=======================================================================
  Use controls
  
! Problem description 
  Character (LEN=80) :: descript(3) ! Control file description
! Input filename bases  
  Character (LEN=80) :: inab_file_base,thermo_file_base
! Output filename bases  
  Character (LEN=80) :: ev_file_base,bin_file_base
! Names of nuclei to be output in ev file
  Character (LEN=5)  :: output_nuc(14)

End Module file_data

Subroutine control_bcast(data_dir)
!=======================================================================
! This routine reads and distrbutes the control file which contains the 
! parameters which control the action of the network.  Control also 
! indicates the relative directory from which the nuclear data should 
! be loaded, as well as the names of the files containing the initial 
! abundances and thermodynamic trajectories.
!=======================================================================
  Use controls
  Use file_data
  Use mpi
! Use neutrino_data                     !NNU
  Character (LEN=80), INTENT(OUT) :: data_dir
  Character (LEN=80) :: control_char(8)
  Integer :: lun_control,control_int(12),i,ierr
  real(8) :: control_real(9)
 
! PE0 ...
  If(myid==0) Then

!-----------------------------------------------------------------------
! The control file contains the parameters which determine the actions 
! of XNet.  
!-----------------------------------------------------------------------
    lun_control=5
    Open(lun_control,FILE='control')                                         

! Read Problem Description
    call find_controls_block(lun_control,'Problem Description',ierr)  
    Read(lun_control,"(a80)") (descript(i),i=1,3) ! text description of the problem.
  
! Read Job Controls
    call find_controls_block(lun_control,'Job Controls',ierr)  
    Read(lun_control,*) szone        ! number of the zone with which to begin
    Read(lun_control,*) nzone        ! total # of zones
    Read(lun_control,*) iweak0       ! controls the treatment of weak reactions 
    Read(lun_control,*) iscrn        ! controls the treatment of nuclear screening
    Read(lun_control,*) iprocess     ! controls the runtime pre-processing of the network data 
    iprocess = 0   ! Always assume network data is pre-processed

! Read Integration Controls
    call find_controls_block(lun_control,'Integration Controls',ierr)  
    Read(lun_control,*) isolv        ! Choice of integrations scheme
    Read(lun_control,*) kstmx        ! max # of timesteps for each zone
    Read(lun_control,*) kitmx        ! max # of iterations before retry
    Read(lun_control,*) ijac         ! rebuild jacobian every ijac iterations after the first
    Read(lun_control,*) iconvc       ! determines which convergence condition is Used
    Read(lun_control,*) changemx     ! allowed abundance change used to set the timestep.
    Read(lun_control,*) yacc         ! abundances > yacc Used for timestep calculation
    Read(lun_control,*) tolm         ! mass conservation convergence criterion
    Read(lun_control,*) tolc         ! convergence limit on the iterative abundance change
    Read(lun_control,*) ymin         ! abundance < ymin is set to 0.0
    Read(lun_control,*) tdel_maxmult ! max factor by which the timestep is changed

! Read Self-heating Controls
    call find_controls_block(lun_control,'Self-heating Controls',ierr)  
    If(ierr/=0) Then
      Read(lun_control,*) iheat      ! controls the coupling of the network to temperature
      Read(lun_control,*) changemxt  ! allowed temperature change used to set the timestep.
      Read(lun_control,*) tolt9      ! convergence limit on the iterative temperature change
    Else
      Write(6,*) 'Using Default Self-heating behavior'
      iheat = 0
      changemxt = 1.0e-2
      tolt9 = 1.0e-4
    EndIf

! Read NSE Initial Abundance Controls
! temperature at which NSE initial conditions are used, t9nse =8 is default.
!   call find_controls_block(lun_control,'NSE Initial Conditions',ierr)   !NSE
!   If(ierr/=0) Then                                                      !NSE
!     Read(lun_control,*) t9nse                                           !NSE
!   Else                                                                  !NSE
!     Write(6,*) 'Using Default NSE behavior'                             !NSE
!     t9nse = 8.0                                                         !NSE
!   Endif                                                                 !NSE

! Read Neutrino Controls
!   call find_controls_block(lun_control,'Neutrinos',ierr)                !NNU
!   If(ierr/=0) Then                                                      !NNU
!     Read(lun_control,*) ineutrino   ! controls neutrino reactions       !NNU
!   Else                                                                  !NNU
!     Write(6,*) 'Using Default behavior'                                 !NNU 
!     ineutrino=0                                                         !NNU
!   Endif                                                                 !NNU

!-------------------------------------------------------------------------------
! XNet output controls include the base of the filenames to which ASCII and 
! binary output are written at the end of each timestep, and a subset of nuclei
! to be included in the per timestep ASCII output file.
!-------------------------------------------------------------------------------
! Read Output Controls
    call find_controls_block(lun_control,'Output Controls',ierr)  
    Read(lun_control,*) idiag        ! sets diagnostic output level
    Read(lun_control,*) itsout       ! sets per timestep output level
    Read(lun_control,"(72x)")
    Read(lun_control,"(a80)") ev_file_base            
    Read(lun_control,"(72x)")
    Read(lun_control,"(a80)") bin_file_base            
    Read(lun_control,"(72x)")
    Read(lun_control,"(14a5)") output_nuc         

!-------------------------------------------------------------------------------
! XNet input controls include the relative directory from which the nuclear data 
! should be loaded, as well as the names of the files containing the initial 
! abundances and thermodynamic trajectories.
!-------------------------------------------------------------------------------
! Read Input Controls
    call find_controls_block(lun_control,'Input Controls',ierr)  
    Read(lun_control,"(72x)")
    Read(lun_control,"(a80)") data_dir
    Read(lun_control,"(72x)")
    Read(lun_control,"(a80)",IOSTAT=ierr) inab_file_base
    Read(lun_control,"(a80)",IOSTAT=ierr) thermo_file_base
    If(ierr > 0) Then
      Write(6,*) 'Problem reading Input Filenames'
    Endif
    Close(lun_control)

!   Load Control passing arrays
    control_int(1)  = nzone ; control_int(2)  = iweak0; control_int(3)  = iscrn
    control_int(4)  = isolv ; control_int(5)  = kstmx ; control_int(6)  = kitmx
    control_int(7)  = ijac  ; control_int(8)  = iconvc; control_int(9)  = idiag
    control_int(10) = itsout; control_int(11) = iheat
!   control_int(12) = ineutrino !NNU
    
    control_real(1) = changemx; control_real(2) = yacc     ; control_real(3) = tolm
    control_real(4) = tolc    ; control_real(5) = ymin     ; control_real(6) = tdel_maxmult
    control_real(7) = t9nse   ; control_real(8) = changemxt; control_real(9) = tolt9
    
    control_char(1:3) = descript      ; control_char(4) = data_dir  
    control_char(5)   = bin_file_base ; control_char(6) = ev_file_base  
    control_char(7)   = inab_file_base; control_char(8) = thermo_file_base
  Endif

! All PE 
  
! Broadcast network control parameters
  call mpi_bcast(control_int,12,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(control_real,9,MPI_REAL8, 0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(control_char,8*80,MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(output_nuc,14*5,MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)

! Unpack network control parameters
  If(myid/=0) Then
    nzone  = control_int(1) ; iweak0 = control_int(2) ; iscrn = control_int(3)
    isolv  = control_int(4) ; kstmx  = control_int(5) ; kitmx = control_int(6)
    ijac   = control_int(7) ; iconvc = control_int(8) ; idiag = control_int(9)
    itsout = control_int(10); iheat  = control_int(11) 
!   ineutrino = control_int(12) !NNU
    
    changemx = control_real(1); yacc      = control_real(2); tolm         = control_real(3)
    tolc     = control_real(4); ymin      = control_real(5); tdel_maxmult = control_real(6)
    t9nse    = control_real(7); changemxt = control_real(8); tolt9        = control_real(9)
    
    descript       = control_char(1:3) ; data_dir         = control_char(4)
    bin_file_base  = control_char(5)   ; ev_file_base     = control_char(6)    
    inab_file_base = control_char(7)   ; thermo_file_base = control_char(8)

  Endif

! Retain the value of iweak
! iweak0=iweak
!$OMP PARALLEL DEFAULT(SHARED)
  iweak=iweak0
  
  Write(lun_diag,*) 'CBcast',idiag,kitmx,kstmx
!$OMP End PARALLEL
End subroutine control_bcast

Subroutine netdata_bcast(data_dir,data_desc)
!=======================================================================
! This routine handles nuclear data I/O for MPI versions, broadcasting 
! the necessary nuclear and reaction data from PE0 to the production PEs
!=======================================================================
  Use controls
  Use file_data
  Use nuclear_data
  Use part_funct_data 
  Use cross_sect_data
  Use ffn_data
! Use nnu_data                  !NNU
! Use neutrino_data             !NNU
  Use reac_rate_data
  Use mpi
  Character (LEN=80), INTENT(IN) :: data_dir  
  Character (LEN=80), INTENT(OUT) :: data_desc  
  Integer :: matrix_shape(5)
  Integer :: nds,ndi,idel,jdel,ijd
  Integer :: nbc,ierr,nr,j

! On PE0 ...
  If(myid==0) Then

! ... read the nuclear and reaction data
    call read_nuclear_data(data_dir,data_desc)
    call read_reaction_data(data_dir)

  Endif

! Share data for nuc_number module
  call mpi_bcast(ny,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

! Share data description
  call mpi_bcast(data_desc,80,MPI_CHARACTER, 0,MPI_COMM_WORLD,ierr)

! Share data for the nuclear_data module
  If(myid/=0) Allocate(nname(0:ny),aa(ny),zz(ny),nn(ny),be(ny),mex(ny))
  call mpi_bcast(aa,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(zz,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(nn,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(be,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(mex,ny,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(izmax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=5*(ny+1)
  call mpi_bcast(nname,nbc,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

! Share data for the part_funct_data module
  If(myid/=0) Then
    Allocate(g(24,ny),angm(0:ny))
!$OMP PARALLEL
    Allocate (gg(0:ny))
    If(iheat>0) Allocate (dlngdt9(0:ny))
!$OMP End PARALLEL
  Endif
  call mpi_bcast(t9i,24,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=24*ny
  call mpi_bcast(g,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=ny+1
  call mpi_bcast(angm,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)

! Share data for the cross_sect_data module
  call mpi_bcast(nreac,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If(myid/=0) Then
    nr=nreac(1)
    Allocate (rc1(7,nr),q1(nr),n1i(4,nr),iwk1(nr),ires1(nr),irev1(nr),iffn(nr),innu(nr))
!$OMP PARALLEL      
    Allocate (csect1(nr),rpf1(nr),h1(nr))
    If(iheat>0) Allocate (dcsect1dt9(nr),dlnrpf1dt9(nr),dh1dt9(nr))
!$OMP End PARALLEL
  Endif
  nbc=7*nreac(1)
  call mpi_bcast(rc1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=4*nreac(1)
  call mpi_bcast(n1i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=nreac(1)
  call mpi_bcast(iffn,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(innu,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iwk1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(ires1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(irev1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  If(myid/=0) Then
    nr=nreac(2)
    Allocate (rc2(7,nr),q2(nr),n2i(5,nr),iwk2(nr),ires2(nr),irev2(nr))
!$OMP PARALLEL      
    Allocate (csect2(nr),rpf2(nr),h2(nr))
    If(iheat>0) Allocate (dcsect2dt9(nr),dlnrpf2dt9(nr),dh2dt9(nr))
!$OMP End PARALLEL
  Endif
  nbc=7*nreac(2)
  call mpi_bcast(rc2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=5*nreac(2)
  call mpi_bcast(n2i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=nreac(2)
  call mpi_bcast(iwk2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(ires2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(irev2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  If(myid/=0) Then
    nr=nreac(3)
    Allocate (rc3(7,nr),q3(nr),n3i(6,nr),iwk3(nr),ires3(nr),irev3(nr))
!$OMP PARALLEL      
    Allocate (csect3(nr),h3(nr))
    If(iheat>0) Allocate (dcsect3dt9(nr),dh3dt9(nr))
!$OMP End PARALLEL
  Endif
  nbc=7*nreac(3)
  call mpi_bcast(rc3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=6*nreac(3)
  call mpi_bcast(n3i,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=nreac(3)
  call mpi_bcast(iwk3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(ires3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 
  call mpi_bcast(irev3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr) 

! Share the data for the reac_rate_data module
  call mpi_bcast(nan,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If(myid/=0) Then
    Allocate (la(3,ny),le(3,ny))
    Allocate (mu1(nan(1)),a1(nan(1)))
    Allocate (mu2(nan(2)),a2(nan(2)))
    Allocate (mu3(nan(3)),a3(nan(3)))
!$OMP PARALLEL
    Allocate (b1(nan(1)),b2(nan(2)),b3(nan(3)))
    If(iheat>0) Allocate (dr1dt9(nan(1)),dr2dt9(nan(2)),dr3dt9(nan(3)))
!$OMP End PARALLEL
  Endif
  nbc=3*ny
  call mpi_bcast(la,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(le,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=nan(1)
  call mpi_bcast(mu1,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(a1,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=nan(2)
  call mpi_bcast(mu2,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(a2,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  nbc=nan(3)
  call mpi_bcast(mu3,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(a3,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  If(myid/=0) Then
    Allocate (n11(nan(1)),n21(nan(2)),n22(nan(2)))
    Allocate (n31(nan(3)),n32(nan(3)),n33(nan(3)))
    Do j=1,nan(1)
      n11(j)=n1i(1,mu1(j))
    Enddo
    Do j=1,nan(2)
      n21(j)=n2i(1,mu2(j))
      n22(j)=n2i(2,mu2(j))
    Enddo
    Do j=1,nan(3)
      n31(j)=n3i(1,mu3(j))
      n32(j)=n3i(2,mu3(j))
      n33(j)=n3i(3,mu3(j))
    Enddo
    If(iscrn>0) Then
      nr=nreac(2)
      Allocate (z21(nr),z22(nr),iz21(nr),iz22(nr))
      z21=zz(n2i(1,:)) ; z22=zz(n2i(2,:))
      iz21=int(z21)    ; iz22=int(z22)
      Allocate (z12(nr),iz12(nr),z12i(nr))
      z12=z21*z22
      iz12=iz21+iz22
      z12i=((z21+z22)**1.86-z21**1.86-z22**1.86)
      nr=nreac(3)
      Allocate (z31(nr),z32(nr),z33(nr),iz31(nr),iz32(nr),iz33(nr))
      z31=zz(n3i(1,:)) ; z32=zz(n3i(2,:)) ; z33=zz(n3i(3,:))
      iz31=int(z31)    ; iz32=int(z32)    ; iz33=int(z33)
      Allocate (z123(nr),iz123(nr))
      z123=z31*z32+z31*z33+z32*z33
      iz123=iz31+iz32+iz33
    EndIf
  Endif

! Share the matrix shape parameters
! If(myid==0) Then
!   matrix_shape(1) = nds ; matrix_shape(2) = ndi; matrix_shape(3) = idel
!   matrix_shape(4) = jdel; matrix_shape(5) = ijd
! Endif
! nbc=5
! call mpi_bcast(matrix_shape,nbc,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
! If(myid/=0) Then
!   nds  = matrix_shape(1); ndi = matrix_shape(2); idel = matrix_shape(3)
!   jdel = matrix_shape(4); ijd = matrix_shape(5)
! Endif

! call mpi_bcast(nnnu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)      !NNU
! If(nnnu>0) Then                                               !NNU
!   If(myid/=0) Then                                            !NNU
!     Allocate (sigmanu(nnnu,7))                                !NNU
!!$OMP PARALLEL DEFAULT(SHARED)                                 !NNU
!     Allocate (rcsnu(nnnu,4))                                  !NNU
!!$OMP End PARALLEL                                             !NNU
!   Endif                                                       !NNU
!   nbc=nnnu*7                                                  !NNU
!   call mpi_bcast(sigmanu,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr) !NNU
! ElseIf(myid/=0) Then                                          !NNU
!   Allocate (sigmanu(1,7))                                     !NNU
!   Allocate(rcsnu(1,7))                                        !NNU
! Endif                                                         !NNU

! Share the data for the ffn_data
  call mpi_bcast(nffn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If(nffn>0) Then
    If(myid/=0) Then
      Allocate (ffnsum(nffn,143),ffnenu(nffn,143))
    Endif
    nbc=nffn*143
    call mpi_bcast(ffnsum,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
    call mpi_bcast(ffnenu,nbc,MPI_REAL8,0,MPI_COMM_WORLD,ierr)
  ElseIf(myid/=0) Then
    Allocate(ffnsum(1,143),ffnenu(1,143))
  Endif
!$OMP PARALLEL DEFAULT(SHARED)
  Write(lun_diag,*) 'NBcast',ny,nffn,nreac
!$OMP End PARALLEL
  Return
  End

  Subroutine match_bcast(data_dir)
!=======================================================================
! This routine reads in the match_data modules on PE 0 and broadcasts 
! the data to other processors.  It also initializes the flux data
!=======================================================================
  Use controls
  Use file_data
  Use nuclear_data
  Use cross_sect_data
  Use match_data
  Use flux_data
  Use mpi
  Character (LEN=80), INTENT(IN) :: data_dir  
  Integer nbc,ierr

! The control PE...
  If(myid==0) Then

! ... reads in the flux matching data
    call read_match_data(data_dir)
  Endif

! Share the match data
  call mpi_bcast(mflx,   1,    MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  If(myid/=0) Then
      Allocate(ifl1(nreac(1)),ifl2(nreac(2)),ifl3(nreac(3)))
      Allocate (nflx(7,mflx),qflx(mflx),iwflx(mflx),descx(mflx))
  Endif
  call mpi_bcast(ifl1,nreac(1),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ifl2,nreac(2),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(ifl3,nreac(3),MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=7*mflx
  call mpi_bcast(nflx,nbc,     MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(qflx,mflx,    MPI_REAL8,  0,MPI_COMM_WORLD,ierr)
  call mpi_bcast(iwflx,mflx,   MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nbc=4*mflx
  call mpi_bcast(descx,nbc,    MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

! The arrays necessary to compute fluxes are initialized
  call flux_init

!$OMP PARALLEL DEFAULT(SHARED)
  Write(lun_diag,*) 'MBcast',mflx
!$OMP End PARALLEL
  Return
  End

