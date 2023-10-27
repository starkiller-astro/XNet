!***************************************************************************************************
! nse_slice.f90 10/18/17
! The program is an example driver to the provided NSE routines.
!
! For a given density, temperature,  and Ye, this program solves for NSE
!***************************************************************************************************

Program nse_slice
  Use, Intrinsic :: iso_fortran_env, Only: iostat_end
  Use nuclear_data, Only: ny, aa, zz, benuc, nname, read_nuclear_data
  Use reaction_data, Only: read_reaction_data
  Use xnet_controls, Only: descript, iconvc, idiag, iheat, iprocess, iscrn, isolv, itsout, iweak0, &
    & changemx, tolm, tolc, yacc, ymin, tdel_maxmult, kstmx, kitmx, bin_file_base, lun_diag, &
    & lun_stdin, lun_stdout, read_controls, myid, nproc
  Use xnet_eos, Only: eos_initialize
  Use xnet_nse, Only: xnse, ynse, i_nn, i_pp, i_he4, i_ni56, knrtot, nse_initialize, nse_solve
  Use xnet_parallel, Only: parallel_initialize, parallel_finalize, parallel_myproc, parallel_nprocs
  Use xnet_preprocess, Only: net_preprocess
  Use xnet_timers, Only: timer_nse, timer_nseinit, timer_nsesolv, timer_nsenrap, timer_nsels, &
    & timer_nseeval, timer_nsescrn
  Use xnet_types, Only: dp
  Use xnet_util, Only: name_ordered, replace_tabs, readnext, xnet_terminate
  Implicit None

  ! Local variables
  Integer, Parameter :: itimer_reset = 0
  Integer, Parameter :: max_line_length = 1024
  Character(80) :: data_dir, data_desc, thermo_desc, abund_desc
  Character(80) :: diag_file_base = 'nse_diag'
  Character(80) :: diag_file, bin_file, inab_file, thermo_file
  Character(max_line_length) :: line
  Real(dp) :: rho, ye, t9, t, tdel, flx, edot, enb, enm, uvec(2)
  Integer :: i, kstep, mflx, ifl_orig, ifl_term, lun_ts, pos, ierr
  Logical :: uvec_guess

  ! Identify number of MPI nodes and ID number for this PE
  Call parallel_initialize()
  myid = parallel_myproc()
  nproc = parallel_nprocs()

  ! Read user-defined controls
  Call read_controls(data_dir)

  ! Open diagnositic output file
  If ( idiag >= 0 ) Then
    diag_file = trim(diag_file_base)
    Call name_ordered(diag_file,0,1)
    Call name_ordered(diag_file,1,1)
    Open(newunit=lun_diag, file=diag_file)
  EndIf

  ! Open the binary time series file
  If ( itsout > 0 ) Then
    bin_file = trim(bin_file_base)
    Call name_ordered(bin_file,1,1)
    Open(newunit=lun_ts, file=bin_file, form='unformatted')

    ! Write Control Parameters to ts file
    Write(lun_ts) (descript(i),i=1,3),data_desc
    Write(lun_ts) kstmx,kitmx,iweak0,iscrn,iconvc,changemx,tolm,tolc,yacc,ymin,tdel_maxmult,iheat,isolv

    ! Write abundance description to ts file
    inab_file = " "
    abund_desc = " "
    Write(lun_ts) inab_file,abund_desc

    ! Write thermo description to ts file
    thermo_file = " "
    thermo_desc = " "
    Write(lun_ts) thermo_file,thermo_desc

    ! Write species identifiers to ts file
    Write(lun_ts) ny,zz,aa

    ! Write flux identifiers to ts file
    mflx = 1
    ifl_orig = 0
    ifl_term = 0
    Write(lun_ts) mflx,ifl_orig,ifl_term
  EndIf

  ! In requested, pre-process the nuclear and reaction data.
  If ( iprocess > 0 ) Call net_preprocess(lun_stdout,data_dir,data_dir)

  ! Read nuclear dataset
  Call read_nuclear_data(data_dir,data_desc)
  Call read_reaction_data(data_dir)

  ! Initialize EoS for screening
  If ( iscrn > 0 ) Call eos_initialize

  ! Initialize NSE
  Call nse_initialize

  kstep = 0
  Do
    kstep = kstep + 1
    uvec_guess = .false.
    line(:) = ' '
    If ( itsout >= 0 ) Write(lun_stdout,"(a)") 'Rho? T9? Ye? Mu_n guess? Mu_p guess?'
    Read(lun_stdin,"(a)",iostat=ierr) line
    If ( ierr == iostat_end .or. len_trim(line) <= 0 ) Then
      Exit
    ElseIf ( ierr /= 0 ) Then
      Call xnet_terminate('Failed while trying to inputs')
    Else
      Call replace_tabs(line)

      pos = 1
      Call readnext(line,pos,rho)
      Call readnext(line,pos,t9)
      Call readnext(line,pos,ye)
      If ( pos == 0 ) Call xnet_terminate('Not enough inputs')

      ! See if guess for chemical potentials provided
      Call readnext(line,pos,uvec(1))
      Call readnext(line,pos,uvec(2))
      If ( pos /= 0 ) uvec_guess = .true.
    End If

    If ( itsout >= 0 )  Then
      Write(lun_stdout,"(a,3es15.7)") 'NSE for Rho, T9, Ye=',rho,t9,ye
      If ( uvec_guess ) Write(lun_stdout,"(2(a,es15.7))") 'Mu_n=',uvec(1),'and Mu_p=',uvec(2)
    EndIf

    ! Calculate NSE abundances
    If ( uvec_guess ) Then
      Call nse_solve(rho,t9,ye,uvec)
    Else
      Call nse_solve(rho,t9,ye)
    EndIf

    If ( idiag >= 0 ) Then
      Write(lun_diag,'(a12,7es12.5)') 'NSE solved',rho,t9,ye,xnse([i_nn,i_pp,i_he4,i_ni56])
      If ( idiag >=1 ) Then
        Write(lun_diag,"(a5,es14.7)") (nname(i),xnse(i),i=1,ny)
      EndIf
      Write(lun_diag,"(a14,a5,3a10)") 'NSE Counters: ','Zone','NR','LineSrch','FuncEval'
      Write(lun_diag,"(14x,i5,3i10)") kstep,(knrtot(i),i=1,3)
      Write(lun_diag,"(a8,7a10)") &
        & 'Timers: ','Total','Solve','NR','LineSrch','FuncEval','Screening','Init'
      Write(lun_diag,"(8x,7es10.3)") &
        & timer_nse,timer_nsesolv,timer_nsenrap,timer_nsels,timer_nseeval,timer_nsescrn,timer_nseinit
    EndIf
    If ( itimer_reset > 0 ) Call reset_timers

    If ( itsout > 0 ) Then
      Call benuc(ynse,enb,enm)
      t = 0.0_dp
      tdel = 1.0_dp
      flx = 0.0_dp
      edot = -enm
      Write(lun_ts) kstep, t, t9, rho, tdel, edot, ynse, flx
    EndIf
  EndDo

  If ( idiag >= 0 ) Close(lun_diag)
  If ( itsout >= 1 ) Close(lun_ts)

  ! Wait for all nodes to finish
  Call parallel_finalize()

End Program nse_slice
