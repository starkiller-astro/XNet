!***************************************************************************************************
! nse_slice.f90 10/18/17
! The program is an example driver to the provided NSE routines.
!
! For a given density, temperature,  and Ye, this program solves for NSE
!***************************************************************************************************

Program nse_slice
  Use nuclear_data, Only: ny, aa, zz, benuc, read_nuclear_data
  Use reaction_data, Only: read_reaction_data
  Use xnet_controls, Only: descript, iconvc, idiag, iheat, iprocess, iscrn, isolv, itsout, iweak0, &
    & changemx, tolm, tolc, yacc, ymin, tdel_maxmult, kstmx, kitmx, bin_file_base, lun_diag, &
    & lun_stdin, lun_stdout, read_controls
  Use xnet_eos, Only: eos_initialize
  Use xnet_nse, Only: ynse, nse_initialize, nse_solve
  Use xnet_preprocess, Only: net_preprocess
  Use xnet_types, Only: dp
  Use xnet_util, Only: name_ordered
  Implicit None

  ! Local variables
  Character(80) :: data_dir, data_desc, thermo_desc, abund_desc
  Character(80) :: diag_file_base = 'nse_diag'
  Character(80) :: diag_file, bin_file, inab_file, thermo_file
  Real(dp) :: rho, ye, t9, t, tdel, flx, edot, enb, enm, ytot, ztot, atot
  Integer :: i, kstep, mflx, ifl_orig, ifl_term, lun_ts

  Write(lun_stdout,"(a)") 'Rho? Ye? T9?'
  Read(lun_stdin,*) rho, ye, t9
  Write(lun_stdout,"(a,f7.4,a,es11.3,a,f7.3)") 'NSE for Ye=',ye,'and Rho=',rho,', stopping at T9=',t9

  ! Read user-defined controls
  Call read_controls(data_dir)

  ! In requested, pre-process the nuclear and reaction data.
  If ( iprocess > 0 ) Call net_preprocess(lun_stdout,data_dir,data_dir)

  ! Open diagnositic output file
  If ( idiag >= 0 ) Then
    diag_file = trim(diag_file_base)
    Call name_ordered(diag_file,0,1)
    Call name_ordered(diag_file,1,1)
    Open(newunit=lun_diag, file=diag_file)
  EndIf

  ! Read nuclear dataset
  Call read_nuclear_data(data_dir,data_desc)
  Call read_reaction_data(data_dir)

  ! Initialize EoS for screening
  If ( iscrn > 0 ) Call eos_initialize

  ! Initialize NSE
  Call nse_initialize

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

  ! Calculate NSE abundances
  Call nse_solve(rho,t9,ye)

  If ( itsout > 0 ) Then
    Call benuc(ynse,enb,enm,ytot,ztot,atot)
    kstep = 0
    t = 0.0_dp
    tdel = 1.0_dp
    flx = 0.0_dp
    edot = -enm
    Write(lun_ts) kstep, t, t9, rho, tdel, edot, ynse, flx
    Close(lun_ts)
  EndIf

  If ( idiag >= 0 ) Close(lun_diag)

End Program nse_slice
