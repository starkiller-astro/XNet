MODULE file_module
  IMPLICIT NONE

  INTEGER            :: lun_input
  CHARACTER(LEN=256) :: input_fname  = 'input.namelist'

  CHARACTER(LEN=256) :: raw_data_dir = './raw_data'
  CHARACTER(LEN=256) :: new_data_dir = './new_data'

  namelist /file_input/ &
    raw_data_dir, &
    new_data_dir

  CONTAINS

  SUBROUTINE build_net_error( message )
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: message

    WRITE(*,'(2a)') 'ERROR: ',TRIM(message)
    STOP 1
  END SUBROUTINE build_net_error

  SUBROUTINE safe_open_old( lun, dir, fname, ierr )
    IMPLICIT NONE

    ! Input variables
    CHARACTER(LEN=*), INTENT(IN) :: dir, fname

    ! Output variables
    INTEGER, INTENT(OUT) :: lun, ierr

    ! Local variables
    CHARACTER(LEN=130) :: tmp_fname

    tmp_fname = TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fname))
    OPEN( NEWUNIT=lun, FILE=TRIM(tmp_fname), STATUS='OLD', IOSTAT=ierr )
    IF ( ierr /= 0 ) WRITE(*,'(a,i3,2a)') 'ERROR(',ierr,'): Cannot open ',tmp_fname

    RETURN
  END SUBROUTINE safe_open_old

  SUBROUTINE safe_open_new( lun, dir, fname, ierr )
    IMPLICIT NONE

    ! Input variables
    CHARACTER(LEN=*), INTENT(IN) :: dir, fname

    ! Output variables
    INTEGER, INTENT(OUT) :: lun, ierr

    ! Local variables
    CHARACTER(LEN=130) :: tmp_fname

    tmp_fname = TRIM(ADJUSTL(dir))//'/'//TRIM(ADJUSTL(fname))
    OPEN( NEWUNIT=lun, FILE=TRIM(tmp_fname), STATUS='NEW', IOSTAT=ierr )
    IF ( ierr /= 0 ) THEN
      OPEN( NEWUNIT=lun, FILE=TRIM(tmp_fname), STATUS='REPLACE', IOSTAT=ierr )
      IF ( ierr /= 0 ) WRITE(*,'(a,i3,2a)') 'ERROR(',ierr,'): Cannot open ',tmp_fname
    END IF

    RETURN
  END SUBROUTINE safe_open_new

  SUBROUTINE file_init
    USE net_module, ONLY: lun_sunet_in, lun_sunet_out, lun_netsu_in, lun_netsu_out, &
    & sunet_fname, netsu_in_fname, netsu_out_fname, netsu_data_dir, net_input, reaclib_ver, no910
    USE partf_module, ONLY: lun_netwinv_in, lun_netwinv_out, &
    & netwinv_in_fname, netwinv_out_fname, lun_ame11, lun_reac1, &
    & lun_ame11extrap, lun_frdm, lun_ame03, lun_ame03extrap, &
    & ame11_fname, reac1_fname, ame11extrap_fname, frdm_fname, &
    & ame03_fname, ame03extrap_fname, netwinv_data_dir, mass_data_dir, partf_input
    USE ffn_module, ONLY: lun_netweak_in, lun_netweak_out, &
    & netweak_in_fname, netweak_out_fname, netweak_flag, netweak_data_dir, ffn_input, &
    & element_list_fname, lun_element_in
    USE nnu_module, ONLY: lun_netneutr_in, lun_netneutr_out, &
    & netneutr_in_fname, netneutr_out_fname, netneutr_flag, netneutr_data_dir, nnu_input
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=256) :: iomsg
    INTEGER :: ierr

    CALL safe_open_old( lun_input, '.', input_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Cannot read '//TRIM(input_fname) )

    READ( lun_input, NML=file_input, IOSTAT=ierr, IOMSG=iomsg )
    IF ( ierr /= 0 ) CALL build_net_error( 'Malformed file_input namelist: '//TRIM(iomsg) )
    READ( lun_input, NML=net_input, IOSTAT=ierr, IOMSG=iomsg )
    IF ( ierr /= 0 ) CALL build_net_error( 'Malformed net_input namelist: '//TRIM(iomsg) )
    READ( lun_input, NML=ffn_input, IOSTAT=ierr, IOMSG=iomsg )
    IF ( ierr /= 0 ) CALL build_net_error( 'Malformed ffn_input namelist: '//TRIM(iomsg) )
    READ( lun_input, NML=partf_input, IOSTAT=ierr, IOMSG=iomsg )
    IF ( ierr /= 0 ) CALL build_net_error( 'Malformed partf_input namelist: '//TRIM(iomsg) )
    READ( lun_input, NML=nnu_input, IOSTAT=ierr, IOMSG=iomsg )
    IF ( ierr /= 0 ) CALL build_net_error( 'Malformed nnu_input namelist: '//TRIM(iomsg) )
    CLOSE( lun_input )

    ! Open raw data files
    CALL safe_open_old( lun_sunet_in, '.', sunet_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing requested-species file' )
    
    CALL safe_open_old( lun_netsu_in, netsu_data_dir, netsu_in_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing REACLIB source' )

    CALL safe_open_old( lun_netwinv_in, netwinv_data_dir, netwinv_in_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing partition-function source' )

    CALL safe_open_old( lun_ame03, mass_data_dir, ame03_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing AME03 mass source' )

    CALL safe_open_old( lun_ame03extrap, mass_data_dir, ame03extrap_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing extrapolated AME03 mass source' )

    CALL safe_open_old( lun_ame11, mass_data_dir, ame11_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing AME11 mass source' )

    CALL safe_open_old( lun_reac1, mass_data_dir, reac1_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing REAC1 mass source' )

    CALL safe_open_old( lun_ame11extrap, mass_data_dir, ame11extrap_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing extrapolated AME11 mass source' )

    CALL safe_open_old( lun_frdm, mass_data_dir, frdm_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Missing FRDM mass source' )

    IF ( netweak_flag ) THEN
      CALL safe_open_old( lun_netweak_in, netweak_data_dir, netweak_in_fname, ierr )
      IF ( ierr /= 0 ) CALL build_net_error( 'Missing enabled weak-rate source' )
    END IF
    
    IF ( netweak_flag ) THEN
      CALL safe_open_old( lun_element_in, netweak_data_dir, element_list_fname, ierr)
      IF ( ierr /= 0 ) CALL build_net_error( 'Missing enabled weak-rate element list' )
    END IF

    IF ( netneutr_flag ) THEN
      CALL safe_open_old( lun_netneutr_in, netneutr_data_dir, netneutr_in_fname, ierr )
      IF ( ierr /= 0 ) netneutr_flag = .false.
    END IF

    ! Open output files for new network
    CALL safe_open_new( lun_sunet_out, new_data_dir, 'sunet', ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Cannot create sunet' )

    CALL safe_open_new( lun_netsu_out, new_data_dir, netsu_out_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Cannot create netsu' )

    CALL safe_open_new( lun_netwinv_out, new_data_dir, netwinv_out_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Cannot create netwinv' )

    CALL safe_open_new( lun_netweak_out, new_data_dir, netweak_out_fname, ierr )
    IF ( ierr /= 0 ) CALL build_net_error( 'Cannot create netweak' )

    IF ( netneutr_flag ) THEN
      CALL safe_open_new( lun_netneutr_out, new_data_dir, netneutr_out_fname, ierr )
      IF ( ierr /= 0 ) CALL build_net_error( 'Cannot create netneutr' )
    END IF

    RETURN
  END SUBROUTINE file_init

  SUBROUTINE file_finalize
    USE net_module, ONLY: lun_sunet_in, lun_sunet_out, lun_netsu_in, lun_netsu_out
    USE partf_module, ONLY: lun_netwinv_in, lun_netwinv_out, lun_ame11, &
    & lun_reac1, lun_ame11extrap, lun_frdm, lun_ame03, lun_ame03extrap
    USE ffn_module, ONLY: lun_netweak_in, lun_netweak_out, netweak_flag, lun_element_in
    USE nnu_module, ONLY: lun_netneutr_in, lun_netneutr_out, netneutr_flag
    IMPLICIT NONE

    CLOSE( lun_sunet_in )
    CLOSE( lun_sunet_out )
    CLOSE( lun_netsu_in )
    CLOSE( lun_netwinv_in )
    CLOSE( lun_ame03 )
    CLOSE( lun_ame03extrap )
    CLOSE( lun_ame11 )
    CLOSE( lun_reac1 )
    CLOSE( lun_ame11extrap )
    CLOSE( lun_frdm )
    CLOSE( lun_netsu_out )
    CLOSE( lun_netwinv_out )
    CLOSE( lun_netweak_out )
    IF ( netweak_flag ) CLOSE( lun_netweak_in )
    IF ( netweak_flag ) CLOSE( lun_element_in )
    IF ( netneutr_flag ) CLOSE( lun_netneutr_in )
    IF ( netneutr_flag ) CLOSE( lun_netneutr_out )

    RETURN
  END SUBROUTINE file_finalize

END MODULE file_module
