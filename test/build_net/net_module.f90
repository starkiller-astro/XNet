MODULE net_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nnet = 7855
  INTEGER            :: nnet

  CHARACTER(LEN=5)   :: nname_net(max_nnet)

  INTEGER            :: lun_sunet_in
  CHARACTER(LEN=256) :: sunet_fname        = 'sunet'

  CHARACTER(LEN=256) :: netsu_data_dir     = './reaclib_data'

  INTEGER            :: lun_netsu_in
  CHARACTER(LEN=256) :: netsu_in_fname     = 'reaclib_JINAv22_ks02'
! CHARACTER(LEN=256) :: netsu_in_fname     = 'reaclib_JINAv21'
! CHARACTER(LEN=256) :: netsu_in_fname     = 'reaclib_jonas'
! CHARACTER(LEN=256) :: netsu_in_fname     = 'reaclib_JINAv11'
! CHARACTER(LEN=256) :: netsu_in_fname     = 'reaclib_JINAv05'
  INTEGER            :: reaclib_ver        = 1
  LOGICAL            :: no910              = .true.

  INTEGER            :: lun_sunet_out

  INTEGER            :: lun_netsu_out
  CHARACTER(LEN=5)   :: netsu_out_fname    = 'netsu'

  namelist /net_input/ &
    sunet_fname, &
    netsu_data_dir, &
    netsu_in_fname, &
    reaclib_ver, &
    no910, &
    netsu_out_fname

  CONTAINS

  SUBROUTINE nuc_rename( nname )
    IMPLICIT NONE

    ! Input variables
    CHARACTER(LEN=5), INTENT(INOUT) :: nname

    ! Local variables
    INTEGER, PARAMETER :: lc_a_ascii=IACHAR('a')
    INTEGER, PARAMETER :: uc_a_ascii=IACHAR('A')
    INTEGER, PARAMETER :: lc_z_ascii=IACHAR('z')
    INTEGER, PARAMETER :: uc_z_ascii=IACHAR('Z')

    CHARACTER(LEN=5) :: amass, letters, cname
    INTEGER :: name_len, in_ascii, jj, kk
    LOGICAL :: digit

    amass    = '     '
    letters  = '     '
    cname    = '     '

    ! Copy to local
    cname    = TRIM(ADJUSTL(nname))

    ! Determine length of name
    name_len = LEN_TRIM(cname)

    ! Take care of special cases
    IF ( cname(1:name_len) == 'H1' .or. cname(1:name_len) == 'h1' .or. &
    &    cname(1:name_len) == 'H'  .or. cname(1:name_len) == 'h' ) THEN
      nname = 'p'
      name_len = 1
    ELSE IF ( cname(1:name_len) == 'H2' .or. cname(1:name_len) == 'h2' ) THEN
      nname = 'd'
      name_len = 1
    ELSE IF ( cname(1:name_len) == 'H3' .or. cname(1:name_len) == 'h3' ) THEN
      nname = 't'
      name_len = 1
    ELSE IF ( cname(1:name_len) == 'N01' .or. cname(1:name_len) == 'N' .or. &
    &         cname(1:name_len) == 'nt1' .or. cname(1:name_len) == 'neut' .or. &
    &         cname(1:name_len) == 'n' ) THEN
      nname = 'n'
      name_len = 1
    ELSE IF ( cname(1:name_len) == 'Al-6' ) THEN
      nname = 'al-6'
      name_len = 4
    ELSE IF ( cname(1:name_len) == 'Al*6' ) THEN
      nname = 'al*6'
      name_len = 4
    ELSE IF ( name_len > 1 ) THEN

      ! Scan name until you find the digit
      digit = .false.
      DO jj = 1, name_len
        SELECT CASE( cname(jj:jj) )
        CASE( '0':'9' )
          digit = .true.
        CASE DEFAULT
          digit = .false.
        END SELECT

        IF ( digit ) THEN
          letters(1:(jj-1)) = cname(1:(jj-1)) ! pick out letters
          letters           = TRIM(ADJUSTL(letters))

          amass = cname(jj:name_len) ! pick out digits
          amass = TRIM(ADJUSTL(amass))

          ! Convert name to lower-case
          DO kk = 1, LEN(letters)
            in_ascii = IACHAR( letters(kk:kk) )
            IF ( in_ascii >= uc_a_ascii .AND. in_ascii <= uc_z_ascii ) THEN
              letters(kk:kk) = ACHAR( in_ascii + (lc_a_ascii - uc_a_ascii) )
            END IF
          END DO

          ! Copy back to species name array
          WRITE(nname,'(a5)') TRIM(ADJUSTL(letters))//TRIM(ADJUSTL(amass))

          EXIT

        END IF ! digit

      END DO ! jj = 1, name_len
    END IF ! name_len > 1

    nname = ADJUSTR( nname )

    RETURN
  END SUBROUTINE nuc_rename

  SUBROUTINE net_index_from_name( nname, inuc )
    IMPLICIT NONE
    
    ! Input variables
    CHARACTER(LEN=5), INTENT(IN) :: nname

    ! Output variables
    INTEGER, INTENT(OUT) :: inuc

    ! Local variables
    CHARACTER(LEN=5) :: cname
    INTEGER :: ii, name_len

    cname    = '     '
    name_len = 0
  
    ! Copy to local
    cname    = TRIM(ADJUSTL(nname))
    name_len = LEN_TRIM(cname)

    inuc = 0
    IF ( name_len > 0 ) THEN
      DO ii = 1, nnet
        IF ( cname == TRIM(ADJUSTL(nname_net(ii))) ) THEN
          inuc = ii
          EXIT
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE net_index_from_name

  SUBROUTINE write_net_rate( lun_out, k, nname, desc, rflag, wflag, q, rc )
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: lun_out
    INTEGER, INTENT(IN) :: k
    CHARACTER(LEN=5), INTENT(IN) :: nname(6)
    CHARACTER(LEN=4), INTENT(IN) :: desc
    CHARACTER(LEN=1), INTENT(IN) :: rflag, wflag
    REAL(8), INTENT(IN) :: q, rc(7)

    ! Local variables
    INTEGER :: jj

    IF ( k < 10 ) THEN
      WRITE(lun_out,'(i1,4x,6a5,8x,a4,a1,a1,3x,1pe12.5)') k,(nname(jj),jj=1,6),desc,rflag,wflag,q
    ELSE
      WRITE(lun_out,'(i2,3x,6a5,8x,a4,a1,a1,3x,1pe12.5)') k,(nname(jj),jj=1,6),desc,rflag,wflag,q
    END IF
    WRITE(lun_out,'(4e13.6)') (rc(jj),jj=1,7)

    RETURN
  END SUBROUTINE write_net_rate

  SUBROUTINE read_sunet
    IMPLICIT NONE

    ! Local variables
    INTEGER :: inuc, ierr

    inuc = 1
    DO
      READ(lun_sunet_in,'(a5)',IOSTAT=ierr) nname_net(inuc)
      IF ( ierr /= 0 ) EXIT
      inuc = inuc + 1
      IF ( inuc > max_nnet ) THEN
        WRITE(*,'(a)') 'ERROR: sunet contains too many species'
        STOP
      END IF
    END DO
    nnet = inuc - 1
!   WRITE(*,'(15a5)') (nname_net(inuc),inuc=1,nnet)
    WRITE(*,'(a,i5)') '# species: ',nnet

    RETURN
  END SUBROUTINE read_sunet

  SUBROUTINE write_sunet
    IMPLICIT NONE

    ! Local variables
    INTEGER :: i, ierr

    DO i = 1, nnet
      WRITE(lun_sunet_out,'(a5)',IOSTAT=ierr) ADJUSTR(nname_net(i))
    END DO

    RETURN
  END SUBROUTINE write_sunet

END MODULE net_module