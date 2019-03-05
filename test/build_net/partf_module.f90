MODULE partf_module
  USE net_module, ONLY: max_nnet
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nwinv = max_nnet
  INTEGER            :: nwinv

  INTEGER            :: inuc_winv(max_nnet)
  CHARACTER(LEN=5)   :: nname_winv(max_nnet)

  CHARACTER(LEN=256) :: netwinv_data_dir = './partf_data'

  INTEGER            :: lun_netwinv_in
! CHARACTER(LEN=14)  :: netwinv_in_fname   = 'winvne_JINAv21'
! CHARACTER(LEN=11)  :: netwinv_in_fname   = 'winvn_jonas'
! CHARACTER(LEN=14)  :: netwinv_in_fname   = 'winvne_JINAv11'
  CHARACTER(LEN=256) :: netwinv_in_fname   = 'winvne_JINAv05'
! CHARACTER(LEN=256) :: netwinv_in_fmt1    = '(a5,f12.3,2i4,f6.1,f10.3)'
! CHARACTER(LEN=256) :: netwinv_in_fmt2    = '(8f9.2)'
  CHARACTER(LEN=256) :: netwinv_in_fmt1    = '(a5,f12.3,2i4,f6.1,f10.3,1x,a11)'
  CHARACTER(LEN=256) :: netwinv_in_fmt2    = '(8es12.5)'

  INTEGER            :: lun_netwinv_out
  CHARACTER(LEN=256) :: netwinv_out_fname  = 'netwinv'
! CHARACTER(LEN=256) :: netwinv_out_fmt1   = '(a5,f12.3,2i4,f6.1,f10.3)'
! CHARACTER(LEN=256) :: netwinv_out_fmt2   = '(8f9.2)'
  CHARACTER(LEN=256) :: netwinv_out_fmt1   = '(a5,f12.3,2i4,f6.1,f15.8)'
  CHARACTER(LEN=256) :: netwinv_out_fmt2   = '(8es12.5)'

  CHARACTER(LEN=256) :: mass_data_dir = './mass_data'

  INTEGER            :: lun_ame03
  CHARACTER(LEN=256) :: ame03_fname = 'mass_ame03.dat'

  INTEGER            :: lun_ame03extrap
  CHARACTER(LEN=256) :: ame03extrap_fname = 'mass_ame03extrap.dat'

  INTEGER            :: lun_ame11
  CHARACTER(LEN=256) :: ame11_fname = 'mass_ame11.dat'

  INTEGER            :: lun_reac1
  CHARACTER(LEN=256) :: reac1_fname = 'mass_reac1.dat'

  INTEGER            :: lun_ame11extrap
  CHARACTER(LEN=256) :: ame11extrap_fname = 'mass_ame11extrap.dat'

  INTEGER            :: lun_frdm
  CHARACTER(LEN=256) :: frdm_fname = 'mass_frdm.dat'

  INTEGER, PARAMETER :: max_iz = 136
  INTEGER, PARAMETER :: max_in = 236

  namelist /partf_input/ &
    netwinv_data_dir, &
    netwinv_in_fname, &
    mass_data_dir, &
    ame03_fname, &
    ame03extrap_fname, &
    ame11_fname, &
    ame11extrap_fname, &
    reac1_fname, &
    frdm_fname, &
    netwinv_out_fname


  CONTAINS

  SUBROUTINE build_netwinv
    USE net_module, ONLY: nnet, nname_net
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=72) :: line_read
    CHARACTER(LEN=5)  :: nname_read
    CHARACTER(LEN=11) :: source_read
    REAL(8) :: aa_read, sp_read, mex_read, g_read(24)
    INTEGER :: iz_read, in_read

    REAL(8) :: ame03_mex(0:max_iz,0:max_in)
    REAL(8) :: ame03extrap_mex(0:max_iz,0:max_in)
    REAL(8) :: ame11_mex(0:max_iz,0:max_in)
    REAL(8) :: reac1_mex(0:max_iz,0:max_in)
    REAL(8) :: ame11extrap_mex(0:max_iz,0:max_in)
    REAL(8) :: frdm_mex(0:max_iz,0:max_in)

    INTEGER :: inuc, jnuc, knuc
    INTEGER :: k, m, ierr

    ! Read mass data
    CALL read_mass( lun_ame03, ame03_mex )
    CALL read_mass( lun_ame03extrap, ame03extrap_mex )
    CALL read_mass( lun_ame11, ame11_mex )
    CALL read_mass( lun_reac1, reac1_mex )
    CALL read_mass( lun_ame11extrap, ame11extrap_mex )
    CALL read_mass( lun_frdm, frdm_mex )

    ! Read old header
    READ(lun_netwinv_in,'(i5)') nwinv
    IF ( nwinv > max_nwinv ) THEN
      WRITE(*,'(a)') 'ERROR: '//TRIM(ADJUSTL(netwinv_in_fname))//' contains too mannet species'
      STOP
    END IF
    READ(lun_netwinv_in,'(a72)') line_read

    ! Write new header
    WRITE(lun_netwinv_out,'(i5)') nnet
    WRITE(lun_netwinv_out,'(a72)') line_read

    jnuc = 1
    DO inuc = 1, max_nwinv
      READ(lun_netwinv_in,'(a5)') nname_read
      nname_winv(inuc) = ADJUSTR( nname_read )
      IF ( inuc > 1 ) THEN
        IF ( nname_winv(inuc) == nname_winv(inuc-1) ) EXIT
      END IF

      IF ( nname_winv(inuc) == nname_net(jnuc) ) THEN
        WRITE(lun_netwinv_out,'(a5)') nname_winv(inuc)
        inuc_winv(jnuc) = inuc
        jnuc = jnuc + 1
      END IF
    END DO
    
    inuc = 1
    DO jnuc = 1, nnet
      knuc = inuc_winv(jnuc) - inuc + 1
      DO k = 1, knuc
        READ(lun_netwinv_in,netwinv_in_fmt1,IOSTAT=ierr) nname_read,aa_read,iz_read,in_read,sp_read,mex_read,source_read
        READ(lun_netwinv_in,*) (g_read(m),m=1,24)
      END DO

      IF ( TRIM(ADJUSTL(source_read)) == 'ame11' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,ame11_mex(iz_read,in_read)
      ELSE IF ( TRIM(ADJUSTL(source_read)) == 'reac1' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,reac1_mex(iz_read,in_read)
      ELSE IF ( TRIM(ADJUSTL(source_read)) == 'ame11extrap' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,ame11extrap_mex(iz_read,in_read)
      ELSE IF ( TRIM(ADJUSTL(source_read)) == 'ame03extrap' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,ame03extrap_mex(iz_read,in_read)
      ELSE IF ( TRIM(ADJUSTL(source_read)) == 'ame03' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,ame03_mex(iz_read,in_read)
      ELSE IF ( TRIM(ADJUSTL(source_read)) == 'frdm' ) THEN
        WRITE(lun_netwinv_out,'(a5,f12.3,2i4,f6.1,f15.8)') &
        & nname_read,aa_read,iz_read,in_read,sp_read,frdm_mex(iz_read,in_read)
      ELSE
        WRITE(lun_netwinv_out,netwinv_out_fmt1) &
        & nname_read,aa_read,iz_read,in_read,sp_read,mex_read
      END IF
      WRITE(lun_netwinv_out,netwinv_out_fmt2) (g_read(m),m=1,24)
      inuc = inuc_winv(jnuc) + 1
    END DO
      
    RETURN
  END SUBROUTINE build_netwinv

  SUBROUTINE read_mass( lun_mass, mex )
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: lun_mass

    ! Output variables
    REAL(8), INTENT(OUT) :: mex(0:max_iz,0:max_in)

    ! Local variables
    INTEGER :: iz_read, ia_read
    CHARACTER(LEN=15) :: eval_read
    REAL(8) :: mex_read

    INTEGER :: i, ierr

    ! Initialize mass excesses to unrealistically large value
    mex(:,:) = 1.0d99

    ! Skip header
    DO i = 1, 15
      READ(lun_mass,*) 
    END DO

    DO 
!     READ(lun_mass,'(2(1x,i3),1x,a15,1x,f12.5)',IOSTAT=ierr) &
      READ(lun_mass,*,IOSTAT=ierr) iz_read, ia_read, eval_read, mex_read
      IF ( ierr < 0 ) EXIT

      mex(iz_read,ia_read-iz_read) = mex_read * 1.0d-3
    END DO

    RETURN
  END SUBROUTINE read_mass

END MODULE partf_module