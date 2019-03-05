MODULE nnu_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nnnu = 3307
  INTEGER            :: nnnu
  LOGICAL            :: netneutr_flag = .true.

  INTEGER            :: inuc_nnu(2,max_nnnu)
  CHARACTER(LEN=5)   :: nname_nnu(2,max_nnnu)
  CHARACTER(LEN=4)   :: desc_nnu(max_nnnu)
  REAL(8)            :: q_nnu(max_nnnu)

  CHARACTER(LEN=256) :: netneutr_data_dir  = './neutrino_data'

  INTEGER            :: lun_netneutr_in
  CHARACTER(LEN=256) :: netneutr_in_fname  = 'neutrino.data'

  INTEGER            :: lun_netneutr_out
  CHARACTER(LEN=256) :: netneutr_out_fname = 'netneutr'

  namelist /nnu_input/ &
    netneutr_flag, &
    netneutr_data_dir, &
    netneutr_in_fname, &
    netneutr_out_fname

  CONTAINS

  SUBROUTINE nnu_index_from_name( nname, innu )
    IMPLICIT NONE
    
    ! Input variables
    CHARACTER(LEN=5), INTENT(IN) :: nname

    ! Output variables
    INTEGER, INTENT(OUT) :: innu(2)

    ! Local variables
    CHARACTER(LEN=5) :: cname
    INTEGER :: ii, jj, name_len

    cname    = '     '
    name_len = 0
  
    ! Copy to local
    cname    = TRIM(ADJUSTL(nname))
    name_len = LEN_TRIM(cname)

    innu(1:2) = 0

    IF ( name_len > 0 .and. nnnu > 0 ) THEN
      jj = 0
      DO ii = 1, nnnu
        IF ( cname == TRIM(ADJUSTL(nname_nnu(1,ii))) ) THEN
          jj = jj + 1
          innu(jj) = ii
          IF ( jj == 2 ) EXIT
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE nnu_index_from_name

  SUBROUTINE write_nnu_rate( innu )
    USE net_module, ONLY: write_net_rate, lun_netsu_out
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: innu

    ! Local variables
    INTEGER :: k
    CHARACTER(LEN=5) :: nname(6)
    CHARACTER(LEN=4) :: desc
    CHARACTER(LEN=1) :: rflag, wflag
    REAL(8) :: q, rc(7)

    IF ( innu /= 0 ) THEN
      k          = 0
      nname(1)   = nname_nnu(1,innu)
      nname(2)   = nname_nnu(2,innu)
      nname(3:6) = '     '
      desc       = desc_nnu(innu)
      rflag      = ' '
      wflag      = ' '
      q          = q_nnu(innu)
      rc(1)      = DBLE( innu )
      rc(2:7)    = 0.0d0
      CALL write_net_rate( lun_netsu_out, k, nname, desc, rflag, wflag, q, rc )
    END IF

    RETURN
  END SUBROUTINE write_nnu_rate

  SUBROUTINE build_netneutr
    USE net_module, ONLY: nuc_rename, net_index_from_name
    IMPLICIT NONE

    ! Local variables
    INTEGER :: k_read
    CHARACTER(LEN=5) :: nname_read(2)
    CHARACTER(LEN=4) :: desc_read
    REAL(8) :: sigma_read(7)

    INTEGER :: inuc(2)
    INTEGER :: ii, jj, ierr
    LOGICAL :: keep_rate

    nnnu = 0
    IF ( .not. netneutr_flag ) RETURN

    q_nnu = 0.0d0

    ! Extract rates needed by the new network
    DO
      READ(lun_netneutr_in,'(i1,4x,2a5,28x,a4)',IOSTAT=ierr) k_read, nname_read(1), nname_read(2), desc_read
      IF ( ierr /= 0 ) EXIT

      ! Make sure all nuclei in the rate are in the network
      keep_rate = .true.
      DO ii = 1, 2

        ! Convert to lower-case and rename some species
        CALL nuc_rename( nname_read(ii) )

        CALL net_index_from_name( nname_read(ii), inuc(ii) )
        IF ( inuc(ii) == 0 ) keep_rate = .false.

      END DO

      ! Read the tabulated rate
      READ(lun_netneutr_in,'(7(f6.2,4x))') (sigma_read(jj),jj=1,7)

      ! Write the rate to the netneutr file
      IF ( keep_rate ) THEN
        nnnu = nnnu + 1
        inuc_nnu(1,nnnu)  = inuc(1)
        inuc_nnu(2,nnnu)  = inuc(2)
        nname_nnu(1,nnnu) = nname_read(1)
        nname_nnu(2,nnnu) = nname_read(2)
        desc_nnu(nnnu)    = desc_read
        WRITE(lun_netneutr_out,'(i1,4x,2a5,28x,a4)',IOSTAT=ierr) k_read, nname_read(1), nname_read(2), desc_read
        WRITE(lun_netneutr_out,'(7(f6.2,4x))') (sigma_read(jj),jj=1,7)
      END IF

    END DO

    RETURN
  END SUBROUTINE build_netneutr

END MODULE nnu_module