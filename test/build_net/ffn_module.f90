MODULE ffn_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nffn = 466
  INTEGER            :: nffn
  LOGICAL            :: netweak_flag = .true.

  INTEGER            :: inuc_ffn(2,max_nffn)
  CHARACTER(LEN=5)   :: nname_ffn(2,max_nffn)
  CHARACTER(LEN=4)   :: desc_ffn = ' ffn'
  REAL(8)            :: q_ffn(max_nffn)

  CHARACTER(LEN=256) :: netweak_data_dir   = './ffn_data'

  INTEGER            :: lun_netweak_in
  CHARACTER(LEN=256) :: netweak_in_fname   = 'lmpffnoda.data'

  INTEGER            :: lun_netweak_out
  CHARACTER(LEN=256) :: netweak_out_fname  = 'netweak'

  namelist /ffn_input/ &
    netweak_flag, &
    netweak_data_dir, &
    netweak_in_fname, &
    netweak_out_fname

  CONTAINS

  SUBROUTINE ffn_index_from_name( nname, iffn )
    IMPLICIT NONE
    
    ! Input variables
    CHARACTER(LEN=5), INTENT(IN) :: nname

    ! Output variables
    INTEGER, INTENT(OUT) :: iffn(2)

    ! Local variables
    CHARACTER(LEN=5) :: cname
    INTEGER :: ii, jj, name_len

    cname    = '     '
    name_len = 0
  
    ! Copy to local
    cname    = TRIM(ADJUSTL(nname))
    name_len = LEN_TRIM(cname)

    iffn(1:2) = 0

    IF ( name_len > 0 .and. nffn > 0 ) THEN
      jj = 0
      DO ii = 1, nffn
        IF ( cname == TRIM(ADJUSTL(nname_ffn(1,ii))) ) THEN
          jj = jj + 1
          iffn(jj) = ii
          IF ( jj == 2 ) EXIT
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE ffn_index_from_name

  SUBROUTINE write_ffn_rate( iffn )
    USE net_module, ONLY: write_net_rate, lun_netsu_out
    IMPLICIT NONE

    ! Input variables
    INTEGER, INTENT(IN) :: iffn

    ! Local variables
    INTEGER :: k
    CHARACTER(LEN=5) :: nname(6)
    CHARACTER(LEN=4) :: desc
    CHARACTER(LEN=1) :: rflag, wflag
    REAL(8) :: q, rc(7)

    IF ( iffn /= 0 ) THEN
      k          = 0
      nname(1)   = nname_ffn(1,iffn)
      nname(2)   = nname_ffn(2,iffn)
      nname(3:6) = '     '
      desc       = desc_ffn
      rflag      = ' '
      wflag      = ' '
      q          = q_ffn(iffn)
      rc(1)      = DBLE( iffn )
      rc(2:7)    = 0.0d0
      CALL write_net_rate( lun_netsu_out, k, nname, desc, rflag, wflag, q, rc )
    END IF

    RETURN
  END SUBROUTINE write_ffn_rate

  SUBROUTINE build_netweak
    USE net_module, ONLY: nuc_rename, net_index_from_name
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=5) :: nname_read(2)
    REAL(8) :: q_read(2), ffnsum_read(2,143), ffnenu_read(2,143)

    REAL(8) :: ffnsum(143,max_nffn), ffnenu(143,max_nffn)
    INTEGER :: inuc(2), iffn(2), iffn_sort(max_nffn)
    INTEGER :: ii, jj, inucmin, inucmax, ierr
    LOGICAL :: keep_rate, is_sorted

    nffn = 0
    IF ( .not. netweak_flag ) RETURN

    ! Extract rates needed by the new network
    DO
      READ(lun_netweak_in,'(14x,a5,25x,f8.4)',IOSTAT=ierr) nname_read(1), q_read(1)
      IF ( ierr /= 0 ) EXIT

      READ(lun_netweak_in,'(14x,a5,25x,f8.4)') nname_read(2), q_read(2)
      READ(lun_netweak_in,*)
      READ(lun_netweak_in,*)

      ! Make sure all nuclei in the rate are in the network
      keep_rate = .true.
      DO ii = 1, 2

        ! Convert to lower-case and rename some species
        CALL nuc_rename( nname_read(ii) )

        CALL net_index_from_name( nname_read(ii), inuc(ii) )
        IF ( inuc(ii) == 0 ) keep_rate = .false.

      END DO

      ! Read the tabulated rate
      DO jj = 1, 143
        READ(lun_netweak_in,'(40x,2(f8.3),18x,2(f9.3))') ffnsum_read(1,jj),ffnenu_read(1,jj),ffnsum_read(2,jj),ffnenu_read(2,jj)
      END DO
      READ(lun_netweak_in,*)

      IF ( keep_rate ) THEN

        ! Pre-sort the rate by the parent nucleus
        inucmin = MINLOC(inuc,1)
        inucmax = MAXLOC(inuc,1)

        nffn = nffn + 1
        inuc_ffn(1,nffn)  = inuc(inucmin)
        inuc_ffn(2,nffn)  = inuc(inucmax)
        nname_ffn(1,nffn) = nname_read(inucmin)
        nname_ffn(2,nffn) = nname_read(inucmax)
        ffnsum(:,nffn)    = ffnsum_read(inucmin,:)
        ffnenu(:,nffn)    = ffnenu_read(inucmin,:)
        q_ffn(nffn)       = q_read(inucmax)

        nffn = nffn + 1
        inuc_ffn(1,nffn)  = inuc(inucmax)
        inuc_ffn(2,nffn)  = inuc(inucmin)
        nname_ffn(1,nffn) = nname_read(inucmax)
        nname_ffn(2,nffn) = nname_read(inucmin)
        ffnsum(:,nffn)    = ffnsum_read(inucmax,:)
        ffnenu(:,nffn)    = ffnenu_read(inucmax,:)
        q_ffn(nffn)       = -q_read(inucmax)
      END IF ! keep_rate

    END DO

    ! Generate a sorted index vector for sorting the reactions
    DO ii = 1, nffn
      iffn_sort(ii) = ii
    END DO

    ! First, sort on reactants
    is_sorted = .false.
    DO WHILE ( .not. is_sorted )
      is_sorted = .true.
      DO ii = 2, nffn
        iffn(1) = iffn_sort(ii-1)
        iffn(2) = iffn_sort(ii)
        IF ( inuc_ffn(1,iffn(2)) < inuc_ffn(1,iffn(1)) ) THEN
          iffn_sort(ii-1) = iffn(2)
          iffn_sort(ii)   = iffn(1)
          is_sorted       = .false.
        END IF
      END DO
    END DO

    ! Next, sort on products
    is_sorted = .false.
    DO WHILE ( .not. is_sorted )
      is_sorted = .true.
      DO ii = 2, nffn
        iffn(1) = iffn_sort(ii-1)
        iffn(2) = iffn_sort(ii)
        IF ( inuc_ffn(1,iffn(2)) == inuc_ffn(1,iffn(1)) .and. &
        &    inuc_ffn(2,iffn(2)) <  inuc_ffn(2,iffn(1)) ) THEN
          iffn_sort(ii-1) = iffn(2)
          iffn_sort(ii)   = iffn(1)
          is_sorted       = .false.
        END IF
      END DO
    END DO

    ! Apply the sorted index vector to the arrays
    inuc_ffn(1,1:nffn)  = inuc_ffn(1,iffn_sort(1:nffn))
    inuc_ffn(2,1:nffn)  = inuc_ffn(2,iffn_sort(1:nffn))
    nname_ffn(1,1:nffn) = nname_ffn(1,iffn_sort(1:nffn))
    nname_ffn(2,1:nffn) = nname_ffn(2,iffn_sort(1:nffn))
    q_ffn(1:nffn)       = q_ffn(iffn_sort(1:nffn))
    DO jj = 1,143
      ffnsum(jj,1:nffn) = ffnsum(jj,iffn_sort(1:nffn))
      ffnenu(jj,1:nffn) = ffnenu(jj,iffn_sort(1:nffn))
    END DO

    ! Write the netweak file
    DO ii = 1, nffn
      WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') nname_ffn(1,ii),nname_ffn(2,ii),'ecr',q_ffn(ii)
      WRITE(lun_netweak_out,'(9(f8.3))') (ffnsum(jj,ii),ffnenu(jj,ii),jj=1,143)
    END DO

    RETURN
  END SUBROUTINE build_netweak

END MODULE ffn_module