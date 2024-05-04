MODULE ffn_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nffn = 466
  INTEGER            :: nffn = 0
  INTEGER            :: nlogft =0 
  LOGICAL            :: netweak_flag = .true.
  LOGICAL            :: read_logft_flag = .false.

  INTEGER            :: inuc_ffn(2,max_nffn)
  INTEGER            :: has_logft(max_nffn)=0
  CHARACTER(LEN=5)   :: nname_ffn(2,max_nffn)
  CHARACTER(LEN=4)   :: desc_ffn = ' ffn'
  REAL(8)            :: q_ffn(max_nffn)

  CHARACTER(LEN=256) :: netweak_data_dir   = './ffn_data'

  INTEGER            :: lun_netweak_in, lun_ffngeff_in
  CHARACTER(LEN=256) :: netweak_in_fname   = 'lmpffnoda.data'
  CHARACTER(LEN=256) :: ffngeff_in_fname   = 'ffngeff.dat'

  INTEGER            :: lun_netweak_out
  CHARACTER(LEN=256) :: netweak_out_fname  = 'netweak'

  namelist /ffn_input/ &
    netweak_flag, &
    netweak_data_dir, &
    netweak_in_fname, &
    netweak_out_fname, &
    read_logft_flag , &
    ffngeff_in_fname

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

  SUBROUTINE read_ffngeff(inuc_ffn,ffnenu,ffn_beta,ffn_ft)
    USE net_module, ONLY: nuc_rename, net_index_from_name
    implicit none
    ! In/Out variables
    INTEGER, intent(in) :: inuc_ffn(2,max_nffn)
    REAL(8), intent(out) :: ffn_beta(143,max_nffn), ffn_ft(143,max_nffn)
    REAL(8), intent(inout) :: ffnenu(143,max_nffn)
    ! Local variables
    REAL(8) :: q_read(2), ffnsum_read(2,143), ffnenu_read(2,143)
    REAL(8) :: ffn_beta_read(2,143), ffn_ft_read(2,143)
    INTEGER :: ierr, ii, jj, inucmax,inucmin
    INTEGER :: inuc(2), nffn
    CHARACTER(LEN=5) :: nname_read(2)
    LOGICAL :: keep_rate


    DO
      READ(lun_ffngeff_in,'(14X,A5,18X,F8.4)',IOSTAT=ierr) nname_read(1), q_read(1)
     ! write(*,*) nname_read(1),q_read(1)
      IF ( ierr /= 0 ) THEN
         EXIT
      ENDIF

      READ(lun_ffngeff_in,'(14X,A5,18X,F8.4)',IOSTAT=ierr) nname_read(2), q_read(2)
      READ(lun_ffngeff_in,*)
      READ(lun_ffngeff_in,*)

      ! Make sure all nuclei in the rate are in the network
      keep_rate = .true.
      DO ii = 1, 2

        ! Convert to lower-case and rename some species
        CALL nuc_rename( nname_read(ii) )

        CALL net_index_from_name( nname_read(ii), inuc(ii) )
        IF ( inuc(ii) == 0 ) THEN
          ! write(*,*) "Not keeping rate",nname_read(ii),inuc(ii)
           keep_rate=.false.
        Endif

      END DO

      ! Read the tabulated rate
      DO jj = 1, 143
      READ(lun_ffngeff_in,'(19X,6(f9.3))') ffn_beta_read(1,jj), ffn_ft_read(1,jj), ffnenu_read(1,jj), &
     &        ffn_beta_read(2,jj), ffn_ft_read(2,jj), ffnenu_read(2,jj) 
      END DO

      READ(lun_ffngeff_in,*)

      IF ( keep_rate ) THEN
        DO jj=1,max_nffn
           IF ( inuc_ffn(1,jj) == inuc(1) .and. inuc_ffn(2,jj)==inuc(2) ) THEN
                   nffn=jj
                   EXIT
           ENDIF
        ENDDO
        IF (nffn<0) THEN
                WRITE(*,*) "Reaction in ffngeff but not in lmpffnoda. Skipping."
                CYCLE
        ENDIF

        nlogft = nlogft + 1  ! count accepted rate
        ffn_beta(:,nffn)    = ffn_beta_read(1,:)
        ffn_ft(:,nffn)    = ffn_ft_read(1,:)
        ffnenu(:,nffn)    = ffnenu_read(1,:)
        q_ffn(nffn)       = -q_read(1)!+0.511d0

        has_logft(nffn)=1
   !     write(*,*)"ft-",inuc(1),inuc(2),nffn

        DO jj=1,max_nffn
           IF ( inuc_ffn(1,jj) == inuc(2) .and. inuc_ffn(2,jj)==inuc(1) ) THEN
                   nffn=jj
                   EXIT
           ENDIF
        ENDDO
        IF (nffn<0) THEN
                WRITE(*,*) "Reaction in ffngeff but not in lmpffnoda. Skipping."
                CYCLE
        ENDIF
 
        nlogft = nlogft + 1  ! count accepted rate
        ffn_beta(:,nffn)    = ffn_beta_read(2,:)
        ffn_ft(:,nffn)    = ffn_ft_read(2,:)
        ffnenu(:,nffn)    = ffnenu_read(2,:)       
        q_ffn(nffn)       = q_read(1)!-0.511d0

        has_logft(nffn)=2
   !     write(*,*)"ft+",inuc(2),inuc(1),nffn

      END IF ! keep_rate

    END DO ! read netweak_in

  END SUBROUTINE read_ffngeff

  SUBROUTINE build_netweak
    USE net_module, ONLY: nuc_rename, net_index_from_name
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=5) :: nname_read(2)
    REAL(8) :: q_read(2), ffnsum_read(2,143), ffnenu_read(2,143)
    REAL(8) :: ffn_beta_read(2,143), ffn_ft_read(2,143)

    REAL(8) :: ffnsum(143,max_nffn), ffnenu(143,max_nffn)
    REAL(8) :: ffn_beta(143,max_nffn), ffn_ft(143,max_nffn)
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
    ! Read alternative values from ffngeff
    if (read_logft_flag) then
        CALL read_ffngeff(inuc_ffn,ffnenu,ffn_beta,ffn_ft)
    endif
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
    has_logft(1:nffn)   = has_logft(iffn_sort(1:nffn))
    DO jj = 1,143
      ffnsum(jj,1:nffn) = ffnsum(jj,iffn_sort(1:nffn))
      ffn_beta(jj,1:nffn) = ffn_beta(jj,iffn_sort(1:nffn))
      ffn_ft(jj,1:nffn) = ffn_ft(jj,iffn_sort(1:nffn))
      ffnenu(jj,1:nffn) = ffnenu(jj,iffn_sort(1:nffn))
    END DO

    ! Write the netweak file
    DO ii = 1, nffn
      IF (has_logft(ii)==1 ) THEN
         WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') nname_ffn(1,ii),nname_ffn(2,ii),'ft-',q_ffn(ii)
         WRITE(lun_netweak_out,'(9(f8.3))') (ffn_beta(jj,ii),ffn_ft(jj,ii),ffnsum(jj,ii),ffnenu(jj,ii),jj=1,143)
      ELSE IF (has_logft(ii)==2 ) THEN
         WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') nname_ffn(1,ii),nname_ffn(2,ii),'ft+',q_ffn(ii)
         WRITE(lun_netweak_out,'(9(f8.3))') (ffn_beta(jj,ii),ffn_ft(jj,ii),ffnsum(jj,ii),ffnenu(jj,ii),jj=1,143)
      ELSE
         WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') nname_ffn(1,ii),nname_ffn(2,ii),'ecr',q_ffn(ii)
         WRITE(lun_netweak_out,'(9(f8.3))') (ffnsum(jj,ii),ffnenu(jj,ii),jj=1,143)
      ENDIF
      
    END DO

    RETURN
  END SUBROUTINE build_netweak

END MODULE ffn_module
