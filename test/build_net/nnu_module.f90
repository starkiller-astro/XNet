MODULE nnu_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nnnu = 900000
  INTEGER            :: nnnu
  LOGICAL            :: netneutr_flag = .true.

  INTEGER            :: inuc_nnu(6,max_nnnu)
  INTEGER            :: pem_nnu(3,max_nnnu) ! Particle emission data (n,p,he4)
  CHARACTER(LEN=5)   :: nname_nnu(6,max_nnnu)
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

  SUBROUTINE nnu_index_from_name_particle_em( nname, innu, krate )
    IMPLICIT NONE

    ! Input variables
    CHARACTER(LEN=5), INTENT(IN) :: nname
    INTEGER, INTENT(IN) :: krate

    ! Output variables
    INTEGER, INTENT(OUT) :: innu(27)

    ! Local variables
    CHARACTER(LEN=5) :: cname
    INTEGER :: ii, jj, name_len, pem
   

    cname    = '     '
    name_len = 0

    ! Copy to local
    cname    = TRIM(ADJUSTL(nname))
    name_len = LEN_TRIM(cname)

    
    innu(1:27) = 0
    IF (krate<=3) THEN
            pem = krate-1
    ELSEIF (krate==11) THEN
            pem =3
    ELSE
            jj = 0
            RETURN
    ENDIF
    IF ( name_len > 0 .and. nnnu > 0 ) THEN
      jj = 0
      DO ii = 1, nnnu
        IF ( cname == TRIM(ADJUSTL(nname_nnu(1,ii))) & 
         & .and. (SUM(pem_nnu(:,ii))==pem)  ) THEN
          jj = jj + 1
          innu(jj) = ii
          IF ( jj == 27 ) EXIT
        END IF
      END DO
    END IF

    RETURN
  END SUBROUTINE nnu_index_from_name_particle_em

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
        IF ( SUM(pem_nnu(:,ii)) > 0 ) CYCLE ! Skip reactions with particle emission
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
      nname   = nname_nnu(:,innu)
!      nname(2)   = nname_nnu(2,innu)
!      nname(3:6) = 
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
    CHARACTER(LEN=5) :: nname_read(6)
    CHARACTER(LEN=4) :: desc_read
    REAL(8) :: sigma_read(7)
    CHARACTER(LEN=200) :: dummy

    INTEGER :: inuc(6), pem(3)    
    INTEGER :: ii, jj, ierr, nprod
    LOGICAL :: keep_rate

    nnnu = 0
    IF ( .not. netneutr_flag ) RETURN

    q_nnu = 0.0d0

    ! Extract rates needed by the new network
    DO
      READ(lun_netneutr_in,'(i1,4x,6a5,8x,a4)',IOSTAT=ierr) & 
      &  k_read, nname_read(1), nname_read(2), nname_read(3), nname_read(4), nname_read(5) , nname_read(6), &
      &  desc_read        
!        write(*,*) "[NU] err",ierr
!      write(*,*) "[NU] read: ",k_read, nname_read, desc_read

      IF ( ierr /= 0 ) EXIT

      ! Make sure all nuclei in the rate are in the network
      keep_rate = .true.
      inuc = 0
      pem = 0
      nprod = 0
      DO ii = 1, 6
        IF (nname_read(ii) == '     ')   CYCLE
        ! Check for particle emission
        nprod = nprod +1
        IF (ii > 1) THEN
             IF (nname_read(ii) == '    n')  THEN                      
                     pem(1) = pem(1) + 1
             ELSEIF (nname_read(ii) == '    p') THEN
                     pem(2) = pem(2) + 1
             ELSEIF (nname_read(ii) == '  he4') THEN
                     pem(3) = pem(3) + 1
             ENDIF
!             write(*,*) ii,nname_read(ii),pem
        ENDIF

        ! Convert to lower-case and rename some species
        CALL nuc_rename( nname_read(ii) )

        CALL net_index_from_name( nname_read(ii), inuc(ii) )
        IF ( inuc(ii) == 0 ) keep_rate = .false.
        

      END DO
      IF ( nprod == 2 ) pem = 0
      ! Read the tabulated rate
!            READ(lun_netneutr_in,'(7(f6.2,4x))') (sigma_read(jj),jj=1,7)
      READ(lun_netneutr_in,'(7(es8.2,2x))') (sigma_read(jj),jj=1,7)

      !Do not include 2+ particle emission channels
      if (sum(pem)>2) keep_rate=.false.

      ! Write the rate to the netneutr file
      IF ( keep_rate ) THEN
        nnnu = nnnu + 1
        inuc_nnu(:,nnnu)  = inuc(:)
        nname_nnu(:,nnnu) = nname_read(:)
        desc_nnu(nnnu)    = desc_read
        pem_nnu(:,nnnu) = pem(:)
        WRITE(lun_netneutr_out,'(i1,4x,6a5,8x,a4)',IOSTAT=ierr) k_read, & 
        &   nname_read(1), nname_read(2), nname_read(3), nname_read(4), nname_read(5), nname_read(6), desc_read
!        WRITE(lun_netneutr_out,'(7(f6.2,4x))') (sigma_read(jj),jj=1,7)
        WRITE(lun_netneutr_out,'(4e13.6)') (sigma_read(jj),jj=1,7)
      END IF

    END DO

    RETURN
  END SUBROUTINE build_netneutr

END MODULE nnu_module
