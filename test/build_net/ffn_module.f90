MODULE ffn_module
  IMPLICIT NONE

  INTEGER, PARAMETER :: max_nffn = 1600
  INTEGER            :: nffn
  LOGICAL            :: netweak_flag = .true.

  INTEGER            :: inuc_ffn(2,max_nffn)
  CHARACTER(LEN=5)   :: nname_ffn(2,max_nffn)
  CHARACTER(LEN=4)   :: desc_ffn = ' ffn'
  REAL(8)            :: q_ffn(max_nffn)

  CHARACTER(LEN=256) :: netweak_data_dir   = './ffn_data'

  INTEGER            :: lun_netweak_in
  CHARACTER(LEN=256) :: netweak_in_fname   = 'updated_rate_table.txt'

  INTEGER            :: lun_netweak_out
  CHARACTER(LEN=256) :: netweak_out_fname  = 'netweak'

  INTEGER            :: lun_sunet_in
  CHARACTER(LEN=256) :: sunet_fname = 'sunet.sn160'

  INTEGER            :: lun_element_in
  CHARACTER(LEN=256) :: element_list_fname = 'element_list.txt'
  
!  INTEGER            :: lun_netwinv_in
!  CHARACTER(LEN=256) :: netwinv_in_fname = 'winvne_JINAv22'

  namelist /ffn_input/ &
    netweak_flag, &
    netweak_data_dir, &
    netweak_in_fname, &
    netweak_out_fname

  namelist /net_input/ &
    sunet_fname

!  namelist /partf_input/ &
!    netwinv_in_fname

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
  
  SUBROUTINE check_network(a_in,z_in,e_list,spec_list,num_spec,rate)
        IMPLICIT NONE
        INTEGER             :: z_in,num_spec
        CHARACTER(LEN=3)    :: a_in
        CHARACTER(LEN=*)    :: e_list(119)
        CHARACTER(LEN=2)    ::element, element_min1
        CHARACTER(LEN=*)    :: spec_list(num_spec)
        CHARACTER(LEN=5)    :: spec1, spec2
        LOGICAL             :: rate

        IF ( a_in == '1  ' .and. z_in == 1 ) THEN
                spec1 = '    p'
                spec2 = '    n'
        !this prevents the code from trying to read data for elements past z=118
        ELSEIF ( z_in > 118 ) THEN
                spec1 = 'no777'
                spec2 = 'no777'
        ELSE
                element = e_list(z_in+1)
                element_min1 = e_list(z_in)
                spec1 = adjustr(element//a_in)
                spec2 = adjustr(element_min1//a_in)
        END IF

        IF (ANY(spec_list == spec1) .AND. ANY(spec_list == spec2)) THEN
                rate = .true.
        ELSE
                rate = .false.
        END IF
  END SUBROUTINE check_network

  SUBROUTINE build_element_list(e_in, e_file_num, zz, e_list)
        IMPLICIT NONE
        CHARACTER(LEN=256)    :: e_in
        INTEGER               :: zz(119)
        CHARACTER(LEN=*)      :: e_list(119)
        INTEGER               :: z_val,ii,e_file_num
        CHARACTER(LEN=2)      :: element

        DO ii = 1,119
                READ(e_file_num,'(i3,1x,a2)') z_val, element
                zz(ii) = z_val
                element = adjustr(element)
                e_list(ii) = element
        END DO

  END SUBROUTINE build_element_list

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
  
  SUBROUTINE calculate_q_value( spec1, spec2, q_value )
    ! calculates the Q-value of a reaction in units of MeV
    USE partf_module, ONLY: mex_list
    USE net_module, ONLY: net_index_from_name
    IMPLICIT NONE
    CHARACTER(5)         :: spec1, spec2
    INTEGER              :: s1_index, s2_index
    REAL(8)              :: mex1, mex2
    REAL(8), INTENT(OUT) ::q_value
    
    !Determine mass excess of each species
    call net_index_from_name( spec1, s1_index )
    call net_index_from_name( spec2, s2_index )
    mex1 = mex_list(s1_index)
    mex2 = mex_list(s2_index)

    q_value = mex1 - mex2
    
  END SUBROUTINE calculate_q_value

  SUBROUTINE build_netweak
    USE net_module, ONLY: nuc_rename, net_index_from_name, nname_net, nnet
    IMPLICIT NONE

    ! Local variables
    CHARACTER(LEN=5) :: nname_read(2)
    REAL(8) :: q_read(2), ffnsum_read(2,143), ffnenu_read(2,143)

    REAL(8) :: ffnsum(143,max_nffn), ffnenu(143,max_nffn)
    INTEGER :: inuc(2), iffn(2), iffn_sort(max_nffn)
    INTEGER :: ii, jj, inucmin, inucmax, ierr
    LOGICAL :: keep_rate, is_sorted

    INTEGER,PARAMETER  :: num_reac = 7677
    CHARACTER(LEN=5)   :: sunet_list(nnet),s1,s2
    INTEGER            :: z_list(119),z
    CHARACTER(LEN=2)   :: element_list(119),el,el_min1
    INTEGER            :: i,j,k,x
    REAL               :: beta_plus, beta_minus, leps_plus, leps_minus, nu, nubar
    REAL(8)            :: q
    !REAL(8)            :: lsum_plus(num_reac,143),lsum_minus(num_reac,143),nu_list(num_reac,143),nubar_list(num_reac,143)
    REAL               :: lsum_plus(143),lsum_minus(143),nu_list(143),nubar_list(143)
    CHARACTER(LEN=3)   :: a
    !REAL               :: temp_sum_plus(num_reac,11,13),temp_sum_minus(num_reac,11,13),temp_nu(num_reac,11,13),temp_nubar(num_reac,11,13)
    REAL               :: temp_sum_plus(11,13),temp_sum_minus(11,13),temp_nu(11,13),temp_nubar(11,13)
    
    nffn = 0
    !skip routine if weak reactions are ignored
    IF ( .not. netweak_flag ) RETURN
    
    IF ( netweak_in_fname == 'updated_rate_table.txt') THEN
      !build element list for use in headers and network check
      CALL build_element_list(element_list_fname, lun_element_in, z_list, element_list)

      READ(lun_netweak_in,*, IOSTAT=ierr)

      DO i = 1, num_reac
        !read first line from each block of data
        READ(lun_netweak_in,'(1x,a3,1x,i3)', advance = 'no') a, z
        READ(lun_netweak_in, '(33x,f8.3,1x,f8.3,2x,f8.3,2x,f8.3,1x,f8.3,2x,f8.3)') beta_plus, leps_minus, nu, &
                                          beta_minus, leps_plus, nubar
        a = trim(adjustl(a))
                
        !determine if all species are in our network
        call check_network(a, z, element_list, nname_net, nnet, keep_rate)
                
        !if all species are in our network, read in relevant data, else go to next reaction
        IF ( keep_rate ) THEN
          
          lsum_plus(1) = log10(10.d0**beta_plus + 10.d0**leps_minus)
          lsum_minus(1) = log10(10.d0**beta_minus + 10.d0**leps_plus)
          nu_list(1) = nu
          nubar_list(1) = nubar
          
          !rename H1 and free neutrons or construct species name
          IF ( z == 1 .and. a == '1  ') THEN
            s1 = '    p'
            s2 = '    n'
          ELSE
            el = element_list(z+1)
            el_min1 = element_list(z)
            s1 = adjustr(el//a)
            s2 = adjustr(el_min1//a)
          END IF
                        
          call calculate_q_value( s1, s2, q ) 

          !count electron capture rate
          nffn = nffn + 1
          nname_ffn(1, nffn) = s1
          nname_ffn(2, nffn) = s2
          q_ffn(nffn) = q

          !count beta decay rate
          nffn = nffn + 1
          nname_ffn(1, nffn) = s2
          nname_ffn(2, nffn) = s1
          q_ffn(nffn) = -q
                        
          IF ( nu_list(1) < -99.999 ) THEN
            nu_list(1) = -99.99
          END IF

          IF ( nubar_list(1) < -99.999 ) THEN
            nubar_list(1) = -99.999
          END IF

          IF ( lsum_plus(1) < -99.999 ) THEN
            lsum_plus(1) = -99.999
          END IF

          IF ( lsum_minus(1) < -99.999 ) THEN
            lsum_minus(1) = -99.999
          END IF
                        
          DO j = 2,143
            READ(lun_netweak_in, '(41x, f8.3,1x,f8.3,2x,f8.3,2x,f8.3,1x,f8.3,2x,f8.3)') beta_plus, leps_minus, nu, &
                    beta_minus, leps_plus, nubar
            
            lsum_plus(j) = log10(10.d0**beta_plus + 10.d0**leps_minus)
            lsum_minus(j) = log10(10.d0**beta_minus + 10.d0**leps_plus)
            nu_list(j) = nu
            nubar_list(j) = nubar

            IF ( nu_list(j) < -99.999 ) THEN  
              nu_list(j) = -99.999
            END IF

            IF ( nubar_list(j) < -99.999 ) THEN  
              nubar_list(j) = -99.999
            END IF

            IF ( lsum_plus(j) < -99.999 ) THEN
              lsum_plus(j) = -99.999
            END IF

            IF ( lsum_minus(j) < -99.999 ) THEN
              lsum_minus(j) = -99.999
            END IF
          END DO
        ELSE
          !if reaction is not in the network, mark it as such
          a = 'NiN'
          z = -1
          
          !skip reaction that is not in our network
          DO j = 2,143
            READ(lun_netweak_in,*)
          END DO
        END IF
        
        !re-index data that is in the network
        IF ( a .ne. 'NiN' ) THEN
          x = 0
          DO j = 1,13
            DO k = 1,11
              x = x+1
              temp_sum_plus(k,j) = lsum_plus(x)
              temp_sum_minus(k,j) = lsum_minus(x)
              temp_nu(k,j) = nu_list(x)
              temp_nubar(k,j) = nubar_list(x)
            END DO
          END DO
          
          x = 0
          DO j = 1,11
            DO k = 1,13
              x = x+1
              lsum_plus(x) = temp_sum_plus(j,k)
              lsum_minus(x) = temp_sum_minus(j,k)
              nu_list(x) = temp_nu(j,k)
              nubar_list(x) = temp_nubar(j,k)
            END DO
          END DO

          !write reaction data to netweak file
          WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') s1,s2,'ecr',q
          WRITE(lun_netweak_out,'(9(f8.3))') (lsum_plus(j),nu_list(j), j=1,143)
          WRITE(lun_netweak_out,'(5x,2a5,28x,a3,6x,1pe12.5)') s2,s1,'ecr',-1*q
          WRITE(lun_netweak_out,'(9(f8.3))') (lsum_minus(j),nubar_list(j), j=1,143)

        END IF  


      END DO
        
    ELSE
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
    END IF
  END SUBROUTINE build_netweak

END MODULE ffn_module
