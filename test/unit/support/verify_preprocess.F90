Program verify_preprocess
  Use, Intrinsic :: iso_fortran_env, Only: iostat_end
  Use nuclear_data, Only: aa, angm, be, g, mex, nname, nn, ny, read_nuclear_data, t9i, zz
  Use reaction_data, Only: a1, a2, a3, a4, ires1, ires2, ires3, ires4, irev1, irev2, irev3, irev4, &
    & iwk1, iwk2, iwk3, iwk4, &
    & mu1, mu2, mu3, mu4, n1i, n2i, n3i, n4i, n10, n11, n20, n21, n22, n30, n31, n32, n33, &
    & n40, n41, n42, n43, n44, nan, nreac, q1, q2, q3, q4, rc1, rc2, rc3, rc4, read_reaction_data
  Use xnet_controls, Only: idiag, iheat, iscrn, iweak0, nzbatchmx, nzevolve, szbatch, tid, zb_hi, zb_lo
  Use xnet_jacobian, Only: cidx, l1s, l2s, l3s, l4s, lval, ns11, ns21, ns22, ns31, ns32, ns33, &
    & ns41, ns42, ns43, ns44, pb, read_jacobian_data, ridx
  Use xnet_match, Only: descx, ifl1, ifl2, ifl3, ifl4, iwflx, mflx, nflx, qflx, read_match_data
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Parameter :: tolerance = 1.0e-12_dp
  Character(*), Parameter :: expected_description = 'preprocess contract fixture'
  Character(80) :: data_desc
  Character(256) :: data_dir, message, summary_file
  Integer :: chapter(11), chapter_end(11), chapter_start(11)
  Integer :: saved_index, saved_integer
  Logical :: ok

  If ( command_argument_count() /= 2 ) Then
    Write(*,*) 'usage: verify_preprocess DATA_DIR SUMMARY_FILE'
    Stop 1
  EndIf
  Call get_command_argument(1,data_dir)
  Call get_command_argument(2,summary_file)

  idiag = -1
  iheat = 0
  iscrn = 0
  iweak0 = 0
  nzbatchmx = 1
  nzevolve = 1
  szbatch = 1
  tid = 1
  zb_lo = 1
  zb_hi = 1

  Call read_nuclear_data(trim(data_dir),data_desc)
  Call require(trim(data_desc) == expected_description,'production net_desc description mismatch')
  Call read_reaction_data(trim(data_dir))
  Call read_match_data(trim(data_dir))
  Call read_jacobian_data(trim(data_dir))
  Call check_generated_artifacts(trim(data_dir),chapter,chapter_start,chapter_end)
  Call validate_loaded(ok,message)
  Call require(ok,'positive semantic validation: '//trim(message))

  saved_integer = mu1(1)
  mu1(1) = nreac(1) + 1
  Call validate_loaded(ok,message)
  Call require(.not. ok,'wrong reaction index was accepted')
  mu1(1) = saved_integer

  saved_integer = cidx(pb(1))
  cidx(pb(1)) = ny + 1
  Call validate_loaded(ok,message)
  Call require(.not. ok,'out-of-range sparse column was accepted')
  cidx(pb(1)) = saved_integer

  saved_index = find_diagonal(1)
  saved_integer = cidx(saved_index)
  cidx(saved_index) = 2
  Call validate_loaded(ok,message)
  Call require(.not. ok,'missing sparse diagonal was accepted')
  cidx(saved_index) = saved_integer

  saved_integer = pb(2)
  pb(2) = pb(1) - 1
  Call validate_loaded(ok,message)
  Call require(.not. ok,'nonmonotone sparse row pointer was accepted')
  pb(2) = saved_integer

  saved_integer = ifl1(2)
  ifl1(2) = -ifl1(2)
  Call validate_loaded(ok,message)
  Call require(.not. ok,'corrupt reverse/match association was accepted')
  ifl1(2) = saved_integer

  Call validate_loaded(ok,message)
  Call require(ok,'state restoration after controlled corruptions: '//trim(message))
  Call write_summary(trim(summary_file),chapter,chapter_start,chapter_end)
  Write(*,*) 'preprocess semantic contract passed'

Contains

  Subroutine require(condition,diagnostic)
    Implicit None

    Logical, Intent(in) :: condition
    Character(*), Intent(in) :: diagnostic

    If ( .not. condition ) Then
      Write(*,*) trim(diagnostic)
      Stop 1
    EndIf

    Return
  End Subroutine require

  Subroutine invalidate(ok_out,message_out,diagnostic)
    Implicit None

    Logical, Intent(out) :: ok_out
    Character(*), Intent(out) :: message_out
    Character(*), Intent(in) :: diagnostic

    ok_out = .False.
    message_out = diagnostic

    Return
  End Subroutine invalidate

  Subroutine validate_loaded(ok_out,message_out)
    Implicit None

    Logical, Intent(out) :: ok_out
    Character(*), Intent(out) :: message_out

    Integer :: i, j, k

    ok_out = .True.
    message_out = ''
    If ( ny /= 6 ) Then
      Call invalidate(ok_out,message_out,'wrong species count')
      Return
    EndIf
    If ( any(nname(1:ny) /= (/ '    n', '    p', '  he4', '  c12', '  o16', ' ne20' /)) ) Then
      Call invalidate(ok_out,message_out,'wrong species ordering')
      Return
    EndIf
    If ( any(nreac /= (/ 3, 1, 1, 0 /)) ) Then
      Call invalidate(ok_out,message_out,'wrong reaction counts; unavailable participant may have been retained')
      Return
    EndIf
    If ( any(ires1 /= 0) .or. any(ires2 /= 0) .or. any(ires3 /= 0) .or. any(ires4 /= 0) ) Then
      Call invalidate(ok_out,message_out,'reaction resonance flag mismatch')
      Return
    EndIf
    If ( any(nan /= (/ 7, 3, 2, 0 /)) ) Then
      Call invalidate(ok_out,message_out,'wrong extended reaction counts')
      Return
    EndIf
    If ( any(n1i(:,1) /= (/ 1, 2, 0, 0, 0 /)) .or. iwk1(1) /= 1 .or. irev1(1) /= 0 ) Then
      Call invalidate(ok_out,message_out,'weak reaction translation mismatch')
      Return
    EndIf
    If ( any(n1i(:,2) /= (/ 5, 3, 4, 0, 0 /)) .or. irev1(2) /= 1 ) Then
      Call invalidate(ok_out,message_out,'reverse alpha-capture translation mismatch')
      Return
    EndIf
    If ( any(n1i(:,3) /= (/ 4, 3, 3, 3, 0 /)) .or. irev1(3) /= 1 ) Then
      Call invalidate(ok_out,message_out,'reverse triple-alpha translation mismatch')
      Return
    EndIf
    If ( any(n2i(:,1) /= (/ 3, 4, 5, 0, 0, 0 /)) .or. irev2(1) /= 0 ) Then
      Call invalidate(ok_out,message_out,'forward alpha-capture translation mismatch')
      Return
    EndIf
    If ( any(n3i(:,1) /= (/ 3, 3, 3, 4, 0, 0 /)) .or. irev3(1) /= 0 ) Then
      Call invalidate(ok_out,message_out,'forward triple-alpha translation mismatch')
      Return
    EndIf
    If ( abs(rc1(1,2)-94.3131_dp) > tolerance .or. abs(rc2(1,1)-69.6526_dp) > tolerance .or. &
      & abs(rc3(1,1)+0.971052_dp) > tolerance ) Then
      Call invalidate(ok_out,message_out,'reaction coefficient mismatch')
      Return
    EndIf
    If ( abs(q1(1)-expected_q(n1i(:,1),1)) > tolerance .or. &
      & abs(q1(2)-expected_q(n1i(:,2),1)) > tolerance .or. &
      & abs(q1(3)-expected_q(n1i(:,3),1)) > tolerance .or. &
      & abs(q2(1)-expected_q(n2i(:,1),2)) > tolerance .or. &
      & abs(q3(1)-expected_q(n3i(:,1),3)) > tolerance .or. &
      & abs(q1(2)+q2(1)) > tolerance .or. abs(q1(3)+q3(1)) > tolerance ) Then
      Call invalidate(ok_out,message_out,'reaction Q was not recomputed from nuclear masses')
      Return
    EndIf
    If ( .not. any(mu1 == 3 .and. abs(a1-3.0_dp) < tolerance) .or. &
      & .not. any(mu3 == 1 .and. abs(a3+0.5_dp) < tolerance) .or. &
      & .not. any(mu3 == 1 .and. abs(a3-1.0_dp/6.0_dp) < tolerance) ) Then
      Call invalidate(ok_out,message_out,'repeated-participant multiplicity mismatch')
      Return
    EndIf

    If ( any(mu1 < 1) .or. any(mu1 > nreac(1)) .or. any(mu2 < 1) .or. any(mu2 > nreac(2)) .or. &
      & any(mu3 < 1) .or. any(mu3 > nreac(3)) ) Then
      Call invalidate(ok_out,message_out,'reaction-to-species map has invalid reaction index')
      Return
    EndIf
    Do j = 1, nan(1)
      If ( n11(j) /= n1i(1,mu1(j)) ) Then
        Call invalidate(ok_out,message_out,'one-reactant extended map mismatch')
        Return
      EndIf
    EndDo
    Do j = 1, nan(2)
      If ( n21(j) /= n2i(1,mu2(j)) .or. n22(j) /= n2i(2,mu2(j)) ) Then
        Call invalidate(ok_out,message_out,'two-reactant extended map mismatch')
        Return
      EndIf
    EndDo
    Do j = 1, nan(3)
      If ( n31(j) /= n3i(1,mu3(j)) .or. n32(j) /= n3i(2,mu3(j)) .or. &
        & n33(j) /= n3i(3,mu3(j)) ) Then
        Call invalidate(ok_out,message_out,'three-reactant extended map mismatch')
        Return
      EndIf
    EndDo

    If ( mflx /= 3 .or. any(abs(ifl1) > mflx) .or. any(abs(ifl2) > mflx) .or. any(abs(ifl3) > mflx) ) Then
      Call invalidate(ok_out,message_out,'match index range mismatch')
      Return
    EndIf
    Call validate_match_group(n1i,q1,iwk1,ifl1,1,ok_out,message_out)
    If ( .not. ok_out ) Return
    Call validate_match_group(n2i,q2,iwk2,ifl2,2,ok_out,message_out)
    If ( .not. ok_out ) Return
    Call validate_match_group(n3i,q3,iwk3,ifl3,3,ok_out,message_out)
    If ( .not. ok_out ) Return

    If ( lval /= size(ridx) .or. lval /= size(cidx) .or. size(pb) /= ny+1 ) Then
      Call invalidate(ok_out,message_out,'incomplete sparse records')
      Return
    EndIf
    If ( pb(1) /= 1 .or. pb(ny+1) /= lval+1 ) Then
      Call invalidate(ok_out,message_out,'wrong sparse terminal pointer')
      Return
    EndIf
    Do i = 1, ny
      If ( pb(i+1) <= pb(i) ) Then
        Call invalidate(ok_out,message_out,'nonmonotone sparse row pointer')
        Return
      EndIf
      Do k = pb(i), pb(i+1)-1
        If ( ridx(k) /= i .or. cidx(k) < 1 .or. cidx(k) > ny ) Then
          Call invalidate(ok_out,message_out,'out-of-range sparse coordinate')
          Return
        EndIf
        If ( k > pb(i) ) Then
          If ( cidx(k) <= cidx(k-1) ) Then
            Call invalidate(ok_out,message_out,'sparse columns are not strictly ordered')
            Return
          EndIf
        EndIf
      EndDo
      If ( find_diagonal(i) == 0 ) Then
        Call invalidate(ok_out,message_out,'missing sparse diagonal')
        Return
      EndIf
    EndDo
    If ( any((/ l1s, l2s, l3s, l4s /) /= nan) ) Then
      Call invalidate(ok_out,message_out,'sparse map sizes disagree with reaction data')
      Return
    EndIf
    Call validate_sparse_maps(ok_out,message_out)

    Return
  End Subroutine validate_loaded

  Subroutine validate_match_group(indices,qvalues,weak,ifl,nreactant,ok_out,message_out)
    Implicit None

    Integer, Intent(in) :: indices(:,:), weak(:), ifl(:), nreactant
    Real(dp), Intent(in) :: qvalues(:)
    Logical, Intent(out) :: ok_out
    Character(*), Intent(out) :: message_out

    Integer :: canonical(8), flux_index, j, participant_count, product_count

    ok_out = .True.
    message_out = ''
    Do j = 1, size(ifl)
      flux_index = abs(ifl(j))
      If ( flux_index < 1 .or. flux_index > mflx ) Then
        Call invalidate(ok_out,message_out,'reaction has invalid flux association')
        Return
      EndIf
      participant_count = count(indices(:,j) /= 0)
      product_count = participant_count - nreactant
      canonical = 0
      If ( qvalues(j) > 0.0_dp ) Then
        canonical(1:nreactant) = indices(1:nreactant,j)
        canonical(5:4+product_count) = indices(nreactant+1:participant_count,j)
        If ( ifl(j) < 0 ) Then
          Call invalidate(ok_out,message_out,'forward reaction has reverse flux sign')
          Return
        EndIf
      Else
        canonical(1:product_count) = indices(nreactant+1:participant_count,j)
        canonical(5:4+nreactant) = indices(1:nreactant,j)
        If ( ifl(j) > 0 ) Then
          Call invalidate(ok_out,message_out,'reverse reaction has forward flux sign')
          Return
        EndIf
      EndIf
      If ( any(canonical /= nflx(:,flux_index)) .or. abs(abs(qvalues(j))-qflx(flux_index)) > tolerance .or. &
        & weak(j) /= iwflx(flux_index) ) Then
        Call invalidate(ok_out,message_out,'reaction and match/flux semantics disagree')
        Return
      EndIf
    EndDo

    Return
  End Subroutine validate_match_group

  Real(dp) Function expected_q(indices,nreactant)
    Implicit None

    Integer, Intent(in) :: indices(:), nreactant

    Integer :: i, participant_count

    participant_count = count(indices /= 0)
    expected_q = 0.0_dp
    Do i = 1, nreactant
      expected_q = expected_q + mex(indices(i))
    EndDo
    Do i = nreactant+1, participant_count
      expected_q = expected_q - mex(indices(i))
    EndDo

    Return
  End Function expected_q

  Subroutine validate_sparse_maps(ok_out,message_out)
    Implicit None

    Logical, Intent(out) :: ok_out
    Character(*), Intent(out) :: message_out

    Integer :: j

    ok_out = .True.
    message_out = ''
    Do j = 1, nan(1)
      If ( .not. sparse_map_matches(ns11(j),n10(j),n11(j)) ) Then
        Call invalidate(ok_out,message_out,'ns11 sparse map mismatch')
        Return
      EndIf
    EndDo
    Do j = 1, nan(2)
      If ( .not. sparse_map_matches(ns21(j),n20(j),n21(j)) .or. &
        & .not. sparse_map_matches(ns22(j),n20(j),n22(j)) ) Then
        Call invalidate(ok_out,message_out,'two-reactant sparse map mismatch')
        Return
      EndIf
    EndDo
    Do j = 1, nan(3)
      If ( .not. sparse_map_matches(ns31(j),n30(j),n31(j)) .or. &
        & .not. sparse_map_matches(ns32(j),n30(j),n32(j)) .or. &
        & .not. sparse_map_matches(ns33(j),n30(j),n33(j)) ) Then
        Call invalidate(ok_out,message_out,'three-reactant sparse map mismatch')
        Return
      EndIf
    EndDo
    Do j = 1, nan(4)
      If ( .not. sparse_map_matches(ns41(j),n40(j),n41(j)) .or. &
        & .not. sparse_map_matches(ns42(j),n40(j),n42(j)) .or. &
        & .not. sparse_map_matches(ns43(j),n40(j),n43(j)) .or. &
        & .not. sparse_map_matches(ns44(j),n40(j),n44(j)) ) Then
        Call invalidate(ok_out,message_out,'four-reactant sparse map mismatch')
        Return
      EndIf
    EndDo

    Return
  End Subroutine validate_sparse_maps

  Logical Function sparse_map_matches(map_index,row_index,column_index)
    Implicit None

    Integer, Intent(in) :: map_index, row_index, column_index

    sparse_map_matches = map_index >= 1 .and. map_index <= lval
    If ( sparse_map_matches ) Then
      sparse_map_matches = ridx(map_index) == row_index .and. cidx(map_index) == column_index
    EndIf

    Return
  End Function sparse_map_matches

  Integer Function find_diagonal(row_index)
    Implicit None

    Integer, Intent(in) :: row_index
    Integer :: k

    find_diagonal = 0
    If ( row_index < 1 .or. row_index > ny ) Return
    If ( pb(row_index) < 1 .or. pb(row_index+1)-1 > lval ) Return
    Do k = pb(row_index), pb(row_index+1)-1
      If ( ridx(k) == row_index .and. cidx(k) == row_index ) Then
        find_diagonal = k
        Exit
      EndIf
    EndDo

    Return
  End Function find_diagonal

  Integer Function find_coordinate(row_index,column_index)
    Implicit None

    Integer, Intent(in) :: row_index, column_index
    Integer :: k

    find_coordinate = 0
    Do k = 1, lval
      If ( ridx(k) == row_index .and. cidx(k) == column_index ) Then
        find_coordinate = k
        Exit
      EndIf
    EndDo

    Return
  End Function find_coordinate

  Subroutine check_generated_artifacts(directory,chapter_out,start_out,end_out)
    Implicit None

    Character(*), Intent(in) :: directory
    Integer, Intent(out) :: chapter_out(11), start_out(11), end_out(11)

    Character(5) :: blank_name(6), expected_names(8), match_arrow, match_desc, match_names(8)
    Character(17) :: label17
    Character(32) :: label32
    Character(256) :: filename, line
    Character(32), Parameter :: artifacts(10) = (/ &
      & 'nuc_data                        ', 'nets3                           ', &
      & 'nets4                           ', 'ab_blank                        ', &
      & 'match_data                      ', 'match_read                      ', &
      & 'sparse_ind                      ', 'matr_shape                      ', &
      & 'net_desc                        ', 'net_diag                        ' /)
    Integer :: column_index, counts(4), i, ierr, ifl_orig, ifl_term, j, lun, match_index, row_index
    Integer :: binary_ny, nonreaclib_counts(2), widths(2)
    Logical :: exists
    Real(dp) :: blank_value(6), match_coordinates(4), match_value, rhs(6), value
    Character(5), Allocatable :: binary_name(:)
    Real(dp), Allocatable :: binary_aa(:), binary_angm(:), binary_be(:), binary_g(:,:), binary_nn(:)
    Real(dp), Allocatable :: binary_t9(:), binary_zz(:), shape_values(:)

    Do i = 1, size(artifacts)
      filename = trim(directory)//'/'//trim(artifacts(i))
      Inquire(file=trim(filename),exist=exists)
      Call require(exists,'missing generated artifact '//trim(artifacts(i)))
    EndDo

    Open(newunit=lun,file=trim(directory)//'/nuc_data',status='old',form='unformatted',action='read')
    Read(lun) binary_ny
    Call require(binary_ny == ny,'nuc_data species count mismatch')
    Allocate(binary_name(ny),binary_aa(ny),binary_angm(ny),binary_be(ny),binary_g(24,ny), &
      & binary_nn(ny),binary_t9(24),binary_zz(ny))
    Read(lun) binary_t9
    Read(lun) (binary_name(i),binary_aa(i),binary_zz(i),binary_nn(i),binary_be(i), &
      & binary_g(:,i),binary_angm(i),i=1,ny)
    Close(lun)
    Call require(all(binary_name == nname(1:ny)) .and. maxval(abs(binary_aa-aa)) < tolerance .and. &
      & maxval(abs(binary_zz-zz)) < tolerance .and. maxval(abs(binary_nn-nn)) < tolerance .and. &
      & maxval(abs(binary_be-be)) < tolerance .and. maxval(abs(binary_t9-t9i)) < tolerance .and. &
      & maxval(abs(binary_g-g)) < tolerance .and. maxval(abs(binary_angm-angm(1:ny))) < tolerance, &
      & 'nuc_data semantic round-trip mismatch')
    Deallocate(binary_name,binary_aa,binary_angm,binary_be,binary_g,binary_nn,binary_t9,binary_zz)

    Open(newunit=lun,file=trim(directory)//'/net_desc',status='old',action='read')
    Read(lun,'(a)') line
    Call require(trim(line) == trim(data_desc),'net_desc description payload mismatch')
    Read(lun,'(a)') line
    Call require(trim(adjustl(line)) == 'Number of Nuclear Species=    6', &
      & 'net_desc species count mismatch')
    Read(lun,'(a)') line
    Call require(trim(line) == 'Reaction Count for the 11 different types', &
      & 'net_desc reaction-count heading mismatch')
    Read(lun,*) (chapter_out(i),start_out(i),end_out(i),i=1,11)
    Call require(all(chapter_out == (/ (i,i=1,11) /)),'net_desc chapter numbering mismatch')
    Call require(all(start_out == (/ 1,2,3,1,2,2,2,1,2,1,4 /)) .and. &
      & all(end_out == (/ 1,2,3,1,1,1,1,1,1,0,3 /)),'net_desc chapter ranges mismatch')
    Read(lun,'(a)') line
    Call require(trim(line) == 'Necessary dimensions','net_desc dimension heading mismatch')
    Read(lun,'(a17,4i8)') label17, counts
    Call require(trim(adjustl(label17)) == 'nreac(1,2,3,4)=' .and. all(counts == nreac), &
      & 'net_desc reaction dimensions mismatch')
    Read(lun,'(a17,4i8)') label17, counts
    Call require(trim(adjustl(label17)) == 'nan(1,2,3,4)=' .and. all(counts == nan), &
      & 'net_desc extended dimensions mismatch')
    Read(lun,'(a)') line
    Call require(trim(line) == 'Reaction Count for non-REACLIB rates', &
      & 'net_desc non-REACLIB heading mismatch')
    Read(lun,'(a7,i8,a7,i8)') label32(1:7), nonreaclib_counts(1), label32(8:14), nonreaclib_counts(2)
    Call require(trim(adjustl(label32(1:7))) == 'nffn=' .and. trim(adjustl(label32(8:14))) == 'nnu=' .and. &
      & all(nonreaclib_counts == 0),'net_desc non-REACLIB counts mismatch')
    Read(lun,'(a)') line
    Call require(trim(line) == 'Matrix Sparseness parameters','net_desc sparse heading mismatch')
    Read(lun,'(a16,2i5)') label32(1:16), widths
    Call require(trim(label32(1:16)) == 'Border Widths' .and. all(widths == ny), &
      & 'net_desc border widths mismatch')
    Read(lun,'(a16,2i5)') label32(1:16), widths
    Call require(trim(label32(1:16)) == 'Diagonal Widths' .and. all(widths == 0), &
      & 'net_desc diagonal widths mismatch')
    Read(lun,'(a)',iostat=ierr) line
    Call require(ierr == iostat_end,'net_desc has unexpected trailing records')
    Close(lun)

    Open(newunit=lun,file=trim(directory)//'/ab_blank',status='old',action='read')
    Read(lun,'(a)') line
    Read(lun,*) (blank_name(i),blank_value(i),i=1,ny)
    Close(lun)
    Call require(trim(line) == 'Abundance Description' .and. maxval(abs(blank_value)) < tolerance, &
      & 'ab_blank semantic mismatch')
    Do i = 1, ny
      Call require(trim(adjustl(blank_name(i))) == trim(adjustl(nname(i))), &
        & 'ab_blank species ordering mismatch')
    EndDo

    Open(newunit=lun,file=trim(directory)//'/match_read',status='old',action='read')
    Do i = 1, mflx
      Read(lun,'(10a5,1es10.3)',iostat=ierr) match_names(1:4), match_arrow, match_names(5:8), &
        & match_desc, match_value
      Call require(ierr == 0,'match_read participant record is unreadable')
      Do j = 1, 8
        expected_names(j) = participant_name(nflx(j,i))
      EndDo
      Call require(all(match_names == expected_names),'match_read participant names mismatch')
      Call require(match_arrow == ' --> ','match_read arrow mismatch')
      Call require(trim(adjustl(match_desc)) == trim(adjustl(descx(i))), &
        & 'match_read participant descriptor mismatch')
      Call require(abs(match_value) < tolerance,'match_read participant value mismatch')
    EndDo
    Do i = 1, mflx
      Read(lun,'(i5,4f6.1,es13.5,a5)',iostat=ierr) match_index, match_coordinates, match_value, match_desc
      Call require(ierr == 0,'match_read coordinate record is unreadable')
      ifl_orig = nflx(count(nflx(1:4,i) /= 0),i)
      ifl_term = nflx(count(nflx(5:8,i) /= 0)+4,i)
      Call require(match_index == i .and. &
        & maxval(abs(match_coordinates-(/ zz(ifl_orig),nn(ifl_orig),zz(ifl_term),nn(ifl_term) /))) < tolerance .and. &
        & abs(match_value-1.0_dp) < tolerance .and. &
        & trim(adjustl(match_desc)) == trim(adjustl(descx(i))), &
        & 'match_read coordinate record mismatch')
    EndDo
    Read(lun,'(a)',iostat=ierr) line
    Close(lun)
    Call require(ierr == iostat_end,'match_read has unexpected trailing records')

    Open(newunit=lun,file=trim(directory)//'/matr_shape',status='old',action='read')
    Read(lun,'(a)') line
    Call require(trim(adjustl(line)) == 'NY=   6','matr_shape dimension mismatch')
    Read(lun,*) rhs
    Call require(maxval(abs(rhs-1.0_dp)) < tolerance,'matr_shape trial right-hand side mismatch')
    Allocate(shape_values(lval))
    shape_values = 0.0_dp
    Do i = 1, ny
      shape_values(find_diagonal(i)) = shape_values(find_diagonal(i)) + 10.0_dp
    EndDo
    Do i = 1, nan(1)
      shape_values(ns11(i)) = shape_values(ns11(i)) + 1.0_dp
    EndDo
    Do i = 1, nan(2)
      shape_values(ns21(i)) = shape_values(ns21(i)) + 1.0_dp
      shape_values(ns22(i)) = shape_values(ns22(i)) + 1.0_dp
    EndDo
    Do i = 1, nan(3)
      shape_values(ns31(i)) = shape_values(ns31(i)) + 1.0_dp
      shape_values(ns32(i)) = shape_values(ns32(i)) + 1.0_dp
      shape_values(ns33(i)) = shape_values(ns33(i)) + 1.0_dp
    EndDo
    Do i = 1, nan(4)
      shape_values(ns41(i)) = shape_values(ns41(i)) + 1.0_dp
      shape_values(ns42(i)) = shape_values(ns42(i)) + 1.0_dp
      shape_values(ns43(i)) = shape_values(ns43(i)) + 1.0_dp
      shape_values(ns44(i)) = shape_values(ns44(i)) + 1.0_dp
    EndDo
    Do i = 1, lval
      Read(lun,*,iostat=ierr) row_index,column_index,value
      Call require(ierr == 0,'matr_shape coordinate is unreadable')
      Call require(row_index == ridx(i) .and. column_index == cidx(i) .and. &
        & abs(value-shape_values(i)) < tolerance,'matr_shape ordered coordinate/value mismatch')
    EndDo
    Read(lun,'(a)',iostat=ierr) line
    Close(lun)
    Call require(ierr == iostat_end,'matr_shape has unexpected trailing coordinates')
    Deallocate(shape_values)

    Call check_net_diagnostics(trim(directory)//'/net_diag')

    Return
  End Subroutine check_generated_artifacts

  Function participant_name(index) Result(name)
    Implicit None

    Integer, Intent(in) :: index
    Character(5) :: name

    If ( index >= lbound(nname,1) .and. index <= ubound(nname,1) ) Then
      name = nname(index)
    Else
      name = '     '
    EndIf

    Return
  End Function participant_name

  Subroutine check_q_diagnostic(lun,descriptor,names,input_q,computed_q)
    Implicit None

    Integer, Intent(in) :: lun
    Character(*), Intent(in) :: descriptor
    Character(5), Intent(in) :: names(6)
    Real(dp), Intent(in) :: input_q, computed_q

    Character(256) :: actual_line, expected_line
    Integer :: i, ierr

    Read(lun,'(a)',iostat=ierr) actual_line
    Call require(ierr == 0,'net_diag Q diagnostic is unreadable')
    Write(expected_line,"(a,7a6,a,es9.2,a,es9.2)") &
      & 'Inconsistent q-value for ',descriptor,(names(i),i=1,6),'  netsu: ',input_q,' netwinv: ',computed_q
    Call require(trim(actual_line) == trim(expected_line),'net_diag Q diagnostic mismatch')

    Return
  End Subroutine check_q_diagnostic

  Subroutine check_net_diagnostics(filename)
    Implicit None

    Character(*), Intent(in) :: filename

    Character(5), Parameter :: blank5 = '     ', arrow = ' --> '
    Character(256) :: actual_line, expected_line
    Integer :: i, ierr, ii, k, lun, row_index

    Open(newunit=lun,file=filename,status='old',action='read')
    Call check_q_diagnostic(lun,descx(1), &
      & (/ '    n', '    p', '     ', '     ', '     ', '     ' /),9.9_dp,q1(1))
    Call check_q_diagnostic(lun,descx(2), &
      & (/ '  o16', '  he4', '  c12', '     ', '     ', '     ' /),-1.0_dp,q1(2))
    Call check_q_diagnostic(lun,descx(3), &
      & (/ '  c12', '  he4', '  he4', '  he4', '     ', '     ' /),-2.0_dp,q1(3))
    Call check_q_diagnostic(lun,descx(2), &
      & (/ '  he4', '  c12', '  o16', '     ', '     ', '     ' /),3.0_dp,q2(1))
    Call check_q_diagnostic(lun,descx(3), &
      & (/ '  he4', '  he4', '  he4', '  c12', '     ', '     ' /),4.0_dp,q3(1))

    Do ii = 1, mflx
      Read(lun,'(a)',iostat=ierr) actual_line
      Call require(ierr == 0 .and. trim(actual_line) == '--','net_diag match-group marker mismatch')
      Do k = 1, nreac(1)
        If ( abs(ifl1(k)) == ii ) Then
          Write(expected_line,"(9a5,2i2,i6,1es12.4)") participant_name(n1i(1,k)),blank5,blank5,arrow, &
            & participant_name(n1i(2,k)),participant_name(n1i(3,k)),participant_name(n1i(4,k)), &
            & blank5,descx(ii),ires1(k),irev1(k),ifl1(k),q1(k)
          Read(lun,'(a)',iostat=ierr) actual_line
          Call require(ierr == 0 .and. trim(actual_line) == trim(expected_line), &
            & 'net_diag one-reactant match record mismatch')
        EndIf
      EndDo
      Do k = 1, nreac(2)
        If ( abs(ifl2(k)) == ii ) Then
          Write(expected_line,"(9a5,2i2,i6,1es12.4)") participant_name(n2i(1,k)), &
            & participant_name(n2i(2,k)),blank5,arrow,participant_name(n2i(3,k)), &
            & participant_name(n2i(4,k)),participant_name(n2i(5,k)),participant_name(n2i(6,k)), &
            & descx(ii),ires2(k),irev2(k),ifl2(k),q2(k)
          Read(lun,'(a)',iostat=ierr) actual_line
          Call require(ierr == 0 .and. trim(actual_line) == trim(expected_line), &
            & 'net_diag two-reactant match record mismatch')
        EndIf
      EndDo
      Do k = 1, nreac(3)
        If ( abs(ifl3(k)) == ii ) Then
          Write(expected_line,"(9a5,2i2,i6,1es12.4)") participant_name(n3i(1,k)), &
            & participant_name(n3i(2,k)),participant_name(n3i(3,k)),arrow,participant_name(n3i(4,k)), &
            & participant_name(n3i(5,k)),blank5,blank5,descx(ii),ires3(k),irev3(k),ifl3(k),q3(k)
          Read(lun,'(a)',iostat=ierr) actual_line
          Call require(ierr == 0 .and. trim(actual_line) == trim(expected_line), &
            & 'net_diag three-reactant match record mismatch')
        EndIf
      EndDo
    EndDo

    Do i = 1, ny
      Read(lun,'(a)',iostat=ierr) actual_line
      Call require(ierr == 0,'net_diag sparse-row index is unreadable')
      Read(actual_line,*,iostat=ierr) row_index
      Call require(ierr == 0 .and. row_index == i,'net_diag sparse-row index mismatch')
      Read(lun,'(a)',iostat=ierr) actual_line
      Call require(ierr == 0,'net_diag sparse-row columns are unreadable')
      Write(expected_line,"(18i4)") cidx(pb(i):pb(i+1)-1)
      Call require(trim(actual_line) == trim(expected_line),'net_diag sparse-row columns mismatch')
    EndDo
    Read(lun,'(a)',iostat=ierr) actual_line
    Call require(ierr == iostat_end,'net_diag has unexpected trailing records')
    Close(lun)

    Return
  End Subroutine check_net_diagnostics

  Subroutine write_summary(filename,chapter_in,start_in,end_in)
    Implicit None

    Character(*), Intent(in) :: filename
    Integer, Intent(in) :: chapter_in(11), start_in(11), end_in(11)

    Integer :: lun

    Open(newunit=lun,file=filename,status='replace',action='write')
    Write(lun,'(a)') trim(data_desc)
    Write(lun,*) ny,nname(1:ny)
    Write(lun,*) aa,zz,nn,be,t9i,g,angm(1:ny)
    Write(lun,*) chapter_in,start_in,end_in
    Write(lun,*) nreac,nan
    Write(lun,*) n1i,iwk1,ires1,irev1,rc1,q1
    Write(lun,*) n2i,iwk2,ires2,irev2,rc2,q2
    Write(lun,*) n3i,iwk3,ires3,irev3,rc3,q3
    Write(lun,*) mu1,a1,mu2,a2,mu3,a3
    Write(lun,*) mflx,ifl1,ifl2,ifl3,nflx,qflx,iwflx,descx
    Write(lun,*) lval,ridx,cidx,pb,l1s,l2s,l3s,l4s
    Write(lun,*) ns11,ns21,ns22,ns31,ns32,ns33
    Close(lun)

    Return
  End Subroutine write_summary

End Program verify_preprocess
