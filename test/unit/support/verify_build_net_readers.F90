Program verify_build_net_readers
  Use nuclear_data, Only: aa, angm, g, mex, nname, nn, ny, read_nuclear_data, zz
  Use reaction_data, Only: ires1, ires2, ires3, ires4, irev1, irev2, iwk1, iwk2, &
    & mu1, mu2, n1i, n2i, n10, n11, n20, n21, n22, nan, nreac, q1, q2, rc1, rc2, &
    & read_reaction_data
  Use xnet_controls, Only: idiag, iheat, iscrn, iweak0, nzevolve, nzbatchmx, szbatch, tid, &
    & zb_hi, zb_lo
  Use xnet_jacobian, Only: cidx, l1s, l2s, l3s, l4s, lval, ns11, ns21, ns22, pb, &
    & read_jacobian_data, ridx
  Use xnet_match, Only: descx, ifl1, ifl2, iwflx, mflx, nflx, qflx, &
    & read_match_data
  Use xnet_types, Only: dp
  Implicit None

  Real(dp), Parameter :: tolerance = 1.0e-8_dp
  Character(80) :: data_desc
  Character(256) :: data_dir
  Integer :: i, j, k

  If ( command_argument_count() /= 1 ) Then
    Write(*,*) 'usage: verify_build_net_readers DATA_DIR'
    Stop 1
  EndIf
  Call get_command_argument(1,data_dir)

  idiag = -1
  iheat = 0
  iscrn = 0
  iweak0 = 1
  nzbatchmx = 1
  nzevolve = 1
  szbatch = 1
  tid = 1
  zb_lo = 1
  zb_hi = 1

  Call read_nuclear_data(trim(data_dir),data_desc)
  Call read_reaction_data(trim(data_dir))
  Call read_match_data(trim(data_dir))
  Call read_jacobian_data(trim(data_dir))

  Call require(trim(data_desc) == 'build_net synthetic fixture', &
    & 'production nuclear reader returned the wrong description')
  Call require(ny == 5,'production nuclear reader returned the wrong species count')
  Call require(all(nname(1:ny) == (/ '    n', '    p', '  he4', '  c12', '  o16' /)), &
    & 'production nuclear reader returned the wrong species order')
  Call require(maxval(abs(aa-(/ 1.0_dp, 1.0_dp, 4.0_dp, 12.0_dp, 16.0_dp /))) < tolerance .and. &
    & maxval(abs(zz-(/ 0.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, 8.0_dp /))) < tolerance .and. &
    & maxval(abs(nn-(/ 1.0_dp, 0.0_dp, 2.0_dp, 6.0_dp, 8.0_dp /))) < tolerance, &
    & 'production nuclear reader returned the wrong identities')
  Call require(maxval(abs(mex-(/ 8.07131710_dp, 7.28897060_dp, 2.42491560_dp, &
    & 0.0_dp, -4.73700140_dp /))) < tolerance, &
    & 'production nuclear reader returned the wrong selected masses')
  Call require(maxval(abs(angm(1:ny)-(/ 2.0_dp, 2.0_dp, 1.0_dp, 1.0_dp, 1.0_dp /))) < tolerance .and. &
    & maxval(abs(g-1.0_dp)) < tolerance, &
    & 'production nuclear reader returned the wrong partition data')

  Call require(all(nreac == (/ 3, 1, 0, 0 /)),'production reaction reader returned the wrong counts')
  Call require(all(nan == (/ 7, 3, 0, 0 /)),'production reaction reader returned the wrong map counts')
  Call require(all(n1i(:,1) == (/ 1, 2, 0, 0, 0 /)) .and. &
    & all(n1i(:,2) == (/ 2, 1, 0, 0, 0 /)) .and. &
    & all(n1i(:,3) == (/ 5, 3, 4, 0, 0 /)) .and. &
    & all(n2i(:,1) == (/ 3, 4, 5, 0, 0, 0 /)), &
    & 'production reaction reader returned the wrong participants')
  Call require(all(iwk1 == (/ 2, 2, 0 /)) .and. all(iwk2 == 0), &
    & 'production reaction reader returned the wrong weak flags')
  Call require(all(irev1 == (/ 0, 0, 1 /)) .and. all(irev2 == 0), &
    & 'production reaction reader returned the wrong reverse flags')
  Call require(all(ires1 == 0) .and. all(ires2 == 0) .and. all(ires3 == 0) .and. &
    & all(ires4 == 0),'production reaction reader returned the wrong resonance flags')
  Call require(abs(q1(1)-expected_q(n1i(:,1),1)) < tolerance .and. &
    & abs(q1(2)-expected_q(n1i(:,2),1)) < tolerance .and. &
    & abs(q1(3)-expected_q(n1i(:,3),1)) < tolerance .and. &
    & abs(q2(1)-expected_q(n2i(:,1),2)) < tolerance .and. &
    & abs(q1(1)+q1(2)) < tolerance .and. abs(q1(3)+q2(1)) < tolerance, &
    & 'production reaction reader did not preserve mass-consistent Q values')
  Call require(abs(rc1(1,1)-2.0_dp) < tolerance .and. abs(rc1(1,2)-1.0_dp) < tolerance .and. &
    & abs(rc1(1,3)+100.0_dp) < tolerance .and. abs(rc2(1,1)+100.0_dp) < tolerance, &
    & 'production reaction reader returned the wrong rate coefficients')

  Call require(mflx == 3,'production match reader returned the wrong rate-match count')
  Call require(all(descx == (/ ' ffn', ' ffn', 'syn1' /)) .and. &
    & all(iwflx == (/ 2, 2, 0 /)),'production match reader returned the wrong metadata')
  Call validate_match_group(n1i,q1,iwk1,ifl1,1)
  Call validate_match_group(n2i,q2,iwk2,ifl2,2)

  Call require(lval == 13 .and. size(ridx) == lval .and. size(cidx) == lval, &
    & 'production sparse reader returned the wrong coordinate count')
  Call require(all(pb == (/ 1, 3, 5, 8, 11, 14 /)), &
    & 'production sparse reader returned the wrong row pointers')
  Call require(all(ridx == (/ 1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5 /)) .and. &
    & all(cidx == (/ 1, 2, 1, 2, 3, 4, 5, 3, 4, 5, 3, 4, 5 /)), &
    & 'production sparse reader returned the wrong ordered coordinates')
  Call require(all((/ l1s, l2s, l3s, l4s /) == nan), &
    & 'production sparse reader map sizes disagree with reactions')

  Do j = 1, nan(1)
    Call require(mu1(j) >= 1 .and. mu1(j) <= nreac(1) .and. &
      & n11(j) == n1i(1,mu1(j)),'one-reactant extended map mismatch')
    Call require(sparse_map_matches(ns11(j),n10(j),n11(j)),'ns11 sparse map mismatch')
  EndDo
  Do j = 1, nan(2)
    Call require(mu2(j) >= 1 .and. mu2(j) <= nreac(2) .and. &
      & n21(j) == n2i(1,mu2(j)) .and. n22(j) == n2i(2,mu2(j)), &
      & 'two-reactant extended map mismatch')
    Call require(sparse_map_matches(ns21(j),n20(j),n21(j)) .and. &
      & sparse_map_matches(ns22(j),n20(j),n22(j)),'two-reactant sparse map mismatch')
  EndDo
  Do i = 1, ny
    k = 0
    Do j = pb(i), pb(i+1)-1
      If ( cidx(j) == i ) k = k + 1
    EndDo
    Call require(k == 1,'production sparse reader row lacks exactly one diagonal')
  EndDo

  Write(*,*) 'build_net production reader semantics passed'

Contains

  Subroutine require(condition,message)
    Implicit None

    Logical, Intent(in) :: condition
    Character(*), Intent(in) :: message

    If ( .not. condition ) Then
      Write(*,*) trim(message)
      Stop 1
    EndIf

    Return
  End Subroutine require

  Real(dp) Function expected_q(indices,nreactant)
    Implicit None

    Integer, Intent(in) :: indices(:), nreactant

    Integer :: index, participant_count

    participant_count = count(indices /= 0)
    expected_q = 0.0_dp
    Do index = 1, nreactant
      expected_q = expected_q + mex(indices(index))
    EndDo
    Do index = nreactant+1, participant_count
      expected_q = expected_q - mex(indices(index))
    EndDo

    Return
  End Function expected_q

  Subroutine validate_match_group(indices,qvalues,weak,ifl,nreactant)
    Implicit None

    Integer, Intent(in) :: indices(:,:), weak(:), ifl(:), nreactant
    Real(dp), Intent(in) :: qvalues(:)

    Integer :: canonical(8), flux_index, index, participant_count, product_count

    Do index = 1, size(ifl)
      flux_index = abs(ifl(index))
      Call require(flux_index >= 1 .and. flux_index <= mflx, &
        & 'reaction has an invalid match association')
      participant_count = count(indices(:,index) /= 0)
      product_count = participant_count - nreactant
      canonical = 0
      If ( qvalues(index) > 0.0_dp ) Then
        canonical(1:nreactant) = indices(1:nreactant,index)
        canonical(5:4+product_count) = indices(nreactant+1:participant_count,index)
        Call require(ifl(index) > 0,'forward reaction has a reverse match sign')
      Else
        canonical(1:product_count) = indices(nreactant+1:participant_count,index)
        canonical(5:4+nreactant) = indices(1:nreactant,index)
        Call require(ifl(index) < 0,'reverse reaction has a forward match sign')
      EndIf
      Call require(all(canonical == nflx(:,flux_index)) .and. &
        & abs(abs(qvalues(index))-qflx(flux_index)) < tolerance .and. &
        & weak(index) == iwflx(flux_index), &
        & 'reaction and match semantics disagree')
    EndDo

    Return
  End Subroutine validate_match_group

  Logical Function sparse_map_matches(map_index,row_index,column_index)
    Implicit None

    Integer, Intent(in) :: map_index, row_index, column_index

    sparse_map_matches = map_index >= 1 .and. map_index <= lval
    If ( sparse_map_matches ) Then
      sparse_map_matches = ridx(map_index) == row_index .and. cidx(map_index) == column_index
    EndIf

    Return
  End Function sparse_map_matches

End Program verify_build_net_readers
