!***************************************************************************************************
! Direct Frontier qualification of XNet's production batched GPU factor/solve path.
!***************************************************************************************************

Program frontier_gpu_linalg_probe
  Use, Intrinsic :: iso_c_binding, Only: C_LOC, C_SIZEOF
  Use, Intrinsic :: iso_fortran_env, Only: output_unit
  Use xnet_controls, Only: myid, nproc, tid, nthread
  Use xnet_gpu, Only: &
    & device_is_present, deviceCount, gpu_finalize, gpu_init, mydevice, on_device
  Use xnet_linalg, Only: LinearSolveBatched
  Use xnet_types, Only: dp
  Implicit None

  Integer, Parameter :: matrix_size = 3
  Integer, Parameter :: batch_count = 2
  Integer, Parameter :: right_hand_sides = 1
  Real(dp), Parameter :: residual_limit = 1.0e-12_dp

  Integer :: ibatch
  Integer :: column_offset
  Integer, Target :: info(batch_count)
  Integer, Target :: ipiv(matrix_size,batch_count)
  Integer :: failure_count
  Real(dp), Target :: matrix(matrix_size,matrix_size*batch_count)
  Real(dp) :: matrix_original(matrix_size,matrix_size*batch_count)
  Real(dp), Target :: rhs(matrix_size,batch_count)
  Real(dp) :: rhs_original(matrix_size,batch_count)
  Real(dp) :: residual(matrix_size)
  Real(dp) :: relative_residual
  Real(dp) :: denominator
  Logical :: offloaded
  Logical :: data_present

  myid = 0
  nproc = 1
  tid = 1
  nthread = 1
  failure_count = 0
  info = -1
  ipiv = 0

  matrix(:,1:3) = Reshape( &
    & [4.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 0.0_dp, 2.0_dp, 0.0_dp, 5.0_dp], &
    & [matrix_size,matrix_size])
  matrix(:,4:6) = Reshape( &
    & [3.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, 4.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], &
    & [matrix_size,matrix_size])
  rhs(:,1) = [4.0_dp, 7.0_dp, -3.0_dp]
  rhs(:,2) = [-6.5_dp, 7.0_dp, 6.5_dp]
  matrix_original = matrix
  rhs_original = rhs

  Call gpu_init()

  offloaded = .false.
  !$omp target map(from:offloaded)
  offloaded = on_device()
  !$omp end target

  Write(output_unit,"(a,i4)") 'XNET_GPU_LINALG device_count ',deviceCount
  Write(output_unit,"(a,i4)") 'XNET_GPU_LINALG device ',mydevice
  Write(output_unit,"(a,l1)") 'XNET_GPU_LINALG offloaded ',offloaded
  If ( deviceCount < 1 .or. .not. offloaded ) failure_count = failure_count + 1

  !$omp target enter data map(to:matrix,rhs) map(alloc:ipiv,info)
  data_present = &
    & device_is_present( &
      & C_LOC(matrix(1,1)), mydevice, Size(matrix) * C_SIZEOF(matrix(1,1))) &
    & .And. device_is_present( &
      & C_LOC(rhs(1,1)), mydevice, Size(rhs) * C_SIZEOF(rhs(1,1))) &
    & .And. device_is_present( &
      & C_LOC(ipiv(1,1)), mydevice, Size(ipiv) * C_SIZEOF(ipiv(1,1))) &
    & .And. device_is_present( &
      & C_LOC(info(1)), mydevice, Size(info) * C_SIZEOF(info(1)))
  Write(output_unit,"(a,l1)") 'XNET_GPU_LINALG data_present ',data_present
  If ( data_present ) Then
    Call LinearSolveBatched( &
      & 'N', matrix_size, right_hand_sides, matrix, matrix_size, ipiv(1,1), rhs, matrix_size, &
      & info(1), batch_count)
  Else
    failure_count = failure_count + 1
  EndIf
  !$omp target update from(rhs,info)
  !$omp target exit data map(release:matrix,rhs,ipiv,info)

  Do ibatch = 1, batch_count
    column_offset = (ibatch - 1) * matrix_size
    residual = Matmul( &
      & matrix_original(:,column_offset+1:column_offset+matrix_size), rhs(:,ibatch)) &
      & - rhs_original(:,ibatch)
    denominator = Maxval( &
      & Sum(Abs(matrix_original(:,column_offset+1:column_offset+matrix_size)), Dim=2)) &
      & * Maxval(Abs(rhs(:,ibatch))) + Maxval(Abs(rhs_original(:,ibatch)))
    relative_residual = Maxval(Abs(residual)) / denominator
    Write(output_unit,"(a,i4,a,i4,a,es24.16)") &
      & 'XNET_GPU_LINALG batch ',ibatch,' info ',info(ibatch), &
      & ' relative_residual ',relative_residual
    If ( info(ibatch) /= 0 .or. relative_residual > residual_limit ) &
      & failure_count = failure_count + 1
  EndDo

  Call gpu_finalize()

  If ( failure_count == 0 ) Then
    Write(output_unit,"(a)") 'XNET_GPU_LINALG status passed'
  Else
    Write(output_unit,"(a,i4)") 'XNET_GPU_LINALG status failed ',failure_count
    Stop 1
  EndIf

End Program frontier_gpu_linalg_probe
