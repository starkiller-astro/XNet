!***************************************************************************************************
! xnet_gpu.f90 10/18/17
! This file contains modules and subroutines to control the GPU execution of XNet.
!***************************************************************************************************

Module xnet_gpu
  Use, Intrinsic :: iso_c_binding
  Use xnet_controls, Only: lun_stdout, myid, tid
  Use xnet_types, Only: dp
  Use xnet_util, Only: xnet_terminate

#if defined(XNET_CUDA)
  Use cudaf, Only: &
    cuda_stream=>stream, &
    cudaGetDeviceCount, &
    cudaSetDevice, &
    cudaStreamCreate, &
    cudaStreamDestroy, &
    cudaStreamSynchronize, &
    cudaDeviceSynchronize, &
    cudaSuccess
  Use cublasf, Only: &
    cublas_handle, &
    cublasCreate_v2, &
    cublasDestroy_v2, &
    cublasGetStream_v2, &
    cublasSetStream_v2
  Use cusolverf, Only: &
    cusolver_handle, &
    cusolverDnCreate, &
    cusolverDnDestroy, &
    cusolverDnSetStream
  !Use cusparsef, Only: &
  !  cusparse_handle, &
  !  cusparseCreate, &
  !  cusparseDestroy, &
  !  cusparseSetStream
#elif defined(XNET_HIP)
  Use hipf, Only: &
    hip_stream=>stream, &
    hipGetDeviceCount, &
    hipSetDevice, &
    hipStreamCreate, &
    hipStreamDestroy &
    hipStreamSynchronize, &
    hipCheck, &
    hipblasCheck, &
    hipsparseCheck, &
    rocblasCheck, &
    rocsparseCheck, &
    rocsolverCheck
  Use rocblasf, Only: &
    rocblas_handle, &
    rocblas_create_handle, &
    rocblas_destroy_handle, &
    rocblas_get_stream, &
    rocblas_set_stream, &
    rocblas_set_pointer_mode, &
    rocblas_pointer_mode_device
  Use rocsolverf, Only: &
    rocsolver_handle
  !Use rocsparsef, Only: &
  !  rocsparse_handle
  Use hipblasf, Only: &
    hipblas_handle, &
    hipblasCreate, &
    hipblasDestroy, &
    hipblasGetStream, &
    hipblasSetStream
  !Use hipsparsef, Only: &
  !  hipsparse_handle, &
  !  hipsparseCreate, &
  !  hipsparseSetStream
#endif

#if defined(XNET_LA_MAGMA)
  Use magmaf, Only: &
    magma_device, &
    magma_queue, &
    magma_getdevice, &
    magma_init, &
#if defined(XNET_CUDA)
    magma_queue_create_from_cuda
#elif defined(XNET_HIP)
    magma_queue_create_from_hip
#endif
#endif

#if defined(XNET_OMP_OL)
  Use openmpf, Only: &
    omp_set_default_device, &
    omp_get_default_device, &
    omp_is_initial_device, &
    omp_target_is_present
  Use omp_lib, Only : &
    omp_get_num_devices
#endif

#if defined(XNET_OACC)
  Use openaccf, Only: &
    acc_set_device_num, &
    acc_get_device_num, &
    acc_on_device, &
    acc_is_present, &
    acc_get_cuda_stream, &
    acc_set_cuda_stream, &
    acc_get_default_async, &
    acc_set_default_async, &
    acc_device_host, &
    acc_device_default, &
    acc_async_sync, &
    acc_async_noval, &
    acc_queue
#endif

  Implicit None
  Private

  Integer(C_INT), Public :: deviceCount
  Integer, Public :: mydevice
  Type(C_PTR), Pointer, Public :: stream

  Interface dev_ptr
    Module Procedure dev_ptr_int
    Module Procedure dev_ptr_dp
    Module Procedure dev_ptr_cptr
  End Interface

  Public :: gpu_init
  Public :: gpu_finalize
  Public :: device_is_present
  Public :: get_device_num
  Public :: on_device
  Public :: stream_sync
  Public :: dev_ptr

Contains


  Subroutine gpu_init
    Implicit None

    ! Local variables
    Integer :: ierr

    ! Initialize GPU
#if defined(XNET_GPU)
#if defined(XNET_CUDA)
    ierr = cudaGetDeviceCount( deviceCount )
    If ( ierr /= cudaSuccess ) Write(lun_stdout,*) "cudaGetDeviceCount, ierr", ierr
#elif defined(XNET_HIP)
    Call hipCheck( hipGetDeviceCount( deviceCount ) )
#elif defined(XNET_OMP_OL)
    deviceCount = omp_get_num_devices()
#endif
    If ( deviceCount > 0 ) Then
      mydevice = mod(myid,deviceCount)
    Else
      Call xnet_terminate('No GPU found',deviceCount)
    EndIf
#if defined(XNET_CUDA)
    ierr = cudaSetDevice( mydevice )
#elif defined(XNET_HIP)
    Call hipCheck( hipSetDevice( mydevice ) )
#elif defined(XNET_OMP_OL)
    Call omp_set_default_device( mydevice )
#endif
#else
    mydevice = -1
    deviceCount = 0
#endif

    !$omp parallel default(shared) private(ierr)

    ! Setup linear algebra library handles
#if defined(XNET_CUDA)
    ierr = cublasCreate_v2( cublas_handle )
    !ierr = cusparseCreate( cusparse_handle )
    ierr = cusolverDnCreate( cusolver_handle )
    stream => cuda_stream
#elif defined(XNET_HIP)
    Call hipblasCheck( hipblasCreate( hipblas_handle ) )
    !Call hipsparseCheck( hipsparseCreate( hipsparse_handle ) )
    Call rocblasCheck( rocblas_create_handle( rocblas_handle ) )
    rocsolver_handle = rocblas_handle
    !rocsparse_handle = rocblas_handle
    stream => hip_stream
#endif

    ! Create a stream and associate with linear algebra libraries
#if defined(XNET_OACC)
    stream = acc_get_cuda_stream( INT( acc_async_noval, KIND=C_LONG_LONG ) )
    Call acc_set_cuda_stream( INT( acc_async_sync, KIND=C_LONG_LONG ), stream )
    !acc_queue = tid
    !stream = acc_get_cuda_stream(acc_queue)
    !ierr = acc_set_cuda_stream(acc_queue, stream)
#elif defined(XNET_CUDA)
    ierr = cudaStreamCreate( stream )
#elif defined(XNET_HIP)
    Call hipCheck( hipStreamCreate( stream ) )
#endif

#if defined(XNET_CUDA)
    ierr = cublasSetStream_v2( cublas_handle, stream )
    !ierr = cusparseSetStream( cusparse_handle, stream )
    ierr = cusolverDnSetStream( cusolver_handle, stream )
#elif defined(XNET_HIP)
    Call hipblasCheck( hipblasSetStream( hipblas_handle, stream ) )
    !Call hipsparseCheck( hipsparseSetStream( hipsparse_handle, stream ) )
    Call rocblasCheck( rocblas_set_stream( rocblas_handle, stream ) )
#endif

    !$omp end parallel

#if defined(XNET_LA_MAGMA)
    Call magma_get_device( magma_device )
    Call magma_init()
    !$omp parallel default(shared) private(ierr)
#if defined(XNET_CUDA)
    Call magma_queue_create_from_cuda &
           ( magma_device, stream, cublas_handle, C_NULL_PTR, magma_queue )
#elif defined(XNET_HIP)
    Call magma_queue_create_from_hip &
           ( magma_device, stream, hipblas_handle, C_NULL_PTR, magma_queue )
#endif
    !$omp end parallel
#endif

    Return
  End Subroutine gpu_init


  Subroutine gpu_finalize
    Implicit None

    ! Local variables
    Integer :: ierr

#if defined(XNET_GPU)
    !$omp parallel default(shared) private(ierr)
#if defined(XNET_LA_MAGMA)
    Call magma_queue_destroy( magma_queue )
#endif
#if defined(XNET_CUDA)
    ierr = cudaStreamDestroy( stream )
    ierr = cublasDestroy_v2( cublas_handle )
    ierr = cusolverDnDestroy( cusolver_handle )
#elif defined(XNET_HIP)
    Call hipCheck( hipStreamDestroy( stream ) )
    Call hipblasCheck( hipblasDestroy( hipblas_handle ) )
    Call rocblasCheck( rocblas_destroy_handle( rocsolver_handle )
#endif
    !$omp end parallel
#endif

    Return
  End Subroutine gpu_finalize


  Logical Function device_is_present( hostptr, device, bytes )
    Type(C_PTR), Intent(in) :: hostptr
    Integer, Intent(in) :: device
    Integer(C_SIZE_T), Intent(in) :: bytes
#if defined(XNET_OMP_OL)
    device_is_present = ( omp_target_is_present( hostptr, device ) > 0 )
#elif defined(XNET_OACC)
    device_is_present = ( acc_is_present( hostptr, bytes ) > 0 )
#else
    device_is_present = .false.
#endif
    Return
  End Function device_is_present


  Integer Function get_device_num()
#if defined(XNET_OMP_OL)
    get_device_num = omp_get_default_device()
#elif defined(XNET_OACC)
    get_device_num = acc_get_device_num( acc_device_default )
#else
    get_device_num = -1
#endif
    Return
  End Function get_device_num


  Logical Function on_device()
#if defined(XNET_OMP_OL)
    !$OMP DECLARE TARGET
    on_device = ( .not. omp_is_initial_device() )
#elif defined(XNET_OACC)
    !$ACC ROUTINE SEQ
    on_device = ( .not. acc_on_device( acc_device_host ) )
#else
    on_device = .false.
#endif
    Return
  End Function on_device


  Subroutine stream_sync( stream )
    Type(C_PTR), Intent(in) :: stream
    Integer :: ierr
#if defined(XNET_CUDA)
    ierr = cudaStreamSynchronize( stream )
#elif defined(XNET_HIP)
    Call hipCheck( hipStreamSynchronize( stream ) )
#endif
    Return
  End Subroutine stream_sync


  Type(C_PTR) Function dev_ptr_int( a )
#if defined(XNET_OMP_OL)
    Integer, Target, Intent(in) :: a
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(XNET_OACC)
    Integer, Target, Intent(in) :: a
    !$ACC HOST_DATA USE_DEVICE( a )
#else
    Integer, Target, Intent(in) :: a
#endif
    dev_ptr_int = C_LOC( a )
#if defined(XNET_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(XNET_OACC)
    !$ACC END HOST_DATA
#endif
  End Function dev_ptr_int

  Type(C_PTR) Function dev_ptr_dp( a )
#if defined(XNET_OMP_OL)
    Real(dp), Target, Intent(in) :: a
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(XNET_OACC)
    Real(dp), Target, Intent(in) :: a
    !$ACC HOST_DATA USE_DEVICE( a )
#else
    Real(dp), Target, Intent(in) :: a
#endif
    dev_ptr_dp = C_LOC( a )
#if defined(XNET_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(XNET_OACC)
    !$ACC END HOST_DATA
#endif
  End Function dev_ptr_dp

  Type(C_PTR) Function dev_ptr_cptr( a )
#if defined(XNET_OMP_OL)
    Type(C_PTR), Target, Intent(in) :: a
    !$OMP TARGET DATA USE_DEVICE_PTR( a )
#elif defined(XNET_OACC)
    Type(C_PTR), Target, Intent(in) :: a
    !$ACC HOST_DATA USE_DEVICE( a )
#else
    Type(C_PTR), Target, Intent(in) :: a
#endif
    dev_ptr_cptr = C_LOC( a )
#if defined(XNET_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(XNET_OACC)
    !$ACC END HOST_DATA
#endif
  End Function dev_ptr_cptr

End Module xnet_gpu
