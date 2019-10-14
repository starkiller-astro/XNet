!***************************************************************************************************
! xnet_gpu.f90 10/18/17
! This file contains modules and subroutines to control the GPU execution of XNet.
!***************************************************************************************************

Module xnet_gpu
  Use, Intrinsic :: iso_c_binding
  Use xnet_controls, Only: lun_stdout, myid, tid
#if defined(XNET_GPU)
  Use cudaf, Only: cudaGetDeviceCount, cudaSetDevice, cudaDeviceSynchronize, cudaSuccess, stream
#if defined(XNET_CUBLAS) || defined(XNET_MAGMA)
  Use cublasf, Only: cublasCreate_v2, cublasSetStream_v2, cublas_handle
#endif
#if defined(XNET_MAGMA)
  Use magmaf, Only: magma_getdevice, magma_init, magma_queue_create_from_cuda, magma_device, &
    & magma_queue
#endif
#if defined(XNET_OACC)
  Use openaccf, Only: acc_get_default_async, acc_set_default_async, acc_get_cuda_stream, &
    & acc_set_cuda_stream, acc_is_present, acc_queue
#endif
#if defined(XNET_OMP_OL)
  !Use openmpf, Only: omp_target_is_present
#endif
#endif

  Implicit None

  Integer(C_INT) :: deviceCount
  Integer :: mydevice

Contains

  Subroutine gpu_init
    Implicit None

    ! Local variables
    Integer :: istat

#if defined(XNET_GPU)
    ! Initialize GPU
    istat = cudaGetDeviceCount(deviceCount)
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaGetDeviceCount, istat", istat
    If ( deviceCount > 0 ) Then
      mydevice = mod(myid,deviceCount)
    Else
      Write(lun_stdout,*) 'No CUDA capable device found'
    EndIf
    istat = cudaSetDevice(mydevice)

    !$omp parallel default(shared) private(istat)

#if defined(XNET_CUBLAS) || defined(XNET_MAGMA)
    ! Create cublas handles
    istat = cublasCreate_v2(cublas_handle)
#endif

#if defined(XNET_OACC)
    ! Associate OpenACC async queue with CUDA stream
    acc_queue = tid
    stream = acc_get_cuda_stream(acc_queue)
#else
    ! Create CUDA streams
    istat = cudaStreamCreate(stream)
#endif

#if defined(XNET_CUBLAS) || defined(XNET_MAGMA)
    ! Associate each cublas cublas_handle with a CUDA stream
    istat = cublasSetStream_v2(cublas_handle, stream)
#endif

    !$omp end parallel

#if defined(XNET_MAGMA)
    call magma_getdevice( magma_device )
    call magma_init()
    !$omp parallel default(shared) private(istat)
    call magma_queue_create_from_cuda &
      ( magma_device, stream, cublas_handle, C_NULL_PTR, magma_queue )
    !$omp end parallel
#endif

    istat = cudaDeviceSynchronize()
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaDeviceSynchronize, istat", istat

#else
    mydevice = -1
    deviceCount = 0
#endif

    Return
  End Subroutine gpu_init

  Subroutine gpu_finalize
#if defined(XNET_GPU)
    Use cudaf, only: cudaStreamDestroy, stream
#if defined(XNET_CUBLAS) || defined(XNET_MAGMA)
    Use cublasf, Only: cublasDestroy_v2, cublas_handle
#endif
#endif
    Implicit None

    ! Local variables
    Integer :: istat

#if defined(XNET_GPU)
    !$omp parallel default(shared) private(istat)
    istat = cudaStreamDestroy(stream)
#if defined(XNET_CUBLAS) || defined(XNET_MAGMA)
    istat = cublasDestroy_v2(cublas_handle)
#endif
    !$omp end parallel
#endif

    Return
  End Subroutine gpu_finalize

  Logical Function device_is_present( hostptr, device, bytes )
    Type(C_PTR), intent(in) :: hostptr
    Integer, intent(in) :: device
    Integer(C_SIZE_T), intent(in) :: bytes
#if defined(XNET_OMP_OL)
    !device_is_present = ( omp_target_is_present( hostptr, device ) > 0 )
#elif defined(XNET_OACC)
    device_is_present = ( acc_is_present( hostptr, bytes ) > 0 )
#else
    device_is_present = .false.
#endif
    Return
  End Function device_is_present

End Module xnet_gpu
