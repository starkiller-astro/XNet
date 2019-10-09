!***************************************************************************************************
! xnet_gpu.f90 10/18/17
! This file contains modules and subroutines to control the GPU execution of XNet.
!***************************************************************************************************

Module xnet_gpu
  Use, Intrinsic :: iso_c_binding, Only: C_INT, C_PTR
  Implicit None

  ! CUDA/CUBLAS management
  Type(C_PTR) :: handle, stream
  !$omp threadprivate(handle,stream)

  Integer(C_INT) :: deviceCount
  Integer :: mydevice

Contains

  Subroutine gpu_init
    Use xnet_controls, Only: lun_stdout, myid, tid
    Use cublasf
    Use cudaf
    Use magmaf
    Use openaccf
    Implicit None

    ! Local variables
    Integer :: istat

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
    
    ! Create cublas handles
    istat = cublasCreate_v2(handle)

    ! Create CUDA streams
    !istat = cudaStreamCreate(stream)

    ! Associate OpenACC async queue with CUDA stream
    !acc_async_default = acc_get_default_async()
    !call acc_set_default_async(acc_async_default)
    acc_queue = tid
    stream = acc_get_cuda_stream(acc_queue)
    !istat = acc_set_cuda_stream(acc_queue, stream)

    ! Associate each cublas handle with a CUDA stream
    istat = cublasSetStream_v2(handle, stream)

    !$omp end parallel

    call magma_getdevice( magma_device )
    call magma_init()

    !$omp parallel default(shared) private(istat)

    call magma_queue_create_from_cuda &
      ( magma_device, stream, handle, C_NULL_PTR, magma_queue )

    !$omp end parallel

    istat = cudaDeviceSynchronize()
    If ( istat /= cudaSuccess ) Write(lun_stdout,*) "cudaDeviceSynchronize, istat", istat

    Return
  End Subroutine gpu_init

  Subroutine gpu_finalize
    Use cublasf
    Use cudaf
    Implicit None

    ! Local variables
    Integer :: istat

    !$omp parallel default(shared) private(istat)

    istat = cudaStreamDestroy(stream)
    istat = cublasDestroy_v2(handle)

    !$omp end parallel

    Return
  End Subroutine gpu_finalize

End Module xnet_gpu
