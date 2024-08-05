#ifdef XNET_DEBUG
#define XNET_DEBUG_LA
#endif
#if defined(XNET_LA_ONEMKL)
include "mkl_omp_offload.f90"
#endif

#include "xnet_macros.fh"

Module xnet_linalg

  Use, Intrinsic :: iso_c_binding
  Use xnet_types, Only: dp
  Use xnet_constants, Only: pi
  Use xnet_gpu, Only: &
    mydevice, &
    device_is_present, &
    dev_ptr, &
    stream_sync, &
    stream

#if defined(XNET_LA_CUBLAS)
  Use cublasf, Only: &
    cublas_handle, &
    cublasDnrm2_v2, &
    cublasDaxpy_v2, &
    cublasDgemm_v2, &
    cublasDgemmStridedBatched, &
    cublasDgetrfBatched, &
    cublasDgetrsBatched, &
    cublasDgemv_v2, &
    cublasDtrsv_v2, &
    cublasDtrsm_v2, &
    cublasDgeam, &
    cublasDdgmm, &
    CUBLAS_OP_N, CUBLAS_OP_T, &
    CUBLAS_SIDE_LEFT, &
    CUBLAS_FILL_MODE_UPPER, &
    CUBLAS_DIAG_NON_UNIT
  Use cusolverf, Only: &
    cusolver_handle, &
    cusolverDnDgeqrf_bufferSize, &
    cusolverDnDgeqrf, &
    cusolverDnDormqr, &
    cusolverDnDgetrf_bufferSize, &
    cusolverDnDgetrf, &
    cusolverDnDgetrs
  Use cusparsef, Only: &
    cusparse_handle, &
    cusparseDgthr, &
    CUSPARSE_INDEX_BASE_ONE
#elif defined(XNET_LA_ROCM)
  Use hipf, Only: &
    hipCheck, &
    hipblasCheck, &
    hipsparseCheck, &
    rocblasCheck, &
    rocsparseCheck, &
    rocsolverCheck
  Use rocblasf, Only: &
    rocblas_handle, &
    rocblas_dnrm2, &
    rocblas_daxpy, &
    rocblas_dgemm, &
    rocblas_dgemm_strided_batched, &
    rocblas_dgemv, &
    rocblas_dtrsv, &
    rocblas_dtrsm, &
    rocblas_dgeam, &
    rocblas_ddgmm, &
    rocblas_ddot_strided_batched, &
    rocblas_operation_none, &
    rocblas_operation_transpose, &
    rocblas_side_left, &
    rocblas_fill_upper, &
    rocblas_diagonal_non_unit
  Use rocsolverf, Only: &
    rocsolver_handle, &
    rocsolver_dgeqrf, &
    rocsolver_dormqr, &
    rocsolver_dgetrf, &
    rocsolver_dgetrf_batched, &
    rocsolver_dgetrs, &
    rocsolver_dgetrs_batched
  Use rocsparsef, Only: &
    rocsparse_handle, &
    rocsparse_dgthr, &
    rocsparse_index_base_one
  Use hipblasf, Only: &
    hipblas_handle, &
    hipblasDnrm2, &
    hipblasDaxpy, &
    hipblasDgemm, &
    hipblasDgemmStridedBatched, &
    hipblasDgetrf, &
    hipblasDgetrfBatched, &
    hipblasDgetrs, &
    hipblasDgetrsBatched, &
    hipblasDgemv, &
    hipblasDtrsv, &
    hipblasDtrsm, &
    hipblasDgeam, &
    hipblasDdgmm, &
    HIPBLAS_OP_N, HIPBLAS_OP_T, &
    HIPBLAS_SIDE_LEFT, &
    HIPBLAS_FILL_MODE_UPPER, &
    HIPBLAS_DIAG_NON_UNIT
  Use hipsparsef, Only: &
    hipsparse_handle, &
    hipsparseDgthr, &
    HIPSPARSE_INDEX_BASE_ONE
#elif defined(XNET_LA_ONEMKL)
  Use onemkl_blas_omp_offload_lp64
#elif defined(XNET_LA_MAGMA)
  Use magmaf, Only: &
    magma_queue, &
    magma_dnrm2, &
    magma_daxpy, &
    magma_dgemm, &
    magmablas_dgemm_batched_strided, &
    magma_dgetrf_native, &
    magma_dgetrf_batched, &
    magma_dgetrs_gpu, &
    magma_dgetrs_batched, &
    magma_dgemv, &
    magma_dgels_gpu, &
    magmablas_dtranspose, &
    magmablas_dlacpy, &
    magmablas_dlascl2, &
    magmablas_dgeadd2, &
    magma_dmalloc, &
    MagmaNoTrans, MagmaTrans, &
    MagmaUpper, MagmaLower, MagmaFull, MagmaGeneral, &
    MagmaUnit, MagmaNonUnit, &
    MagmaLeft, MagmaRight
#endif

  Implicit None
  Private

  Public :: MatrixMatrixAdd
  Public :: MatrixMatrixMultiply
  Public :: MatrixMatrixMultiplyBatched
  Public :: MatrixVectorMultiply
  Public :: MatrixDiagScale
  Public :: VectorDotProductBatched
  Public :: VectorNorm2
  Public :: VectorNorm2_Kernel
  Public :: VectorVectorAdd
  Public :: LinearLeastSquares_LWORK
  Public :: LinearLeastSquares
  Public :: EigenvaluesSymmetric3

  Public :: LinearSolve
  Public :: LinearSolve_CPU
  Public :: LinearSolve_GPU
  Public :: LinearSolveBatched
  Public :: LinearSolveBatched_CPU
  Public :: LinearSolveBatched_GPU

  Public :: LUDecomp_LWORK
  Public :: LUDecomp
  Public :: LUDecomp_CPU
  Public :: LUDecomp_GPU
  !Public :: LUDecompBatched
  Public :: LUDecompBatched_CPU
  Public :: LUDecompBatched_GPU

  Public :: LUBksub
  Public :: LUBksub_CPU
  Public :: LUBksub_GPU
  !Public :: LUBksubBatched
  Public :: LUBksubBatched_CPU
  Public :: LUBksubBatched_GPU


Contains

  Integer Function itrans_from_char( ctrans )
    Character, Intent(in) :: ctrans
    itrans_from_char = 0
#if defined(XNET_LA_CUBLAS)
    If ( ctrans == 'T' ) Then
      itrans_from_char = CUBLAS_OP_T
    Else
      itrans_from_char = CUBLAS_OP_N
    End If
#elif defined(XNET_LA_ROCM)
    If ( ctrans == 'T' ) Then
      !itrans_from_char = rocblas_operation_transpose
      itrans_from_char = HIPBLAS_OP_T
    Else
      !itrans_from_char = rocblas_operation_none
      itrans_from_char = HIPBLAS_OP_N
    End If
#elif defined(XNET_LA_ONEMKL)
#elif defined(XNET_LA_MAGMA)
    If ( ctrans == 'T' ) Then
      itrans_from_char = MagmaTrans
    Else
      itrans_from_char = MagmaNoTrans
    End If
#endif
    Return
  End Function itrans_from_char


  Subroutine MatrixMatrixAdd( transa, transb, m, n, alpha, a, lda, beta, b, ldb, c, ldc )

    Character                          :: transa, transb
    Integer                            :: m, n, lda, ldb, ldc
    Real(dp)                           :: alpha, beta
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Real(dp), Dimension(ldc,*), Target :: c

    Integer                            :: i, j, ierr, info
    Integer(C_INT)                     :: itransa, itransb
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c, mn
    Real(dp), Dimension(:,:), Pointer  :: pa, pb, pc
    Type(C_PTR)                        :: ha, hb, hc
    Type(C_PTR)                        :: da, db, dc
    Type(C_PTR)                        :: dat, dbt
    Integer                            :: ka, kb
    Logical                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = m * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)
    mn = m * n

    If ( transa == 'N' ) Then
      ka = n
    Else
      ka = m
    End If

    If ( transb == 'N' ) Then
      kb = n
    Else
      kb = m
    End If

    pa => a(:,1:ka)
    pb => b(:,1:kb)
    pc => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDgeam &
             ( cublas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc )
#elif defined(XNET_LA_ROCM)
      !Call rocblasCheck( rocblas_dgeam &
      !       ( rocblas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc ) )
      Call hipblasCheck( hipblasDgeam &
             ( hipblas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc ) )
#elif defined(XNET_LA_ONEMKL)
      !!$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      !Call DGEAM &
      !       ( transa, transb, m, n, alpha, a, lda, beta, b, ldb, c, ldc )
      !!$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
      If ( transb  == 'N' ) Then
        Call magmablas_dlacpy &
               ( MagmaFull, m, n, db, ldb, dc, ldc, magma_queue )
      Else
        Call magmablas_dtranspose &
               ( n, m, db, ldb, dc, ldc, magma_queue )
      End If
      If ( transa == 'N' ) Then
        Call magmablas_dgeadd2 &
               ( m, n, alpha, da, lda, beta, dc, ldc, magma_queue )
      Else
        ierr = magma_dmalloc( dat, mn )
        Call magmablas_dtranspose &
               ( n, m, da, lda, dat, m, magma_queue )
        Call magmablas_dgeadd2 &
               ( m, n, alpha, dat, m, beta, dc, ldc, magma_queue )
      End If
#endif
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixAdd] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixAdd]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixAdd]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixAdd]   C missing'
#endif
#endif

      If ( alpha == 0.0_DP .AND. beta == 0.0_DP ) Then

        Do j = 1, n
          Do i = 1, m
            c(i,j) = 0.0_DP
          End Do
        End Do

      Else If ( alpha == 0.0_DP ) Then

        If ( transb == 'N' ) Then
          Do j = 1, n
            Do i = 1, m
              c(i,j) = + beta * b(i,j)
            End Do
          End Do
        Else
          Do j = 1, n
            Do i = 1, m
              c(i,j) = + beta * b(j,i)
            End Do
          End Do
        End If

      Else If ( beta == 0.0_DP ) Then

        If ( transa == 'N' ) Then
          Do j = 1, n
            Do i = 1, m
              c(i,j) = alpha * a(i,j)
            End Do
          End Do
        Else
          Do j = 1, n
            Do i = 1, m
              c(i,j) = alpha * a(j,i)
            End Do
          End Do
        End If

      Else

        If ( transa == 'N' ) Then
          If ( transb == 'N' ) Then
            Do j = 1, n
              Do i = 1, m
                c(i,j) = alpha * a(i,j) + beta * b(i,j)
              End Do
            End Do
          Else
            Do j = 1, n
              Do i = 1, m
                c(i,j) = alpha * a(i,j) + beta * b(j,i)
              End Do
            End Do
          End If
        Else
          If ( transb == 'N' ) Then
            Do j = 1, n
              Do i = 1, m
                c(i,j) = alpha * a(j,i) + beta * b(i,j)
              End Do
            End Do
          Else
            Do j = 1, n
              Do i = 1, m
                c(i,j) = alpha * a(j,i) + beta * b(j,i)
              End Do
            End Do
          End If
        End If

      End If

    End If

  End Subroutine MatrixMatrixAdd


  Subroutine MatrixMatrixMultiply( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    Character                          :: transa, transb
    Integer                            :: m, n, k, lda, ldb, ldc
    Real(dp)                           :: alpha, beta
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Real(dp), Dimension(ldc,*), Target :: c

    Integer                            :: ierr
    Integer(C_INT)                     :: itransa, itransb
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    Real(dp), Dimension(:,:), Pointer  :: pa, pb, pc
    Type(C_PTR)                        :: ha, hb, hc
    Type(C_PTR)                        :: da, db, dc
    Integer                            :: ka, kb
    Logical                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * k * c_sizeof(0.0_DP)
    sizeof_b = k * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)

    If ( transa == 'N' ) Then
      ka = k
    Else
      ka = m
    End If

    If ( transb == 'N' ) Then
      kb = n
    Else
      kb = k
    End If

    pa => a(:,1:ka)
    pb => b(:,1:kb)
    pc => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDgemm_v2 &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc )
#elif defined(XNET_LA_ROCM)
      !Call rocblasCheck( rocblas_dgemm &
      !       ( rocblas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc ) )
      Call hipblasCheck( hipblasDgemm &
             ( hipblas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc ) )
#elif defined(XNET_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      Call DGEMM &
             ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
      Call magma_dgemm &
             ( itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, magma_queue )
#endif
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixMultiply] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   C missing'
#endif
#endif

      Call DGEMM &
             ( transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc )

    End If

  End Subroutine MatrixMatrixMultiply


  Subroutine MatrixMatrixMultiplyBatched( transa, transb, m, n, k, alpha, a, lda, stridea, &
                                          b, ldb, strideb, beta, c, ldc, stridec, batchcount )

    Character                          :: transa, transb
    Integer                            :: m, n, k, lda, ldb, ldc, stridea, strideb, stridec, batchcount
    Real(dp)                           :: alpha, beta
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Real(dp), Dimension(ldc,*), Target :: c

    Integer                            :: ierr, i
    Integer(C_INT)                     :: itransa, itransb
    Integer(C_INT64_T)                 :: stridea_64, strideb_64, stridec_64
    Integer(C_LONG_LONG)               :: stridea_l, strideb_l, stridec_l
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    Real(dp), Dimension(:,:), Pointer  :: pa, pb, pc
    Type(C_PTR)                        :: ha, hb, hc
    Type(C_PTR)                        :: da, db, dc
    Integer                            :: ka, kb
    Integer                            :: osa, osb, osc
    Integer                            :: ia, ib, ic
    Integer                            :: ja, jb, jc
    Logical                            :: data_on_device

    data_on_device = .false.
    sizeof_a = m * k * batchcount * c_sizeof(0.0_DP)
    sizeof_b = k * n * batchcount * c_sizeof(0.0_DP)
    sizeof_c = m * n * batchcount * c_sizeof(0.0_DP)

    If ( transa == 'N' ) Then
      ka = k
    Else
      ka = m
    End If

    If ( transb == 'N' ) Then
      kb = n
    Else
      kb = k
    End If

    pa => a(:,1:ka*batchcount)
    pb => b(:,1:kb*batchcount)
    pc => c(:,1:n *batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dc = dev_ptr( pc(1,1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDgemmStridedBatched &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount )
#elif defined(XNET_LA_ROCM)
      !stridea_64 = stridea
      !strideb_64 = strideb
      !stridec_64 = stridec
      !Call rocblasCheck( rocblas_dgemm_strided_batched &
      !       ( rocblas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea_64, &
      !         db, ldb, strideb_64, beta, dc, ldc, stridec_64, batchcount ) )
      stridea_l = stridea
      strideb_l = strideb
      stridec_l = stridec
      Call hipblasCheck( hipblasDgemmStridedBatched &
             ( hipblas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea_l, &
               db, ldb, strideb_l, beta, dc, ldc, stridec_l, batchcount ) )
#elif defined(XNET_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, c )
      Call DGEMM_BATCH_STRIDED &
             ( transa, transb, m, n, k, alpha, a, lda, stridea, b, ldb, strideb, beta, c, ldc, stridec, batchcount )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
      Call magmablas_dgemm_batched_strided &
             ( itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount, magma_queue )
#endif
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixMultiplyBatched] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   C missing'
#endif
#endif

      Do i = 1, batchcount
        osa = (i-1) * stridea + 1
        osb = (i-1) * strideb + 1
        osc = (i-1) * stridec + 1
        ia = mod( (osa-1), lda ) + 1
        ib = mod( (osb-1), ldb ) + 1
        ic = mod( (osc-1), ldc ) + 1
        ja = 1
        jb = 1
        jc = 1
        If ( stridea /= 0 ) ja = mod( (osa-1)/lda, stridea ) + 1
        If ( strideb /= 0 ) jb = mod( (osb-1)/ldb, strideb ) + 1
        If ( stridec /= 0 ) jc = mod( (osc-1)/ldc, stridec ) + 1
        Call DGEMM &
               ( transa, transb, m, n, k, alpha, a(ia,ja), lda, b(ib,jb), ldb, beta, c(ic,jc), ldc )
      End Do

    End If

  End Subroutine MatrixMatrixMultiplyBatched


  Subroutine MatrixVectorMultiply( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    Character                          :: trans
    Integer                            :: m, n, lda, incx, incy
    Real(dp)                           :: alpha, beta
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(*)    , Target :: x
    Real(dp), Dimension(*)    , Target :: y

    Integer                            :: ierr, lenx, leny
    Integer(C_INT)                     :: itrans
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_x, sizeof_y
    Real(dp), Dimension(:,:), Pointer  :: pa
    Real(dp), Dimension(:)  , Pointer  :: px, py
    Type(C_PTR)                        :: ha, hx, hy
    Type(C_PTR)                        :: da, dx, dy
    Logical                            :: data_on_device

    data_on_device = .false.

    If ( trans == 'T' ) Then
      lenx = m
      leny = n
    Else
      lenx = n
      leny = m
    End If

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_x = lenx * c_sizeof(0.0_DP)
    sizeof_y = leny * c_sizeof(0.0_DP)

    pa(1:lda,1:n) => a(:,1:n)
    px(1:lenx) => x(1:lenx)
    py(1:leny) => y(1:leny)

    ha = C_LOC( pa )
    hx = C_LOC( px )
    hy = C_LOC( py )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_y )

    If ( data_on_device ) Then

      itrans = itrans_from_char( trans )

      da = dev_ptr( pa(1,1) )
      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDgemv_v2 &
             ( cublas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy )
#elif defined(XNET_LA_MAGMA)
      Call magma_dgemv &
             ( itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy, magma_queue )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[MatrixVectorMultiply] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixVectorMultiply]   A missing'
      If ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        Write(*,*) '[MatrixVectorMultiply]   x missing'
      If ( .not. device_is_present( hy, mydevice, sizeof_y ) ) &
        Write(*,*) '[MatrixVectorMultiply]   y missing'
#endif

      Call DGEMV &
             ( trans, m, n, alpha, a, lda, x, incx, beta, y, incy )

    End If

  End Subroutine MatrixVectorMultiply


  Subroutine VectorDotProductBatched( n, x, incx, stridex, y, incy, stridey, batchcount, xy )

    Integer                         :: n, incx, incy, stridex, stridey, batchcount
    Real(dp), Dimension(*), Target  :: x, y, xy

    Integer                         :: ierr, i
    Integer(C_INT64_T)              :: stridex_64, stridey_64
    Integer(C_SIZE_T)               :: sizeof_x, sizeof_y, sizeof_xy
    Real(dp), Dimension(:), Pointer :: px, py, pxy
    Type(C_PTR)                     :: hx, hy, hxy
    Type(C_PTR)                     :: dx, dy, dxy
    Integer                         :: osx, osy
    Logical                         :: data_on_device
    Real(dp), EXTERNAL              :: DDOT

    data_on_device = .false.
    sizeof_x  = n * batchcount * c_sizeof(0.0_DP)
    sizeof_y  = n * batchcount * c_sizeof(0.0_DP)
    sizeof_xy = batchcount * c_sizeof(0.0_DP)

    px => x(1:n*batchcount)
    py => y(1:n*batchcount)
    pxy => xy(1:batchcount)

    hx = C_LOC( px )
    hy = C_LOC( py )
    hxy = C_LOC( pxy )

    data_on_device = device_is_present( hx,  mydevice, sizeof_x  ) &
               .AND. device_is_present( hy,  mydevice, sizeof_y  ) &
               .AND. device_is_present( hxy, mydevice, sizeof_xy )

    If ( data_on_device ) Then

      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )
      dxy = dev_ptr( pxy(1) )

#if defined(XNET_LA_CUBLAS)
      ! Currently unavailable
      !ierr = cublasDdot_v2( cublas_handle, n, dx, incx, xnorm )
#elif defined(XNET_LA_ROCM)
      ! Currently unavailable
      stridex_64 = stridex
      stridey_64 = stridey
      Call rocblasCheck( rocblas_ddot_strided_batched &
             ( rocblas_handle, n, dx, incx, stridex_64, dy, incy, stridey_64, batchcount, hxy ) )
#elif defined(XNET_LA_ONEMKL)
#elif defined(XNET_LA_MAGMA)
      ! Currently unavailable
      !xnorm = magma_ddot( n, dx, incx, magma_queue )
#endif
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[VectorDotProductBatched] Data not present on device'
      If ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        Write(*,*) '[VectorDotProductBatched]   x missing'
      If ( .not. device_is_present( hy, mydevice, sizeof_y ) ) &
        Write(*,*) '[VectorDotProductBatched]   y missing'
      If ( .not. device_is_present( hxy, mydevice, sizeof_xy ) ) &
        Write(*,*) '[VectorDotProductBatched]  xy missing'
#endif
#endif

      Do i = 1, batchcount
        osx = (i-1) * n + 1
        osy = (i-1) * n + 1
        xy(i) = DDOT( n, x(osx), incx, y(osy), incy )
      End Do

    End If

  End Subroutine VectorDotProductBatched


  Subroutine MatrixDiagScale( m, n, a, lda, x, incx, c, ldc )

    Character                          :: trans
    Integer                            :: m, n, lda, incx, ldc
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldc,*), Target :: c
    Real(dp), Dimension(*)    , Target :: x

    Integer                            :: ierr, info, i, j, ix
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_c, sizeof_x
    Real(dp), Dimension(:,:), Pointer  :: pa, pc
    Real(dp), Dimension(:)  , Pointer  :: px
    Type(C_PTR)                        :: ha, hx, hc
    Type(C_PTR)                        :: da, dx, dc
    Logical                            :: data_on_device

    data_on_device = .false.

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_c = m * n * c_sizeof(0.0_DP)
    sizeof_x = m * c_sizeof(0.0_DP)

    pa(1:lda,1:n) => a(:,1:n)
    pc(1:ldc,1:n) => c(:,1:n)
    px(1:m) => x(1:m)

    ha = C_LOC( pa )
    hc = C_LOC( pc )
    hx = C_LOC( px )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hc, mydevice, sizeof_c ) &
               .AND. device_is_present( hx, mydevice, sizeof_x )

    If ( data_on_device ) Then

      da = dev_ptr( pa(1,1) )
      dc = dev_ptr( pc(1,1) )
      dx = dev_ptr( px(1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDdgmm &
             ( cublas_handle, CUBLAS_SIDE_LEFT, m, n, da, lda, dx, incx, dc, ldc )
#elif defined(XNET_LA_MAGMA)
      Call magmablas_dlacpy &
             ( MagmaGeneral, m, n, da, lda, dc, ldc, magma_queue )
      Call magmablas_dlascl2 &
             ( MagmaGeneral, m, n, dx, dc, ldc, magma_queue, info )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[MatrixDiagScale] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixDiagScale]   A missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixDiagScale]   C missing'
      If ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        Write(*,*) '[MatrixDiagScale]   x missing'
#endif

      If ( incx == 1 ) Then
        Do j = 1, n
          Do i = 1, m
            c(i,j) = a(i,j) * x(i)
          End Do
        End Do
      Else
        Do j = 1, n
          ix = 1
          If ( incx < 0 ) ix = (-m+1)*incx + 1
          Do i = 1, m
            c(i,j) = a(i,j) * x(ix)
            ix = ix + incx
          End Do
        End Do
      End If

    End If

  End Subroutine MatrixDiagScale


  Subroutine VectorNorm2( n, x, incx, xnorm )

    Integer                         :: n, incx
    Real(dp), Dimension(*), Target  :: x
    Real(dp)                        :: xnorm

    Integer                         :: ierr
    Integer(C_SIZE_T)               :: sizeof_x
    Real(dp), Dimension(:), Pointer :: px
    Type(C_PTR)                     :: hx
    Type(C_PTR)                     :: dx
    Logical                         :: data_on_device
    Real(dp), External              :: DNRM2

    data_on_device = .false.
    sizeof_x = n * c_sizeof(0.0_DP)

    px(1:n) => x(1:n)

    hx = C_LOC( px )

    data_on_device = device_is_present( hx, mydevice, sizeof_x )

    If ( data_on_device ) Then

      dx = dev_ptr( px(1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDnrm2_v2( cublas_handle, n, dx, incx, xnorm )
#elif defined(XNET_LA_MAGMA)
      xnorm = magma_dnrm2( n, dx, incx, magma_queue )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[VectorNorm2] Data not present on device'
      If ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        Write(*,*) '[VectorNorm2]   x missing'
#endif

      xnorm = DNRM2( n, x, incx )

    End If

  End Subroutine VectorNorm2


  Subroutine VectorNorm2_Kernel( n, x, incx, xnorm )
    !__dir_routine_seq

    Integer                         :: n, incx
    Real(dp), Dimension(*), Target  :: x
    Real(dp)                        :: xnorm

    Integer                         :: ix
    Real(dp)                        :: xscale, xssq, absxi

    If ( n < 1 .OR. incx < 1 ) Then
      xnorm = 0.0d0
    Else If ( n == 1 ) Then
      xnorm = ABS( x(1) )
    Else
      xscale = 0.0d0
      xssq = 1.0d0
      Do ix = 1, 1 + (n-1)*incx, incx
        If ( x(ix) /= 0.0d0 ) Then
          absxi = ABS( x(ix) )
          If ( xscale < absxi ) Then
            xssq = 1.0d0 + xssq * (xscale/absxi)**2
            xscale = absxi
          Else
            xssq = xssq + (absxi/xscale)**2
          End If
        End If
      End Do
      xnorm = xscale * SQRT(xssq)
    End If

  End Subroutine VectorNorm2_Kernel


  Subroutine VectorVectorAdd( n, alpha, x, incx, y, incy )

    Integer                         :: n, incx, incy
    Real(dp)                        :: alpha
    Real(dp), Dimension(*), Target  :: x, y

    Integer                         :: ierr
    Integer(C_SIZE_T)               :: sizeof_x, sizeof_y
    Real(dp), Dimension(:), Pointer :: px, py
    Type(C_PTR)                     :: hx, hy
    Type(C_PTR)                     :: dx, dy
    Logical                         :: data_on_device

    data_on_device = .false.

    sizeof_x = n * c_sizeof(0.0_DP)
    sizeof_y = n * c_sizeof(0.0_DP)

    px(1:n) => x(1:n)
    py(1:n) => y(1:n)

    hx = C_LOC( px )
    hy = C_LOC( py )

    data_on_device = device_is_present( hx, mydevice, sizeof_x ) &
               .AND. device_is_present( hy, mydevice, sizeof_y )

    If ( data_on_device ) Then

      dx = dev_ptr( px(1) )
      dy = dev_ptr( py(1) )

#if defined(XNET_LA_CUBLAS)
      ierr = cublasDaxpy_v2( cublas_handle, n, alpha, dx, incx, dy, incy )
#elif defined(XNET_LA_MAGMA)
      Call magma_daxpy( n, alpha, dx, incx, dy, incy, magma_queue )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[VectorVectorAdd] Data not present on device'
      If ( .not. device_is_present( hx, mydevice, sizeof_x ) ) &
        Write(*,*) '[VectorVectorAdd]   x missing'
#endif

      Call DAXPY( n, alpha, x, incx, y, incy )

    End If

  End Subroutine VectorVectorAdd


  Subroutine LinearLeastSquares_LWORK( trans, m, n, nrhs, a, lda, b, ldb, work, lwork )

    Character                          :: trans
    Integer                            :: m, n, nrhs, lda, ldb, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Real(dp), Dimension(*)    , Target :: work

    Integer                            :: ierr, info, max_mn
    Integer(C_INT)                     :: itrans
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_work
    Real(dp), Dimension(:,:), Pointer  :: pa, pb
    Real(dp), Dimension(:)  , Pointer  :: pwork
    Type(C_PTR)                        :: ha, hb, hwork
    Type(C_PTR)                        :: da, db

    lwork = -1

    max_mn = MAX(m,n)

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = max_mn * nrhs * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    pwork => work(1:lwork)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hwork = C_LOC( pwork )

    da = dev_ptr( pa(1,1) )
    db = dev_ptr( pb(1,1) )

    itrans = itrans_from_char( trans )

#if defined(XNET_LA_CUBLAS)
    ierr = cusolverDnDgeqrf_bufferSize &
           ( cusolver_handle, m, n, da, lda, lwork )
#elif defined(XNET_LA_ROCM)
#elif defined(XNET_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, work )
    Call DGELS &
           ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
    !$OMP END TARGET VARIANT DISPATCH
    !$OMP TARGET UPDATE FROM( work(1) )
    lwork = INT( work(1) )
#elif defined(XNET_LA_MAGMA)
    Call magma_dgels_gpu &
           ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info )
    lwork = INT( work(1) )
#else
    Call DGELS &
           ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
    lwork = INT( work(1) )
#endif

  End Subroutine LinearLeastSquares_LWORK


  Subroutine LinearLeastSquares( trans, m, n, nrhs, a, lda, b, ldb, tau, work, lwork, info )

    Character                          :: trans
    Integer                            :: m, n, nrhs, lda, ldb, lwork
    Integer                   , Target :: info
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Real(dp), Dimension(*)    , Target :: tau
    Real(dp), Dimension(*)    , Target :: work

    Integer                            :: ierr, max_mn, min_mn
    Integer(C_INT)                     :: itrans
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_tau, sizeof_work, sizeof_info
    Real(dp), Dimension(:,:), Pointer  :: pa, pb
    Real(dp), Dimension(:)  , Pointer  :: ptau, pwork
    Integer                 , Pointer  :: pinfo
    Type(C_PTR)                        :: ha, hb, htau, hwork, hinfo
    Type(C_PTR)                        :: da, db, dtau, dwork, dinfo
    Logical                            :: data_on_device

    data_on_device = .false.
    max_mn = MAX(m,n)
    min_mn = MIN(m,n)

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_b = max_mn * nrhs * c_sizeof(0.0_DP)
    sizeof_tau = min_mn * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)
    sizeof_info = c_sizeof(info)

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    ptau => tau(1:min_mn)
    pwork => work(1:lwork)
    pinfo => info

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    htau = C_LOC( ptau )
    hwork = C_LOC( pwork )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( htau,  mydevice, sizeof_tau  ) &
               .AND. device_is_present( hwork, mydevice, sizeof_work ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

      itrans = itrans_from_char( trans )

      da = dev_ptr( pa(1,1) )
      db = dev_ptr( pb(1,1) )
      dtau = dev_ptr( ptau(1) )
      dwork = dev_ptr( pwork(1) )
      dinfo = dev_ptr( pinfo )

#if defined(XNET_LA_CUBLAS)
      ierr = cusolverDnDgeqrf &
             ( cusolver_handle, m, n, da, lda, dtau, dwork, lwork, dinfo )
      ierr = cusolverDnDormqr &
             ( cusolver_handle, &
               CUBLAS_SIDE_LEFT, CUBLAS_OP_T, &
               m, nrhs, n, da, lda, dtau, db, ldb, dwork, lwork, dinfo )

      If ( nrhs == 1 ) Then

        ierr = cublasDtrsv_v2 &
               ( cublas_handle, &
                 CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, &
                 n, da, lda, db, 1 )

      Else

        ierr = cublasDtrsm_v2 &
               ( cublas_handle, &
                 CUBLAS_SIDE_LEFT, CUBLAS_FILL_MODE_UPPER, CUBLAS_OP_N, CUBLAS_DIAG_NON_UNIT, &
                 n, nrhs, 1.0_dp, da, lda, db, ldb )

      End If
#elif defined(XNET_LA_ROCM)
      Call rocsolverCheck( rocsolver_dgeqrf &
             ( rocsolver_handle, m, n, da, lda, dtau ) )
      Call rocsolverCheck( rocsolver_dormqr &
             ( rocsolver_handle, &
               rocblas_side_left, rocblas_operation_transpose, &
               m, nrhs, n, da, lda, dtau, db, ldb ) )

      If ( nrhs == 1 ) Then

        !Call rocblasCheck( rocblas_dtrsv &
        !       ( rocblas_handle, &
        !         rocblas_fill_upper, rocblas_operation_none, rocblas_diagonal_non_unit, &
        !         n, da, lda, db, 1 ) )
        Call hipblasCheck( hipblasDtrsv &
               ( hipblas_handle, &
                 HIPBLAS_FILL_MODE_UPPER, HIPBLAS_OP_N, HIPBLAS_DIAG_NON_UNIT, &
                 n, da, lda, db, 1 ) )

      Else

        !Call rocblasCheck( rocblas_dtrsm &
        !       ( rocblas_handle, &
        !         rocblas_side_left, rocblas_fill_upper, &
        !         rocblas_operation_none, rocblas_diagonal_non_unit, &
        !         n, nrhs, 1.0_dp, da, lda, db, ldb ) )
        Call hipblasCheck( hipblasDtrsm &
               ( hipblas_handle, &
                 HIPBLAS_SIDE_LEFT, HIPBLAS_FILL_MODE_UPPER, HIPBLAS_OP_N, HIPBLAS_DIAG_NON_UNIT, &
                 n, nrhs, 1.0_dp, da, lda, db, ldb ) )

      End If
#elif defined(XNET_LA_ONEMKL)
      !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, work )
      Call DGELS &
             ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )
      !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
      Call magma_dgels_gpu &
             ( itrans, m, n, nrhs, da, lda, db, ldb, hwork, lwork, info )
#endif
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[LinearLeastSquares] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[LinearLeastSquares]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[LinearLeastSquares]   b missing'
      If ( .not. device_is_present( htau, mydevice, sizeof_tau ) ) &
        Write(*,*) '[LinearLeastSquares]   tau missing'
      If ( .not. device_is_present( hwork, mydevice, sizeof_work ) ) &
        Write(*,*) '[LinearLeastSquares]   work missing'
      If ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        Write(*,*) '[LinearLeastSquares]   info missing'
#endif
#endif

      Call DGELS &
             ( trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info )

    End If

  End Subroutine LinearLeastSquares


  Subroutine EigenvaluesSymmetric3( A, Lambda )
    !__dir_routine_seq

    Real(dp), Intent(in)  :: A(3,3)
    Real(dp), Intent(out) :: Lambda(3)

    Real(dp) :: B11, B22, B33
    Real(dp) :: B12, B13, B21, B23, B31, B32
    Real(dp) :: P1, P2, P, Q, R, PHI, DETB

    P1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2

    If ( P1 == 0.0_DP ) Then

      Lambda(1) = A(1,1)
      Lambda(2) = A(2,2)
      Lambda(3) = A(3,3)

    Else

      Q = ( A(1,1) + A(2,2) + A(3,3) ) / 3.0_DP
      P2 = 2.0_DP * P1 &
           + ( A(1,1) - Q )**2 &
           + ( A(2,2) - Q )**2 &
           + ( A(3,3) - Q )**2
      P = SQRT( P2 / 6.0_DP )

      B11 = ( A(1,1) - Q ) / P
      B22 = ( A(2,2) - Q ) / P
      B33 = ( A(3,3) - Q ) / P
      B12 = A(1,2) / P ; B21 = B12
      B13 = A(1,3) / P ; B31 = B13
      B23 = A(2,3) / P ; B32 = B23
      DETB =   B11 * B22 * B33  &
             - B11 * B23 * B32  &
             - B12 * B21 * B33  &
             + B12 * B23 * B31  &
             + B13 * B21 * B32  &
             - B13 * B22 * B31
      R = DETB * 0.5_DP
      If ( R <= - 1.0_DP ) Then
        PHI = pi
      Else If ( R >= 1.0_DP ) Then
        PHI = 0.0_DP
      Else
        PHI = ACOS( R ) / 3.0_DP
      End If

      Lambda(1) = Q + 2.0_DP * P * COS( PHI )
      Lambda(3) = Q + 2.0_DP * P * COS( PHI + ( 2.0_DP * pi / 3.0_DP ) )
      Lambda(2) = 3.0_DP * Q - Lambda(1) - Lambda(3)

    End If

  End Subroutine EigenvaluesSymmetric3


  Subroutine LUDecomp_LWORK( m, n, a, lda, work, lwork )

    Integer                            :: m, n, lda, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(*)    , Target :: work

    Integer                           :: ierr, info
    Integer(C_SIZE_T)                 :: sizeof_a, sizeof_work
    Real(dp), Dimension(:,:), Pointer :: pa
    Real(dp), Dimension(:)  , Pointer :: pwork
    Type(C_PTR)                       :: ha, hwork
    Type(C_PTR)                       :: da

    lwork = -1

    sizeof_a = m * n * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)

    pa => a(:,1:n)
    pwork => work(1:lwork)

    ha = C_LOC( pa )
    hwork = C_LOC( pwork )

    da = dev_ptr( pa(1,1) )

#if defined(XNET_LA_CUBLAS)
    ierr = cusolverDnDgetrf_bufferSize &
           ( cusolver_handle, m, n, da, lda, lwork )
#elif defined(XNET_LA_ROCM)
#elif defined(XNET_LA_ONEMKL)
#elif defined(XNET_LA_MAGMA)
#else
#endif

  End Subroutine LUDecomp_LWORK


  Subroutine LUDecomp( m, n, a, lda, work, lwork, ipiv, info )

    Integer                            :: m, n, lda, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(*)    , Target :: work
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Integer                           :: min_mn
    Integer(C_SIZE_T)                 :: sizeof_a, sizeof_work, sizeof_ipiv, sizeof_info
    Real(dp), Dimension(:,:), Pointer :: pa
    Real(dp), Dimension(:)  , Pointer :: pwork
    Integer,  Dimension(:),   Pointer :: pipiv
    Integer,                  Pointer :: pinfo
    Type(C_PTR)                       :: ha, hwork, hipiv, hinfo
    Logical                           :: data_on_device

    data_on_device = .false.
    min_mn = MIN(m,n)

    sizeof_a    = m * n * c_sizeof(0.0_DP)
    sizeof_work = lwork * c_sizeof(0.0_DP)
    sizeof_ipiv = min_mn * c_sizeof(0)
    sizeof_info = c_sizeof(0)

    pa => a(:,1:n)
    pwork => work(1:lwork)
    pipiv => ipiv(1:min_mn)
    pinfo => info

    ha = C_LOC( pa )
    hwork = C_LOC( pwork )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hwork, mydevice, sizeof_work ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

      Call LUDecomp_GPU( m, n, a, lda, work, lwork, ipiv, info )
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[LUDecomp] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[LUDecomp]   A missing'
      If ( .not. device_is_present( hwork, mydevice, sizeof_work ) ) &
        Write(*,*) '[LUDecomp]   work missing'
      If ( .not. device_is_present( hipiv, mydevice, sizeof_ipiv ) ) &
        Write(*,*) '[LUDecomp]   ipiv missing'
      If ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        Write(*,*) '[LUDecomp]   info missing'
#endif
#endif

      Call LUDecomp_CPU( m, n, a, lda, ipiv, info )

    End If

  End Subroutine LUDecomp


  Subroutine LUDecomp_CPU( m, n, a, lda, ipiv, info )

    Integer                            :: m, n, lda
    Real(dp), Dimension(lda,*), Target :: a
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Call DGETRF &
           ( m, n, a, lda, ipiv, info )

  End Subroutine LUDecomp_CPU


  Subroutine LUDecomp_GPU( m, n, a, lda, work, lwork, ipiv, info )

    Integer                            :: m, n, lda, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(*)    , Target :: work
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Integer                           :: ierr, min_mn
    Real(dp), Dimension(:,:), Pointer :: pa
    Real(dp), Dimension(:)  , Pointer :: pwork
    Integer,  Dimension(:)  , Pointer :: pipiv
    Integer                 , Pointer :: pinfo
    Type(C_PTR)                       :: da, dwork, dipiv, dinfo

    min_mn = MIN(m,n)

    pa => a(:,1:n)
    pipiv => ipiv(1:min_mn)
    pwork => work(1:lwork)
    pinfo => info

    da = dev_ptr( pa(1,1) )
    dwork = dev_ptr( pwork(1) )
    dipiv = dev_ptr( pipiv(1) )
    dinfo = dev_ptr( pinfo )

#if defined(XNET_LA_CUBLAS)
    ierr = cusolverDnDgetrf &
           ( cusolver_handle, m, n, da, lda, dwork, dipiv, dinfo )
#elif defined(XNET_LA_ROCM)
    Call hipblasCheck( hipblasDgetrf &
           ( hipblas_handle, n, da, lda, dipiv, dinfo ) )
#elif defined(XNET_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, ipiv, info )
    Call DGETRF &
           ( m, n, a, lda, ipiv, info )
    !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
    Call magma_dgetrf_native &
           ( m, n, da, lda, ipiv, info )
#endif

  End Subroutine LUDecomp_GPU


  Subroutine LUBksub( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Integer(C_SIZE_T)                 :: sizeof_a, sizeof_b, sizeof_ipiv, sizeof_info
    Real(dp), Dimension(:,:), Pointer :: pa, pb
    Integer,  Dimension(:),   Pointer :: pipiv
    Integer,                  Pointer :: pinfo
    Type(C_PTR)                       :: ha, hb, hipiv, hinfo
    Logical                           :: data_on_device

    data_on_device = .false.
    sizeof_a    = n * n * c_sizeof(0.0_DP)
    sizeof_b    = n * nrhs * c_sizeof(0.0_DP)
    sizeof_ipiv = n * c_sizeof(0)
    sizeof_info = c_sizeof(0)

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    pipiv => ipiv(1:n)
    pinfo => info

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

      Call LUBksub_GPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[LUBksub] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[LUBksub]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[LUBksub]   B missing'
      If ( .not. device_is_present( hipiv, mydevice, sizeof_ipiv ) ) &
        Write(*,*) '[LUBksub]   ipiv missing'
      If ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        Write(*,*) '[LUBksub]   info missing'
#endif
#endif

      Call LUBksub_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    End If

  End Subroutine LUBksub


  Subroutine LUBksub_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Call DGETRS &
           ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

  End Subroutine LUBksub_CPU


  Subroutine LUBksub_GPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Integer                           :: ierr
    Integer(C_INT)                    :: itrans
    Real(dp), Dimension(:,:), Pointer :: pa, pb
    Integer,  Dimension(:)  , Pointer :: pipiv
    Integer,                  Pointer :: pinfo
    Type(C_PTR)                       :: hinfo
    Type(C_PTR)                       :: da, db, dipiv, dinfo

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    pipiv => ipiv(1:n)
    pinfo => info

    hinfo = C_LOC( pinfo )

    da = dev_ptr( pa(1,1) )
    db = dev_ptr( pb(1,1) )
    dipiv = dev_ptr( pipiv(1) )
    dinfo = dev_ptr( pinfo )

    itrans = itrans_from_char( trans )

#if defined(XNET_LA_CUBLAS)
    ierr = cusolverDnDgetrs &
           ( cusolver_handle, itrans, n, nrhs, da, lda, dipiv, db, ldb, dinfo )
#elif defined(XNET_LA_ROCM)
    Call hipblasCheck( hipblasDgetrs &
           ( hipblas_handle, itrans, n, nrhs, da, lda, dipiv, db, ldb, hinfo ) )
#elif defined(XNET_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, ipiv )
    Call DGETRS &
           ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )
    !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
    Call magma_dgetrs_gpu &
           ( itrans, n, nrhs, da, lda, ipiv, db, ldb, hinfo )
#endif

  End Subroutine LUBksub_GPU


  Subroutine LinearSolve( trans, n, nrhs, a, lda, ipiv, b, ldb, info, work, lwork )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info
    Real(dp), Dimension(*)    , Target :: work

    Integer(C_SIZE_T)                 :: sizeof_a, sizeof_b, sizeof_ipiv, sizeof_info
    Real(dp), Dimension(:,:), Pointer :: pa, pb
    Integer,  Dimension(:),   Pointer :: pipiv
    Integer,                  Pointer :: pinfo
    Type(C_PTR)                       :: ha, hb, hipiv, hinfo
    Logical                           :: data_on_device

    data_on_device = .false.
    sizeof_a    = n * n * c_sizeof(0.0_DP)
    sizeof_b    = n * nrhs * c_sizeof(0.0_DP)
    sizeof_ipiv = n * c_sizeof(0)
    sizeof_info = c_sizeof(0)

    pa => a(:,1:n)
    pb => b(:,1:nrhs)
    pipiv => ipiv(1:n)
    pinfo => info

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

      Call LinearSolve_GPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, work, lwork )
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[LinearSolve] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[LinearSolve]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[LinearSolve]   B missing'
      If ( .not. device_is_present( hipiv, mydevice, sizeof_ipiv ) ) &
        Write(*,*) '[LinearSolve]   ipiv missing'
      If ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        Write(*,*) '[LinearSolve]   info missing'
#endif
#endif

      Call LinearSolve_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    End If

  End Subroutine LinearSolve


  Subroutine LinearSolve_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info

    Call LUDecomp_CPU &
      & ( n, n, a, lda, ipiv, info )
    Call LUBksub_CPU &
      & ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

  End Subroutine LinearSolve_CPU


  Subroutine LinearSolve_GPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, work, lwork )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, lwork
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,                    Target :: info
    Real(dp), Dimension(*)    , Target :: work

    Call LUDecomp_GPU &
      & ( n, n, a, lda, work, lwork, ipiv, info )
    Call LUBksub_GPU &
      & ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

  End Subroutine LinearSolve_GPU


  Subroutine LinearSolveBatched( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info

    Integer                                    :: i
    Integer(C_SIZE_T)                          :: sizeof_a, sizeof_b, sizeof_ipiv, sizeof_info
    Real(dp), Dimension(:,:), Pointer          :: pa, pb
    Integer,  Dimension(:),   Pointer          :: pipiv, pinfo
    Type(C_PTR)                                :: ha, hb, hipiv, hinfo
    Type(C_PTR), Dimension(batchcount), Target :: da, db, dipiv
    Integer                                    :: osa, osb
    Logical                                    :: data_on_device

    data_on_device = .false.
    sizeof_a    = n * n * batchcount * c_sizeof(0.0_DP)
    sizeof_b    = n * nrhs * batchcount * c_sizeof(0.0_DP)
    sizeof_ipiv = n * batchcount * c_sizeof(0)
    sizeof_info = batchcount * c_sizeof(0)

    pa => a(:,1:n*batchcount)
    pb => b(:,1:nrhs*batchcount)
    pipiv => ipiv(1:n*batchcount)
    pinfo => info(1:batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

      !__dir_enter_data &
      !__dir_create(da,db,dipiv)
      Do i = 1, batchcount
        osa = (i-1) * n + 1
        osb = (i-1) * nrhs + 1
        da(i) = dev_ptr( pa(1,osa) )
        db(i) = dev_ptr( pb(1,osb) )
        dipiv(i) = dev_ptr( pipiv(osa) )
      End Do
      !__dir_update_gpu(da,db,dipiv)

      Call LinearSolveBatched_GPU &
        &  ( trans, n, nrhs, a, da(1), lda, ipiv, dipiv(1), b, db(1), ldb, info, batchcount )
#if defined(XNET_OMP_OL)
      Call stream_sync( stream )
#endif

      !__dir_exit_data &
      !__dir_delete(da,db,dipiv)

    Else

#if defined(XNET_DEBUG_LA)
#if defined(XNET_GPU)
      Write(*,*) '[LinearSolveBatched] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[LinearSolveBatched]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[LinearSolveBatched]   B missing'
      If ( .not. device_is_present( hipiv, mydevice, sizeof_ipiv ) ) &
        Write(*,*) '[LinearSolveBatched]   ipiv missing'
      If ( .not. device_is_present( hinfo, mydevice, sizeof_info ) ) &
        Write(*,*) '[LinearSolveBatched]   info missing'
#endif
#endif

      Call LinearSolveBatched_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    End If

  End Subroutine LinearSolveBatched


  Subroutine LinearSolveBatched_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info

    Call LUDecompBatched_CPU &
      & ( n, n, a, lda, ipiv, info, batchcount )
    Call LUBksubBatched_CPU &
      & ( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

  End Subroutine LinearSolveBatched_CPU


  Subroutine LinearSolveBatched_GPU( trans, n, nrhs, a, da, lda, ipiv, dipiv, b, db, ldb, info, batchcount )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info
    Type(C_PTR), Dimension(*),  Target :: da, dipiv, db

    Call LUDecompBatched_GPU &
      & ( n, n, a, da(1), lda, ipiv, dipiv(1), info, batchcount )
    Call LUBksubBatched_GPU &
      & ( trans, n, nrhs, a, da(1), lda, ipiv, dipiv(1), b, db(1), ldb, info, batchcount )

  End Subroutine LinearSolveBatched_GPU


  Subroutine LUDecompBatched_CPU( m, n, a, lda, ipiv, info, batchcount )

    Integer                            :: m, n, lda, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info

    Integer :: i
    Integer :: osa

    Do i = 1, batchcount
      osa = (i-1) * n + 1
      Call DGETRF &
             ( m, n, a(1,osa), lda, ipiv(osa), info(i) )
    End Do

  End Subroutine LUDecompBatched_CPU


  Subroutine LUDecompBatched_GPU( m, n, a, da, lda, ipiv, dipiv, info, batchcount )

    Integer                            :: m, n, lda, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info
    Type(C_PTR), Dimension(*),  Target :: da, dipiv

    Integer                         :: ierr, i, stridea, strideipiv
    Integer(C_INT64_T)              :: strideP_64
    Integer,  Dimension(:), Pointer :: pinfo
    Type(C_PTR)                     :: da_array, dipiv_array, dinfo

    stridea    = n * n
    strideipiv = n

    pinfo => info(1:batchcount)

    da_array = dev_ptr( da(1) )
    dipiv_array = dev_ptr( dipiv(1) )
    dinfo = dev_ptr( pinfo(1) )

#if defined(XNET_LA_CUBLAS)
    ierr = cublasDgetrfBatched &
           ( cublas_handle, n, da_array, lda, dipiv(1), dinfo, batchcount )
#elif defined(XNET_LA_ROCM)
    !strideP_64 = n
    !Call rocsolverCheck( rocsolver_dgetrf_batched &
    !       ( rocsolver_handle, n, n, da_array, lda, dipiv(1), strideP_64, dinfo, batchcount ) )
    Call hipblasCheck( hipblasDgetrfBatched &
           ( hipblas_handle, n, da_array, lda, dipiv(1), dinfo, batchcount ) )
#elif defined(XNET_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, ipiv )
    Call DGETRF_BATCH_STRIDED &
           ( n, n, a, lda, stridea, ipiv, strideipiv, batchcount, info )
    !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
    Call magma_dgetrf_batched &
           ( n, n, da_array, lda, dipiv_array, dinfo, batchcount, magma_queue )
#endif

  End Subroutine LUDecompBatched_GPU


  Subroutine LUBksubBatched_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info

    Integer :: i
    Integer :: osa, osb

    Do i = 1, batchcount
      osa = (i-1) * n + 1
      osb = (i-1) * nrhs + 1
      Call DGETRS &
             ( trans, n, nrhs, a(1,osa), lda, ipiv(osa), b(1,osb), ldb, info(i) )
    End Do

  End Subroutine LUBksubBatched_CPU


  Subroutine LUBksubBatched_GPU( trans, n, nrhs, a, da, lda, ipiv, dipiv, b, db, ldb, info, batchcount )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info
    Type(C_PTR), Dimension(*),  Target :: da, dipiv, db

    Integer                         :: ierr, i, stridea, strideb, strideipiv
    Integer(C_INT)                  :: itrans
    Integer(C_INT64_T)              :: strideP_64
    Integer,  Dimension(:), Pointer :: pinfo
    Type(C_PTR)                     :: hinfo
    Type(C_PTR)                     :: da_array, db_array, dipiv_array, dinfo

    stridea    = n * n
    strideb    = n * nrhs
    strideipiv = n

    pinfo => info(1:batchcount)

    hinfo = C_LOC( pinfo )

    da_array = dev_ptr( da(1) )
    db_array = dev_ptr( db(1) )
    dipiv_array = dev_ptr( dipiv(1) )
    dinfo = dev_ptr( pinfo(1) )

    itrans = itrans_from_char( trans )

#if defined(XNET_LA_CUBLAS)
    ierr = cublasDgetrsBatched &
           ( cublas_handle, itrans, n, nrhs, da_array, lda, dipiv(1), db_array, ldb, hinfo, batchcount )
#elif defined(XNET_LA_ROCM)
    !strideP_64 = n
    !Call rocsolverCheck( rocsolver_dgetrs_batched &
    !       ( rocsolver_handle, itrans, n, nrhs, da_array, lda, dipiv(1), strideP_64, db_array, ldb, batchcount ) )
    Call hipblasCheck( hipblasDgetrsBatched &
           ( hipblas_handle, itrans, n, nrhs, da_array, lda, dipiv(1), db_array, ldb, hinfo, batchcount ) )
#elif defined(XNET_LA_ONEMKL)
    !$OMP TARGET VARIANT DISPATCH USE_DEVICE_PTR( a, b, ipiv )
    Call dgetrs_batch_strided &
           ( trans, n, nrhs, a, lda, stridea, ipiv, strideipiv, b, ldb, strideb, info )
    !$OMP END TARGET VARIANT DISPATCH
#elif defined(XNET_LA_MAGMA)
    Call magma_dgetrs_batched &
           ( itrans, n, nrhs, da_array, lda, dipiv_array, db_array, ldb, batchcount, magma_queue )
#endif

  End Subroutine LUBksubBatched_GPU


End Module xnet_linalg
