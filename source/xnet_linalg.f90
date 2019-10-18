MODULE xnet_linalg

  Use, Intrinsic :: iso_c_binding
  Use xnet_gpu, Only: &
    mydevice, &
    device_is_present
  Use xnet_types, Only: dp

#if defined(XNET_CUBLAS)
  Use cudaf, Only: &
    stream
  Use cublasf, Only: &
    cublas_handle, &
    cublasDnrm2_v2, &
    cublasDaxpy_v2, &
    cublasDgemm_v2, &
    cublasDgemmStridedBatched, &
    cublasDgetrfBatched, &
    cublasDgetrsBatched, &
    cublasDgemv_v2, &
    cublasDgeam, &
    cublasDdgmm, &
    CUBLAS_OP_N, CUBLAS_OP_T, &
    CUBLAS_SIDE_LEFT, &
    CUBLAS_FILL_MODE_UPPER, &
    CUBLAS_DIAG_NON_UNIT
#endif

#if defined(XNET_MAGMA)
  Use magmaf, Only: &
    magma_queue, &
    magma_dnrm2, &
    magma_daxpy, &
    magma_dgemm, &
    magmablas_dgemm_batched_strided, &
    magma_dgetrf_batched, &
    magma_dgetrf_nopiv_batched, &
    magma_dgetrs_batched, &
    magma_dlaswp_rowserial_batched, &
    magmablas_dtrsv_outofplace_batched, &
    magma_dgemv, &
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
  Public :: VectorNorm2
  Public :: VectorNorm2_Kernel
  Public :: VectorVectorAdd
  Public :: LinearSolve_CPU
  Public :: LinearSolveBatched
  Public :: LinearSolveBatched_CPU
  Public :: LinearSolveBatched_GPU
  Public :: LUDecomp_CPU
  Public :: LUDecompBatched_CPU
  Public :: LUDecompBatched_GPU
  Public :: LUBksub_CPU
  Public :: LUBksubBatched_CPU
  Public :: LUBksubBatched_GPU


Contains

  Integer Function itrans_from_char( ctrans )
    Character, Intent(in) :: ctrans
#if defined(XNET_CUBLAS)
    If ( ctrans == 'T' ) Then
      itrans_from_char = CUBLAS_OP_T
    Else
      itrans_from_char = CUBLAS_OP_N
    End If
#elif defined(XNET_MAGMA)
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

    pa(1:lda,1:ka) => a(:,1:ka)
    pb(1:ldb,1:kb) => b(:,1:kb)
    pc(1:ldc,1:n ) => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, pc )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, pc )
#endif
      da = C_LOC( pa )
      db = C_LOC( pb )
      dc = C_LOC( pc )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDgeam &
             ( cublas_handle, itransa, itransb, m, n, alpha, da, lda, beta, db, ldb, dc, ldc )
#elif defined(XNET_MAGMA)
      If ( transb  == 'N' ) Then
        Call magmablas_dlacpy &
               ( MagmaGeneral, m, n, db, ldb, dc, ldc, magma_queue )
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

    Else

#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixAdd] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixAdd]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixAdd]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixAdd]   C missing'
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

    pa(1:lda,1:ka) => a(:,1:ka)
    pb(1:ldb,1:kb) => b(:,1:kb)
    pc(1:ldc,1:n ) => c(:,1:n )

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, pc )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, pc )
#endif
      da = C_LOC( pa )
      db = C_LOC( pb )
      dc = C_LOC( pc )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDgemm_v2 &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc )
#elif defined(XNET_MAGMA)
      Call magma_dgemm &
             ( itransa, itransb, m, n, k, alpha, da, lda, db, ldb, beta, dc, ldc, magma_queue )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixMultiply] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixMultiply]   C missing'
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
    Integer(C_SIZE_T)                  :: sizeof_a, sizeof_b, sizeof_c
    Real(dp), Dimension(:,:), Pointer  :: pa, pb, pc
    Type(C_PTR)                        :: ha, hb, hc
    Type(C_PTR)                        :: da, db, dc
    Integer                            :: ka, kb
    Integer                            :: osa, osb, osc
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

    pa(1:lda,1:ka*batchcount) => a(:,1:ka*batchcount)
    pb(1:ldb,1:kb*batchcount) => b(:,1:kb*batchcount)
    pc(1:ldc,1:n *batchcount) => c(:,1:n *batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hc = C_LOC( pc )

    data_on_device = device_is_present( ha, mydevice, sizeof_a ) &
               .AND. device_is_present( hb, mydevice, sizeof_b ) &
               .AND. device_is_present( hc, mydevice, sizeof_c )

    If ( data_on_device ) Then

      itransa = itrans_from_char( transa )
      itransb = itrans_from_char( transb )

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, pc )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, pc )
#endif
      da = C_LOC( pa )
      db = C_LOC( pb )
      dc = C_LOC( pc )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDgemmStridedBatched &
             ( cublas_handle, itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount )
#elif defined(XNET_MAGMA)
      Call magmablas_dgemm_batched_strided &
             ( itransa, itransb, m, n, k, alpha, da, lda, stridea, &
               db, ldb, strideb, beta, dc, ldc, stridec, batchcount, magma_queue )
#endif

    Else

#if defined(XNET_GPU)
      Write(*,*) '[MatrixMatrixMultiplyBatched] Data not present on device'
      If ( .not. device_is_present( ha, mydevice, sizeof_a ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   A missing'
      If ( .not. device_is_present( hb, mydevice, sizeof_b ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   B missing'
      If ( .not. device_is_present( hc, mydevice, sizeof_c ) ) &
        Write(*,*) '[MatrixMatrixMultiplyBatched]   C missing'
#endif

      Do i = 1, batchcount
        osa = (i-1) * ka + 1
        osb = (i-1) * kb + 1
        osc = (i-1) * n  + 1
        Call DGEMM &
               ( transa, transb, m, n, k, alpha, a(1,osa), lda, b(1,osb), ldb, beta, c(1,osc), ldc )
      End Do

    End If

  End Subroutine MatrixMatrixMultiplyBatched


  Subroutine LinearSolve_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                  :: trans
    Integer                    :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*) :: a
    Real(dp), Dimension(ldb,*) :: b
    Integer,  Dimension(*)     :: ipiv
    Integer                    :: info

    Call LUDecomp_CPU &
      & ( n, n, a, lda, ipiv, info )
    Call LUBksub_CPU &
      & ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

  End Subroutine LinearSolve_CPU


  Subroutine LinearSolveBatched( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount, pivot )

    Character                          :: trans
    Integer                            :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*), Target :: a
    Real(dp), Dimension(ldb,*), Target :: b
    Integer,  Dimension(*),     Target :: ipiv
    Integer,  Dimension(*),     Target :: info
    Logical,  Optional                 :: pivot

    Integer,  Dimension(n,batchcount),      Target :: ipivinfo
    Real(dp), Dimension(n*nrhs,batchcount), Target :: work

    Integer                                    :: ierr, i
    Integer(C_SIZE_T)                          :: sizeof_a, sizeof_b, sizeof_ipiv, sizeof_info
    Real(dp), Dimension(:,:), Pointer          :: pa, pb
    Integer,  Dimension(:),   Pointer          :: pipiv, pinfo
    Type(C_PTR)                                :: ha, hb, hipiv, hinfo
    Type(C_PTR), Dimension(batchcount), Target :: da, db, dipiv, dipivinfo, dwork
    Integer                                    :: osa, osb, oswork
    Logical                                    :: data_on_device

    data_on_device = .false.
    sizeof_a    = n * n * batchcount * c_sizeof(0.0_DP)
    sizeof_b    = n * nrhs * batchcount * c_sizeof(0.0_DP)
    sizeof_ipiv = n * batchcount * c_sizeof(0)
    sizeof_info = batchcount * c_sizeof(0)

    pa(1:lda,1:n*batchcount) => a(:,1:n*batchcount)
    pb(1:ldb,1:nrhs*batchcount) => b(:,1:nrhs*batchcount)
    pipiv(1:n*batchcount) => ipiv(1:n*batchcount)
    pinfo(1:batchcount) => info(1:batchcount)

    ha = C_LOC( pa )
    hb = C_LOC( pb )
    hipiv = C_LOC( pipiv )
    hinfo = C_LOC( pinfo )

    data_on_device = device_is_present( ha,    mydevice, sizeof_a    ) &
               .AND. device_is_present( hb,    mydevice, sizeof_b    ) &
               .AND. device_is_present( hipiv, mydevice, sizeof_ipiv ) &
               .AND. device_is_present( hinfo, mydevice, sizeof_info )

    If ( data_on_device ) Then

#if defined(XNET_OMP_OL)
      !$OMP TARGET ENTER DATA &
      !$OMP MAP( alloc: da, db, dipiv, ipivinfo, work, dipivinfo, dwork )
#elif defined(XNET_OACC)
      !$ACC ENTER DATA &
      !$ACC CREATE( da, db, dipiv, ipivinfo, work, dipivinfo, dwork )
#endif

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pb, pipiv, pinfo, ipivinfo, work )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pb, pipiv, pinfo, ipivinfo, work )
#endif
      Do i = 1, batchcount
        osa = (i-1) * n + 1
        osb = (i-1) * nrhs + 1
        oswork = (i-1) * n * nrhs + 1
        da(i) = C_LOC( pa(1,osa) )
        db(i) = C_LOC( pb(1,osb) )
        dipiv(i) = C_LOC( pipiv(osa) )
        dipivinfo(i) = C_LOC( ipivinfo(1,osa) )
        dwork(i) = C_LOC( work(1,oswork) )
      End Do
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
      !$OMP TARGET UPDATE TO( da, db, dipiv, dipivinfo, dwork )
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
      !$ACC UPDATE DEVICE( da, db, dipiv, dipivinfo, dwork )
#endif

      Call LinearSolveBatched_GPU( trans, n, nrhs, da, lda, dipiv, dipivinfo, db, ldb, dwork, info, batchcount, pivot )

#if defined(XNET_OMP_OL)
      !$OMP TARGET EXIT DATA &
      !$OMP MAP( release: da, db, dipiv, dipivinfo, dwork )
#elif defined(XNET_OACC)
      !$ACC EXIT DATA &
      !$ACC DELETE( da, db, dipiv, dipivinfo, dwork )
#endif

    Else

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


  Subroutine LinearSolveBatched_GPU( trans, n, nrhs, da, lda, dipiv, dipivinfo, db, ldb, dwork, info, batchcount, pivot )

    Character                         :: trans
    Integer                           :: n, nrhs, lda, ldb, batchcount
    Type(C_PTR), Dimension(*), Target :: da, dipiv, dipivinfo, db, dwork
    Integer,     Dimension(*), Target :: info
    Logical,     Optional             :: pivot

    Call LUDecompBatched_GPU &
      & ( n, n, da(1), lda, dipiv(1), dipivinfo(1), info, batchcount, pivot )
    Call LUBksubBatched_GPU &
      & ( trans, n, nrhs, da(1), lda, dipiv(1), db(1), ldb, dwork(1), info, batchcount, pivot )

  End Subroutine LinearSolveBatched_GPU


  Subroutine LUDecomp_CPU( m, n, a, lda, ipiv, info )

    Integer                    :: m, n, lda
    Real(dp), Dimension(lda,*) :: a
    Integer,  Dimension(*)     :: ipiv
    Integer                    :: info

    Call DGETRF &
      & ( m, n, a, lda, ipiv, info )

  End Subroutine LUDecomp_CPU


  Subroutine LUDecompBatched_CPU( m, n, a, lda, ipiv, info, batchcount )

    Integer                    :: m, n, lda, batchcount
    Real(dp), Dimension(lda,*) :: a
    Integer,  Dimension(*)     :: ipiv
    Integer,  Dimension(*)     :: info

    Integer                    :: ierr, i
    Integer                    :: osa

    Do i = 1, batchcount
      osa = (i-1) * n + 1
      Call DGETRF &
             ( m, n, a(1,osa), lda, ipiv(osa), info(i) )
    End Do

  End Subroutine LUDecompBatched_CPU


  Subroutine LUDecompBatched_GPU( m, n, da, lda, dipiv, dipivinfo, info, batchcount, pivot )

    Integer                           :: m, n, lda, batchcount
    Type(C_PTR), Dimension(*), Target :: da, dipiv, dipivinfo
    Integer,     Dimension(*), Target :: info
    Logical,     Optional             :: pivot

    Logical                           :: lpiv
    Integer                           :: ierr, i
    Type(C_PTR)                       :: da_array, dipiv_array, dipivinfo_array, dinfo

    If ( present(pivot) ) Then
      lpiv = pivot
    Else
      lpiv = .true.
    EndIf

#if defined(XNET_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( da, dipiv, dipivinfo, info )
#elif defined(XNET_OACC)
    !$ACC HOST_DATA USE_DEVICE( da, dipiv, dipivinfo, info )
#endif
    da_array = C_LOC( da(1) )
    dipiv_array = C_LOC( dipiv(1) )
    dipivinfo_array = C_LOC( dipivinfo(1) )
    dinfo = C_LOC( info(1) )
#if defined(XNET_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(XNET_OACC)
    !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
    ierr = cublasDgetrfBatched &
      & ( cublas_handle, m, da_array, lda, dipiv(1), dinfo, batchcount )
#elif defined(XNET_MAGMA)
    If ( lpiv ) Then
      Call magma_dgetrf_batched &
        & ( m, n, da_array, lda, dipiv_array, dipivinfo_array, dinfo, batchcount, magma_queue )
    Else
      Call magma_dgetrf_nopiv_batched &
        & ( m, n, da_array, lda, dinfo, batchcount, magma_queue )
    EndIf
#endif

  End Subroutine LUDecompBatched_GPU


  Subroutine LUBksub_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

    Character                  :: trans
    Integer                    :: n, nrhs, lda, ldb
    Real(dp), Dimension(lda,*) :: a
    Real(dp), Dimension(ldb,*) :: b
    Integer,  Dimension(*)     :: ipiv
    Integer                    :: info

    Call DGETRS &
      & ( trans, n, nrhs, a, lda, ipiv, b, ldb, info )

  End Subroutine LUBksub_CPU


  Subroutine LUBksubBatched_CPU( trans, n, nrhs, a, lda, ipiv, b, ldb, info, batchcount )

    Character                  :: trans
    Integer                    :: n, nrhs, lda, ldb, batchcount
    Real(dp), Dimension(lda,*) :: a
    Real(dp), Dimension(ldb,*) :: b
    Integer,  Dimension(*)     :: ipiv
    Integer,  Dimension(*)     :: info

    Integer                    :: ierr, i
    Integer                    :: osa, osb

    Do i = 1, batchcount
      osa = (i-1) * n + 1
      osb = (i-1) * nrhs + 1
      Call DGETRS &
             ( trans, n, nrhs, a(1,osa), lda, ipiv(osa), b(1,osb), ldb, info(i) )
    End Do

  End Subroutine LUBksubBatched_CPU


  Subroutine LUBksubBatched_GPU( trans, n, nrhs, da, lda, dipiv, db, ldb, dwork, info, batchcount, pivot )

    Character                         :: trans
    Integer                           :: n, nrhs, lda, ldb, batchcount
    Type(C_PTR), Dimension(*), Target :: da, dipiv, db, dwork
    Integer,     Dimension(*), Target :: info
    Logical,     Optional             :: pivot

    Logical                           :: lpiv
    Integer                           :: ierr, i
    Integer(C_INT)                    :: itrans
    Type(C_PTR)                       :: da_array, db_array, dipiv_array, dwork_array
    Type(C_PTR)                       :: dinfo, hinfo

    If ( present(pivot) ) Then
      lpiv = pivot
    Else
      lpiv = .true.
    EndIf

    itrans = itrans_from_char( trans )

#if defined(XNET_OMP_OL)
    !$OMP TARGET DATA USE_DEVICE_PTR( da, db, dipiv, dwork, info )
#elif defined(XNET_OACC)
    !$ACC HOST_DATA USE_DEVICE( da, db, dipiv, dwork, info )
#endif
    da_array = C_LOC( da(1) )
    db_array = C_LOC( db(1) )
    dipiv_array = C_LOC( dipiv(1) )
    dwork_array = C_LOC( dwork(1) )
    dinfo = C_LOC( info(1) )
#if defined(XNET_OMP_OL)
    !$OMP END TARGET DATA
#elif defined(XNET_OACC)
    !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
    hinfo = C_LOC( info(1) )
    ierr = cublasDgetrsBatched &
      & ( cublas_handle, itrans, n, nrhs, da_array, lda, dipiv(1), db_array, ldb, hinfo, batchcount )
#elif defined(XNET_MAGMA)
    !Call magma_dgetrs_batched &
    !  & ( itrans, n, nrhs, da_array, lda, dipiv_array, db_array, ldb, batchcount, magma_queue )
    !Call magma_dgetrs_nopiv_batched &
    !  & ( itrans, n, nrhs, da_array, lda, db_array, ldb, dinfo, batchcount, magma_queue )
    If ( trans == 'N' ) Then
      If ( lpiv ) Then
        Call magma_dlaswp_rowserial_batched &
          & ( nrhs, db_array, ldb, 1, n, dipiv_array, batchcount, magma_queue )
      EndIf
      If ( nrhs == 1 ) Then
        Call magmablas_dtrsv_outofplace_batched &
          & ( MagmaLower, MagmaNoTrans, MagmaUnit, n, da_array, lda, db_array, 1, dwork_array, batchcount, magma_queue, 0 )
        Call magmablas_dtrsv_outofplace_batched &
          & ( MagmaUpper, MagmaNoTrans, MagmaNonUnit, n, da_array, lda, dwork_array, 1, db_array, batchcount, magma_queue, 0 )
      !Else
      !  Call magmablas_dtrsm_batched &
      !    & ( MagmaLeft, MagmaLower, MagmaNoTrans, MagmaUnit, n, nrhs, 1.0d0, da_array, lda, db_array, ldb, batchcount, magma_queue )
      !  Call magmablas_dtrsm_batched &
      !    & ( MagmaLeft, MagmaUpper, MagmaNoTrans, MagmaNonUnit, n, nrhs, 1.0d0, da_array, lda, db_array, ldb, batchcount, magma_queue )
      EndIf
    Else
      If ( nrhs == 1 ) Then
        Call magmablas_dtrsv_outofplace_batched &
          & ( MagmaUpper, itrans, MagmaUnit, n, da_array, lda, db_array, 1, dwork_array, batchcount, magma_queue, 0 )
        Call magmablas_dtrsv_outofplace_batched &
          & ( MagmaLower, itrans, MagmaNonUnit, n, da_array, lda, dwork_array, 1, db_array, batchcount, magma_queue, 0 )
      !Else
      !  Call magmablas_dtrsm_batched &
      !    & ( MagmaLeft, MagmaUpper, itrans, MagmaUnit, n, nrhs, 1.0d0, da_array, lda, db_array, ldb, batchcount, magma_queue )
      !  Call magmablas_dtrsm_batched &
      !    & ( MagmaLeft, MagmaLower, itrans, MagmaNonUnit, n, nrhs, 1.0d0, da_array, lda, db_array, ldb, batchcount, magma_queue )
      EndIf
      If ( lpiv ) Then
        Call magma_dlaswp_rowserial_batched &
          & ( nrhs, db_array, ldb, 1, n, dipiv_array, batchcount, magma_queue )
      EndIf
    EndIf
#endif

  End Subroutine LUBksubBatched_GPU


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

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, px, py )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, px, py )
#endif
      da = C_LOC( pa )
      dx = C_LOC( px )
      dy = C_LOC( py )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDgemv_v2 &
             ( cublas_handle, itrans, m, n, alpha, da, lda, dx, incx, beta, dy, incy )
#elif defined(XNET_MAGMA)
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

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( pa, pc, px )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( pa, pc, px )
#endif
      da = C_LOC( pa )
      dc = C_LOC( pc )
      dx = C_LOC( px )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDdgmm &
             ( cublas_handle, CUBLAS_SIDE_LEFT, m, n, da, lda, dx, incx, dc, ldc )
#elif defined(XNET_MAGMA)
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

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( px )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( px )
#endif
      dx = C_LOC( px )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDnrm2_v2( cublas_handle, n, dx, incx, xnorm )
#elif defined(XNET_MAGMA)
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
#if defined(XNET_OMP_OL)
    !$OMP DECLARE Target
#elif defined(XNET_OACC)
    !$ACC ROUTINE SEQ
#endif

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

#if defined(XNET_OMP_OL)
      !$OMP TARGET DATA USE_DEVICE_PTR( px, py )
#elif defined(XNET_OACC)
      !$ACC HOST_DATA USE_DEVICE( px, py )
#endif
      dx = C_LOC( px )
      dy = C_LOC( py )
#if defined(XNET_OMP_OL)
      !$OMP END TARGET DATA
#elif defined(XNET_OACC)
      !$ACC END HOST_DATA
#endif

#if defined(XNET_CUBLAS)
      ierr = cublasDaxpy_v2( cublas_handle, n, alpha, dx, incx, dy, incy )
#elif defined(XNET_MAGMA)
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


End Module xnet_linalg
