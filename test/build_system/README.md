# Production build-system checks

Run the focused build-system checks from the repository root:

```bash
python3 test/build_system/test_build_system.py
```

The script invokes the repository-root Makefile and builds the three programs
in temporary GNU build directories. It
checks generic and Perlmutter defaults, explicit machine selection, separate
parallel builds, automatic and explicit build-directory names, incremental
rebuilds, configuration-reuse protection,
configuration-specific cleaning, solver selection, CUDA path selection, and
early rejection of incompatible options. The CUDA checks exercise Make
selection only; they do not qualify NVIDIA hardware or a CUDA toolchain.
