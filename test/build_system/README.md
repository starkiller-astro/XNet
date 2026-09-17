# Production build-system checks

Run the focused build-system checks from the repository root:

```bash
python3 test/build_system/test_build_system.py
```

The script builds the three direct products in a temporary GNU build directory,
checks generic, Perlmutter, and retired Summit/Cori host selection, verifies
configuration-reuse protection and configuration-local cleaning, and checks
CUDA selector resolution and early rejection. The CUDA checks establish Make
selector behavior only; they do not qualify NVIDIA hardware or a CUDA
toolchain.
