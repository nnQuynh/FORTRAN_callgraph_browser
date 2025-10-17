**** How to compile PHITS with KURBUC function ****
1. Compile PHITS as usual
2. Copy replace kurbuc.o or kurbuc.obj by one of the following files suitable for your environment
3. Save a source file except for kurbuc.f (e.g. main.f), and recompile PHITS

**** List of object files *****
kurbuc_WinIFX.obj    : Windows Intel Fortran (IFX) Single or MPI
kurbuc_WinIFX_OMP.obj: Windows Intel Fortran (IFX) OpenMP or Hybrid
kurbuc_LinIFX.o      : Linux Intel Fortran (IFX) Single or MPI
kurbuc_LinIFX_OMP.o  : Linux Intel Fortran (IFX) OpenMP or Hybrid
kurbuc_MacGfort_arm64.o      : Mac gfortran (14.2.0) Single or MPI for arm64
kurbuc_MacGfort_OMP_arm64.o      : Mac gfortran (14.2.0) OpenMP or Hybrid for arm64
kurbuc_MacGfort_x86_64.o      : Mac gfortran (14.2.0) Single or MPI for x86_64
kurbuc_MacGfort_OMP_x86_64.o      : Mac gfortran (14.2.0) OpenMP or Hybrid for x86_64

If you get "Segmentation fault" or "Access violation" error during
the executation of PHITS with KURBUC mode,
please change the version of your compiler equivalent to ours (as listed above)
