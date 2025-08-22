# numerical-methods-physical-chemistry

This repository contains implementations of numerical methods developed as part of the "Metodi Numerici per la Chimica Fisica" course.

## Implemented Algorithms 

- Cholesky Decomposition with linear systems solver. 
- Jacobi method.
- Conjugate gradient with preconditioning.
- Davidson method eigensolver.

## Building and Running 

Compiled with GNU Fortran (`gfortran`) compiler 

`BLAS` and `LAPACK` libraries needed 

```
# Navigate to desired algorithm directory
cd Cholesky/  # or Jacobi/ or Conjugate_gradient/ or Davidson/

# Compile using the provided Makefile
make

# Run the executablek
./main.exe

```
## Test Matrix 
All implementations use a common test matrix bulid with:

- Dominant diagonal elements defined as  A_(ii)=i+1, 
- and off-diagonal elements given by A_(ij)=1/(i+j). 

This matrix structure mimics the properties usually found in quantum chemistry problems such as being symmetric, positive definite, diagonally dominant, and sparse.    
