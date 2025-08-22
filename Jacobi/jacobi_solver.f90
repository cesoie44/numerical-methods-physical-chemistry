
! ==============================================================================
! Module      : jacobi_solver
! Purpose     : Solve linear systems using the Jacobi iterative method
! Author      : [Your Name]
! Date        : [Date]
! Dependencies: utils (general utilities)
! Description : 
!     - Provides subroutines for diagonal and off-diagonal matrix operations
!     - Implements the Jacobi iteration method
!     - Checks convergence based on residual norms
! ==============================================================================

module jacobi_solver

        use utils
        implicit none

contains

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        !---------------------------------------------------------------------------
        ! Subroutine to divide vector v element-wise by the diagonal elements of A
        ! Diagonal elements are assumed to be (i+1) for row i
        !---------------------------------------------------------------------------
        
        subroutine diagonal_A(n,x,v)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: x(n)
        real(dp) , intent(in) :: v(n)

        integer :: i

        do i=1,n
                x(i)=v(i)/real(i+1,dp)
        end do

        end subroutine diagonal_A


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        !---------------------------------------------------------------------------
        ! Subroutine to compute the off-diagonal contribution for the Jacobi method
        ! Off-diagonal elements are assumed to be 1/(i+j) for (i,j)
        !---------------------------------------------------------------------------

        subroutine fuori_diag_A(n,z,x_old)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: z(n)
        real(dp) , intent(in) :: x_old(n)

        integer :: i,j

        z=0.0_dp

        do i=1,n
                do j=1,n
                        if (i .NE. j) then
                                z(i)=z(i)+x_old(j)/real(i+j,dp)
                        end if
                end do
        end do

        end subroutine fuori_diag_A


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        

        !---------------------------------------------------------------------------
        ! Subroutine to solve a linear system using the Jacobi iterative method
        !---------------------------------------------------------------------------

        subroutine jacobi(n,x_old,b)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: x_old(n) , b(n)

        real(dp) :: tau
        integer :: k_max , k
        real(dp) , allocatable :: z(:) , v(:) , x(:)
        logical :: converged

        k_max=200               ! Maximum number of iterations
        tau=1.0e-10_dp          ! Convergence tolerance
        k=0
        allocate(z(n),v(n),x(n))
        
        
        do while (k .LE. k_max)
                k=k+1
                print *
                write(6,*) 'Iteration: ', k

                ! Compute the off-diagonal part: z = A_off_diag * x_old
                call fuori_diag_A(n,z,x_old)
                
                ! Compute vector v = b - z
                v=b-z

                ! Compute new approximation: x = D^(-1) * v
                call diagonal_A(n,x,v)

                ! Check for convergence
                call residuo(n,tau,x,x_old,converged)
                if (converged) exit
        end do

        deallocate(z,v,x)

        end subroutine jacobi


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        

        !---------------------------------------------------------------------------
        ! Subroutine to compute the residual and check for convergence
        ! Convergence is achieved if the 2-norm and infinity norm are below tau
        !---------------------------------------------------------------------------

        subroutine residuo(n,tau,x,x_old,converged)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(in) :: tau
        real(dp) , intent(in) :: x(n)
        real(dp) , intent(inout) :: x_old(n)
        logical , intent(out) :: converged

        real(dp) , allocatable :: r(:)
        real(dp) :: norma,norma_infty,dnrm2

        allocate(r(n))
        r=x-x_old
        norma=dnrm2(n, r, 1)            ! Euclidean (L2) norm of the residual
        norma_infty = maxval(abs(r))    ! Infinity norm (max absolute value)
        write(6,*) 'The 2-norm of the residual is: ', norma
        write(6,*) 'The infinity norm of the residual is: ', norma_infty
        
        ! Check if the solution has converged
        if (norma .le. tau .and. norma_infty .le. tau) then
                print *
                write(6,*) 'The solution is x: ', x(:)
                converged=.TRUE.
        else
                x_old=x
                converged=.FALSE.
        end if
        deallocate(r)

        end subroutine residuo


end module jacobi_solver
