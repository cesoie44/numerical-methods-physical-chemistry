
! ==============================================================================
! Module     : conjugate_gradient_solver
! Purpose    : Implements the Conjugate Gradient method for solving linear systems
! Author     : [Pietro Gozzoli, Alessandro Michieletti]
! Description: 
!     - Provides subroutines for the Conjugate Gradient method
!     - Includes utility subroutines for matrix-vector multiplication (Matvec) and
!       preconditioning (Mr)
!     - Solves A*x = b using Conjugate Gradient method
! ==============================================================================


module conjugate_gradient_solver

        use utils
        implicit none

contains

        ! -----------------------------------------------------------------------------
        ! Subroutine : Matvec
        ! Purpose    : Performs matrix-vector multiplication A*x = b
        ! Arguments  : n - size of the system
        !             x_guess - input vector
        !             x - result vector
        ! -----------------------------------------------------------------------------

        subroutine Matvec(n,x_guess,x)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(in) :: x_guess(n)
        real(dp) , intent(inout) :: x(n)

        integer :: i,j

        do i=1,n
                x(i)=0.0_dp
                do j=1,n
                        if (i .NE. j) then
                                x(i)=x(i)+x_guess(j)/real(i+j,dp)
                        else
                                x(i)=x(i)+x_guess(i)*real(i+1,dp)
                        end if
                end do
        end do

        end subroutine Matvec


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        

        ! -----------------------------------------------------------------------------
        ! Subroutine : Mr
        ! Purpose    : Preconditioner for the Conjugate Gradient method
        ! Arguments  : n - size of the system
        !             r - residual vector
        !             z - preconditioned residual vector
        ! ----------------------------------------------------------------------------


        subroutine Mr(n,r,z)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: z(n)
        real(dp) , intent(in) :: r(n)

        integer :: i

        z=0.0_dp

        do i=1,n
                z(i)=r(i)/real(i+1,dp)
        end do

        end subroutine Mr


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        
        ! -----------------------------------------------------------------------------
        ! Subroutine : conj_grad
        ! Purpose    : Solves A*x = b using the Conjugate Gradient method
        ! Arguments  : n - size of the system
        !             x_guess - initial guess for the solution
        !             b - right-hand side vector
        ! -----------------------------------------------------------------------------


        subroutine conj_grad(n,x_guess,b)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: x_guess(n)
        real(dp) , intent(in) :: b(n)

        integer :: k_max,k
        real(dp) , allocatable :: x(:) , r(:) , z(:) , p(:) , h(:)
        real(dp) :: gamma_old,norm_x,tau,f,alpha,gamma_new,beta,dnrm2,ddot
        logical :: converged

        ! Initialize variables
        
        k=0
        tau=1.0e-10_dp
        k_max=200
        allocate(x(n))
        call Matvec(n,x_guess,x)
        norm_x=dnrm2(n, x, 1)
        
        ! Calculate initial residual
        
        allocate(r(n))
        if (norm_x .EQ. 0.0_dp) then
                r=b
        else
                r=b-x
        end if
        deallocate(x)
        
        ! Apply preconditioner
        
        allocate(z(n))
        call Mr(n,r,z)
        allocate(p(n))
        p=z
        gamma_old=ddot(n,r,1,z,1)
        
        ! Begin iterative loop
        
        do while (k .LE. k_max)

                k=k+1
                write(6,*) 'Iteration: ', k
                allocate(h(n))
                call Matvec(n,p,h)

                ! Compute step size and update solution
                
                f=ddot(n,h,1,p,1)
                alpha=gamma_old/f
                x_guess=x_guess+alpha*p
                r=r-alpha*h
                deallocate(h)
                
                ! Check convergence
                
                call residuo(n,tau,r,x_guess,converged)
                call Mr(n,r,z)
                gamma_new=ddot(n,r,1,z,1)
                beta=gamma_new/gamma_old
                p=z+beta*p
                gamma_old=gamma_new
                if(converged) exit

        end do
        deallocate(z,p,r)

        end subroutine conj_grad


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        
        ! -----------------------------------------------------------------------------
        ! Subroutine : residuo
        ! Purpose    : Checks the residuals and convergence of the solution
        ! Arguments  : n - size of the system
        !             tau - tolerance for convergence
        !             r - residual vector
        !             x_guess - current guess of the solution
        !             converged - logical flag for convergence
        ! -----------------------------------------------------------------------------


        subroutine residuo(n,tau,r,x_guess,converged)

        integer , intent(in) :: n
        real(dp) , intent(in) :: tau
        real(dp) ,  intent(inout) :: x_guess(n) , r(n)
        logical ,  intent(out) :: converged

        real(dp) :: norma_r,norma_infty_r,dnrm2

        ! Compute residual norms

        norma_r=dnrm2(n, r, 1)
        norma_infty_r = maxval(abs(r))
        write(6,*) 'The 2-norm of the residual is: ' , norma_r
        write(6,*) 'The infinity norm of the residual is: ', norma_infty_r
        
        ! Check for convergence

        if (norma_r .LE. tau .AND. norma_infty_r .LE. tau) then
                print *
                write(6,*) 'The solution of the linear system is x: ', x_guess(:)
                converged=.TRUE.
        else
                converged=.FALSE.
                write(6,*) 'Convergence not reached.'
                print *
        end if

        end subroutine residuo

end module conjugate_gradient_solver
