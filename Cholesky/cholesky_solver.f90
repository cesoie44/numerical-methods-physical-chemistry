
! ==============================================================================
! Module     : cholesky_solver
! Purpose    : Cholesky decomposition and solving linear systems
! Description: 
!     - Provides a subroutine to perform Cholesky decomposition A = L*L^T
!     - Provides a subroutine to solve linear systems using the decomposition
! ==============================================================================

module cholesky_solver
        
        use utils
        implicit none

        contains

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        subroutine decomposizione_L(n,L,A_loc,success)
        
        !---------------------------------------------------------------------------
        ! Subroutine to perform Cholesky decomposition of a positive definite matrix A
        ! A = L * L^T, where L is lower triangular
        !---------------------------------------------------------------------------  

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: L(n,n)
        real(dp) , intent(in) :: A_loc(n,n)
        logical , intent(out) :: success

        integer :: i,j,k
        real(dp) :: s
        
        L=0.0_dp
        success=.TRUE.

        do j = 1, n
                s = 0.0_dp
                do k = 1, j-1
                        s = s + L(j,k)**2
                end do
                if (A(j,j) - s .LE. 0.0_dp) then
                        success=.FALSE.
                        return
                end if
                L(j,j) = sqrt(A(j,j) - s)
    
                do i = j+1, n
                        s = 0.0_dp
                        do k = 1, j-1
                                s = s + L(i,k)*L(j,k)
                        end do
                        L(i,j) = (A(i,j) - s)/L(j,j)
                end do
        end do

        end subroutine decomposizione_L
        
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        subroutine sistema_lineare(n,b,L)

        !---------------------------------------------------------------------------
        ! Subroutine to solve A*x = b using the Cholesky decomposition A = L * L^T
        ! It performs forward and backward substitution
        !---------------------------------------------------------------------------

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(inout) :: b(n) , L(n,n)

        integer :: i,j
        real(dp) , allocatable :: x(:) , y(:)
        real(dp) :: s

        allocate(y(n))
        do i=1,n
                s=0.0_dp
                        do j=1,i-1
                                s=s+L(i,j)*y(j)
                        end do
                        y(i)=(b(i)-s)/L(i,i)
        end do
        allocate(x(n))
        do i=n,1,-1
                s=0.0_dp
                        do j=i+1,n
                                s=s+L(j,i)*x(j)
                        end do
                        x(i)=(y(i)-s)/L(i,i)
        end do
        write(6,*) 'The solution of the linear system is: '
        print *
        write(6,*) x
        deallocate(x,y)

        end subroutine sistema_lineare

end module cholesky_solver

