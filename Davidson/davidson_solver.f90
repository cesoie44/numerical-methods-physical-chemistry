
! ==============================================================================
! Module     : davidson
! Purpose    : Solve for eigenvalues using the Davidson algorithm
! Description:
!     - Iterative solver for a few eigenvalues of large symmetric matrices.
!     - Expands a subspace and projects the eigenproblem.
!     - Handles convergence and restarting of the subspace.
! ==============================================================================


module davidson_solver

        use utils
        use ortho
contains

!-----------------------------------------------------------------------------------------------------------------------------------
        
                !---------------------------------------------------------------------------
                ! Subroutine to solve for the lowest eigenvalues of a symmetric matrix A
                ! using the Davidson iterative method.
                !
                ! Arguments:
                !     n      : dimension of the matrix A
                !     n_eig  : number of eigenvalues to compute
                !---------------------------------------------------------------------------
                
                subroutine davidson_algorithm(n,n_eig)

                integer , intent(in) :: n,n_eig
                
                real(dp) :: tau,denominator
                integer :: ldu ,i,m_max,m,k_max,lda,lwork,ind,flag
                real(dp) , allocatable :: V(:,:) , W(:,:) , E(:,:) , U(:,:) , X(:,:) , R(:,:) , Y(:,:), YV(:,:) , eigenvalues(:) , work(:)
                logical :: logic
                
                
                flag=0                          ! Check if n_eig > n
                if (n_eig .gt. n) then
                        flag=1
                        write(6,*) 'The number of eigenvalues exceeds the matrix size.'
                        stop
                end if
                
                
                k=0                             ! Initialize parameters
                m=0
                k_max=200
                m_max=20
                tau=1.0D-8

                lda=m_max*n_eig
                lwork=lda**2
               

                allocate(V(n,lda))              ! Allocate matrices and vectors
                allocate(W(n,lda))  
                allocate(E(lda,lda))
                allocate(U(lda,lda))
                allocate(X(n,n_eig))
                allocate(R(n,n_eig))
                allocate(Y(n,n_eig))
                allocate(YV(lda,n_eig))
                allocate(eigenvalues(lda))
                allocate(work(lwork))
                

                V=0.0_dp                        ! Initialize V
                do i=1,n_eig
                        V(i,i)=1.0_dp
                end do

                do while (k .LE. k_max)         ! Davidson main iteration loop
                        !write(6,*) 'Iteration: ', k                       
                        k=k+1
                        m=m+1
                        ldu=m*n_eig
                        logic=.FALSE.        

                        
                        ! Compute W = A * V
                        do i=1, n_eig
                                ind=(m-1)*n_eig+i
                                call calcAV(n, V(1,ind), W(1,ind))
                        end do
                       
                        call dgemm("T", "N", ldu, ldu, n, 1.d0, V, n, W, n, 0.d0, E, lda)       ! Compute E = V^T * W
                        
                        
                        U=E
                        call dsyev("V", "L", ldu, U, lda, eigenvalues, work, lwork, info)       ! Diagonalize E
                     


                        if (info .NE. 0) then
                                write(6,*) 'Error during diagonalization: info = ', info
                                return
                        end if
                        

                        call dgemm("N", "N", n, n_eig, ldu, 1.d0, V, n, U, lda, 0.d0, X, n)     ! Compute approximate eigenvectors X and residuals R
                        call dgemm("N", "N", n, n_eig, ldu, 1.d0, W, n, U, lda, 0.d0, R, n) 

                        do i = 1, n_eig
                                R(:, i) = R(:, i) - eigenvalues(i) * X(:, i)
                        end do
                        

                        call convergenza(n,R,tau,logic,n_eig)                                   ! Check convergence
                

                        if (logic) then
                                write(6,*) 'Convergence reached.The eigenvalues are: '
                                print *
                                do i=1,n_eig
                                        write(6,*) eigenvalues(i)
                                end do
                                exit
                        end if

                        do i=1,n                                                                ! Build correction vectors Y
                                do j=1,n_eig
                                        denominator = real(i+1,dp)-eigenvalues(1)
                                        if (abs(denominator) < 1.0D-10) denominator = sign(1.0D-10, denominator)
                                        Y(i,j) = R(i,j)/denominator
                                end do
                        end do
                        
                        call ortho_vs_x(n, ldu, n_eig, V, Y)                                    ! Orthogonalize Y against V
                        

                        if (m.lt.m_max) then                                                    ! Update the basis V
                        do i=1,n_eig
                                V(:,ldu+i)=Y(:,i)
                        end do
                        else
                                m=1
                                do i=1,n_eig
                                        V(:,i)=X(:,i)
                                        V(:,n_eig+i)=Y(:,i)
                                end do
                        end if
                end do
                
                deallocate(V,W,E,U,X,R,Y,YV,eigenvalues,work) 
        
        
        
        end subroutine davidson_algorithm   

end module davidson_solver
