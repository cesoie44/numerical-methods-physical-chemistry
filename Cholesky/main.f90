
! ==============================================================================
! Program    : main
! Purpose    : Solve a linear system A*x = b using Cholesky decomposition
! Author     : [Pietro Gozzoli, Alessandro Michieletti]
! Modules    : cholesky_solver (Cholesky decomposition and solver), utils (general utilities)
! Description: 
!     - Reads the system size from user input
!     - Defines a symmetric positive definite matrix A
!     - Performs Cholesky decomposition A = L*L^T
!     - Generates a random right-hand side vector b
!     - Solves the system L*L^T*x = b
!     - Prints the execution time
! ==============================================================================


program main

        use cholesky_solver                                     ! Module for Cholesky decomposition and linear system solving
        use utils                                               ! Module containing general utilities
        implicit none
        real(dp) , allocatable :: L(:,:) , b(:)
        integer :: n,i
        real(dp) :: start,finish
        logical :: success

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call cpu_time(start)                                    ! Start CPU time measurement
        write(6,*) 'Enter the system dimension: '
        read(5,*) n
        print *

        call define_a(n)                                       
        write(6,*) 'Matrix A is: '
        print *
        do i=1,n
                write(6,*) A(i,:)                              
        end do
        print*

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(L(n,n))                                        ! Allocate memory for lower triangular matrix L
        call decomposizione_L(n,L,A,success)                    ! Perform Cholesky decomposition
        if (.NOT. success) then
                write(6,*) 'The matrix is not positive definite'
                deallocate(L,A)
                stop
        end if
        write(6,*) 'The lower triangular matrix L is: '
        print *
        do i=1,n
                write(6,*) L(i,:)                              
        end do
        print *

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(b(n))                                         
        call random_number(b)                                   ! Fill vector b with random numbers
        write(6,*) 'Vector b is: ', b
        print *

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call sistema_lineare(n,b,L)                             ! Solve the linear system L*L^T*x = b
        deallocate(L,b)                                       

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call free_a                                             ! Deallocate matrix A
        call cpu_time(finish)                                   ! End CPU time measurement
        print *
        write(6,*) 'CPU time: ', finish-start


end program main

