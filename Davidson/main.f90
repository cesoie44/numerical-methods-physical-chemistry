
! ==============================================================================
! Program    : main
! Purpose    : Compute a few eigenvalues of a symmetric matrix using the Davidson algorithm
! Author     : [Pietro Gozzoli, Alessandro Michieletti]
! Modules    : utils (matrix utilities), davidson_solver (Davidson method), ortho (orthonormalization utilities)
! Description:
!     - Reads the system dimension and the number of eigenvalues from user input
!     - Defines a symmetric matrix A
!     - Applies the Davidson algorithm to compute the lowest eigenvalues
!     - Deallocates the matrix A
! ==============================================================================


program main

        use utils                                       ! Module containing general utilities (define_a, free_a)
        use davidson_solver                             ! Module implementing the Davidson eigenvalue solver
        use ortho                                       ! Module for orthonormalization routines

        implicit none
        integer :: n,n_eig
        real(dp) :: start,finish
       
!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        
        call cpu_time(start) 
        write(6,*) 'Enter the system dimension: '
        read(5,*) n
        print *        

        write(6,*) 'Enter the number of eigenvalues to compute: '
        read(5,*) n_eig
        print *

!--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

        call davidson_algorithm(n,n_eig)                ! Solve the eigenvalue problem for the n_eig lowest eigenvalues
        call cpu_time(finish)
        print *
        write(6,*) 'CPU time: ', finish-start , 'seconds'
        print *

end program main
