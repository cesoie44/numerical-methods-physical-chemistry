
! ==============================================================================
! Program     : main
! Purpose     : Solve a linear system using the Jacobi iterative method
! Author      : [Pietro Gozzoli, Alessandro Michieletti]
! Modules     : jacobi_solver (for Jacobi method), utils (for general utilities)
! Description : 
!     - Generates a random vector b
!     - Solves the system using Jacobi method
!     - Measures and prints CPU execution time
! ==============================================================================

program main

        use jacobi_solver                       ! Module for solving systems with Jacobi method
        use utils                               ! Module containing general utilities
        implicit none
        real(dp) , allocatable :: x_old(:) , b(:)
        integer :: n,i
        real(dp) :: start,finish

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call cpu_time(start)                    ! Start CPU time measurement
        write(6,*) 'Enter the system dimension: '
        read(5,*) n
        print *

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(b(n))
        call random_number(b)                   ! Fill vector b with random numbers
        write(6,*) 'The vector b is: ', b
        print *
        print*

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(x_old(n))
        x_old=0.0_dp                            ! Initialize x_old to zeros
        
        call jacobi(n,x_old,b)                  ! Solve the system using Jacobi method
        deallocate(b,x_old)
        
        call cpu_time(finish)
        print *
        write(6,*) 'The execution time is: ', finish - start, 'seconds'
        print *

end program main

