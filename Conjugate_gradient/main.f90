! ==============================================================================
! Program    : main
! Purpose    : Solve a linear system A*x = b using the Conjugate Gradient method
! Author     : [Pietro Gozzoli, Alessandro Michieletti]
! Modules    : conjugate_gradient_solver, utils (general utilities)
! Description: 
!     - Reads system size from user
!     - Generates random right-hand side vector b
!     - Generates random initial guess x_guess
!     - Solves the system with the conjugate gradient method
!     - Prints the execution time
! ==============================================================================

program main

        use conjugate_gradient_solver                   ! Module with Conjugate Gradient solver
        use utils                                       ! Module with general utilities (e.g., precision definition)
        implicit none
        integer :: n,i
        real(dp) , allocatable :: b(:) , x_guess(:)
        real(dp) :: start,finish

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call cpu_time(start)                            ! Record CPU start time
        write(6,*) 'Enter the system dimension: '            ! Ask user to input the system size
        read(5,*) n
        print *

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(b(n))                                  ! Allocate and generate random vector b
        call random_number(b)
        write(6,*) 'The vector b is: ', b
        print *

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        allocate(x_guess(n))                            ! Allocate and generate random initial guess x_guess
        call random_number(x_guess)
        print *
        write(6,*) 'The initial guess vector is: ', x_guess
        print *

!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------        
        call conj_grad(n,x_guess,b)                     ! Solve A*x = b using Conjugate Gradient method from the solver module
        deallocate(x_guess,b)
        call cpu_time(finish)
        print *
        write(6,*) 'The execution time is: ', finish-start, 'seconds'
        print *

end program main

