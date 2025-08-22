module utils
  implicit none
  integer,  parameter :: dp=kind(1.d0)
!
  real(dp), allocatable :: a(:,:)
!
  contains
!
  subroutine define_a(n)
    implicit none
!
!   n is input only: it is therefore declared as "intent in"
!
    integer, intent(in) :: n
!
    integer  :: istat
    integer  :: i, j
    real(dp) :: fac
!
    allocate (a(n,n), stat=istat)
    if (istat.ne.0) then
      write(6,*) 'allocation error in define_a.'
      stop
    end if
!
    do j = 1, n
!
!     real(i,kind) converts i into a real number of the desired kind.
!
      a(j,j) = real(j,dp) + 1.0_dp
      do i = 1, j - 1
        fac = 1.0_dp/real(i+j,dp)
        a(j,i) = fac
        a(i,j) = fac
      end do
    end do
    return
!
  end subroutine define_a
!
  subroutine free_a
    implicit none
    integer  :: istat
!
    deallocate (a, stat=istat)
    if (istat.ne.0) then
      write(6,*) 'deallocation error in free_a.'
      stop
    end if
    return
  end subroutine free_a
!
end module utils

