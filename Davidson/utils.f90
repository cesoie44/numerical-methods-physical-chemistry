module utils
  implicit none
  integer, parameter :: dp = 8
  
contains


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        

!The subroutine calcAV(n, V, W) computes the product of a vector V with an implicit matrix A, storing the result in W, without ever forming the matrix A explicitly. The matrix is symmetric, with dominant diagonal elements defined as  A_(ii)=i+1, and off-diagonal elements given by A_(ij)=1/(i+j) within a bandwidth of Delta. This makes A symmetric, banded, and diagonally dominant. The use of an implicit matrix significantly reduces memory and computational cost, which is especially useful for large-scale problems like those tackled by the Davidson method.
        
subroutine calcAV(n,V,W)

        implicit none
        integer , intent(in) :: n
        real(dp) , intent(in) :: V(n)
        real(dp) , intent(out) :: W(n)

        integer :: i,j,delta
        real(dp) :: sum


        delta=1000

        w=0.0_dp

        do i=1, n
                sum = V(i) * real(i+1, dp)
                do j = max(1, i - delta), min(n, i + delta)
                        if (i.ne.j) then
                                sum = sum + V(j) / real(i + j, dp)
                        end if
                end do

                W(i) = sum
        end do
        
end subroutine calcAV


!---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


subroutine convergenza(n,R,tau,logic,n_eig)
        
         implicit none               
         integer , intent(in) :: n,n_eig
         real(dp) , intent(in) :: R(n,n_eig)
         real(dp) :: tau,dnrm2,normr,norminf
         logical, intent(out) :: logic

         logical :: converged
         integer :: i

         converged=.TRUE.

         do i=1,n_eig       
                normr=dnrm2(n, R(:,i), 1)
                write(6,*) 'The 2-norm of column',i,'is', normr
                norminf=maxval(abs(R(:,i)))
                write(6,*) 'The infinity norm of column',i,'is', norminf
                if (.NOT.(normr .LE. tau .AND. norminf .LE. tau)) then
                        converged=.FALSE.
                end if
         end do
         logic=converged
  end subroutine convergenza       


end module utils
