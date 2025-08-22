module ortho
  use utils, only : dp
  implicit none
  real(dp), parameter :: zero = 0.0_dp, one = 1.0_dp, two = 2.0_dp
  real(dp), parameter :: tol_ortho = two * epsilon(one)
!
  contains
!
  subroutine ortho_vs_x(n,m,k,x,u)
    implicit none
!
!   given two sets x(n,m) and u(n,k) of vectors, where x 
!   is assumed to be orthogonal, orthogonalize u against x.
!
!   furthermore, orthonormalize u.
!
!   this routine performs the u vs x orthogonalization and the
!   subsequent orthonormalization of u iteratively, until the
!   overlap between x and the orthogonalized u is smaller than
!   a (tight) threshold. 
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m, k
    real(dp),  dimension(n,m), intent(in)    :: x
    real(dp),  dimension(n,k), intent(inout) :: u
!
!   local variables:
!   ================
!
    logical                :: done, ok
    integer                :: it, info
    real(dp)               :: xu_norm, growth
    real(dp),  allocatable :: xu(:,:)
!
!   external functions:
!   ===================
!
    real(dp)               :: dnrm2
    external               :: dnrm2, dgemm
!   
    integer, parameter     :: maxit = 10
!
!   allocate space for the overlap between x and u.
!
    ok = .false.
    allocate (xu(m,k))
    done = .false.
    it   = 0
!
!   start with an initial orthogonalization to improve conditioning.
!
    call ortho_cd(n,k,u,growth,ok)
!
!   iteratively orthogonalize u against x, and then orthonormalize u.
!
    do while (.not. done)
      it = it + 1
!
!     u = u - x (x^t u)
!
      call dgemm('t','n',m,k,n,one,x,n,u,n,zero,xu,m)
      call dgemm('n','n',n,k,m,-one,x,n,xu,m,one,u,n)
!
!     now, orthonormalize u.
!
      call ortho_cd(n,k,u,growth,ok)
!
!     the orthogonalization has introduced an error that makes the new
!     vector no longer fully orthogonal to x. assuming that u was 
!     orthogonal to x to machine precision before, we estimate the 
!     error with growth * eps, where growth is the product of the norms
!     of all the linear transformations applied to u.
!     if ortho_cd has failed, we just compute the overlap and its norm.
!
      xu_norm = growth * epsilon(one)
      done    = xu_norm.lt.tol_ortho
!
!     if things went really wrong, abort.
!
      if (it.gt.maxit) stop ' catastrophic failure of ortho_vs_x'
    end do
!
    deallocate(xu)
!
    return
  end subroutine ortho_vs_x
!
  subroutine ortho_cd(n,m,u,growth,ok)
    implicit none
!
!   orthogonalize m vectors of lenght n using the Cholesky factorization
!   of their overlap. 
!   this is done by metric = U^t U and then by computing its cholesky 
!   decompositoin metric = L L^t. The orthogonal vectors are obtained then
!   by solving the triangular linear system
!
!     U(ortho)L^T = U
!
!   as cholesky decomposition is not the most stable way of orthogonalizing
!   a set of vectors, the orthogonalization is refined iteratively. 
!   a conservative estimate of the orthogonalization error is used to 
!   assess convergence. 
!
!   this routine returns a growth factor, which can be used in (b_)ortho_vs_x
!   to estimate the orthogonality error introduced by ortho_cd.
!
!   while it is very unlikely to do so, this routine can fail. 
!   a logical flag is then set to false, so that the calling program can 
!   call a more robust orthogonalization routine without aborting.
!
!   arguments:
!   ==========
!
    integer,                   intent(in)    :: n, m
    real(dp),  dimension(n,m), intent(inout) :: u
    real(dp),                  intent(inout) :: growth
    logical,                   intent(inout) :: ok
!
!   local variables
!   ===============
!
    integer               :: it, it_micro, info
    real(dp)              :: error, dnrm2, alpha, unorm, shift
    real(dp)              :: rcond, l_norm, linv_norm
    logical               :: macro_done, micro_done
    real(dp), parameter   :: tol_ortho_cd = two * epsilon(one)
    integer,  parameter   :: maxit = 10
!
!   local scratch
!   =============
!
    real(dp), allocatable :: metric(:,:), msave(:,:)
!
!   external functions:
!   ===================
!
    external              :: dgemm, dpotrf, dtrsm, dtrmm, dtrtri
!
!   get memory for the metric.
!
    allocate (metric(m,m), msave(m,m))
    metric = zero
    macro_done = .false.
!
!   assemble the metric
!
    it = 0
    growth = one
    do while(.not. macro_done)
      it = it + 1
      if (it .gt. maxit) then
!
!       ortho_cd failed. return with an error message
!
        ok = .false.
        write(6,100) ' maximum number of iterations reached.'
        return
      end if
      call dgemm('t','n',m,m,n,one,u,n,u,n,zero,metric,m)
      msave = metric
!
!   compute the cholesky factorization of the metric.
!
      call dpotrf('l',m,metric,m,info)
!
!     if dpotrf failed, try a second time, after level-shifting the diagonal of the metric.
!
      if (info.ne.0) then
!
        alpha      = 100.0_dp
        unorm      = dnrm2(n*m,u,1)
        it_micro   = 0
        micro_done = .false.
!
!       add larger and larger shifts to the diagonal until dpotrf manages to factorize it.
!
        do while (.not. micro_done)
          it_micro = it_micro + 1
          if (it_micro.gt.maxit) then
!
!           something went very wrong. return with an error status, the orthogonalization
!           will be carried out using a different algorithm.
!
            ok = .false.
            write(6,100) ' maximum number of iterations for factorization reached.'
            stop
            return
          end if
!
          shift = max(epsilon(one)*alpha*unorm,tol_ortho)
          metric = msave
          call diag_shift(m,shift,metric)
          call dpotrf('l',m,metric,m,info)
          alpha = alpha * 10.0_dp
          micro_done = info.eq.0
        end do
!
      end if
!
!     we assume that the error on the orthogonality is of order k(l)^2 * eps,
!     where eps is the machine precision. 
!     the condition number k(l) is estimated by computing 
!
!     k(l) ||l|| ||l^-1||,
!
!     where the norm used is the following (see norm_estimate):
!
!     || A || = || D + O || <= || D ||_inf + || O ||_2
!
!     compute l^-1, using msave to store the inverse cholesky factor
!
      msave = metric
      call dtrtri('l','n',m,msave,m,info)
!
!     compute the norm of l, l^-1 and the condition number:
!
      l_norm    = norm_est(m,metric)
      linv_norm = norm_est(m,msave)
      rcond     = l_norm * linv_norm
!
!     in each iteration of ortho_cd, we apply l^-t to u, which introduces
!     a numerical error of order ||l^-1||. 
!     this error is saved in growth and used in ortho_vs_x to check how much
!     ortho_cd spoiled the previously computed orthogonality to x. 
!
      growth    = growth * linv_norm
!
!     orthogonalize u by applying l^(-t)
!    
      call dtrmm('r','l','t','n',n,m,one,msave,m,u,n)
!
!     check the error:
!
      error      = epsilon(one) * rcond*rcond
      macro_done = error .lt. tol_ortho_cd
    end do
!
    100 format(t3,'ortho_cd failed with the following error:',a)
!
    ok = .true.
!
    deallocate (metric, msave)
    return
  end subroutine ortho_cd
!
  real(dp) function norm_est(m,a)
!
!   compute a cheap estimate of the norm of a lower triangular matrix.
!   let a = d + o, where d = diag(a). as
!   
!   || a || <= || d || + || o ||
!
!   we compute || d || as max_i |d(i)| and || o || as its frobenius norm.
!   this is tight enough, and goes to 1 when a approaches the identity.
!
    implicit none
    integer,                  intent(in) :: m
    real(dp), dimension(m,m), intent(in) :: a
!
    integer  :: i, j
    real(dp) :: diag_norm, od_norm
!
    diag_norm = zero
    do i = 1, m
      diag_norm = max(diag_norm,abs(a(i,i)))
    end do
!
    od_norm = zero
    do i = 1, m
      do j = 1, i - 1
        od_norm = od_norm + a(i,j)**2
      end do
    end do
    od_norm = sqrt(od_norm)
!
    norm_est = diag_norm + od_norm
    return
  end function norm_est
!
  subroutine diag_shift(n,shift,a)
    implicit none
!  
!   add shift to the diagonal elements of the matric a
!  
    integer,                  intent(in)    :: n
    real(dp),                 intent(in)    :: shift
    real(dp), dimension(n,n), intent(inout) :: a
!  
    integer :: i
!  
    do i = 1, n
      a(i,i) = a(i,i) + shift
    end do
!  
    return
  end subroutine diag_shift
!
end module ortho

