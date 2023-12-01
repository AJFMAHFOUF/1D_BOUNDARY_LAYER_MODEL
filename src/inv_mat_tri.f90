subroutine inv_mat_tri (n,a,b,c,d,x)
 implicit none
 integer :: i
 integer, intent (in) :: n
 real :: w
 real, dimension (n), intent(inout) :: a,b,c,d
 real, dimension (n), intent(out) :: x
 a(1) = 0.0
 c(n) = 0.0
 do i = 2, n
   w = a(i)/b(i-1)
   b(i) = b(i) - w*c(i-1)
   d(i) = d(i) - w*d(i-1)
 enddo
 x(n) = d(n)/b(n)
 do i = n-1,1,-1
   x(i) = (d(i) - c(i)*x(i+1))/b(i)
 enddo
 return
end subroutine inv_mat_tri

