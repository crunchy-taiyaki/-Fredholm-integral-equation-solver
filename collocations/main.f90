program main
    use collocations
    use input_functions
    use precision
implicit none
real(mp) :: bottom_border,upper_border !interval borders
integer(mp) :: n
real(mp) :: alpha
real(mp), allocatable :: x(:), y(:), residuals(:)
integer :: i ! counter
bottom_border = 0.0_mp
upper_border = 1.0_mp
n = 10
alpha = 0.001
allocate(x(0:n-1))
allocate(y(0:n-1))
allocate(residuals(0:n-1))
call collocations_solver(bottom_border,upper_border,n,alpha,x,y,residuals)
do i=0,n-1
  write(*,*) x(i), y(i), residuals(i)
enddo
deallocate(x,y,residuals)   
end program