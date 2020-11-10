module collocations
    use precision
    use constants
    use input_functions
    use gauss_module
    use SLAU_solver
    implicit none
    contains 
    
subroutine legendre_coord_grid(a,b,n,x)
    real(mp), intent(in) :: a,b !interval borders
    integer(mp), intent(in) :: n !number of points
    real(mp), intent(out) :: x(0:n-1) !chebyshev nodes or grid
    integer(mp) :: k !counter
    call solution_legendre(n,x)
    x = x*(b-a)/2.0_mp + (a+b)/2.0_mp
end subroutine legendre_coord_grid

subroutine chebyshev_coord_grid(a,b,n,x)
    real(mp), intent(in) :: a,b !interval borders
    integer(mp), intent(in) :: n !number of points
    real(mp), intent(out) :: x(0:n-1) !chebyshev nodes or grid
    integer(mp) :: k !counter
    do k=1,n
        x(k-1)=(a+b)/2.0_mp + (b-a)*cos((2.0_mp*k-1.0_mp)*pi/2.0_mp/n)/2.0_mp
    enddo
end subroutine chebyshev_coord_grid
    
recursive function chebyshev_pol(n,x) result(res)
     !Chebyshev polynomials of the first kind     
     integer(mp), intent(in) :: n !number of points
     real(mp), intent(in) :: x !chebyshev node
     real(mp) :: res
     real(mp) :: t_n, t_n_prev !two previous chebyshev polynomials
     if (n.eq.0) then
         res = 1.0_mp
     elseif (n.eq.1) then
         res = x
     else
         t_n_prev = chebyshev_pol(n-2,x)
         t_n = chebyshev_pol(n-1,x)
         res = 2.0_mp*x*t_n - t_n_prev
    endif
end function chebyshev_pol
     
function KK(e) !K(e,x0)*K(e,t0)
    real(mp), intent(in) :: e
    real(mp) :: KK
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    KK = K(e,x0)*K(e,t0)
end function KK
     
function func_T(x,t)
    real(mp), intent(in) :: x,t
    real(mp) :: func_T
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    real(mp) :: a, b
    common /borders/ a, b
    x0 = x; t0=t
    func_T = gauss_integrator(a,b,KK) !integrate [K(e,x0)*K(e,t0)de] from a to b
end function func_T

function T_x_cheb(t) !T(x0,t)*legendre_polinom(j0,t)
    real(mp), intent(in) :: t
    real(mp) :: T_x_cheb
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    T_x_cheb = func_T(x0,t)*chebyshev_pol(j0,t)
end function T_x_cheb

function T_x_legendre(t) !T(x0,t)*legendre_polinom(j0,t)
    real(mp), intent(in) :: t
    real(mp) :: T_x_legendre
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    T_x_legendre = func_T(x0,t)*legendre_polynomial(j0,t)
end function T_x_legendre

function Kf(e) !K(e,x)*f(e)
    real(mp), intent(in) :: e
    real(mp) :: Kf
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    Kf = K(e,x0)*f(e)
end function Kf

function f1(x)
    real(mp), intent(in) :: x
    real(mp) :: f1
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    real(mp) :: a, b
    common /borders/ a, b
    x0 = x
    f1 = gauss_integrator(a,b,Kf) !integrate [K(e,x0)*f(e)de] from a to b
end function f1

function solution_y(x)
    real(mp), intent(in) :: x
    real(mp) :: solution_y
    integer(mp) :: j
    integer(mp) :: polynomial_number
    real(mp), pointer :: legendre_series_coeff_pointer(:)
    common /legendre_polinomial_series/ polynomial_number, legendre_series_coeff_pointer
    solution_y = 0.0_mp
      do j=1,polynomial_number ! for each coordinate function
        solution_y = solution_y + legendre_series_coeff_pointer(j)*legendre_polynomial(j,x)
      enddo
end function solution_y

function K_x_solution_y(t)
    real(mp), intent(in) :: t
    real(mp) :: K_x_solution_y
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
     K_x_solution_y = K(x0,t)*solution_y(t)
end function K_x_solution_y
        
subroutine collocations_solver(bottom_border,upper_border,n,alpha,x,y,residuals)
    ! it would be better to sort x,y,residuals by x in ascending order        
    real(mp), intent(in) :: bottom_border,upper_border !interval borders
    integer(mp), intent(in) :: n !number of grid points
    real(mp), intent(in) :: alpha !regularization parameter
    real(mp), intent(out) :: x(0:n-1), y(0:n-1), residuals(0:n-1) !solution
    real(mp) :: L(0:n-1,1:n) !IE operator
    real(mp) :: matrix_coef(1:n+1,0:n-1)
    real(mp), target :: legendre_series_coeff(1:n)
    integer(mp) :: i,j !counters
    
    !common blocks    
    integer(mp) :: j0
    real(mp) :: x0, t0
    common /grid_nodes/ j0, x0, t0
    real(mp) :: a, b
    common /borders/ a, b
    integer(mp) :: polynomial_number
    real(mp), pointer :: legendre_series_coeff_pointer(:)
    common /legendre_polinomial_series/ polynomial_number, legendre_series_coeff_pointer
    legendre_series_coeff_pointer => legendre_series_coeff
    
    call gauss_coefficents !generate file 'gauss_coef.dat' with coefficents for gauss integrating algorithm
    a = bottom_border; b = upper_border
    call legendre_coord_grid(a,b,n,x)
      do j=1,n ! for each coordinate function
        do i=0,n-1 !for each grid point
        j0 = j; x0=x(i)
        matrix_coef(j,i) = legendre_polynomial(j0,x0) + gauss_integrator(a,b,T_x_legendre)
      enddo
      enddo
    do i=0,n-1 !for each grid point
        matrix_coef(n+1,i) = f1(x(i))
    enddo 
    call lead(matrix_coef,legendre_series_coeff)
    polynomial_number = n
    
    do i=0,n-1 !for each grid point
        x0 = x(i)
        y(i) = solution_y(x0)
        residuals(i) = abs(gauss_integrator(a,b,K_x_solution_y) - f(x0))
    enddo
end subroutine collocations_solver


end module collocations






