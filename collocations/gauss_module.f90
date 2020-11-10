module gauss_module
    use precision
    use gauss_parameters
    use legendre_module
    use SLAU_solver
    implicit none
    contains

    subroutine gauss_coefficents
    ! *** Программа вычисляет коэфф. квадратурной формулы Гаусса
    integer(mp) :: n
    integer(mp) :: i
    real(mp), allocatable :: coef(:), nodes(:), M(:,:)
    n = gauss_n
    allocate(coef(1:n),nodes(n),M(1:n+1,0:n-1))
    call solution_legendre(n,nodes)
    forall (i=0:n-1) M(1:n,i)=nodes**(i)
    M(n+1,1:n-1:2)=0
    forall (i=0:n-1:2) M(n+1,i)=2.0/(i+1)
    call lead(M,coef)

    open(1,file='gauss_coef.dat')
        do i=1,n
            write(1,*) coef(i), nodes(i)
        enddo
    close(1)
    end subroutine gauss_coefficents
    
    function gauss_integrator(a,b,f)
    interface
        function f(x)
        integer, parameter :: mp=8
        real(mp), intent(in) :: x
        real(mp) :: f
        end function
    end interface
    real(mp), intent(in) :: a, b
    real(mp) :: gauss_integrator
    real(mp) :: coef(1:gauss_n), nodes(1:gauss_n)
    integer :: i
    
    open(2,file='gauss_coef.dat')
        do i=1,gauss_n
            read(2,*) coef(i), nodes(i)
        enddo
    close(2)
    !рассчет интеграла по формуле Гаусса
    gauss_integrator = 0.0_mp
    do i=1,gauss_n
        gauss_integrator = gauss_integrator+coef(i)*(b-a)/2.0_mp*f(nodes(i)*(b-a)/2.0_mp+(a+b)/2.0_mp)
    enddo
    end function gauss_integrator

end module gauss_module
