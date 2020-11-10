module legendre_module
    use bernulli_module
    use precision
    contains
    
    function legendre_polynomial(n,x)
    integer(mp), intent(in) :: n
    real(mp), intent(in) :: x
    real(mp) :: legendre_polynomial
    real(mp) :: coef(0:n)
    integer :: i ! counter
    call legendre_koef(n,coef)
    legendre_polynomial = 0.0_mp
    do i=0,n
        legendre_polynomial = legendre_polynomial + coef(n-i)*x**i
    enddo
    end function legendre_polynomial

    recursive subroutine legendre_koef(n,coef)
    !процедура вычисляет коэфф. полинома Лежандра степени n
    implicit none
    integer(mp), intent(in) :: n
    real(mp), dimension(0:n), intent(out) :: coef
    real(mp), dimension(0:n-1) :: P1
    real(mp), dimension(0:n-2) :: P2
    integer :: i
    !---------------------------------------
    select case(n)
    case(0); coef(0)=1             ! P1(x)=1
    case(1); coef(0:1)=(/1.0,0.0/) ! P2(x)=x
    case default                ! полиномы степени больше 2-х ищем по реккурентной формуле
        call legendre_koef(n-1,P1)
        call legendre_koef(n-2,P2)
        coef(0)=(2*n-1.0)/n*P1(0)
        coef(1)=(2*n-1.0)/n*P1(1)
        forall (i=2:n-1) coef(i)=(2*n-1.0)/n*P1(i)-(n-1.0)/n*P2(i-2)
        coef(n)=-(n-1.0)/n*P2(n-2)
    end select

    end subroutine legendre_koef

    subroutine solution_legendre(n,X)
    !процедура вычисляет корни полинома Лежандра степени n
    implicit none
    integer(mp), intent(in) :: n
    real(mp), dimension(1:n), intent(out) :: X
    real(mp), dimension(0:n) :: coef

    select case(n)
    !корни полинома Лежандра степеней меньше 4
    case(1); X=0.0
    case(2); X=(/0.5773503,-0.5773503/)
    case(3); X=(/0.7745967,-0.7745967,0.0/)
    case default
        call legendre_koef(n,coef)
        call bernulli(coef,X)
    end select

    end subroutine solution_legendre


    end module legendre_module