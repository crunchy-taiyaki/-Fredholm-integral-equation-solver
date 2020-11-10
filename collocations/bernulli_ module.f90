module bernulli_module
    use precision
    contains

    subroutine gorner(A,x0,B)
    !процедура реализует деление многочлена(А - массив коэффициентов) на (x-x0),
    !в результате имеем массив коэф. В, где последний элемент - остаток от деления
    implicit none
    real(mp), dimension(0:), intent(in) :: A
    real(mp), intent(in) :: x0 ! x0 - уже известный корень
    real(mp), dimension(0:size(A)-1), intent(out) :: B
    integer :: n, i

    n=size(A)-1

    B(0)=A(0)
    do i=1,n
    B(i)=x0*B(i-1)+A(i)
    enddo

    end subroutine gorner

    subroutine bernulli(A0,X)
    !процедура возвращает массив корней(Х) полинома(А - массив коэффициентов) методом Бернулли
    implicit none
    real(mp), dimension(0:), intent(in) :: A0
    real(mp), dimension(0:(size(A0)-1)/2) :: A, B ! Массивы коэффициентов полиномов от t=x^2
    real(mp), dimension(1:size(A0)-1), intent(out) :: X
    real(mp), dimension(1:size(A0)-1) :: Y ! Массив из чисел, пределы которых ищутся в методе Бернулли
    real(mp) :: x0,random ! x0 - корень многочлена
    integer :: n, i, j

    n=size(A0)-1


    if (mod(n,2)/=0) then
        X(n)=0.0
    endif

    A(0:n/2)=A0(0:n:2)
    i=0
    do while (i<n/2-2)
        call random_number(Y)
        call solution(A(0:n/2-i),Y(1:n/2-i),x0)
        X(2*i+1)=sqrt(x0); X(2*i+2)=-sqrt(x0)
        call gorner(A(0:n/2-i),x0,B(0:n/2-i))
        A(0:n/2-i-1)=B(0:n/2-i-1)
        i=i+1
    enddo
    !решаем квадратное уравнение
    X(2*i+1)=sqrt((-A(1)+sqrt(A(1)**2-4*A(0)*A(2)))/2/A(0)); X(2*i+2)=-X(2*i+1)
    X(2*i+3)=sqrt((-A(1)-sqrt(A(1)**2-4*A(0)*A(2)))/2/A(0)); X(2*i+4)=-X(2*i+3)

    contains

    subroutine solution(A,Y,x0)
    !ищем наибольший по модулю корень многочлена
    implicit none
    real(mp), dimension(0:), intent(in) :: A
    real(mp) :: x0 ! x0 - искомый корень
    real(mp), dimension(1:size(A)) :: Y
    integer :: i, n
    real(mp) :: eps=0.1**mp*0.1

    n=size(A)-1

    do while (abs(Y(n)/Y(n-1)-Y(n-1)/Y(n-2)) > eps)
        Y(n+1)=sum((A(1:n)/(-A(0)))*Y(n:1:-1))
        forall (i=1:n-1) Y(i)=Y(i+1)
        Y(n)=Y(n+1)
    enddo

    x0=Y(n)/Y(n-1)

    end subroutine solution

    end subroutine bernulli

    end module bernulli_module
