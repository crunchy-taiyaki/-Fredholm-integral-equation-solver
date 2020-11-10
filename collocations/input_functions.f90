module input_functions
    use precision
    implicit none
    contains
    
    real function K(x,t) !core
    real(mp), intent(in) :: x,t
    K = log(2.0_mp+x*t)
    end function K
    
    real function f(x)
    real(mp), intent(in) :: x
    f = (2.0_mp+x)*log(4.0_mp/(2.0_mp+x))/(2.0_mp*(2.0_mp-x))
    end function f
    
    !real function K(x,t) !core
    !real(mp), intent(in) :: x,t
    !K = x*t*1.0_mp
    !end function K
    !
    !real function f(x)
    !real(mp), intent(in) :: x
    !f = x*1.0_mp/3.0_mp
    !end function f  

    
end module input_functions