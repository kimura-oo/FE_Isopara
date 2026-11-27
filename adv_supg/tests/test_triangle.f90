program test_triangle
    use iso_fortran_env, only: real64
    use triangle_quadrature, only: triangle_rule
    use triangle_elements, only: triangle_shape
    implicit none
    integer, parameter :: sizes(3) = [3,7,16], exact_degree(3) = [2,5,6]
    real(real64), allocatable :: xi(:), eta(:), w(:), values(:), deriv(:,:)
    real(real64) :: exact, integral
    integer :: degree, p, q, k, n

    do degree = 1,3
        n = sizes(degree)
        allocate(xi(n),eta(n),w(n))
        allocate(values((degree+1)*(degree+2)/2),deriv(2,(degree+1)*(degree+2)/2))
        call triangle_rule(n,xi,eta,w)
        if (any(w <= 0d0)) error stop 'quadrature weights must be positive'
        do p = 0,exact_degree(degree)
            do q = 0,exact_degree(degree)-p
                exact = gamma(dble(p+1))*gamma(dble(q+1))/gamma(dble(p+q+3))
                integral = sum(w*xi**p*eta**q)
                if (abs(exact-integral) > 5d-14) error stop 'incorrect polynomial integral'
            end do
        end do
        do k = 1,n
            call triangle_shape(degree,xi(k),eta(k),values,deriv)
            if (abs(sum(values)-1d0) > 5d-14) error stop 'shape functions do not sum to one'
            if (maxval(abs(sum(deriv,dim=2))) > 5d-14) error stop 'shape gradients do not sum to zero'
        end do
        deallocate(xi,eta,w,values,deriv)
    end do
end program
