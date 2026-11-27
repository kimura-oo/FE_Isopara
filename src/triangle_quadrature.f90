module triangle_quadrature
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: triangle_rule
contains
    ! Positive rules on (0,0), (1,0), (0,1); weights include the area 1/2.
    subroutine triangle_rule(npoints, xi, eta, weight)
        integer, intent(in) :: npoints
        real(real64), intent(out) :: xi(npoints), eta(npoints), weight(npoints)
        real(real64), parameter :: a = (6d0 + sqrt(15d0))/21d0
        real(real64), parameter :: b = (6d0 - sqrt(15d0))/21d0
        real(real64), parameter :: gauss(4) = [-0.8611363115940526d0, -0.3399810435848563d0, &
                                                0.3399810435848563d0,  0.8611363115940526d0]
        real(real64), parameter :: weights(4) = [0.3478548451374539d0, 0.6521451548625461d0, &
                                                 0.6521451548625461d0, 0.3478548451374539d0]
        integer :: i, j, k
        select case (npoints)
        case (3)
            xi = [1d0, 4d0, 1d0]/6d0
            eta = [1d0, 1d0, 4d0]/6d0
            weight = 1d0/6d0
        case (7)
            xi  = [1d0/3d0, a, 1d0-2*a, a, b, 1d0-2*b, b]
            eta = [1d0/3d0, a, a, 1d0-2*a, b, b, 1d0-2*b]
            weight(1) = 9d0/80d0
            weight(2:4) = (155d0 + sqrt(15d0))/2400d0
            weight(5:7) = (155d0 - sqrt(15d0))/2400d0
        case (16)
            ! A 4 x 4 Gauss rule on the square mapped by eta = (1-xi) v.
            do i = 1,4
                do j = 1,4
                    k = 4*(i-1) + j
                    xi(k) = (gauss(i)+1d0)/2d0
                    eta(k) = (1d0-xi(k))*(gauss(j)+1d0)/2d0
                    weight(k) = weights(i)*weights(j)*(1d0-xi(k))/4d0
                end do
            end do
        case default
            error stop 'triangle_rule: expected 3, 7 or 16 points'
        end select
    end subroutine
end module
