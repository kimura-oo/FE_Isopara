module triangle_elements
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: triangle_shape, edge_nodes
contains
    ! Barycentric multi-indices in mesh order: three vertices, nodes along
    ! edges 1->2, 2->3, 3->1, then the interior node for cubic triangles.
    pure subroutine triangle_shape(degree, xi, eta, values, deriv)
        integer, intent(in) :: degree
        real(real64), intent(in) :: xi, eta
        real(real64), intent(out) :: values(:), deriv(:,:)
        integer, parameter :: linear(3,3) = reshape([1,0,0, 0,1,0, 0,0,1], [3,3])
        integer, parameter :: quadratic(3,6) = reshape([2,0,0, 0,2,0, 0,0,2, 1,1,0, 0,1,1, 1,0,1], [3,6])
        integer, parameter :: cubic(3,10) = reshape([3,0,0, 0,3,0, 0,0,3, 2,1,0, 1,2,0, &
                                                     0,2,1, 0,1,2, 1,0,2, 2,0,1, 1,1,1], [3,10])
        integer :: i, k, axes(3)
        real(real64) :: bary(3), polynomial(0:3,3), derivative(0:3,3), factor(3), grad(3)
        if (degree < 1 .or. degree > 3) error stop 'unsupported triangle degree'
        if (size(values) /= (degree+1)*(degree+2)/2 .or. any(shape(deriv) /= [2,size(values)])) &
            error stop 'triangle shape array dimensions do not match degree'
        bary = [1.0_real64-xi-eta, xi, eta]
        polynomial(0,:) = 1
        derivative(0,:) = 0
        do k = 1, degree
            factor = (degree*bary - (k-1))/k
            polynomial(k,:) = polynomial(k-1,:)*factor
            derivative(k,:) = derivative(k-1,:)*factor + polynomial(k-1,:)*real(degree,real64)/k
        end do
        do i = 1, size(values)
            select case (degree)
            case (1)
                axes = linear(:,i)
            case (2)
                axes = quadratic(:,i)
            case (3)
                axes = cubic(:,i)
            end select
            factor = [polynomial(axes(1),1), polynomial(axes(2),2), polynomial(axes(3),3)]
            values(i) = product(factor)
            grad(1) = derivative(axes(1),1)*factor(2)*factor(3)
            grad(2) = derivative(axes(2),2)*factor(1)*factor(3)
            grad(3) = derivative(axes(3),3)*factor(1)*factor(2)
            deriv(:,i) = grad(2:3)-grad(1)
        end do
    end subroutine

    pure function edge_nodes(degree, side) result(nodes)
        integer, intent(in) :: degree, side
        integer :: nodes(degree+1), k
        nodes(1:2) = [side, mod(side,3)+1]
        do k = 1, degree-1
            nodes(2+k) = 3 + (side-1)*(degree-1) + k
        end do
    end function
end module
