module supg_numerics
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use triangle_elements, only: triangle_shape
    implicit none
    private
    public :: map_triangles, assemble_transport, scale_diagonal, solve_bicgstab
    public :: integrated_mass, apply_element_matrix
contains
    subroutine map_triangles(degree, coords, nodes, xi, eta, shape, dx, dy, det)
        integer, intent(in) :: degree, nodes(:,:)
        real(real64), intent(in) :: coords(:,:), xi(:), eta(:)
        real(real64), intent(out) :: shape(:,:,:), dx(:,:,:), dy(:,:,:), det(:,:)
        real(real64) :: values(size(nodes,1)), deriv(2,size(nodes,1)), j(2,2)
        integer :: m, q
        if (any(nodes < 1) .or. any(nodes > size(coords,2))) error stop 'invalid mesh connectivity'
        do q = 1,size(xi)
            call triangle_shape(degree, xi(q), eta(q), values, deriv)
            do m = 1,size(nodes,2)
                shape(:,q,m) = values
                j = matmul(coords(:,nodes(:,m)), transpose(deriv))
                det(q,m) = j(1,1)*j(2,2) - j(1,2)*j(2,1)
                if (.not. (det(q,m) > 0d0)) then
                    write(*,*) 'Invalid Jacobian at element and integration point', m, q
                    error stop 'inverted or degenerate element'
                end if
                dx(:,q,m) = (j(2,2)*deriv(1,:) - j(2,1)*deriv(2,:))/det(q,m)
                dy(:,q,m) = (j(1,1)*deriv(2,:) - j(1,2)*deriv(1,:))/det(q,m)
            end do
        end do
    end subroutine

    ! Assemble both sides of the theta scheme from the same stabilized weak form.
    subroutine assemble_transport(nodes, u, v, dt, theta, tau, diffusivity, shape, dx, dy, det, weights, &
                                  matrix, diagonal, previous)
        integer, intent(in) :: nodes(:,:)
        real(real64), intent(in) :: u(:), v(:), dt, theta, shape(:,:,:), dx(:,:,:), dy(:,:,:), det(:,:), weights(:)
        real(real64), intent(in) :: tau(:,:), diffusivity
        real(real64), intent(out) :: matrix(:,:,:), diagonal(:), previous(:,:,:)
        real(real64) :: adv(size(nodes,1)), test(size(nodes,1)), uq, vq, jw, mass, stiffness
        integer :: m, q, i, j, ni
        diagonal = 0d0
        matrix = 0d0
        previous = 0d0
        do m = 1,size(nodes,2)
            do q = 1,size(weights)
                uq = dot_product(shape(:,q,m), u(nodes(:,m)))
                vq = dot_product(shape(:,q,m), v(nodes(:,m)))
                adv = uq*dx(:,q,m) + vq*dy(:,q,m)
                test = shape(:,q,m) + tau(q,m)*adv
                jw = det(q,m)*weights(q)
                do j = 1,size(nodes,1)
                    do i = 1,size(nodes,1)
                        mass = test(i)*shape(j,q,m)*jw/dt
                        stiffness = (test(i)*adv(j) + diffusivity*(dx(i,q,m)*dx(j,q,m) + dy(i,q,m)*dy(j,q,m)))*jw
                        matrix(i,j,m) = matrix(i,j,m) + mass + theta*stiffness
                        previous(i,j,m) = previous(i,j,m) + mass - (1d0-theta)*stiffness
                    end do
                end do
            end do
            do i = 1,size(nodes,1)
                ni = nodes(i,m)
                diagonal(ni) = diagonal(ni) + matrix(i,i,m)
            end do
        end do
    end subroutine

    subroutine scale_diagonal(nodes, diagonal, matrix, square_root)
        integer, intent(in) :: nodes(:,:)
        real(real64), intent(inout) :: diagonal(:), matrix(:,:,:)
        logical, intent(in) :: square_root
        integer :: i, j, m
        if (square_root) then
            if (any(diagonal <= 0d0)) error stop 'nonpositive diagonal in symmetric scaling'
            diagonal = 1d0/sqrt(diagonal)
        else
            if (any(diagonal == 0d0)) error stop 'zero diagonal in scaling'
            diagonal = 1d0/diagonal
        end if
        do m = 1,size(nodes,2)
            do j = 1,size(nodes,1)
                do i = 1,size(nodes,1)
                    matrix(i,j,m) = matrix(i,j,m)*diagonal(nodes(i,m))*diagonal(nodes(j,m))
                end do
            end do
        end do
    end subroutine

    subroutine apply_element_matrix(nodes, matrix, vector, fixed, product)
        integer, intent(in) :: nodes(:,:), fixed(:)
        real(real64), intent(in) :: matrix(:,:,:), vector(:)
        real(real64), intent(out) :: product(:)
        integer :: m
        product = 0d0
        do m = 1,size(nodes,2)
            product(nodes(:,m)) = product(nodes(:,m)) + matmul(matrix(:,:,m), vector(nodes(:,m)))
        end do
        product(fixed) = 0d0
    end subroutine

    ! The legacy eps is a tolerance on the squared, scaled residual.
    subroutine solve_bicgstab(eps, limit, diagonal, iterations, solution, rhs, nodes, matrix, fixed)
        real(real64), intent(in) :: eps, diagonal(:), rhs(:), matrix(:,:,:)
        integer, intent(in) :: limit, nodes(:,:), fixed(:)
        integer, intent(out) :: iterations
        real(real64), intent(inout) :: solution(:)
        real(real64), dimension(size(rhs)) :: x, residual, shadow, direction, image, t, at
        real(real64) :: tolerance, residual_sq, rho, rho_old, alpha, omega, denom
        integer :: k
        x = solution/diagonal
        residual = rhs*diagonal
        residual(fixed) = 0d0
        tolerance = eps*dot_product(residual, residual)
        call apply_element_matrix(nodes, matrix, x, fixed, image)
        residual = residual - image
        residual_sq = dot_product(residual, residual)
        if (tolerance == 0d0) tolerance = eps*residual_sq
        shadow = residual
        direction = residual
        rho = residual_sq
        iterations = 0
        do k = 1,limit
            if (residual_sq <= tolerance) exit
            call apply_element_matrix(nodes, matrix, direction, fixed, image)
            denom = dot_product(shadow,image)
            if (denom == 0d0 .or. rho == 0d0) error stop 'BiCGSTAB breakdown in alpha'
            alpha = rho/denom
            t = residual - alpha*image
            x = x + alpha*direction
            iterations = k
            residual_sq = dot_product(t,t)
            if (residual_sq <= tolerance) exit
            call apply_element_matrix(nodes, matrix, t, fixed, at)
            denom = dot_product(at,at)
            if (denom == 0d0) error stop 'BiCGSTAB breakdown in omega'
            omega = dot_product(at,t)/denom
            if (omega == 0d0) error stop 'BiCGSTAB zero omega'
            x = x + omega*t
            residual = t - omega*at
            residual_sq = dot_product(residual,residual)
            rho_old = rho
            rho = dot_product(shadow,residual)
            direction = residual + (alpha/omega)*(rho/rho_old)*(direction-omega*image)
        end do
        if (.not.ieee_is_finite(residual_sq) .or. residual_sq > tolerance) error stop 'BiCGSTAB did not converge'
        solution = x*diagonal
    end subroutine

    subroutine integrated_mass(nodes, field, shape, det, weights, mass, absolute_mass)
        integer, intent(in) :: nodes(:,:)
        real(real64), intent(in) :: field(:), shape(:,:,:), det(:,:), weights(:)
        real(real64), intent(out) :: mass, absolute_mass
        integer :: m, q
        real(real64) :: contribution
        mass = 0d0
        absolute_mass = 0d0
        do m = 1,size(nodes,2)
            do q = 1,size(weights)
                contribution = dot_product(shape(:,q,m), field(nodes(:,m)))*det(q,m)*weights(q)
                mass = mass + contribution
                absolute_mass = absolute_mass + abs(contribution)
            end do
        end do
    end subroutine

end module
