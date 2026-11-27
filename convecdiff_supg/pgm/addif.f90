program convection_diffusion
    use iso_fortran_env, only: real64
    use triangle_quadrature, only: triangle_rule
    use supg_io, only: supg_case, read_supg_case, write_solution
    use supg_numerics, only: map_triangles, assemble_transport, scale_diagonal, solve_bicgstab, apply_element_matrix
    implicit none
    integer, parameter :: degree = 1, nq = 3, nlocal = 3
    type(supg_case) :: c
    real(real64) :: xi(nq), eta(nq), weights(nq), old_time, new_time
    real(real64), allocatable :: shape(:,:,:), dx(:,:,:), dy(:,:,:), det(:,:), tau(:,:)
    real(real64), allocatable :: matrix(:,:,:), previous(:,:,:), diagonal(:), rhs(:), field(:), cosine_load(:), sine_load(:)
    integer :: nnode, nelm, step, result, iterations

    call read_supg_case(degree,.true.,c)
    nelm = size(c%nodes,2)
    nnode = size(c%initial)
    allocate(shape(nlocal,nq,nelm),dx(nlocal,nq,nelm),dy(nlocal,nq,nelm),det(nq,nelm),tau(nq,nelm))
    allocate(matrix(nlocal,nlocal,nelm),previous(nlocal,nlocal,nelm),diagonal(nnode),rhs(nnode))
    allocate(cosine_load(nnode),sine_load(nnode))
    call triangle_rule(nq,xi,eta,weights)
    call map_triangles(degree,c%coords,c%nodes,xi,eta,shape,dx,dy,det)
    call diffusion_stabilization()
    call assemble_transport(c%nodes,c%u,c%v,c%dt,0.5d0,tau,c%diffusion,shape,dx,dy,det,weights,matrix,diagonal,previous)
    call scale_diagonal(c%nodes,diagonal,matrix,square_root=.false.)
    call manufactured_load()
    field = c%initial
    open(newunit=result,file=c%result_file,status='replace',action='write')
    if (c%first_step == 1) call write_solution(result,c,0,field,.false.)
    do step = c%first_step,c%last_step
        old_time = (step-1)*c%dt
        new_time = step*c%dt
        call apply_element_matrix(c%nodes,previous,field,c%fixed,rhs)
        rhs = rhs + 0.5d0*(cos(old_time)+cos(new_time))*cosine_load + 0.5d0*(sin(old_time)+sin(new_time))*sine_load
        call solve_bicgstab(c%eps,c%max_iterations,diagonal,iterations,field,rhs,c%nodes,matrix,c%fixed)
        field(c%fixed) = c%initial(c%fixed)
        call write_solution(result,c,step,field,.false.)
    end do
    close(result)

contains
    subroutine diffusion_stabilization()
        real(real64) :: u, v, speed, inverse_time, h
        integer :: m, q
        do m = 1,nelm
            u = sum(c%u(c%nodes(:,m)))/nlocal
            v = sum(c%v(c%nodes(:,m)))/nlocal
            speed = hypot(u,v)
            do q = 1,nq
                inverse_time = sum(abs(u*dx(:,q,m)+v*dy(:,q,m)))
                if (inverse_time == 0d0) then
                    tau(q,m) = 0d0
                else
                    h = 2d0*speed/inverse_time
                    tau(q,m) = 1d0/hypot(inverse_time,4d0*c%diffusion/h**2)
                end if
            end do
        end do
    end subroutine

    ! Manufactured solution sin(t) sin(pi*x) sin(pi*y), with a prescribed velocity.
    subroutine manufactured_load()
        real(real64), parameter :: pi = acos(-1d0)
        real(real64) :: x, y, vx, vy, s, adv_s, test(nlocal), wj
        integer :: m,q,ids(nlocal)
        cosine_load = 0d0
        sine_load = 0d0
        do m = 1,nelm
            ids = c%nodes(:,m)
            do q = 1,nq
                x = dot_product(shape(:,q,m),c%coords(1,ids))
                y = dot_product(shape(:,q,m),c%coords(2,ids))
                vx = dot_product(shape(:,q,m),c%u(ids))
                vy = dot_product(shape(:,q,m),c%v(ids))
                s = sin(pi*x)*sin(pi*y)
                adv_s = pi*(vx*cos(pi*x)*sin(pi*y) + vy*sin(pi*x)*cos(pi*y))
                test = shape(:,q,m) + tau(q,m)*(vx*dx(:,q,m)+vy*dy(:,q,m))
                wj = det(q,m)*weights(q)
                cosine_load(ids) = cosine_load(ids) + test*s*wj
                sine_load(ids) = sine_load(ids) + test*(2d0*c%diffusion*pi**2*s + adv_s)*wj
            end do
        end do
    end subroutine
end program
