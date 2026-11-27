module supg_advection
    use iso_fortran_env, only: real64
    use triangle_quadrature, only: triangle_rule
    use supg_io, only: supg_case, read_supg_case, write_solution, write_section
    use supg_numerics, only: map_triangles, assemble_transport, scale_diagonal, solve_bicgstab, &
                             apply_element_matrix, integrated_mass
    implicit none
    private
    public :: run_advection
contains
    subroutine run_advection(degree)
        integer, intent(in) :: degree
        type(supg_case) :: c
        integer, parameter :: rule_sizes(3) = [3,7,16]
        integer :: nq, nlocal, nelm, nnode, step, iterations, result, mass_file, section
        real(real64), allocatable :: xi(:),eta(:),weights(:),shape(:,:,:),dx(:,:,:),dy(:,:,:),det(:,:),tau(:,:)
        real(real64), allocatable :: matrix(:,:,:), previous(:,:,:), diagonal(:), field(:), rhs(:)
        real(real64) :: initial_mass, absolute_mass, mass
        logical :: relative_mass

        call read_supg_case(degree,.false.,c)
        nq = rule_sizes(degree)
        nlocal = size(c%nodes,1)
        nelm = size(c%nodes,2)
        nnode = size(c%initial)
        allocate(xi(nq),eta(nq),weights(nq),shape(nlocal,nq,nelm),dx(nlocal,nq,nelm),dy(nlocal,nq,nelm),det(nq,nelm))
        allocate(matrix(nlocal,nlocal,nelm),previous(nlocal,nlocal,nelm),diagonal(nnode),rhs(nnode))
        call triangle_rule(nq,xi,eta,weights)
        call map_triangles(degree,c%coords,c%nodes,xi,eta,shape,dx,dy,det)
        allocate(tau(nq,nelm),source=c%dt/2d0)
        call assemble_transport(c%nodes,c%u,c%v,c%dt,0.5d0,tau,0d0,shape,dx,dy,det,weights,matrix,diagonal,previous)
        ! Retain each driver's original diagonal preconditioner.
        call scale_diagonal(c%nodes,diagonal,matrix,square_root=degree /= 3)

        field = c%initial
        call integrated_mass(c%nodes,field,shape,det,weights,initial_mass,absolute_mass)
        relative_mass = abs(initial_mass) > 64d0*epsilon(initial_mass)*absolute_mass
        open(newunit=result,file=c%result_file,status='replace',action='write')
        open(newunit=mass_file,file=c%mass_file,status='replace',action='write')
        open(newunit=section,file=c%sample_file,status='replace',action='write')
        if (.not.relative_mass) write(mass_file,'(a)') '# step  mass(t) - mass(0) (zero initial integral)'
        print *, 'Initial mass integral', initial_mass
        if (c%first_step == 1) call write_solution(result,c,0,field,.true.)
        do step = c%first_step,c%last_step
            call apply_element_matrix(c%nodes,previous,field,c%fixed,rhs)
            call solve_bicgstab(c%eps,c%max_iterations,diagonal,iterations,field,rhs,c%nodes,matrix,c%fixed)
            ! Reinstate prescribed values exactly after scaling roundoff.
            field(c%fixed) = c%initial(c%fixed)
            call integrated_mass(c%nodes,field,shape,det,weights,mass,absolute_mass)
            if (relative_mass) then
                write(mass_file,*) step, 100d0*(mass/initial_mass)
            else
                write(mass_file,*) step, mass-initial_mass
            end if
            if (step == c%sample_step) call write_section(section,c%coords,field)
            call write_solution(result,c,step,field,.true.)
        end do
        close(section)
        close(mass_file)
        close(result)
    end subroutine
end module
