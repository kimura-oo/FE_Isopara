module supg_io
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: supg_case, read_supg_case, write_solution, write_section
    type supg_case
        integer :: first_step, last_step, max_iterations, sample_step
        real(real64) :: dt, eps, diffusion = 0d0
        character(:), allocatable :: result_file, vtk_prefix, sample_file, mass_file
        real(real64), allocatable :: coords(:,:), u(:), v(:), initial(:)
        integer, allocatable :: nodes(:,:), fixed(:)
    end type
contains
    function next_path(unit, directory) result(path)
        integer, intent(in) :: unit
        character(*), intent(in) :: directory
        character(:), allocatable :: path
        character(4096) :: line
        integer :: ios
        line = ''
        read(unit,'(a)',iostat=ios) line
        if (ios /= 0 .or. len_trim(line) == 0) error stop 'missing path in SUPG file list'
        path = trim(adjustl(line))
        if (path(1:1) /= '/') path = directory // path
    end function

    subroutine read_supg_case(degree, with_diffusion, c)
        integer, intent(in) :: degree
        logical, intent(in) :: with_diffusion
        type(supg_case), intent(out) :: c
        character(4096) :: manifest, key, value
        character(:), allocatable :: dir, input, mesh, velocity, boundary, initial
        real(real64) :: val
        integer :: unit, i, j, id, nnode, nelm, nfixed, nlocal, ios, steps
        integer :: record, node, extra
        logical, allocatable :: seen(:)
        real(real64), allocatable :: boundary_values(:)

        manifest = 'file.txt'
        if (mod(command_argument_count(),2) /= 0) error stop 'SUPG options require option/value pairs'
        do i = 1,command_argument_count(),2
            call get_command_argument(i, key)
            if (key == '--config') call get_command_argument(i+1,manifest)
        end do
        j = index(manifest,'/',back=.true.)
        dir = manifest(:j)
        open(newunit=unit,file=trim(manifest),status='old',action='read')
        input = next_path(unit,dir)
        mesh = next_path(unit,dir)
        velocity = next_path(unit,dir)
        boundary = next_path(unit,dir)
        initial = next_path(unit,dir)
        c%result_file = next_path(unit,dir)
        c%vtk_prefix = next_path(unit,dir)
        if (.not.with_diffusion) then
            c%sample_file = next_path(unit,dir)
            c%mass_file = next_path(unit,dir)
        end if
        close(unit)

        open(newunit=unit,file=input,status='old',action='read')
        read(unit,*) c%first_step, c%last_step
        read(unit,*) c%dt
        read(unit,*) c%eps, c%max_iterations
        c%sample_step = 50
        if (degree == 1) c%sample_step = 1000
        if (with_diffusion) then
            read(unit,*) c%diffusion
        else
            ! Legacy files may have an unused second case after a blank line.
            read(unit,'(a)',iostat=ios) value
            if (ios == 0) then
                value = adjustl(value)
                if (index(value,'sample_step=') == 1) then
                    read(value(13:),*,iostat=ios) c%sample_step
                    if (ios /= 0) error stop 'invalid sample_step in SUPG input'
                end if
            end if
        end if
        close(unit)
        if (c%first_step == -1) c%first_step = 1
        do i = 1,command_argument_count(),2
            call get_command_argument(i,key)
            call get_command_argument(i+1,value)
            ios = 0
            select case(trim(key))
            case ('--config')
                cycle
            case ('--steps')
                read(value,*,iostat=ios) steps
                if (ios == 0) c%last_step = c%first_step + steps - 1
            case ('--dt')
                read(value,*,iostat=ios) c%dt
            case ('--eps')
                read(value,*,iostat=ios) c%eps
            case ('--max-iterations')
                read(value,*,iostat=ios) c%max_iterations
            case ('--diffusion')
                if (.not.with_diffusion) error stop 'diffusion requires convection-diffusion solver'
                read(value,*,iostat=ios) c%diffusion
            case ('--sample-step')
                if (with_diffusion) error stop 'cross sections require advection solver'
                read(value,*,iostat=ios) c%sample_step
            case default
                write(*,*) 'Unknown SUPG option: ',trim(key)
                error stop 1
            end select
            if (ios /= 0) error stop 'invalid numeric SUPG option'
        end do
        if (c%sample_step < 0) error stop 'sample_step must be nonnegative'
        if (c%first_step < 1 .or. c%last_step < c%first_step .or. c%max_iterations < 1) error stop 'invalid step limits'
        if (.not.ieee_is_finite(c%dt) .or. c%dt <= 0d0) error stop 'invalid dt'
        if (.not.ieee_is_finite(c%eps) .or. c%eps <= 0d0) error stop 'invalid solver tolerance'
        if (.not.ieee_is_finite(c%diffusion) .or. c%diffusion < 0d0) error stop 'invalid diffusion'
        print *, 'steps', c%first_step, c%last_step, 'dt', c%dt, 'eps', c%eps, 'iteration limit', c%max_iterations
        if (with_diffusion) print *, 'diffusion', c%diffusion
        if (.not.with_diffusion) print *, 'cross section step', c%sample_step

        nlocal = (degree+1)*(degree+2)/2
        open(newunit=unit,file=mesh,status='old',action='read')
        read(unit,*) nnode, nelm
        if (min(nnode,nelm) < 1) error stop 'empty mesh'
        allocate(c%coords(2,nnode), c%nodes(nlocal,nelm), c%initial(nnode), c%u(nnode), c%v(nnode))
        allocate(seen(nnode))
        seen = .false.
        do i = 1,nnode
            block
                real(real64) :: xy(2)
                read(unit,*) id, xy
                call check_node(id,seen)
                c%coords(:,id) = xy
            end block
        end do
        do i = 1,nelm
            read(unit,'(a)') value
            read(value,*,iostat=ios) id,c%nodes(:,i)
            if (ios /= 0) error stop 'wrong number of nodes in triangle record'
            read(value,*,iostat=ios) id,c%nodes(:,i),extra
            if (ios >= 0) error stop 'unexpected data after triangle connectivity'
        end do
        close(unit)

        open(newunit=unit,file=velocity,status='old',action='read')
        seen = .false.
        do i = 1,nnode
            block
                real(real64) :: uv(2)
                read(unit,*) id,uv
                call check_node(id,seen)
                c%u(id) = uv(1)
                c%v(id) = uv(2)
            end block
        end do
        close(unit)
        open(newunit=unit,file=initial,status='old',action='read')
        seen = .false.
        do i = 1,nnode
            read(unit,*) id,val
            call check_node(id,seen)
            c%initial(id) = val
        end do
        close(unit)
        open(newunit=unit,file=boundary,status='old',action='read')
        read(unit,*) nfixed
        if (nfixed < 0 .or. nfixed > nnode) error stop 'invalid boundary count'
        allocate(c%fixed(nfixed), boundary_values(nfixed))
        seen = .false.
        do i = 1,nfixed
            read(unit,*) record, node, boundary_values(i)
            call check_node(node,seen)
            c%fixed(i) = node
        end do
        close(unit)
        c%initial(c%fixed) = boundary_values
    end subroutine

    subroutine check_node(id,seen)
        integer, intent(in) :: id
        logical, intent(inout) :: seen(:)
        if (id < 1 .or. id > size(seen)) error stop 'node ID outside mesh'
        if (seen(id)) error stop 'duplicated node ID'
        seen(id) = .true.
    end subroutine

    subroutine write_solution(unit, c, step, field, velocity_output)
        integer, intent(in) :: unit, step
        type(supg_case), intent(in) :: c
        real(real64), intent(in) :: field(:)
        logical, intent(in) :: velocity_output
        integer :: vtk, i, m, nnode, nelm, nlocal, cell_type
        character(32) :: num
        nnode = size(field)
        nlocal = size(c%nodes,1)
        nelm = size(c%nodes,2)
        ! Existing .res readers expect step number and step size.
        write(unit,*) step, c%dt
        do i = 1,nnode
            write(unit,*) i, field(i)
        end do
        write(num,'(i0.4)') step
        open(newunit=vtk,file=c%vtk_prefix//trim(num)//'.vtk',status='replace',action='write')
        write(vtk,'(a)') '# vtk DataFile Version 3.0','Advection result','ASCII','DATASET UNSTRUCTURED_GRID'
        write(vtk,'(a,i0,a)') 'POINTS ', nnode, ' double'
        do i = 1,nnode
            write(vtk,'(3f16.8)') c%coords(:,i), 0d0
        end do
        write(vtk,'(a,i0,1x,i0)') 'CELLS ', nelm, nelm*(nlocal+1)
        do m = 1,nelm
            write(vtk,'(*(i8))') nlocal, c%nodes(:,m)-1
        end do
        select case(nlocal)
        case(3)
            cell_type = 5
        case(6)
            cell_type = 22
        case(10)
            cell_type = 69
        case default
            error stop 'unsupported triangle type for VTK'
        end select
        write(vtk,'(a,i0)') 'CELL_TYPES ',nelm
        do m = 1,nelm
            write(vtk,'(i8)') cell_type
        end do
        write(vtk,'(a,i0)') 'POINT_DATA ',nnode
        write(vtk,'(a)') 'SCALARS phi double 1','LOOKUP_TABLE default'
        write(vtk,'(f18.9)') field
        if (velocity_output) then
            write(vtk,'(a)') 'VECTORS velocity double'
            do i = 1,nnode
                write(vtk,'(3f18.9)') c%u(i), c%v(i), 0d0
            end do
        end if
        close(vtk)
    end subroutine

    subroutine write_section(unit,coords,field)
        integer, intent(in) :: unit
        real(real64), intent(in) :: coords(:,:), field(:)
        integer, allocatable :: ids(:)
        integer :: i,j,id
        ! Sorted x profile on the historical y=0.5 cut.
        ids = pack([(i,i=1,size(field))], abs(coords(2,:)-0.5d0) < 1d-9)
        do i = 2,size(ids)
            id = ids(i)
            j = i-1
            do while (j > 0)
                if (coords(1,ids(j)) <= coords(1,id)) exit
                ids(j+1) = ids(j)
                j = j-1
            end do
            ids(j+1) = id
        end do
        do i = 1,size(ids)
            write(unit,'(2(1x,g20.9))') coords(1,ids(i)), field(ids(i))
        end do
    end subroutine
end module
