module poisson_solver
    use iso_fortran_env, only: real64, iostat_end, error_unit
    use ieee_arithmetic, only: ieee_is_finite
    use triangle_elements, only: triangle_shape, edge_nodes
    use triangle_quadrature, only: triangle_rule
    implicit none
    private
    public :: solve_poisson

    type mesh
        real(real64), allocatable :: xy(:,:)
        integer, allocatable :: cells(:,:)
    end type
contains
    subroutine solve_poisson(degree, default_source)
        integer, intent(in) :: degree
        character(*), intent(in) :: default_source
        type(mesh) :: grid
        character(1024) :: paths(6)
        character(32) :: source
        integer :: nlocal, nelem, nnode, nquad, bandwidth, element, node, q
        real(real64), allocatable :: basis(:,:), grad(:,:,:,:), weights(:,:), velocity(:,:)
        real(real64), allocatable :: matrix(:,:), rhs(:), prescribed(:)
        logical, allocatable :: fixed(:)

        call read_config(default_source, paths, source)
        call read_mesh(paths(1), degree, grid)
        nlocal = size(grid%cells,1)
        nelem = size(grid%cells,2)
        nnode = size(grid%xy,2)
        nquad = merge(7,3,degree == 3)
        allocate(basis(nlocal,nquad),grad(2,nlocal,nquad,nelem),weights(nquad,nelem))
        call prepare_geometry(degree,grid,basis,grad,weights)

        ! Symmetric upper band. matrix(j-i,i) stores K(i,j) for i <= j.
        bandwidth = maxval(maxval(grid%cells,dim=1) - minval(grid%cells,dim=1))
        allocate(matrix(0:bandwidth,nnode), rhs(nnode), prescribed(nnode), fixed(nnode))
        call assemble(grid,source,basis,grad,weights,matrix,rhs)
        call read_boundary(paths(2),degree,grid,rhs,prescribed,fixed)
        if (.not. any(fixed)) error stop 'Poisson needs a Dirichlet value to fix the potential'
        do node = 1, nnode
            if (fixed(node)) call apply_dirichlet(matrix,rhs,node,prescribed(node))
        end do
        call solve_band(matrix,rhs)

        allocate(velocity(2,nelem))
        velocity = 0
        do element = 1, nelem
            do q = 1, nquad
                velocity(:,element) = velocity(:,element) - weights(q,element)* &
                    matmul(grad(:,:,q,element),rhs(grid%cells(:,element)))
            end do
            velocity(:,element) = velocity(:,element) / sum(weights(:,element))
        end do
        call write_results(paths(3:6),degree,grid,rhs,velocity)
        print '(a,i0,2a)', 'Poisson P',degree,' source=',trim(source)
    end subroutine

    subroutine read_config(default_source, paths, source)
        character(*), intent(in) :: default_source
        character(*), intent(out) :: paths(:), source
        character(1024) :: arg, file_list
        integer :: unit, i, status
        logical :: exists
        namelist /poisson/ source

        source = default_source
        inquire(file='poisson.nml',exist=exists)
        if (exists) then
            open(newunit=unit,file='poisson.nml',status='old',action='read')
            read(unit,nml=poisson,iostat=status)
            if (status /= 0) error stop 'invalid poisson.nml'
            close(unit)
        end if
        file_list = 'file_FEM.dat'
        i = 1
        do while (i <= command_argument_count())
            call get_command_argument(i,arg)
            if (i == command_argument_count()) error stop 'missing Poisson option value'
            select case (arg)
            case ('--source')
                call get_command_argument(i+1,source,status=status)
            case ('--files')
                call get_command_argument(i+1,file_list,status=status)
            case default
                error stop 'expected --source NAME or --files PATH'
            end select
            if (status /= 0) error stop 'invalid Poisson option value'
            i = i+2
        end do
        select case (source)
        case ('zero','sine','quadratic')
        case default
            error stop 'source must be zero, sine or quadratic'
        end select
        open(newunit=unit,file=file_list,status='old',action='read')
        do i = 1, size(paths)
            read(unit,'(a)',iostat=status) paths(i)
            if (status /= 0 .or. len_trim(paths(i)) == 0) error stop 'expected six paths in Poisson file list'
        end do
        close(unit)
    end subroutine

    subroutine read_mesh(path, degree, grid)
        character(*), intent(in) :: path
        integer, intent(in) :: degree
        type(mesh), intent(out) :: grid
        integer :: unit, nnode, nelem, row, id, status, k, conn((degree+1)*(degree+2)/2)
        real(real64) :: xy(2)
        character(2048) :: line
        logical, allocatable :: seen(:)
        open(newunit=unit,file=path,status='old',action='read')
        read(unit,*) nnode, nelem
        if (nnode < size(conn) .or. nelem < 1) error stop 'invalid mesh dimensions'
        allocate(grid%xy(2,nnode),grid%cells(size(conn),nelem),seen(max(nnode,nelem)))
        seen = .false.
        do row = 1, nnode
            read(unit,'(a)') line
            read(line,*,iostat=status) id, xy
            if (status /= 0) error stop 'malformed node record'
            if (id < 1 .or. id > nnode) error stop 'node ID out of range'
            if (seen(id) .or. .not. all(ieee_is_finite(xy))) error stop 'duplicate or invalid mesh node'
            grid%xy(:,id) = xy
            seen(id) = .true.
        end do
        seen = .false.
        do row = 1, nelem
            read(unit,'(a)') line
            read(line,*,iostat=status) id, conn
            if (status /= 0) error stop 'malformed element connectivity'
            ! Catch a higher-order mesh used by a lower-order entry point.
            read(line,*,iostat=status) id, conn, k
            if (status /= iostat_end) error stop 'element record does not match degree'
            if (id < 1 .or. id > nelem) error stop 'element ID out of range'
            if (seen(id) .or. minval(conn) < 1 .or. maxval(conn) > nnode) error stop 'invalid element connectivity'
            do k = 2, size(conn)
                if (any(conn(:k-1) == conn(k))) error stop 'duplicate node in element'
            end do
            grid%cells(:,id) = conn
            seen(id) = .true.
        end do
        close(unit)
    end subroutine

    subroutine prepare_geometry(degree, grid, basis, grad, weights)
        integer, intent(in) :: degree
        type(mesh), intent(in) :: grid
        real(real64), intent(out) :: basis(:,:), grad(:,:,:,:), weights(:,:)
        real(real64) :: xi(size(basis,2)), eta(size(basis,2)), w(size(basis,2))
        real(real64) :: deriv(2,size(basis,1)), jac(2,2), det
        integer :: q, e
        call triangle_rule(size(w),xi,eta,w)
        do q = 1, size(w)
            call triangle_shape(degree,xi(q),eta(q),basis(:,q),deriv)
            do e = 1, size(grid%cells,2)
                jac = matmul(grid%xy(:,grid%cells(:,e)),transpose(deriv))
                det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
                if (.not. ieee_is_finite(det) .or. det <= 0.0_real64) then
                    write(error_unit,*) 'Invalid Jacobian in element',e,' at quadrature point',q
                    error stop 1
                end if
                weights(q,e) = det*w(q)
                grad(1,:,q,e) = (jac(2,2)*deriv(1,:) - jac(2,1)*deriv(2,:))/det
                grad(2,:,q,e) = (jac(1,1)*deriv(2,:) - jac(1,2)*deriv(1,:))/det
            end do
        end do
    end subroutine

    pure real(real64) function source_at(name, xy) result(value)
        character(*), intent(in) :: name
        real(real64), intent(in) :: xy(2)
        real(real64), parameter :: pi = acos(-1.0_real64)
        select case (name)
        case ('zero')
            value = 0
        case ('quadratic')
            value = -4
        case ('sine')
            value = 2*pi**2*sin(pi*xy(1))*sin(pi*xy(2))
        case default
            error stop 'unknown source'
        end select
    end function

    subroutine assemble(grid, source, basis, grad, weights, matrix, rhs)
        type(mesh), intent(in) :: grid
        character(*), intent(in) :: source
        real(real64), intent(in) :: basis(:,:), grad(:,:,:,:), weights(:,:)
        real(real64), intent(out) :: matrix(0:,:), rhs(:)
        integer :: e, q, i, j, ni, nj
        real(real64) :: force
        matrix = 0
        rhs = 0
        do e = 1, size(grid%cells,2)
            do q = 1, size(basis,2)
                force = source_at(source,matmul(grid%xy(:,grid%cells(:,e)),basis(:,q)))
                do i = 1, size(basis,1)
                    ni = grid%cells(i,e)
                    rhs(ni) = rhs(ni) + weights(q,e)*basis(i,q)*force
                    do j = 1, size(basis,1)
                        nj = grid%cells(j,e)
                        if (nj < ni) cycle
                        matrix(nj-ni,ni) = matrix(nj-ni,ni) + weights(q,e)*dot_product(grad(:,i,q,e),grad(:,j,q,e))
                    end do
                end do
            end do
        end do
    end subroutine

    subroutine read_boundary(path, degree, grid, rhs, prescribed, fixed)
        character(*), intent(in) :: path
        integer, intent(in) :: degree
        type(mesh), intent(in) :: grid
        real(real64), intent(inout) :: rhs(:)
        real(real64), intent(out) :: prescribed(:)
        logical, intent(out) :: fixed(:)
        integer :: unit, status, count, row, id, nodes(degree+1), nids, i, n, section
        character(1) :: extra
        character(2048) :: line
        character(40) :: label
        real(real64) :: value
        logical :: is_flux

        open(newunit=unit,file=path,status='old',action='read')
        fixed = .false.
        prescribed = 0
        section = 0
        do
            call boundary_record(unit,line,status)
            if (status == iostat_end) exit
            section = section+1
            read(line,*,iostat=status) count
            if (status /= 0) error stop 'malformed boundary section header'
            if (count < 0) error stop 'negative boundary count'
            label = ''
            read(line,*,iostat=status) count,label
            if (status /= iostat_end .and. status /= 0) error stop 'malformed boundary label'
            do i = 1, len_trim(label)
                n = iachar(label(i:i))
                if (n >= 97 .and. n <= 122) label(i:i) = achar(n-32)
            end do
            is_flux = .false.
            select case (label)
            case ('')
                is_flux = section > 1
            case ('DIRICHLET','POTENTIAL')
            case ('NEUMANN','VELOCITY_IN','VELOCITY_INLET','VELOCITY_WALL')
                is_flux = .true.
            case default
                error stop 'unknown Poisson boundary section'
            end select
            nids = 1
            if (is_flux .or. label == 'POTENTIAL') nids = degree+1
            do row = 1, count
                call boundary_record(unit,line,status)
                if (status /= 0) error stop 'truncated boundary section'
                read(line,*,iostat=status) id,nodes(1:nids),value
                if (status /= 0) error stop 'malformed boundary record'
                read(line,*,iostat=status) id,nodes(1:nids),value,extra
                if (status /= iostat_end) error stop 'boundary record does not match element degree or section type'
                if (any(nodes(1:nids) < 1) .or. any(nodes(1:nids) > size(rhs))) error stop 'boundary node out of range'
                if (.not. ieee_is_finite(value)) error stop 'non-finite boundary value'
                if (is_flux) then
                    call add_edge_flux(degree,grid,nodes,value,rhs)
                else
                    do i = 1, nids
                        n = nodes(i)
                        if (fixed(n) .and. prescribed(n) /= value) error stop 'conflicting prescribed potentials'
                        prescribed(n) = value
                        fixed(n) = .true.
                    end do
                end if
            end do
        end do
        close(unit)
    end subroutine

    subroutine boundary_record(unit, line, status)
        integer, intent(in) :: unit
        character(*), intent(out) :: line
        integer, intent(out) :: status
        integer :: comment
        do
            read(unit,'(a)',iostat=status) line
            if (status == iostat_end) return
            if (status /= 0) error stop 'could not read boundary record'
            comment = scan(line,'!#')
            if (comment > 0) line(comment:) = ''
            if (len_trim(line) > 0) return
        end do
    end subroutine

    subroutine add_edge_flux(degree, grid, nodes, flux, rhs)
        integer, intent(in) :: degree, nodes(:)
        type(mesh), intent(in) :: grid
        real(real64), intent(in) :: flux
        real(real64), intent(inout) :: rhs(:)
        ! Three-point Gauss on [0,1]. The trace and tangent use exactly
        ! the same isoparametric shape functions as the volume integral.
        real(real64), parameter :: offset = sqrt(0.6_real64)/2
        real(real64), parameter :: point(3) = [0.5_real64-offset,0.5_real64,0.5_real64+offset]
        real(real64), parameter :: weight(3) = [5.0_real64,8.0_real64,5.0_real64]/18
        integer :: e, side, matches, q, ids(degree+1), edge(degree+1)
        real(real64) :: values(size(grid%cells,1)), deriv(2,size(grid%cells,1)), ds
        matches = 0
        do e = 1, size(grid%cells,2)
            do side = 1, 3
                edge = grid%cells(edge_nodes(degree,side),e)
                if (all(edge(1:2) == nodes(2:1:-1))) then
                    edge(1:2) = edge(2:1:-1)
                    edge(3:) = edge(size(edge):3:-1)
                end if
                if (all(edge == nodes)) matches = matches+1
            end do
        end do
        if (matches /= 1) error stop 'Neumann nodes must describe one exterior edge with local node ordering'
        ids = edge_nodes(degree,1)
        do q = 1, size(point)
            call triangle_shape(degree,point(q),0.0_real64,values,deriv)
            ds = norm2(matmul(grid%xy(:,nodes),deriv(1,ids)))
            rhs(nodes) = rhs(nodes) + flux*weight(q)*ds*values(ids)
        end do
    end subroutine

    subroutine apply_dirichlet(matrix, rhs, node, value)
        real(real64), intent(inout) :: matrix(0:,:), rhs(:)
        integer, intent(in) :: node
        real(real64), intent(in) :: value
        integer :: i, width
        width = ubound(matrix,1)
        do i = max(1,node-width), node-1
            rhs(i) = rhs(i) - matrix(node-i,i)*value
            matrix(node-i,i) = 0
        end do
        do i = node+1, min(size(rhs),node+width)
            rhs(i) = rhs(i) - matrix(i-node,node)*value
            matrix(i-node,node) = 0
        end do
        matrix(0,node) = 1
        rhs(node) = value
    end subroutine

    subroutine solve_band(matrix,rhs)
        real(real64), intent(inout) :: matrix(0:,:), rhs(:)
        integer :: i, j, last, width
        real(real64) :: factor, threshold
        width = ubound(matrix,1)
        threshold = epsilon(factor)*maxval(matrix(0,:))
        do i = 1, size(rhs)
            last = min(size(rhs),i+width)
            if (matrix(0,i) <= threshold) error stop 'singular Poisson stiffness matrix'
            do j = i+1, last
                if (matrix(j-i,i) == 0) cycle
                factor = matrix(j-i,i)/matrix(0,i)
                matrix(0:last-j,j) = matrix(0:last-j,j) - factor*matrix(j-i:last-i,i)
                rhs(j) = rhs(j) - factor*rhs(i)
            end do
        end do
        do i = size(rhs), 1, -1
            last = min(size(rhs),i+width)
            rhs(i) = (rhs(i) - dot_product(matrix(1:last-i,i),rhs(i+1:last)))/matrix(0,i)
        end do
        if (.not. all(ieee_is_finite(rhs))) error stop 'non-finite Poisson solution'
    end subroutine

    subroutine write_results(paths, degree, grid, phi, velocity)
        character(*), intent(in) :: paths(:)
        integer, intent(in) :: degree
        type(mesh), intent(in) :: grid
        real(real64), intent(in) :: phi(:),velocity(:,:)
        integer, parameter :: vtk_cell(3) = [5,22,69]
        integer :: unit, i, nnode, nelem, nlocal, axis
        real(real64) :: tolerance
        nnode = size(phi)
        nelem = size(grid%cells,2)
        nlocal = size(grid%cells,1)
        open(newunit=unit,file=paths(1),status='replace',action='write')
        do i = 1, nnode
            if (degree == 2) then
                write(unit,'(i7,d15.6)') i,phi(i)
            else
                write(unit,'(i7,3d15.6)') i,grid%xy(:,i),phi(i)
            end if
        end do
        do i = 1, nelem
            write(unit,'(i7,2d15.6)') i,velocity(:,i)
        end do
        close(unit)

        open(newunit=unit,file=paths(2),status='replace',action='write')
        write(unit,'(a)') '# vtk DataFile Version 3.0','2D FEM Poisson solution','ASCII','DATASET UNSTRUCTURED_GRID'
        write(unit,'(a,i0,a)') 'POINTS ',nnode,' double'
        do i = 1, nnode
            write(unit,'(3(1x,es24.16))') grid%xy(:,i),0.0_real64
        end do
        write(unit,'(a,i0,1x,i0)') 'CELLS ',nelem,nelem*(nlocal+1)
        do i = 1, nelem
            write(unit,'(*(i0,1x))') nlocal,grid%cells(:,i)-1
        end do
        write(unit,'(a,i0)') 'CELL_TYPES ',nelem
        do i = 1, nelem
            write(unit,'(i0)') vtk_cell(degree)
        end do
        write(unit,'(a,i0)') 'POINT_DATA ',nnode
        if (degree == 1) then
            write(unit,'(a)') 'SCALARS Potential double 1'
        else
            write(unit,'(a)') 'SCALARS potential double 1'
        end if
        write(unit,'(a)') 'LOOKUP_TABLE default'
        write(unit,'(es24.16)') phi
        write(unit,'(a,i0)') 'CELL_DATA ',nelem
        if (degree == 1) then
            write(unit,'(a)') 'VECTORS Velocity double'
        else
            write(unit,'(a)') 'VECTORS velocity double'
        end if
        do i = 1, nelem
            write(unit,'(3(1x,es24.16))') velocity(:,i),0.0_real64
        end do
        close(unit)

        tolerance = merge(0.1_real64,0.001_real64,degree == 3)
        do axis = 1, 2
            open(newunit=unit,file=paths(2+axis),status='replace',action='write')
            do i = 1, nnode
                if (abs(grid%xy(axis,i)-0.5_real64) < tolerance) write(unit,'(i6,2f12.6)') i,grid%xy(3-axis,i),phi(i)
            end do
            close(unit)
        end do
    end subroutine
end module
