module mesh_elevation
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: convert_mesh
contains
    subroutine convert_mesh(order)
        integer, intent(in) :: order
        character(1024) :: infile, outfile
        integer :: io, unit, nvertex, nelem, nnode, local_size
        integer :: cell, side, i, j, id, a, b, low, high, edge, nedge, first, slot
        integer, allocatable :: conn(:,:), head(:), endpoint(:), next(:), first_node(:)
        real(real64), allocatable :: xy(:,:)

        select case (command_argument_count())
        case (0)
            open(newunit=unit, file='file.txt', status='old', action='read', iostat=io)
            if (io == 0) then
                read(unit,'(a)') infile
                read(unit,'(a)') outfile
                close(unit)
            else if (order == 3) then
                infile = 'input.txt'
                outfile = 'output.p3.txt'
            else
                error stop 'Expected file.txt with input and output filenames'
            end if
        case (2)
            call get_command_argument(1, infile)
            call get_command_argument(2, outfile)
        case default
            error stop 'Usage: mesh converter [input.mes output.mes]'
        end select
        open(newunit=unit, file=infile, status='old', action='read')
        read(unit,*) nvertex, nelem
        if (nvertex < 3 .or. nelem < 1) error stop 'Empty triangle mesh'
        local_size = (order+1)*(order+2)/2
        allocate(xy(2,nvertex+(local_size-3)*nelem), conn(local_size,nelem))
        do i = 1, nvertex
            read(unit,*) id, xy(:,i)
            if (id /= i) error stop 'Expected consecutive vertex IDs in input order'
        end do
        do cell = 1, nelem
            read(unit,*) id, conn(1:3,cell)
            if (any(conn(1:3,cell)<1) .or. any(conn(1:3,cell)>nvertex)) error stop 'Invalid vertex ID'
        end do
        close(unit)

        ! Index each unique edge by its smaller vertex. Number nodes in first-visit
        ! order to preserve IDs in the meshes and boundary files from earlier runs.
        allocate(head(nvertex), endpoint(3*nelem), next(3*nelem), first_node(3*nelem))
        head = 0
        nedge = 0
        nnode = nvertex
        do cell = 1, nelem
            do side = 1, 3
                a = conn(side,cell)
                b = conn(mod(side,3)+1,cell)
                if (a == b) error stop 'Repeated triangle vertex'
                low = min(a,b)
                high = max(a,b)
                edge = head(low)
                do while (edge /= 0)
                    if (endpoint(edge) == high) exit
                    edge = next(edge)
                end do
                if (edge == 0) then
                    nedge = nedge + 1
                    edge = nedge
                    next(edge) = head(low)
                    head(low) = edge
                    endpoint(edge) = high
                    first_node(edge) = nnode+1
                    do j = 1, order-1
                        nnode = nnode+1
                        xy(:,nnode) = ((order-j)*xy(:,low) + j*xy(:,high))/order
                    end do
                end if
                first = first_node(edge)
                slot = 3+(side-1)*(order-1)
                do j = 1, order-1
                    if (a == low) then
                        conn(slot+j,cell) = first+j-1
                    else
                        conn(slot+j,cell) = first+order-j-1
                    end if
                end do
            end do
        end do
        if (order == 3) then
            do cell = 1, nelem
                nnode = nnode+1
                xy(:,nnode) = sum(xy(:,conn(1:3,cell)),dim=2)/3
                conn(10,cell) = nnode
            end do
        end if

        open(newunit=unit, file=outfile, status='replace', action='write')
        write(unit,'(2i10)') nnode, nelem
        do i = 1, nnode
            write(unit,'(i8,2f14.6)') i, xy(:,i)
        end do
        do cell = 1, nelem
            write(unit,'(*(i8))') cell, conn(:,cell)
        end do
        close(unit)
    end subroutine
end module
