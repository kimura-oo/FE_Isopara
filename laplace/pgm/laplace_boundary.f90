module laplace_boundary
    use iso_fortran_env, only: iostat_end
    use ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: read_edge_conditions
contains
    subroutine read_edge_conditions(unit, edge_size, node_count, potential_nodes, potential_values, flux_nodes, flux_values)
        integer, intent(in) :: unit, edge_size, node_count
        integer, allocatable, intent(out) :: potential_nodes(:,:), flux_nodes(:,:)
        double precision, allocatable, intent(out) :: potential_values(:), flux_values(:)
        integer, allocatable :: nodes(:,:), combined(:,:)
        double precision, allocatable :: values(:)
        character(1024) :: line
        character(32) :: label
        integer :: status, count, row, id, i, k, nflux
        allocate(flux_nodes(edge_size,0), flux_values(0))
        do
            call next_line(unit,line,status)
            if (status == iostat_end) exit
            read(line,*,iostat=status) count,label
            if (status /= 0) error stop 'malformed Laplace boundary header'
            if (count < 0) error stop 'negative Laplace boundary count'
            allocate(nodes(edge_size,count),values(count))
            do row = 1, count
                call next_line(unit,line,status)
                if (status == iostat_end) error stop 'truncated Laplace boundary section'
                read(line,*,iostat=status) id,nodes(:,row),values(row)
                if (status /= 0) error stop 'malformed Laplace boundary edge'
                if (any(nodes(:,row)<1) .or. any(nodes(:,row)>node_count)) error stop 'Laplace boundary node out of range'
                if (.not. ieee_is_finite(values(row))) error stop 'non-finite Laplace boundary value'
            end do
            do i = 1, len_trim(label)
                k = iachar(label(i:i))
                if (k >= 97 .and. k <= 122) label(i:i) = achar(k-32)
            end do
            select case (label)
            case ('POTENTIAL')
                if (allocated(potential_nodes)) error stop 'multiple POTENTIAL sections'
                call move_alloc(nodes,potential_nodes)
                call move_alloc(values,potential_values)
            case ('VELOCITY_IN','VELOCITY_INLET','VELOCITY_WALL','NEUMANN')
                nflux = size(flux_values)
                allocate(combined(edge_size,nflux+count))
                combined(:,1:nflux) = flux_nodes
                combined(:,nflux+1:) = nodes
                call move_alloc(combined,flux_nodes)
                flux_values = [flux_values,values]
                deallocate(nodes,values)
            case default
                error stop 'unknown Laplace boundary section'
            end select
        end do
        if (.not. allocated(potential_nodes)) error stop 'missing Laplace POTENTIAL section'
        if (size(potential_values) == 0) error stop 'Laplace needs a prescribed potential'
    end subroutine

    subroutine next_line(unit,line,status)
        integer, intent(in) :: unit
        character(*), intent(out) :: line
        integer, intent(out) :: status
        do
            read(unit,'(a)',iostat=status) line
            if (status == iostat_end) return
            if (status /= 0) error stop 'failed to read Laplace boundary data'
            if (len_trim(line)>0) return
        end do
    end subroutine
end module
