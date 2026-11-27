module supg_fields
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: prepare_fields
contains
    subroutine prepare_fields(default_velocity, default_initial, default_boundary)
        character(*), intent(in) :: default_velocity, default_initial, default_boundary
        character(4096) :: config, key, value, paths(4), line
        character(:), allocatable :: directory, velocity, initial, boundary
        real(real64), allocatable :: xy(:,:), uv(:,:), phi(:)
        logical, allocatable :: on_boundary(:)
        integer :: i, j, file, ios

        config = 'file.dat'
        velocity = default_velocity
        initial = default_initial
        boundary = default_boundary
        if (mod(command_argument_count(),2) /= 0) error stop 'field generator expects option/value pairs'
        do i = 1,command_argument_count(),2
            call get_command_argument(i,key)
            if (key == '--config') call get_command_argument(i+1,config)
        end do
        directory = config(:index(config,'/',back=.true.))
        open(newunit=file, file=trim(config),status='old',action='read')
        do i = 1,size(paths)
            read(file,'(a)') line
            line = adjustl(line)
            if (line == '') error stop 'empty path in field generator configuration'
            if (line(1:1) == '/') then
                paths(i) = line
            else
                paths(i) = directory//trim(line)
            end if
        end do
        ! The first blank line starts the legacy list of alternative cases.
        do
            read(file,'(a)',iostat=ios) line
            if (ios < 0) exit
            if (ios > 0) error stop 'cannot read field generator options'
            if (line == '') exit
            j = index(line,'=')
            if (j < 2) error stop 'expected velocity=, initial=, or boundary= after the four paths'
            call set_option(trim(adjustl(line(:j-1))),trim(adjustl(line(j+1:))))
        end do
        close(file)
        do i = 1,command_argument_count(),2
            call get_command_argument(i,key)
            call get_command_argument(i+1,value)
            if (key == '--config') cycle
            if (key(1:2) /= '--') error stop 'expected long field generator option'
            call set_option(trim(key(3:)),trim(value))
        end do

        call read_coordinates(paths(1),xy)
        allocate(uv(2,size(xy,2)),phi(size(xy,2)),on_boundary(size(xy,2)))
        call evaluate_fields(xy,velocity,initial,boundary,uv,phi,on_boundary)
        call write_fields(paths(2:),uv,phi,on_boundary)
        print *, 'velocity=',velocity,' initial=',initial,' boundary=',boundary

    contains
        subroutine set_option(name,setting)
            character(*), intent(in) :: name, setting
            select case(name)
            case('velocity')
                velocity = setting
            case('initial')
                initial = setting
            case('boundary')
                boundary = setting
            case default
                error stop 'unknown SUPG field generator option'
            end select
        end subroutine
    end subroutine

    subroutine read_coordinates(filename,xy)
        character(*), intent(in) :: filename
        real(real64), allocatable, intent(out) :: xy(:,:)
        logical, allocatable :: seen(:)
        real(real64) :: point(2)
        integer :: file, nnode, nelm, i, id
        open(newunit=file,file=trim(filename),status='old',action='read')
        read(file,*) nnode, nelm
        if (min(nnode,nelm) < 1) error stop 'empty mesh'
        allocate(xy(2,nnode),seen(nnode))
        seen = .false.
        do i = 1,nnode
            read(file,*) id,point
            if (id < 1 .or. id > nnode) error stop 'mesh node ID outside range'
            if (seen(id)) error stop 'duplicated mesh node'
            seen(id) = .true.
            xy(:,id) = point
        end do
        close(file)
    end subroutine

    subroutine evaluate_fields(xy,velocity,initial,boundary,uv,phi,on_boundary)
        real(real64), intent(in) :: xy(:,:)
        character(*), intent(in) :: velocity,initial,boundary
        real(real64), intent(out) :: uv(:,:),phi(:)
        logical, intent(out) :: on_boundary(:)
        real(real64), parameter :: pi = acos(-1d0), tol = 1d-10
        real(real64) :: radius(size(phi))
        select case(velocity)
        case('down')
            uv(1,:) = 0d0
            uv(2,:) = -0.5d0
        case('right')
            uv(1,:) = 0.5d0
            uv(2,:) = 0d0
        case('rotation')
            uv(1,:) = -xy(2,:)
            uv(2,:) = xy(1,:)
        case default
            error stop 'velocity must be down, right or rotation'
        end select
        select case(initial)
        case('box')
            phi = merge(1d0,0d0,abs(xy(1,:)) <= 0.25d0 .and. xy(2,:) >= 0.25d0 .and. xy(2,:) <= 0.75d0)
        case('cosine')
            radius = hypot(xy(1,:),xy(2,:)-0.5d0)
            phi = merge(0.5d0*(1d0+cos(2d0*pi*radius)),0d0,radius <= 0.5d0)
        case('zero')
            phi = 0d0
        case default
            error stop 'initial must be box, cosine or zero'
        end select
        on_boundary = abs(xy(1,:)-minval(xy(1,:))) <= tol .or. abs(xy(2,:)-minval(xy(2,:))) <= tol &
                      .or. abs(xy(2,:)-maxval(xy(2,:))) <= tol
        select case(boundary)
        case('all')
            on_boundary = on_boundary .or. abs(xy(1,:)-maxval(xy(1,:))) <= tol
        case('open_right')
            ! The original convection-diffusion example leaves the right edge open.
        case default
            error stop 'boundary must be all or open_right'
        end select
    end subroutine

    subroutine write_fields(paths,uv,phi,on_boundary)
        character(*), intent(in) :: paths(:)
        real(real64), intent(in) :: uv(:,:),phi(:)
        logical, intent(in) :: on_boundary(:)
        integer :: vel,initial,boundary,i,j
        open(newunit=vel,file=trim(paths(1)),status='replace',action='write')
        open(newunit=initial,file=trim(paths(2)),status='replace',action='write')
        open(newunit=boundary,file=trim(paths(3)),status='replace',action='write')
        write(boundary,'(i9)') count(on_boundary)
        j = 0
        do i = 1,size(phi)
            write(vel,'(i9,2e15.6)') i,uv(:,i)
            write(initial,'(i9,e15.6)') i,phi(i)
            if (on_boundary(i)) then
                j = j+1
                write(boundary,'(2i9,e15.6)') j,i,0d0
            end if
        end do
        close(vel)
        close(initial)
        close(boundary)
    end subroutine
end module
