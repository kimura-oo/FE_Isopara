program advection_p1
    use supg_advection, only: run_advection
    implicit none
    call run_advection(1)
end program
