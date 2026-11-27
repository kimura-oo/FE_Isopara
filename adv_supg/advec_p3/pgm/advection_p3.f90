program advection_p3
    use supg_advection, only: run_advection
    implicit none
    call run_advection(3)
end program
