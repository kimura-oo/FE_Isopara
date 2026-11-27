program advection_p2
    use supg_advection, only: run_advection
    implicit none
    call run_advection(2)
end program
