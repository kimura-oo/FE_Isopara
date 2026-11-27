program poisson_p1
    use poisson_solver, only: solve_poisson
    implicit none
    call solve_poisson(1, 'sine')
end program
