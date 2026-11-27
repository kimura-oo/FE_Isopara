program poisson_p3
    use poisson_solver, only: solve_poisson
    implicit none
    call solve_poisson(3, 'sine')
end program
