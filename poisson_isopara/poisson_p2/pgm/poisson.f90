program poisson_p2
    use poisson_solver, only: solve_poisson
    implicit none
    call solve_poisson(2, 'quadratic')
end program
