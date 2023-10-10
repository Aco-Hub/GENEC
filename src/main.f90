program main

    use genec

    implicit none

    call initialise_genec()
    call read_parameters()
    call initialise_star()
    call evolve()
    call finalise()
end program main
