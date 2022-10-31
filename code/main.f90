program main

    use genec
    use inputparam, only: nwseq

    implicit none

    call initialise_genec()
    call read_parameters()
    call initialise_star()
    call evolve()
    call finalise()
end program main
