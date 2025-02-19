program main

    use genec
    use WriteSaveClose, only: write_compilation_informations

    implicit none

    call initialise_genec()
    call write_compilation_informations()
    call read_parameters()
    call initialise_star()
    call evolve()
    call finalise()

end program main
