#ifndef GIT_COMMIT
#define GIT_COMMIT "GIT_COMMIT NOT DEFINED"
#endif

#ifndef COMPILATION_DATE
#define COMPILATION_DATE "COMPILATION_DATE NOT DEFINED"
#endif

program main

    use genec

    implicit none

    call initialise_genec()
    call read_parameters()
    call initialise_star()
    call evolve()
    call finalise()
end program main
