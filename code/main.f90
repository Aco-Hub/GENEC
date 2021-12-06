!>   STELLAR EVOLUTION PROGRAM OF THE GENEVA GROUP
!!
!!  @author A. Maeder, G. Meynet, D. Schaerer, R. Hirschi, S. Ekstrom, C. Georgy
!!  @version  283
!!  @date     mars 2013
!!
!!  @brief Kippenhahn program modified for the effects of rotation, advanced phases, ...
! --------------------------------------------------------------------------
program main

    use genec

    implicit none

    call initialise_genec()
    call read_parameters()
    call initialise_star()
    call evolve()
    call finalise()
end program main
