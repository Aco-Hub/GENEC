! This is a version of GENEC that mimics the AMUSE interface
! It should run exactly as a "default" run with AMUSE
! - Steven Rieder
program main

    use genec
    use helpers, only: set_defaults,makeini,mstar,zini
    use evol, only: input_dir
    use inputparam, only: libgenec
    use inputparam,only: nzmod,end_at_phase,end_at_model
    use timestep, only: alter
    implicit none

    libgenec = .true.
    input_dir = "./src/GENEC/code"
    call initialise_genec()
    call set_defaults()
    starname = "AmuseDefaultStar"
    mstar = 7.0
    zini = 0.014

    ! Second, commit_parameters

    ! Third, commit_particles
    call makeini()
    call initialise_star()

    ! Fourth, evolve_for some time (default here: huge number)
    nzmod = 1
    n_snap = 0
    end_at_phase=4
    end_at_model=0

    alter = 0.d0

    do while (alter < huge(1.0d0))
        write(*,*) "Current time: ", alter, " current model: ", nwseq, " current nzmod: ", nzmod
        call evolve()
        call finalise()
    end do

end program main
