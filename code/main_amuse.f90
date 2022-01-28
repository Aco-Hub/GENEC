! This is a version of GENEC that mimics the AMUSE interface
! It should run exactly as a "default" run with AMUSE
! - Steven Rieder
program main

    use genec, only: initialise_genec, writetofiles, read_parameters, modell
    use amuse_helpers, only: set_defaults
    use inputparam, only: nzmod, iprn, nwseq, modanf, vwant, amuseinterface, phase
    use caramodele, only: nwmd

    use amuse_helpers, only: mstar, zini, starname

    !use amuse_helpers, only: initialise_star, makeini
    use amuse_helpers, only: makeini
    use genec, only: initialise_star
    use WriteSaveClose, only: OpenAll,CloseAll,quitafterclosing
    use timestep, only: TimestepControle

    use genec, only: evolve
    use const, only: year
    use timestep, only: alter
    use amuse_helpers, only: nfseq

    !use Sequence, only: SequenceFinish
    use genec, only: finalise
    implicit none

    double precision:: alter_max

    ! First, initialize_code from interface
    amuseinterface = .true.
    writetofiles = .true.
    quitafterclosing = .false.

    call initialise_genec()
    call set_defaults()
    !call read_parameters()
    write(*,*) ' Enter the star name:'
    !read(*,*) starname
    starname = "AmuseDefaultStar"
    write(*,*) 'Enter the desired mass and metallicity:'
    !read(*,*) mstar, zini
    mstar = 7.0
    zini = 0.014d0
    write(*,*) ' Which rotation velocity on the ZAMS?'
    !read(*,*) vwant
    vwant = 0.0

    ! Second, commit_parameters

    ! Third, commit_particles
    call makeini()
    call OpenAll()
    call initialise_star()

    !! Fourth, evolve_for some time (default here: huge number)
    nzmod = 10
    iprn = 10
    nwseq = 1
    modanf = 0

    ! Optionally, set this to end after a specific time (in years)
    alter_max = huge(1.0d0)

    alter = 0.d0

    do while (alter < alter_max)
        write(*,*) "Current time: ", alter, " current model: ", nwseq, " current nzmod: ", nzmod
        write(*,*) "Current phase: ",phase
        write(*,*) "Final model: ", nfseq, " nwmd: ", nwmd
        modell = 1
        call evolve()
        !write(*,*) "*evolve done*************************************************"
        !if (nwseq >= n)
        call finalise()
        !write(*,*) "*finalise done***********************************************"
        !call CloseAll()
        call OpenAll()
        !write(*,*) "*OpenAll done************************************************"
        !call read_parameters()
        !write(*,*) "*read_parameters done*********************************************"
        call initialise_star()
        !write(*,*) "*initialise_star done****************************************"
        !call SequenceFinish()
        !Write(*,*) 'nwseq = ', nwseq, ' modanf = ', modanf
        !Modanf = modanf + 1
        !Nwseq = nwseq + nzmod
        !call INPUTS_Change(xm,ym,c12m,ne20m,O16m,rapom2,m,nzmodini,nzmodnew)
        !call TimestepControle(xcprev,xclast,xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,nzmodini,xcnwant)
    end do

    ! Last, stop
    quitafterclosing = .true.
    call finalise()
    !call SequenceFinish()

end program main
