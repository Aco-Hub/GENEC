module storage
    use evol, only: kindreal,ldi,npondcouche,mbelx

    implicit none

    type genec_star
        ! genec_star type contains all variables that are saved.
        ! It should store everything needed to restore a GENEC simulation to its current state,
        ! without any changes to the result.
        ! All values are initialised to their default value and should be changed afterwards.
        ! Current variable names are taken from GENEC, but should slowly be updated to names
        ! that are more clear.
        ! Values are copied to/from GENEC by means of the functions in the 'helpers' module.

        ! **** AMUSE specific
        integer :: &
                modell
        logical :: &
                initialised=.false.,&  ! .true. if the star has been synchronised to at any time
                synchronised=.false.  ! .true. if the star is up-to-date, set to .false. when the model changes
        integer :: &
                index_of_the_star  ! AMUSE-specific, do not change here!
        real(kindreal) :: &
                radius
        logical :: &
                stopped=.false.  ! a STOP will set this to .true. and then return
        character(256) :: &
                stop_message=""

        ! Variables that would go in the .input file
        ! **** Model characteristics
        character(256) :: &
                star_name=""
        integer :: &
                nwmd,&
                nwseq=1,&
                modanf=0,&
                nzmod,&
                end_at_phase=4,&
                end_at_model=0

        ! **** Physical inputs
        logical :: &
                var_rates=.false.,&
                bintide=.false.,&
                const_per=.true.
        integer :: &
                irot,&
                isol,&
                imagn=0,&
                ialflu,&
                ianiso=0,&
                ipop3=0,&
                ibasnet=0,&
                phase=1,&
                iprezams=1
        real(kindreal) :: &
                binm2=0.d0,&
                periodini=0.d0
        ! **** Chemical composition
        integer :: &
                iopac=3,&
                ikappa=5
        real(kindreal) :: &
                initial_metallicity,&
                zsol=1.40d-2,&
                z
        ! **** Rotation-linked parameters
        integer :: &
                idiff=0,&
                iadvec=0,&
                istati=0,&
                icoeff=11,&
                igamma=0,&
                idialo=0,&
                idialu=0,&
                n_mag=1,&
                nsmooth=1
        real(kindreal) :: &
                fenerg=1.0d0,&
                richac=1.0d0,&
                frein=0.0d0,&
                K_Kawaler=0.0d0,&
                Omega_saturation=14.d0,&
                rapcrilim,&
                xfom=1.0d0,&
                omega,&
                xdial=0.0,&
                B_initial=0.d0,&
                add_diff=0.0d0,&
                alpha_F=1.d0
        real(kindreal) :: &
                zams_velocity=0.0d0 ! vwant
        logical :: &
                Add_Flux=.true.,&
                diff_only=.false.,&
                qminsmooth=.false.
        ! **** Winds parameters
        integer :: &
                OB_Mdot,&
                RSG_Mdot,&
                WR_Mdot,&
                Fallback_Mdot
        real(kindreal) :: &
                fmlos=0.85d0,&
                Be_mdotfrac=0.0d0,&
                start_mdot=0.80d0,&
                Z_dep=0.85d0,&
                Xs_WR=0.3d0,&
                D_clump=10.d0
        logical :: &
                SupraEddMdot=.true.,&
                hardJump=.true.,&
                print_winds=.false.,&
                prezams_winds_not_applied=.true.,&
                winds_not_applied=.false.
        ! **** Surface parameters
        integer :: &
                ifitm,&
                nndr=1
        real(kindreal) :: &
                fitm,&
                fitmi,&
                fitmi_default,&
                deltal=0.02d0,&
                deltat=0.02d0
        ! **** Convection-linked parameters
        integer :: &
                iledou=0,&
                idifcon=0,&
                my,&
                iover=1,&
                iunder=0
        real(kindreal):: &
                elph,&
                dovhp,&
                dunder=0.0d0
        ! **** Convergence-linked parameters
        integer :: &
                nbchx=200,&
                nrband=1
        real(kindreal) :: &
                gkorm=9.d0,&
                alph=0.3d0,&
                agdr=1.d-5,&
                faktor=1.d0,&
                dgrp=0.01d0,&
                dgrl=0.01d0,&
                dgry=0.003d0,&
                dgrc=0.01d0,&
                dgro=0.010d0,&
                dgr20=0.010d0
        ! **** Timestep controle
        integer :: &
                islow=2,&
                icncst=0,&
                tauH_fit=1
        real(kindreal) :: &
                xcn=1.d0
        ! **** Other controls
        integer :: &
                iauto,&
                iprn=10,&
                iout=0,&
                itmin=5,&
                idebug=0,&
                itests=0,&
                n_snap=10
        logical :: &
                display_plot,&
                xyfiles=.false.,&
                verbose=.false.,&
                stop_deg=.false.

        ! bfile stuff
        integer :: &
                m
        real(kindreal) :: &
                gms,&
                alter,&
                gls,&
                teff,&
                glsv,&
                teffv,&
                dzeitj,&
                dzeit,&
                dzeitv,&
                xmini,&
                xini,&
                summas,&
                ab,&
                dm_lost

        !real(kindreal), dimension(15,ldi) :: mainnam  ! FIXME to replace x, y3,y,...?
        real(kindreal), dimension(ldi) :: &
                q,&
                p,&
                t,&
                r,&
                s,&
                x,&
                y3,&
                y,&
                xc12,&
                xc13,&
                xn14,&
                xn15,&
                xo16,&
                xo17,&
                xo18,&
                xne20,&
                xne22,&
                xmg24,&
                xmg25,&
                xmg26,&
                xf19,&
                xne21,&
                xna23,&
                xal27,&
                xsi28,&
                xc14,&
                xf18,&
                xal26,&
                xneut,&
                xprot,&
                omegi,&
                xbid,&
                xbid1,&
                vp,&
                vt,&
                vr,&
                vs,&
                vx,&
                vy,&
                vy3,&
                vxc12,&
                vxc13,&
                vxn14,&
                vxn15,&
                vxo16,&
                vxo17,&
                vxo18,&
                vxne20,&
                vxne22,&
                vxmg24,&
                vxmg25,&
                vxmg26,&
                vxf19,&
                vxne21,&
                vxna23,&
                vxal27,&
                vxsi28,&
                vxc14,&
                vxf18,&
                vxal26,&
                vxneut,&
                vxprot,&
                vomegi,&
                vxbid,&
                vxbid1
        real(kindreal), dimension(3) :: &
                drl,&
                drte,&
                drp,&
                drt,&
                drr
        real(kindreal) :: &
                dk,&
                rlp,&
                rlt,&
                rlc,&
                rrp,&
                rrt,&
                rrc,&
                rtp,&
                rtt,&
                rtc,&
                tdiff,&
                suminenv,&
                vsuminenv
        real(kindreal), dimension(npondcouche) :: &
                CorrOmega
        real(kindreal) :: &
                xLtotbeg,&
                dlelexprev,&
                zams_radius

        ! netalu stuff
        real(kindreal), dimension(5) :: &
                xnetalu

        ! netdef stuff
        integer :: &
                mbelx  ! note that mbelx is imported, may have to sync
        real(kindreal) :: &
                xlostneu
        integer, dimension (mbelx) :: &
                nbzel,&
                nbael
        real(kindreal), dimension (mbelx) :: &
                abels

        ! Stuff for plotting
        real(kindreal), dimension(ldi) :: &
                Nabla_rad,&
                Nabla_ad,&
                Nabla_mu
        real(kindreal), dimension(ldi) :: &
                eps,&
                epsy,&
                eps_c_adv,&
                eps_ne_adv,&
                eps_o_adv,&
                eps_si_adv,&
                eps_grav,&
                eps_nu

        ! not for plotting
        integer :: &
                inum,&
                nsugi,&
                imloss,&
                id1
        real(kindreal) :: &
                period,&
                r_core,&
                vna,&
                vnr
        real(kindreal) :: &
                xtefflast,&
                xllast,&
                xrholast,&
                xclast,&
                xtclast

        real(kindreal), dimension (mbelx,ldi) :: &
                abelx,&
                vabelx
        logical :: &
                veryFirst

        ! stuff from the initial setup of the star
        real(kindreal) :: &
                initial_mass
        integer :: &
                idefaut ! use_default
        integer :: &
                ipoly ! use_polytrope
        real(kindreal) :: &
                n ! polytropic index
        integer :: &
                source,&
                alpha,& ! alpha_enhanced
                formatx ! output_format
    end type

    ! Initialise a genec_star
    !type(genec_star_ini) :: InitialGenecStar
    type(genec_star) :: GenecStar
end module storage
