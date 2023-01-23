module storage
    use evol, only: kindreal,ldi,npondcouche
    use abundmod, only: mbelx

    implicit none
    type genec_star_ini
        ! this type carries all the information needed to initialise a star
        ! i.e. it has all the input that is normally read by the makeini program
        ! calling makeini with this 'genec_star_ini' type will then set up a 'genec_star'
        integer :: index_of_the_star  ! AMUSE property
        character(256) :: starname
        real(kindreal) :: mstar, zini, vwant
        integer :: idefaut
        integer :: ipoly
        real(kindreal) :: n
        integer :: source,alpha,formatx
    end type

    type genec_network_ini
        ! contains netdef.in and netalu.dat
        !netdef
        real(kindreal) :: xlostneu
        integer, dimension(mbelx) :: nbzel,nbael
        real(kindreal), dimension(mbelx) :: abels

        !netalu
        real(kindreal), dimension(5) :: xnetalu
    end type

    type genec_star
        ! genec_star type contains all variables needed to initialise a star
        ! it can be used to save information about a star and if needed roll back

        ! **** AMUSE specific
        logical :: initialised=.false.  ! .true. if the star has been synchronised to at any time
        logical :: synchronised  ! .true. if the star is up-to-date, set to .false. when the model changes
        integer :: index_of_the_star  ! AMUSE-specific, do not change here!
        real(kindreal) :: radius
        logical :: stopped  ! a STOP will set this to .true. and then return
        character(256) :: stop_message

        ! Variables that would go in the .input file
        ! **** Model characteristics
        character(256) :: starname
        integer :: nwmd
        integer :: nwseq,modanf,nzmod,end_at_phase,end_at_model

        ! **** Physical inputs
        logical :: var_rates,bintide,const_per
        integer :: irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,iprezams
        real(kindreal) :: binm2,periodini
        ! **** Chemical composition
        integer :: iopac,ikappa
        real(kindreal) :: zinit,zsol,z
        ! **** Rotation-linked parameters    
        integer :: idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth
        real(kindreal) :: &
                fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,&
                omega,xdial,B_initial,add_diff,alpha_F
        logical :: Add_Flux,diff_only,qminsmooth
        ! **** Surface parameters
        integer :: imloss,ifitm,nndr,RSG_Mdot
        real(kindreal) :: fmlos,fitm,fitmi,fitmi_default,deltal,deltat,Be_mdotfrac,start_mdot
        logical :: SupraEddMdot
        ! **** Convection-linked parameters
        integer :: iledou,idifcon,my,iover,iunder
        real(kindreal):: elph,dovhp,dunder
        ! **** Convergence-linked parameters
        integer :: nbchx,nrband
        real(kindreal) :: gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20
        ! **** Timestep controle
        integer :: islow,icncst,tauH_fit
        real(kindreal) :: xcn
        ! **** Other controles
        integer :: iauto,iprn,iout,itmin,idebug,itests,n_snap
        logical :: display_plot,xyfiles,verbose,stop_deg

        ! bfile stuff
        integer :: m
        real(kindreal) :: &
                gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,summas,ab,&
                dm_lost

        !real(kindreal), dimension(15,ldi) :: mainnam  ! FIXME to replace x, y3,y,...?
        real(kindreal), dimension(ldi) :: &
                q,&
                p,t,r,s,&
                x,y3,y,xc12,xc13,xn14,xn15,xo16,xo17,xo18,&
                xne20,xne22,xmg24,xmg25,xmg26,xf19,xne21,xna23,&
                xal27,xsi28,&
                xc14,xf18,xal26,&
                xneut,xprot,&
                omegi,&
                xbid,xbid1,&
                vp,vt,vr,vs,&
                vx,vy,vy3,vxc12,vxc13,vxn14,vxn15,vxo16,vxo17,vxo18,&
                vxne20,vxne22,vxmg24,vxmg25,vxmg26,vxf19,vxne21,vxna23,&
                vxal27,vxsi28,&
                vxc14,vxf18,vxal26,&
                vxneut,vxprot,&
                vomegi,&
                vxbid,vxbid1
        real(kindreal), dimension(3) :: &
                drl,drte,drp,drt,drr
        real(kindreal) :: &
                dk,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,tdiff,suminenv,vsuminenv
        real(kindreal), dimension(npondcouche) :: &
                CorrOmega
        real(kindreal) :: &
                xLtotbeg,dlelexprev,zams_radius

        ! netalu stuff
        real(kindreal), dimension(5) :: &
                xnetalu

        ! netdef stuff
        integer :: mbelx  ! note that mbelx is imported, may have to sync
        real(kindreal) :: &
                xlostneu
        integer, dimension (mbelx) :: &
                nbzel,nbael
        real(kindreal), dimension (mbelx) :: &
                abels

        ! Stuff for plotting
        real(kindreal), dimension(ldi) :: Nabla_rad,Nabla_ad,Nabla_mu
        real(kindreal), dimension(ldi) :: &
                eps,epsy,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv,eps_grav,eps_nu

        integer :: &
                inum,nsugi
        real(kindreal) :: &
                period,r_core,vna,vnr
        real(kindreal) :: &
                xteffprev,xtefflast,xlprev,xllast,xrhoprev,xrholast,&
                xcprev,xclast,xtcprev,xtclast

        real(kindreal), dimension (mbelx,ldi) :: abelx,vabelx

    end type

    type(genec_star_ini) :: InitialGenecStar
    type(genec_network_ini) :: InitialNetwork
    type(genec_star) :: GenecStar
end module storage
