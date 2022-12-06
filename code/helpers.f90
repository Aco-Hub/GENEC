module helpers
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
    type genec_star
        ! genec_star type contains all variables needed to initialise a star
        ! it can be used to save information about a star and if needed roll back

        ! **** AMUSE specific
        logical :: initialised
        integer :: index_of_the_star

        ! Variables that would go in the .input file
        ! **** Model characteristics
        character(256) :: starname
        integer :: nwseq=1,modanf=0,nzmod=1,end_at_phase=4,end_at_model=0

        ! **** Physical inputs
        integer :: irot=0,isol=0,imagn,ialflu=0,ianiso,ipop3,ibasnet,phase=1,iprezams
        real(kindreal) :: binm2,periodini
        ! **** Chemical composition
        integer :: iopac,ikappa
        real(kindreal) :: zinit,zsol,z
        ! **** Rotation-linked parameters    
        integer :: idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth
        real(kindreal) :: fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,omega,xdial,B_initial,add_diff,alpha_F
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
                gms,alter,gls,teff,glsv,teffv,dzeitj,dzeit,dzeitv,xmini,ab,&
                dm_lost

        real(kindreal), dimension(ldi) :: &
                q,p,t,r,s,x,y,xc12,vp,vt,vr,vs,xo16,vx,vy,vxc12,vxo16,&
                y3,xc13,xn14,xn15,xo17,xo18,vy3,vxc13,vxn14,vxn15,vxo17,&
                vxo18,xne20,xne22,xmg24,xmg25,xmg26,vxne20,vxne22,vxmg24,&
                vxmg25,vxmg26,omegi,vomegi,&
                xf19,xne21,xna23,xal26,xal27,xsi28,vxf19,&
                vxne21,vxna23,vxal26g,vxal27,vxsi28,xneut,xprot,&
                xc14,xf18,xbid,xbid1,vxneut,vxprot,vxc14,vxf18,&
                vxbid,vxbid1
        real(kindreal) :: &
                drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,tdiff,suminenv
        real(kindreal), dimension(npondcouche) :: &
                CorrOmega
        real(kindreal) :: &
                xltotbeg,dlelexprev,zams_radius

        ! netalu stuff
        real(kindreal), dimension(5) :: &
                xnetalu

        ! netdef stuff
        real(kindreal) :: &
                xlostneu
        integer, dimension (mbelx) :: &
                nbzel,nbael
        real(kindreal), dimension (mbelx) :: &
                abels

    end type

    type(genec_star_ini) :: InitialGenecStar
    type(genec_star) :: GenecStar

contains

    subroutine input_ini
        implicit none
        character(len=1):: answer
        character(len=4), dimension(4):: sourceid !< source identifier used in output file name

        
        write(*,*)'Enter the star name:'
        read(5,*) InitialGenecStar%starname
        write(*,*)'Enter the desired mass and metallicity:'
        read(5,*) InitialGenecStar%mstar,InitialGenecStar%zini
        write(*,*) 'Which rotation velocity on the ZAMS?'
        read(5,*) InitialGenecStar%vwant

        !!!!! inichemmod bit
        answer = ''
    
        do while (answer /= 'y' .and. answer /= 'n' .and. answer /= '1' .and. answer /= '0')
            write(*,*)'Default settings are:'
            write(*,*)'         - Asplund-Cunha abundances'
            write(*,*)'         - scaled solar'
            write(*,*)'         - Geneva format'
            write(*,*)'         - structure from pre-calculated model'
            write(*,*)'Is it ok? (y)es (n)o'
            read(5,*) answer
            if (answer /= 'y' .and. answer /= 'n'  .and. answer /= '1' .and. answer /= '0') then
                write(*,*) 'Please type y or n...'
            endif
        enddo
        if (answer == 'y' .or. answer == 'Y' .or. answer == '1') then
            InitialGenecStar%idefaut = 1
        elseif (answer == 'n' .or. answer == 'N' .or. answer == '0') then
            InitialGenecStar%idefaut = 0
        endif

        ! source identifier
        sourceid=(/'AG89','GN93','As05','As09'/)

        if (InitialGenecStar%idefaut == 0) then
        ! choose input file
            write(*,*) 'Choose source'
            write(*,*)'Anders & Grevesse 1989: 1 / Grevesse and Noels 1993: 2 / Asplund 2005: 3 / Asplund 2009: 4'
            read(5,*) InitialGenecStar%source
            write(*,*) sourceid(int(InitialGenecStar%source))
        ! choose alpha enhanced or solar scaled composition
            write(*,*)'scaled solar abundances (0) or alpha-enhanced (1)?'
            read(5,*) InitialGenecStar%alpha
            write(*,*) 'alpha-enhanced:',InitialGenecStar%alpha
            if(InitialGenecStar%alpha/=0 .and. InitialGenecStar%alpha/=1) then
                write(*,*) InitialGenecStar
                stop 'wrong choice - program stopped!'
            endif
        ! choose output format
            write(*,*)'choose format! [1=basnet/2=GENEC/3=PPN input format]'
            read(5,*) InitialGenecStar%formatx
            if(InitialGenecStar%formatx/=1.and.InitialGenecStar%formatx/=2.and.InitialGenecStar%formatx/=3) then
                stop 'wrong format choice'
            endif
        endif  ! idefaut

        !!!!! end inichemmod bit
        select case (InitialGenecStar%idefaut)
            case (0)
                write(*,*)'Do you want to compute a polytropic structure or use ',&
                        'a pre-computed structure? Polytrope:1 - structure:0'
                read(5,*) InitialGenecStar%ipoly
            case (1)
                write(*,*)'Enter the polytropic index (recommended: 2.5):'
                read(5,*) InitialGenecStar%n
        end select
    end subroutine input_ini

end module helpers
