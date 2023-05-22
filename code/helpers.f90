module helpers
    use evol, only: kindreal,ldi,npondcouche
    use abundmod, only: mbelx
    use storage, only: InitialGenecStar,GenecStar,genec_star

    use strucmod, only: &
            m,&
            q,p,t,r,s,&
            vp,vt,vr,vs,&
            drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,&
            Nabla_rad,Nabla_ad,Nabla_mu,&
            vna,vnr
    use diffadvmod, only: tdiff
    use rotmod, only: &
            omegi,vomegi,suminenv,vsuminenv,CorrOmega,dlelexprev
    use bintidemod, only: period
    use abundmod, only: &
            x,y,xc12,&
            xo16,vx,vy,vxc12,vxo16,&
            y3,xc13,xn14,xn15,xo17,xo18,vy3,vxc13,vxn14,vxn15,vxo17,&
            vxo18,xne20,xne22,xmg24,xmg25,xmg26,vxne20,vxne22,vxmg24,&
            vxmg25,vxmg26,&
            xf19,xne21,xna23,xal26,xal27,xsi28,vxf19,&
            vxne21,vxna23,vxal26g,vxal27,vxsi28,xneut,xprot,&
            xc14,xf18,xbid,xbid1,vxneut,vxprot,vxc14,vxf18,&
            vxbid,vxbid1,&
            xlostneu,nbzel,nbael,abels,&
            eps,epsy,eps_c_adv,eps_ne_adv,eps_o_adv,eps_si_adv,eps_grav,eps_nu,&
            abelx,vabelx
    use caramodele, only: &
            gms,gls,teff,glsv,teffv,dm_lost,xmini,ab,xLtotbeg,zams_radius,radius,&
            xtefflast,xllast,xrholast,xclast,xtclast,&
            inum

    use timestep, only: alter,dzeitj,dzeit,dzeitv
    use genec, only: xnetalu,summas,nwmd,veryFirst

    use inputparam, only: &
            starname,&
            nwseq,modanf,nzmod,end_at_phase,end_at_model,&
            var_rates,bintide,const_per,&
            irot,isol,imagn,ialflu,ianiso,ipop3,ibasnet,phase,iprezams,&
            binm2,periodini,&
            iopac,ikappa,&
            zinit,zsol,z,&
            idiff,iadvec,istati,icoeff,igamma,idialo,idialu,n_mag,nsmooth,&
            fenerg,richac,frein,K_Kawaler,Omega_saturation,rapcrilim,vwant,xfom,omega,xdial,B_initial,add_diff,alpha_F,&
            Add_Flux,diff_only,qminsmooth,&
            imloss,ifitm,nndr,RSG_Mdot,&
            fmlos,fitm,fitmi,fitmi_default,deltal,deltat,Be_mdotfrac,start_mdot,&
            SupraEddMdot,&
            iledou,idifcon,my,iover,iunder,&
            elph,dovhp,dunder,&
            nbchx,nrband,&
            gkorm,alph,agdr,faktor,dgrp,dgrl,dgry,dgrc,dgro,dgr20,&
            islow,icncst,tauH_fit,&
            xcn,&
            iauto,iprn,iout,itmin,idebug,itests,n_snap,&
            display_plot,xyfiles,verbose,stop_deg

    use henyey_solver, only: nsugi
    use convection, only: r_core

    implicit none

contains

    !subroutine copy_network_to_star(Network, Star)
    !    implicit none
    !    type(genec_star), intent(inout) :: Star
    !    type(genec_network_ini), intent(in) :: Network
    !    Star%xnetalu = Network%xnetalu
    !    Star%xlostneu = Network%xlostneu
    !    Star%nbzel = Network%nbzel
    !    Star%nbael = Network%nbael
    !    Star%abels = Network%abels
    !end subroutine copy_network_to_star
    !

    subroutine copy_to_genec_star(Star)
        ! copy all values to Star

        implicit none
        type(genec_star), intent(inout) :: Star

        ! Characteristics

        Star%starname         = starname
        Star%nwmd             = nwmd
        Star%nwseq            = nwseq            
        Star%modanf           = modanf
        Star%nzmod            = nzmod
        Star%end_at_phase     = end_at_phase
        Star%end_at_model     = end_at_model

        !Physics
        Star%irot             = irot
        Star%isol             = isol
        Star%imagn            = imagn
        Star%ialflu           = ialflu
        Star%ianiso           = ianiso
        Star%ipop3            = ipop3
        Star%ibasnet          = ibasnet
        Star%phase            = phase
        Star%iprezams         = iprezams
        Star%var_rates        = var_rates
        Star%bintide          = bintide
        Star%binm2            = binm2
        Star%periodini        = periodini
        Star%const_per        = const_per

        ! Composition
        Star%zinit            = zinit
        Star%zsol             = zsol
        Star%z                = z
        Star%iopac            = iopac
        Star%ikappa           = ikappa

        ! Rotation
        Star%idiff            = idiff
        Star%iadvec           = iadvec
        Star%istati           = istati
        Star%icoeff           = icoeff
        Star%igamma           = igamma
        Star%idialo           = idialo
        Star%idialu           = idialu
        Star%n_mag            = n_mag
        Star%nsmooth          = nsmooth
        Star%fenerg           = fenerg
        Star%richac           = richac
        Star%frein            = frein
        Star%K_Kawaler        = K_Kawaler
        Star%Omega_saturation = Omega_saturation
        Star%rapcrilim        = rapcrilim
        Star%vwant            = vwant
        Star%xfom             = xfom
        Star%omega            = omega
        Star%xdial            = xdial
        Star%B_initial        = B_initial
        Star%add_diff         = add_diff
        Star%alpha_F          = alpha_F
        Star%Add_Flux         = Add_Flux
        Star%diff_only        = diff_only
        Star%qminsmooth       = qminsmooth

        ! Surface
        Star%imloss           = imloss
        Star%ifitm            = ifitm
        Star%nndr             = nndr
        Star%RSG_Mdot         = RSG_Mdot
        Star%fmlos            = fmlos
        Star%fitm             = fitm
        Star%fitmi            = fitmi
        Star%fitmi_default    = fitmi_default
        Star%deltal           = deltal
        Star%deltat           = deltat
        Star%Be_mdotfrac      = Be_mdotfrac
        Star%start_mdot       = start_mdot
        Star%SupraEddMdot     = SupraEddMdot

        ! Convection
        Star%iledou           = iledou
        Star%idifcon          = idifcon
        Star%my               = my
        Star%iover            = iover
        Star%iunder           = iunder
        Star%elph             = elph
        Star%dovhp            = dovhp
        Star%dunder           = dunder

        ! Convergence
        Star%gkorm            = gkorm
        Star%alph             = alph
        Star%agdr             = agdr
        Star%faktor           = faktor
        Star%dgrp             = dgrp
        Star%dgrl             = dgrl
        Star%dgry             = dgry
        Star%dgrc             = dgrc
        Star%dgro             = dgro
        Star%dgr20            = dgr20
        Star%nbchx            = nbchx
        Star%nrband           = nrband

        ! Time
        Star%xcn              = xcn
        Star%islow            = islow
        Star%icncst           = icncst
        Star%tauH_fit         = tauH_fit

        ! Various
        Star%iauto            = iauto
        Star%iprn             = iprn
        Star%iout             = iout
        Star%itmin            = itmin
        Star%idebug           = idebug
        Star%itests           = itests
        Star%n_snap           = n_snap
        Star%display_plot     = display_plot
        Star%xyfiles          = xyfiles
        Star%verbose          = verbose
        Star%stop_deg         = stop_deg

        ! Properties
        Star%m           = m
        Star%gms         = gms
        Star%alter       = alter
        Star%gls         = gls
        Star%teff        = teff
        Star%glsv        = glsv
        Star%teffv       = teffv
        Star%dzeitj      = dzeitj
        Star%dzeit       = dzeit
        Star%dzeitv      = dzeitv
        Star%summas      = summas
        Star%xmini       = xmini
        Star%ab          = ab
        Star%dm_lost     = dm_lost

        ! Extra
        Star%xtefflast   = xtefflast
        Star%xllast      = xllast
        Star%xrholast    = xrholast
        Star%xclast      = xclast
        Star%xtclast     = xtclast
        Star%inum        = inum
        Star%nsugi       = nsugi
        Star%period      = period
        Star%r_core      = r_core
        Star%vna         = vna
        Star%vnr         = vnr

        ! Structure
        Star%q           = q
        Star%p           = p
        Star%t           = t
        Star%r           = r
        Star%s           = s
        Star%x           = x
        Star%y3          = y3
        Star%y           = y
        Star%xc12        = xc12
        Star%xc13        = xc13
        Star%xc14        = xc14
        Star%xn14        = xn14
        Star%xn15        = xn15
        Star%xo16        = xo16
        Star%xo17        = xo17
        Star%xo18        = xo18
        Star%xf18        = xf18
        Star%xf19        = xf19
        Star%xne20       = xne20
        Star%xne21       = xne21
        Star%xne22       = xne22
        Star%xna23       = xna23
        Star%xmg24       = xmg24
        Star%xmg25       = xmg25
        Star%xmg26       = xmg26
        Star%xal26       = xal26
        Star%xal27       = xal27
        Star%xsi28       = xsi28
        Star%xneut       = xneut
        Star%xprot       = xprot
        Star%xbid        = xbid
        Star%xbid1       = xbid1
        Star%omegi       = omegi
        Star%vp          = vp
        Star%vt          = vt
        Star%vr          = vr
        Star%vs          = vs
        Star%vx          = vx
        Star%vy3         = vy3
        Star%vy          = vy
        Star%vxc12       = vxc12
        Star%vxc13       = vxc13
        Star%vxc14       = vxc14
        Star%vxn14       = vxn14
        Star%vxn15       = vxn15
        Star%vxo16       = vxo16
        Star%vxo17       = vxo17
        Star%vxo18       = vxo18
        Star%vxf18       = vxf18
        Star%vxf19       = vxf19
        Star%vxne20      = vxne20
        Star%vxne21      = vxne21
        Star%vxne22      = vxne22
        Star%vxna23      = vxna23
        Star%vxmg24      = vxmg24
        Star%vxmg25      = vxmg25
        Star%vxmg26      = vxmg26
        Star%vxal26      = vxal26g
        Star%vxal27      = vxal27
        Star%vxsi28      = vxsi28
        Star%vxneut      = vxneut
        Star%vxprot      = vxprot
        Star%vxbid       = vxbid
        Star%vxbid1      = vxbid1
        Star%vomegi      = vomegi

        ! Structure extra
        Star%mbelx       = mbelx
        Star%abelx       = abelx
        Star%vabelx      = vabelx

        ! network
        Star%xlostneu    = xlostneu
        Star%nbzel       = nbzel
        Star%nbael       = nbael
        Star%abels       = abels

        ! no need to restore (?)
        Star%drl         = drl
        Star%drte        = drte
        Star%dk          = dk
        Star%drp         = drp
        Star%drt         = drt
        Star%drr         = drr
        Star%rlp         = rlp
        Star%rlt         = rlt
        Star%rlc         = rlc
        Star%rrp         = rrp
        Star%rrt         = rrt
        Star%rrc         = rrc
        Star%rtp         = rtp
        Star%rtt         = rtt
        Star%rtc         = rtc
        Star%tdiff       = tdiff
        Star%suminenv    = suminenv
        Star%CorrOmega   = CorrOmega
        Star%xLtotbeg    = xLtotbeg
        Star%dlelexprev  = dlelexprev
        Star%zams_radius = zams_radius
        !Star%xnetalu     = xnetalu

        Star%radius      = radius

        Star%Nabla_rad   = Nabla_rad
        Star%Nabla_ad    = Nabla_ad
        Star%Nabla_mu    = Nabla_mu
        Star%eps         = eps
        Star%epsy        = epsy
        Star%eps_c_adv   = eps_c_adv
        Star%eps_ne_adv  = eps_ne_adv
        Star%eps_o_adv   = eps_o_adv
        Star%eps_si_adv  = eps_si_adv
        Star%eps_grav    = eps_grav
        Star%eps_nu      = eps_nu



        Star%veryFirst = veryFirst
        Star%synchronised = .true.

        !write(*,*) Star
        !write(*,*) "After sync: name/mass: ", Star%starname, Star%gms, gms
    end subroutine copy_to_genec_star

    subroutine copy_namelists_from_genec_star(Star)
        implicit none
        type(genec_star), intent(in) :: Star

        !CharacteristicsParams
        starname = Star%starname
        nwseq = Star%nwseq
        nwmd = Star%nwmd
        modanf = Star%modanf
        nzmod = Star%nzmod
        end_at_phase = Star%end_at_phase
        end_at_model = Star%end_at_model

        !PhysicsParams
        irot = Star%irot
        isol = Star%isol
        imagn = Star%imagn
        ialflu = Star%ialflu
        ianiso = Star%ianiso
        ipop3 = Star%ipop3
        ibasnet = Star%ibasnet
        phase = Star%phase
        var_rates = Star%var_rates
        bintide = Star%bintide
        binm2 = Star%binm2
        periodini = Star%periodini
        const_per = Star%const_per
        iprezams = Star%iprezams

        !CompositionParams
        zinit = Star%zinit
        zsol = Star%zsol
        z = Star%z
        iopac = Star%iopac
        ikappa = Star%ikappa

        !RotationParams
        idiff = Star%idiff
        iadvec = Star%iadvec
        istati = Star%istati
        icoeff = Star%icoeff
        fenerg = Star%fenerg
        richac = Star%richac
        igamma = Star%igamma
        frein = Star%frein
        K_Kawaler = Star%K_Kawaler
        Omega_saturation = Star%Omega_saturation
        rapcrilim = Star%rapcrilim
        vwant = Star%vwant
        xfom = Star%xfom
        omega = Star%omega
        xdial = Star%xdial
        idialo = Star%idialo
        idialu = Star%idialu
        Add_Flux = Star%Add_Flux
        diff_only = Star%diff_only
        B_initial = Star%B_initial
        add_diff = Star%add_diff
        n_mag = Star%n_mag
        alpha_F = Star%alpha_F
        nsmooth = Star%nsmooth
        qminsmooth = Star%qminsmooth

        !SurfaceParams
        imloss = Star%imloss
        fmlos = Star%fmlos
        ifitm = Star%ifitm
        fitm = Star%fitm
        fitmi = Star%fitmi
        deltal = Star%deltal
        deltat = Star%deltat
        nndr = Star%nndr
        RSG_Mdot = Star%RSG_Mdot
        SupraEddMdot = Star%SupraEddMdot
        Be_mdotfrac = Star%Be_mdotfrac
        start_mdot = Star%start_mdot

        !ConvectionParams
        iledou = Star%iledou
        idifcon = Star%idifcon
        iover = Star%iover
        elph = Star%elph
        my = Star%my
        dovhp = Star%dovhp
        iunder = Star%iunder
        dunder = Star%dunder

        !ConvergenceParams
        gkorm = Star%gkorm
        alph = Star%alph
        agdr = Star%agdr
        faktor = Star%faktor
        dgrp = Star%dgrp
        dgrl = Star%dgrl
        dgry = Star%dgry
        dgrc = Star%dgrc
        dgro = Star%dgro
        dgr20 = Star%dgr20
        nbchx = Star%nbchx
        nrband = Star%nrband

        !TimeControle
        xcn = Star%xcn
        islow = Star%islow
        icncst = Star%icncst
        tauH_fit = Star%tauH_fit

        !VariousSettings
        display_plot = Star%display_plot
        iauto = Star%iauto
        iprn = Star%iprn
        iout = Star%iout
        itmin = Star%itmin
        xyfiles = Star%xyfiles
        idebug = Star%idebug
        itests = Star%itests
        verbose = Star%verbose
        stop_deg = Star%stop_deg
        n_snap = Star%n_snap
    end subroutine copy_namelists_from_genec_star

    subroutine copy_initial_structure_from_genec_star(Star)
        implicit none
        type(genec_star), intent(in) :: Star
        gms = Star%gms
        alter = Star%alter
        gls = Star%gls
        teff = Star%teff
        glsv = Star%glsv
        teffv = Star%teffv
        dzeitj = Star%dzeitj
        dzeit = Star%dzeit
        dzeitv = Star%dzeitv
        summas = Star%summas
        ab = Star%ab
        m = Star%m
        q = Star%q
        p = Star%p
        t = Star%t
        r = Star%r
        s = Star%s
        vp = Star%vp
        vt = Star%vt
        vr = Star%vr
        vs = Star%vs
        x = Star%x
        y3 = Star%y3
        y = Star%y
        xc12 = Star%xc12
        xc13 = Star%xc13
        xn14 = Star%xn14
        xn15 = Star%xn15
        xo16 = Star%xo16
        xo17 = Star%xo17
        xo18 = Star%xo18
        xne20 = Star%xne20
        xne22 = Star%xne22
        xmg24 = Star%xmg24
        xmg25 = Star%xmg25
        xmg26 = Star%xmg26
        omegi = Star%omegi

    end subroutine copy_initial_structure_from_genec_star

    subroutine copy_structure_from_genec_star(Star)
        implicit none
        type(genec_star), intent(in) :: Star

        m = Star%m
        gms = Star%gms
        write(*,*) "alter:", alter, "copying from Star%alter:", Star%alter
        alter = Star%alter
        gls = Star%gls
        teff = Star%teff
        glsv = Star%glsv
        teffv = Star%teffv
        dzeitj = Star%dzeitj
        dzeit = Star%dzeit
        dzeitv = Star%dzeitv
        xmini = Star%xmini
        summas = Star%summas
        ab = Star%ab
        dm_lost = Star%dm_lost
        !mbelx = Star%mbelx
        xtefflast = Star%xtefflast
        xllast = Star%xllast
        xrholast = Star%xrholast
        xclast = Star%xclast
        xtclast = Star%xtclast
        inum = Star%inum
        nsugi = Star%nsugi
        period = Star%period
        r_core = Star%r_core
        vna = Star%vna
        vnr = Star%vnr
        q = Star%q
        p = Star%p
        t = Star%t
        r = Star%r
        s = Star%s
        x = Star%x
        y3 = Star%y3
        y = Star%y
        xc12 = Star%xc12
        xc13 = Star%xc13
        xn14 = Star%xn14
        xn15 = Star%xn15
        xo16 = Star%xo16
        xo17 = Star%xo17
        xo18 = Star%xo18
        xne20 = Star%xne20
        xne22 = Star%xne22
        xmg24 = Star%xmg24
        xmg25 = Star%xmg25
        xmg26 = Star%xmg26
        xf19 = Star%xf19
        xne21 = Star%xne21
        xna23 = Star%xna23
        xal27 = Star%xal27
        xsi28 = Star%xsi28
        xc14 = Star%xc14
        xf18 = Star%xf18
        xal26 = Star%xal26
        xneut = Star%xneut
        xprot = Star%xprot
        omegi = Star%omegi
        xbid = Star%xbid
        xbid1 = Star%xbid1
        vp = Star%vp
        vt = Star%vt
        vr = Star%vr
        vs = Star%vs
        vx = Star%vx
        vy = Star%vy
        vy3 = Star%vy3
        vxc12 = Star%vxc12
        vxc13 = Star%vxc13
        vxn14 = Star%vxn14
        vxn15 = Star%vxn15
        vxo16 = Star%vxo16
        vxo17 = Star%vxo17
        vxo18 = Star%vxo18
        vxne20 = Star%vxne20
        vxne22 = Star%vxne22
        vxmg24 = Star%vxmg24
        vxmg25 = Star%vxmg25
        vxmg26 = Star%vxmg26
        vxf19 = Star%vxf19
        vxne21 = Star%vxne21
        vxna23 = Star%vxna23
        vxal27 = Star%vxal27
        vxsi28 = Star%vxsi28
        vxc14 = Star%vxc14
        vxf18 = Star%vxf18
        vxal26g = Star%vxal26
        vxneut = Star%vxneut
        vxprot = Star%vxprot
        vomegi = Star%vomegi
        vxbid = Star%vxbid
        vxbid1 = Star%vxbid1
        abelx = Star%abelx
        vabelx = Star%vabelx
        veryFirst = Star%veryFirst
        write(*,*) 'age set to ', alter
        flush(6)

    end subroutine copy_structure_from_genec_star

    subroutine copy_netdef_from_genec_star(Star)
        implicit none
        type(genec_star), intent(in) :: Star
        xnetalu     = Star%xnetalu    
        xlostneu    = Star%xlostneu   
        nbzel       = Star%nbzel      
        nbael       = Star%nbael      
        abels       = Star%abels      
    end subroutine copy_netdef_from_genec_star

    subroutine copy_from_genec_star(Star)
        implicit none
        type(genec_star), intent(in) :: Star
        call copy_namelists_from_genec_star(Star)
        call copy_netdef_from_genec_star(Star)
        call copy_structure_from_genec_star(Star)
    end subroutine copy_from_genec_star
end module helpers
