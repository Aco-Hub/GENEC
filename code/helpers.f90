module helpers
    use evol, only: kindreal,ldi,npondcouche
    use abundmod, only: mbelx
    use storage, only: InitialGenecStar,InitialNetwork,GenecStar,genec_star,genec_network_ini

    use strucmod, only: &
            m,&
            q,p,t,r,s,&
            vp,vt,vr,vs,&
            drl,drte,dk,drp,drt,drr,rlp,rlt,rlc,rrp,rrt,rrc,rtp,rtt,rtc,&
            Nabla_rad,Nabla_ad,Nabla_mu
    use diffadvmod, only: tdiff
    use rotmod, only: &
            omegi,vomegi,suminenv,CorrOmega,dlelexprev
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
            xlostneu,nbzel,nbael,abels
    use caramodele, only: &
            gms,gls,teff,glsv,teffv,dm_lost,xmini,ab,xLtotbeg,zams_radius,radius
    use timestep, only: alter,dzeitj,dzeit,dzeitv
    use genec, only: xnetalu,summas,nwmd

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

        Star%initialised = .true.

        Star%starname         = starname
        Star%nwmd             = nwmd
        Star%nwseq            = nwseq            
        Star%modanf           = modanf
        Star%nzmod            = nzmod
        Star%end_at_phase     = end_at_phase
        Star%end_at_model     = end_at_model
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
        Star%iopac            = iopac
        Star%ikappa           = ikappa
        Star%zinit            = zinit
        Star%zsol             = zsol
        Star%z                = z
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
        Star%iledou           = iledou
        Star%idifcon          = idifcon
        Star%my               = my
        Star%iover            = iover
        Star%iunder           = iunder
        Star%elph             = elph
        Star%dovhp            = dovhp
        Star%dunder           = dunder
        Star%nbchx            = nbchx
        Star%nrband           = nrband
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
        Star%islow            = islow
        Star%icncst           = icncst
        Star%tauH_fit         = tauH_fit
        Star%xcn              = xcn
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
        Star%q           = q
        Star%p           = p
        Star%t           = t
        Star%r           = r
        Star%s           = s
        Star%x           = x
        Star%y           = y
        Star%xc12        = xc12
        Star%vp          = vp
        Star%vt          = vt
        Star%vr          = vr
        Star%vs          = vs
        Star%xo16        = xo16
        Star%vx          = vx
        Star%vy          = vy
        Star%vxc12       = vxc12
        Star%vxo16       = vxo16
        Star%y3          = y3
        Star%xc13        = xc13
        Star%xn14        = xn14
        Star%xn15        = xn15
        Star%xo17        = xo17
        Star%xo18        = xo18
        Star%vy3         = vy3
        Star%vxc13       = vxc13
        Star%vxn14       = vxn14
        Star%vxn15       = vxn15
        Star%vxo17       = vxo17
        Star%vxo18       = vxo18
        Star%xne20       = xne20
        Star%xne22       = xne22
        Star%xmg24       = xmg24
        Star%xmg25       = xmg25
        Star%xmg26       = xmg26
        Star%vxne20      = vxne20
        Star%vxne22      = vxne22
        Star%vxmg24      = vxmg24
        Star%vxmg25      = vxmg25
        Star%vxmg26      = vxmg26
        Star%omegi       = omegi
        Star%vomegi      = vomegi
        Star%xf19        = xf19
        Star%xne21       = xne21
        Star%xna23       = xna23
        Star%xal26       = xal26
        Star%xal27       = xal27
        Star%xsi28       = xsi28
        Star%vxf19       = vxf19
        Star%vxne21      = vxne21
        Star%vxna23      = vxna23
        Star%vxal26g     = vxal26g
        Star%vxal27      = vxal27
        Star%vxsi28      = vxsi28
        Star%xneut       = xneut
        Star%xprot       = xprot
        Star%xc14        = xc14
        Star%xf18        = xf18
        Star%xbid        = xbid
        Star%xbid1       = xbid1
        Star%vxneut      = vxneut
        Star%vxprot      = vxprot
        Star%vxc14       = vxc14
        Star%vxf18       = vxf18
        Star%vxbid       = vxbid
        Star%vxbid1      = vxbid1
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
        !Star%xlostneu    = xlostneu
        !Star%nbzel       = nbzel
        !Star%nbael       = nbael
        !Star%abels       = abels

        Star%radius      = radius

        Star%Nabla_rad   = Nabla_rad
        Star%Nabla_ad    = Nabla_ad
        Star%Nabla_mu    = Nabla_mu

        Star%synchronised = .true.

        !write(*,*) Star
        !write(*,*) "After sync: name/mass: ", Star%starname, Star%gms, gms
    end subroutine copy_to_genec_star

    subroutine copy_from_genec_star(Star)
        implicit none
        type(genec_star), intent(inout) :: Star

        starname         = Star%starname         
        nwmd             = Star%nwmd
        nwseq            = Star%nwseq            
        modanf           = Star%modanf           
        nzmod            = Star%nzmod            
        end_at_phase     = Star%end_at_phase     
        end_at_model     = Star%end_at_model     
        irot             = Star%irot             
        isol             = Star%isol             
        imagn            = Star%imagn            
        ialflu           = Star%ialflu           
        ianiso           = Star%ianiso           
        ipop3            = Star%ipop3            
        ibasnet          = Star%ibasnet          
        phase            = Star%phase            
        iprezams         = Star%iprezams
        var_rates        = Star%var_rates
        bintide          = Star%bintide
        binm2            = Star%binm2            
        periodini        = Star%periodini
        const_per        = Star%const_per
        iopac            = Star%iopac            
        ikappa           = Star%ikappa           
        zinit            = Star%zinit            
        zsol             = Star%zsol             
        z                = Star%z                
        idiff            = Star%idiff            
        iadvec           = Star%iadvec           
        istati           = Star%istati           
        icoeff           = Star%icoeff           
        igamma           = Star%igamma           
        idialo           = Star%idialo           
        idialu           = Star%idialu           
        n_mag            = Star%n_mag            
        nsmooth          = Star%nsmooth          
        fenerg           = Star%fenerg           
        richac           = Star%richac           
        frein            = Star%frein            
        K_Kawaler        = Star%K_Kawaler        
        Omega_saturation = Star%Omega_saturation 
        rapcrilim        = Star%rapcrilim        
        vwant            = Star%vwant            
        xfom             = Star%xfom             
        omega            = Star%omega            
        xdial            = Star%xdial            
        B_initial        = Star%B_initial        
        add_diff         = Star%add_diff         
        alpha_F          = Star%alpha_F          
        Add_Flux         = Star%Add_Flux         
        diff_only        = Star%diff_only        
        qminsmooth       = Star%qminsmooth       
        imloss           = Star%imloss           
        ifitm            = Star%ifitm            
        nndr             = Star%nndr             
        RSG_Mdot         = Star%RSG_Mdot         
        fmlos            = Star%fmlos            
        fitm             = Star%fitm             
        fitmi            = Star%fitmi            
        fitmi_default    = Star%fitmi_default    
        deltal           = Star%deltal           
        deltat           = Star%deltat           
        Be_mdotfrac      = Star%Be_mdotfrac      
        start_mdot       = Star%start_mdot       
        SupraEddMdot     = Star%SupraEddMdot     
        iledou           = Star%iledou           
        idifcon          = Star%idifcon          
        my               = Star%my               
        iover            = Star%iover            
        iunder           = Star%iunder           
        elph             = Star%elph             
        dovhp            = Star%dovhp            
        dunder           = Star%dunder           
        nbchx            = Star%nbchx            
        nrband           = Star%nrband           
        gkorm            = Star%gkorm            
        alph             = Star%alph             
        agdr             = Star%agdr             
        faktor           = Star%faktor           
        dgrp             = Star%dgrp             
        dgrl             = Star%dgrl             
        dgry             = Star%dgry             
        dgrc             = Star%dgrc             
        dgro             = Star%dgro             
        dgr20            = Star%dgr20            
        islow            = Star%islow            
        icncst           = Star%icncst           
        tauH_fit         = Star%tauH_fit         
        xcn              = Star%xcn              
        iauto            = Star%iauto            
        iprn             = Star%iprn             
        iout             = Star%iout             
        itmin            = Star%itmin            
        idebug           = Star%idebug           
        itests           = Star%itests           
        n_snap           = Star%n_snap           
        display_plot     = Star%display_plot     
        xyfiles          = Star%xyfiles          
        verbose          = Star%verbose          
        stop_deg         = Star%stop_deg         

        m           = Star%m          
        gms         = Star%gms        
        alter       = Star%alter      
        gls         = Star%gls        
        teff        = Star%teff       
        glsv        = Star%glsv       
        teffv       = Star%teffv      
        dzeitj      = Star%dzeitj     
        dzeit       = Star%dzeit      
        dzeitv      = Star%dzeitv     
        xmini       = Star%xmini      
        summas      = Star%summas
        ab          = Star%ab         
        dm_lost     = Star%dm_lost    
        q           = Star%q          
        p           = Star%p          
        t           = Star%t          
        r           = Star%r          
        s           = Star%s          
        x           = Star%x          
        y           = Star%y          
        xc12        = Star%xc12       
        vp          = Star%vp         
        vt          = Star%vt         
        vr          = Star%vr         
        vs          = Star%vs         
        xo16        = Star%xo16       
        vx          = Star%vx         
        vy          = Star%vy         
        vxc12       = Star%vxc12      
        vxo16       = Star%vxo16      
        y3          = Star%y3         
        xc13        = Star%xc13       
        xn14        = Star%xn14       
        xn15        = Star%xn15       
        xo17        = Star%xo17       
        xo18        = Star%xo18       
        vy3         = Star%vy3        
        vxc13       = Star%vxc13      
        vxn14       = Star%vxn14      
        vxn15       = Star%vxn15      
        vxo17       = Star%vxo17      
        vxo18       = Star%vxo18      
        xne20       = Star%xne20      
        xne22       = Star%xne22      
        xmg24       = Star%xmg24      
        xmg25       = Star%xmg25      
        xmg26       = Star%xmg26      
        vxne20      = Star%vxne20     
        vxne22      = Star%vxne22     
        vxmg24      = Star%vxmg24     
        vxmg25      = Star%vxmg25     
        vxmg26      = Star%vxmg26     
        omegi       = Star%omegi      
        vomegi      = Star%vomegi     
        xf19        = Star%xf19       
        xne21       = Star%xne21      
        xna23       = Star%xna23      
        xal26       = Star%xal26      
        xal27       = Star%xal27      
        xsi28       = Star%xsi28      
        vxf19       = Star%vxf19      
        vxne21      = Star%vxne21     
        vxna23      = Star%vxna23     
        vxal26g     = Star%vxal26g    
        vxal27      = Star%vxal27     
        vxsi28      = Star%vxsi28     
        xneut       = Star%xneut      
        xprot       = Star%xprot      
        xc14        = Star%xc14       
        xf18        = Star%xf18       
        xbid        = Star%xbid       
        xbid1       = Star%xbid1      
        vxneut      = Star%vxneut     
        vxprot      = Star%vxprot     
        vxc14       = Star%vxc14      
        vxf18       = Star%vxf18      
        vxbid       = Star%vxbid      
        vxbid1      = Star%vxbid1     
        drl         = Star%drl        
        drte        = Star%drte       
        dk          = Star%dk         
        drp         = Star%drp        
        drt         = Star%drt        
        drr         = Star%drr        
        rlp         = Star%rlp        
        rlt         = Star%rlt        
        rlc         = Star%rlc        
        rrp         = Star%rrp        
        rrt         = Star%rrt        
        rrc         = Star%rrc        
        rtp         = Star%rtp        
        rtt         = Star%rtt        
        rtc         = Star%rtc        
        tdiff       = Star%tdiff      
        suminenv    = Star%suminenv   
        CorrOmega   = Star%CorrOmega  
        xLtotbeg    = Star%xLtotbeg   
        dlelexprev  = Star%dlelexprev 
        zams_radius = Star%zams_radius
        !xnetalu     = Star%xnetalu    
        !xlostneu    = Star%xlostneu   
        !nbzel       = Star%nbzel      
        !nbael       = Star%nbael      
        !abels       = Star%abels      
        
        radius      = Star%radius

    end subroutine copy_from_genec_star

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
                InitialGenecStar%ipoly = 0
        end select

        if (InitialGenecStar%ipoly == 1) then
            write(*,*)'Enter the polytropic index (recommended: 2.5):'
            read(5,*) InitialGenecStar%n
        endif
        
    end subroutine input_ini

end module helpers
