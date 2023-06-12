program main
        use storage, only:InitialGenecStar
        use makeini

        implicit none
        call input_ini()
        call make_initial_star()

contains
    subroutine input_ini
        implicit none
        character(len=1):: answer
        character(len=4), dimension(4):: sourceid !< source identifier used in output file name

        
        write(*,*)'Enter the star name:'
        read(5,*) InitialGenecStar%star_name
        write(*,*)'Enter the desired mass and metallicity:'
        read(5,*) InitialGenecStar%initial_mass,InitialGenecStar%initial_metallicity
        write(*,*) 'Which rotation velocity on the ZAMS?'
        read(5,*) InitialGenecStar%zams_velocity

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
end program main
