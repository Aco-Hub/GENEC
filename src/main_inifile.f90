program main
        use storage, only:GenecStar
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
        read(5,*) &
                GenecStar%star_name
        write(*,*)'Enter the desired mass and metallicity:'
        read(5,*) &
                GenecStar%initial_mass,&
                GenecStar%initial_metallicity
        write(*,*) 'Which rotation velocity on the ZAMS?'
        read(5,*) &
                GenecStar%zams_velocity

        !!!!! inichemmod bit
        answer = ''

        do while (answer /= 'y' .and. answer /= 'n' .and. answer /= '1' .and. answer /= '0')
            write(*,*)'Default settings are:'
            write(*,*)'         - Asplund-Cunha abundances'
            write(*,*)'         - scaled solar'
            write(*,*)'         - Geneva format'
            write(*,*)'         - small network'
            write(*,*)'         - structure from pre-calculated model'
            write(*,*)'Is it ok? (y)es (n)o'
            read(5,*) answer
            if (answer /= 'y' .and. answer /= 'n'  .and. answer /= '1' .and. answer /= '0') then
                write(*,*) 'Please type y or n...'
            endif
        enddo
        if (answer == 'y' .or. answer == 'Y' .or. answer == '1') then
            GenecStar%idefaut = 1
        elseif (answer == 'n' .or. answer == 'N' .or. answer == '0') then
            GenecStar%idefaut = 0
        endif

        ! source identifier
        sourceid=(/'AG89','GN93','As05','As09'/)

        if (GenecStar%idefaut == 0) then
        ! choose input file
            write(*,*) 'Choose source'
            write(*,*)'Anders & Grevesse 1989: 1 / Grevesse and Noels 1993: 2 / Asplund 2005: 3 / Asplund 2009: 4'
            read(5,*) &
                    GenecStar%source
            write(*,*) sourceid(int(GenecStar%source))
        ! choose alpha enhanced or solar scaled composition
            write(*,*)'scaled solar abundances (0) or alpha-enhanced (1)?'
            read(5,*) &
                    GenecStar%alpha
            write(*,*) 'alpha-enhanced:', GenecStar%alpha
            if(GenecStar%alpha/=0 .and. GenecStar%alpha/=1) then
                write(*,*) &
                        GenecStar
                stop 'wrong choice - program stopped!'
            endif
        ! choose output format
            write(*,*)'choose format! [1=basnet/2=GENEC/3=PPN input format]'
            read(5,*) &
                    GenecStar%formatx
            if(GenecStar%formatx/=1.and.GenecStar%formatx/=2.and.GenecStar%formatx/=3) then
                stop 'wrong format choice'
            endif
            write(*,*) 'choose network size:'
            write(*,*) '0: default GENEC - 1: 23 network - 2: 48 network - 3: free network'
            read(5,*) &
                     GenecStar%inetwork
            write(*,*) 'inetwork:',GenecStar%inetwork
        endif  ! idefaut

        !!!!! end inichemmod bit
        select case (GenecStar%idefaut)
            case (0)
                write(*,*)'Do you want to compute a polytropic structure or use ',&
                        'a pre-computed structure? Polytrope:1 - structure:0'
                read(5,*) &
                        GenecStar%ipoly
            case (1)
                GenecStar%ipoly = 0
        end select

        if (GenecStar%ipoly == 1) then
            write(*,*)'Enter the polytropic index (recommended: 2.5):'
            read(5,*) &
                    GenecStar%n
        endif

    end subroutine input_ini
end program main
