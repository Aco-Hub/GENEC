program main
        use storage, only:GenecStar
        use makeini

        implicit none
        call input_ini()
        call make_initial_star()

contains
    subroutine input_ini
        implicit none
        integer:: Input_change, Temp_Var_Int
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

! initialisation to default values
        GenecStar%source=3
        GenecStar%alpha=0
        GenecStar%formatx=2
        GenecStar%inetwork = 0
        GenecStar%ipoly = 0

        ! source identifier
        sourceid=(/'AG89','GN93','As05','As09'/)

        Input_change = 99
        do while (Input_change /= 0)
          write(*,*) '|------------------------------------------------------|'
          write(*,*) '| Default parameters you can change:                   |'
          write(*,*) '|------------------------------------------------------|'
          write(*,*) '|  1: Abundances (default Asplund+2005)                |'
          write(*,*) '|  2: Scaled solar (default) or alpha-enhanced         |'
          write(*,*) '|  3: Chemical network (default inetwork=0)            |'
          write(*,*) '|  4: Initial structure (default pre-calculated model) |'
          write(*,*) '|  5: File format (default Geneva)                     |'
          write(*,*) '|------------------------------------------------------|'
          write(*,*) 'Enter the parameter number you want to change (0 to skip or exit):'
          read(5,*) Input_change
          select case(Input_change)
          case (0)
            write(*,*) 'No (more) changes...'
          case (1)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 4)
              write(*,*) 'Choose the abundance mixture:'
              write(*,*) '------------------------------'
              write(*,*) ' 1: Anders & Grevesse 1989'
              write(*,*) ' 2: Grevesse and Noels 1993'
              write(*,*) ' 3: Asplund 2005 (default)'
              write(*,*) ' 4: Asplund 2009'
              write(*,*) '------------------------------'
              write(*,*) 'Enter your choice:'
              read(5,*) Temp_Var_Int
            enddo
            if (Temp_Var_Int == 0) then
              GenecStar%source = 3
            else
              GenecStar%source = Temp_Var_Int
            endif
            write(*,*) 'Abundance mixture set to ',sourceid(int(GenecStar%source))
          case (2)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 1)
              write(*,*) 'Choose the alpha content:'
              write(*,*) '------------------------------'
              write(*,*) ' 0: solar-scaled (default)'
              write(*,*) ' 1: alpha-enhanced'
              write(*,*) '------------------------------'
              write(*,*) 'Enter your choice:'
              read(5,*) Temp_Var_Int
            enddo
            GenecStar%alpha = Temp_Var_Int
            write(*,*) 'alpha-enhanced:', GenecStar%alpha
          case (3)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 3)
              write(*,*) 'Choose network size:'
              write(*,*) '------------------------------'
              write(*,*) ' 0: default GENEC'
              write(*,*) ' 1: 23 network'
              write(*,*) ' 2: 48 network'
              write(*,*) ' 3: free network'
              write(*,*) '------------------------------'
              write(*,*) 'Enter your choice:'
              read(5,*) Temp_Var_Int
            enddo
            GenecStar%inetwork = Temp_Var_Int
            write(*,*) 'inetwork:',GenecStar%inetwork
          case (4)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 1)
              write(*,*)'Choose the type of initial structure:'
              write(*,*) '------------------------------'
              write(*,*) ' 0: Pre-calculated structure (default)'
              write(*,*) ' 1: Polytrope'
              write(*,*) '------------------------------'
              write(*,*) 'Enter your choice:'
              read(5,*) Temp_Var_Int
            enddo
            GenecStar%ipoly = Temp_Var_Int
            if (GenecStar%ipoly == 1) then
              write(*,*) 'Enter the polytropic index (recommended 2.5)'
              read(5,*) GenecStar%n_poly
            endif
          case (5)
            Temp_Var_Int = 99
            do while (Temp_Var_Int > 3)
              write(*,*)'Choose the file format:'
              write(*,*) '------------------------------'
              write(*,*) ' 1: basnet'
              write(*,*) ' 2: GENEC (default)'
              write(*,*) ' 3: PPN input format'
              write(*,*) '------------------------------'
              write(*,*) 'Enter your choice:'
              read(5,*) Temp_Var_Int
            enddo
            if (Temp_Var_Int == 0) then
              GenecStar%formatx = 2
            else
              GenecStar%formatx = Temp_Var_Int
            endif
          case default
            write(*,*) 'wrong number: should be between 0 (exit) and 5'
          end select
        enddo

    end subroutine input_ini
end program main
