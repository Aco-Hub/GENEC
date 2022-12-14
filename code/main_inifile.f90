program main
        use helpers, only: input_ini
        use makeini

        implicit none
        call input_ini()
        call make_initial_star()
end program main
