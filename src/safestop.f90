module safestop
    use evol, only: libgenec
    use storage, only: GenecStar
    contains

    subroutine safe_stop(message)
        implicit none
        character(*), intent(in) :: message
        if (libgenec) then
            GenecStar%stopped = .true.
            GenecStar%stop_message = message
        else
            stop message
        endif
    end subroutine

end module safestop
