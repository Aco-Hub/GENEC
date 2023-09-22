module stopping
  implicit none
  character(*):: continue_point="start"
  character(*):: stop_reason

  contains
    subroutine stop_with_reason(reason)
      use inputparam, only: libgenec
      character(*), intent(in):: reason
      stop_reason = reason
      continue_point = reason
      if (.not.libgenec) then
        stop
      endif
        

end module stopping
