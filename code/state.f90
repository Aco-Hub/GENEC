module State
  use inputparam, only: amuseinterface
  implicit none
  character(99), save:: stopping_condition=""

  contains

  subroutine conditioned_stop()
    implicit none
    if (amuseinterface) then
      write(*,*) "Stopping: ", stopping_condition
    else
      stop stopping_condition
    endif
  end subroutine
end module
