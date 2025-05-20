module evol
  use ISO_FORTRAN_ENV
  implicit none

  integer, parameter:: kindreal=real64
  integer, parameter:: ldi=5000,mmax=4999,npondcouche=200,npondcoucheAdv=50
  integer, parameter:: mbelx=26

  character (len=256), save:: input_dir

  ! if libgenec is set to .true., no input will be asked.
  logical,save:: &
          libgenec=.false.
  public
end module evol
