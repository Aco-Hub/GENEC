program bench_fortran
  use interpolation
  use SmallFunc
  use ionisation
  use const
  use strucmod, only: vmion, vmol, beta_env, chem, ychem, x_env, cp, vna, vmionp, vmiont
  implicit none

  character(len=64) :: func_name, nstr
  integer :: n, i, j, k
  real(8) :: elapsed, r1
  real(8), allocatable :: log_p(:), log_t(:), results(:)
  integer :: count_start, count_end, count_rate

  call get_command_argument(1, func_name)
  call get_command_argument(2, nstr)
  read(nstr, *) n

  allocate(log_p(n), log_t(n), results(n))

  ! Generate realistic stellar inputs
  call random_seed()
  do i = 1, n
    call random_number(r1)
    log_t(i) = 3.7d0 + r1 * 3.45d0       ! log(T) in [3.7, 7.15]
    call random_number(r1)
    log_p(i) = 0.5d0 + r1 * 16.8d0        ! log(P) in [0.5, 17.3]
  end do

  call system_clock(count_start, count_rate)

  select case (trim(func_name))

  case ('ionpart')
    ! Setup composition (solar-like)
    abond(1) = 0.72d0    ! H
    abond(2) = 0.266d0   ! He
    abond(3) = 2.56d-3   ! C
    abond(4) = 6.42d-3   ! O
    abond(5) = 1.65d-3   ! Ne
    abond(6) = 5.13d-4   ! Mg

    ! Setup list and vnu
    do j = 1, iatoms
      list(j) = j
    end do
    vmol = 0.d0
    do j = 1, iatoms
      vmol = vmol + abond(j) / a_ion(j)
    end do
    vmol = vmol + 0.5d0 * (1.d0 - (abond(1)+abond(2)+abond(3)+abond(4)+abond(5)+abond(6)))
    vmol = 1.d0 / vmol
    do j = 1, iatoms
      vnu(j) = vmol * abond(j) / a_ion(j)
    end do
    ionized = 0
    vmion = 0.d0
    chem = abond(3) + abond(4) + abond(5) + abond(6)
    ychem = abond(2)

    do i = 1, n
      ionized = 0
      vmion = 0.d0
      call ionpart(log_p(i), log_t(i))
      results(i) = vmion
    end do

  case ('eos')
    ! Simplified ideal gas + radiation EOS (same as Python eos_ideal)
    block
      real(8) :: P, T, P_rad, P_gas, beta, rho, mu_val
      real(8), parameter :: cst_a_local = 7.56576738d-15
      real(8), parameter :: cst_k_local = 1.3806504d-16
      real(8), parameter :: cst_mh_local = 1.67372346d-24
      mu_val = 0.62d0  ! typical fully-ionized solar
      do i = 1, n
        P = 10.d0**log_p(i)
        T = 10.d0**log_t(i)
        P_rad = (cst_a_local / 3.d0) * T**4
        P_gas = max(P - P_rad, P * 1.d-10)
        beta = P_gas / P
        rho = P_gas * mu_val * cst_mh_local / (cst_k_local * T)
        results(i) = rho
      end do
    end block

  case ('energy')
    ! Simplified PP-chain + CNO (same formulas as Python)
    block
      real(8) :: T9, T913, eps_pp, eps_cno, X_H, X_CNO, rho_val
      X_H = 0.72d0
      X_CNO = 3.3d-3
      do i = 1, n
        T9 = 10.d0**log_t(i) / 1.d9
        rho_val = 10.d0**(log_p(i) * 0.1d0)  ! approximate density
        if (10.d0**log_t(i) < 4.d6) then
          results(i) = 0.d0
          cycle
        endif
        T913 = T9**(1.d0/3.d0)
        eps_pp = 2.38d6 * rho_val * X_H**2 * T9**(-2.d0/3.d0) * exp(-3.381d0 / T913)
        eps_cno = 8.24d25 * rho_val * X_H * X_CNO * T9**(-2.d0/3.d0) * exp(-15.231d0 / T913)
        results(i) = eps_pp + eps_cno
      end do
    end block

  case default
    write(*,*) 'Unknown function: ', trim(func_name)
    stop 1

  end select

  call system_clock(count_end, count_rate)
  elapsed = dble(count_end - count_start) / dble(count_rate)

  write(*,'(ES20.10)') elapsed

  deallocate(log_p, log_t, results)

end program bench_fortran
