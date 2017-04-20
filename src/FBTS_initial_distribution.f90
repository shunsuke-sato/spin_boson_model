!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS_initial_distribution
  use global_variables
  implicit none
  real(8) :: xx,pp

  select case(calc_mode)
  case('FBTS')
    call gaussian_random_number(xx,pp)
    x_m(1) = sqrt(0.5d0)*xx; p_m(1) = sqrt(0.5d0)*pp
    call gaussian_random_number(xx,pp)
    x_m(2) = sqrt(0.5d0)*xx; p_m(2) = sqrt(0.5d0)*pp

    call gaussian_random_number(xx,pp)
    x_n(1) = sqrt(0.5d0)*xx; p_n(1) = sqrt(0.5d0)*pp
    call gaussian_random_number(xx,pp)
    x_n(2) = sqrt(0.5d0)*xx; p_n(2) = sqrt(0.5d0)*pp

  case('FBTS_mod')
    call correlated_gaussian_random_number(xx,pp)
    x_m(1) = xx; x_n(1) = pp
    call correlated_gaussian_random_number(xx,pp)
    x_m(2) = xx; x_n(2) = pp

    call correlated_gaussian_random_number(xx,pp)
    p_m(1) = xx; p_n(1) = pp
    call correlated_gaussian_random_number(xx,pp)
    p_m(2) = xx; p_n(2) = pp

  case('FBTS_approx')
    call gaussian_random_number(xx,pp)
    x_m(1) = sqrt(0.5d0)*xx; p_m(1) = sqrt(0.5d0)*pp
    call gaussian_random_number(xx,pp)
    x_m(2) = sqrt(0.5d0)*xx; p_m(2) = sqrt(0.5d0)*pp

    x_n = x_m; p_n = p_m

  case default
    call err_finalize('Invalid calc_mode')
  end select

end subroutine FBTS_initial_distribution
