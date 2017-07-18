!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME_initial_distribution
  use global_variables
  use random_number_module
  implicit none
  real(8) :: xx,pp

  call gaussian_random_number(xx,pp)
  x_m(1) = sqrt(0.5d0)*xx; p_m(1) = sqrt(0.5d0)*pp
  call gaussian_random_number(xx,pp)
  x_m(2) = sqrt(0.5d0)*xx; p_m(2) = sqrt(0.5d0)*pp


end subroutine PBME_initial_distribution
