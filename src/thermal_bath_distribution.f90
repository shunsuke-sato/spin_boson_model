!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine thermal_bath_distribution
  use global_variables
  implicit none
  real(8) :: alpha,x1,x2
  integer :: j

  do j = 1,Num_HO
    alpha = 2d0*tanh(beta_kB*Omega_HO(j)/2d0)/Omega_HO(j)
    call gaussian_random_number(x1,x2)
    X_HO(j) = x1/sqrt(alpha*M_HO*Omega_HO(j)**2)
    V_HO(j) = x2/sqrt(M_HO*alpha)
  end do

end subroutine thermal_bath_distribution
