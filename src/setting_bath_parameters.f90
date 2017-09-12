!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine setting_bath_parameters
  use global_variables
  implicit none
  real(8) :: w0,r,ss
  integer :: j
  real(8) :: f1

  w0 = Omega_c*(1d0-exp(-Omega_max/Omega_c))/dble(Num_HO)
  
  do j = 1,Num_HO
    Omega_HO(j) = -Omega_c*log(1d0-dble(j)*w0/Omega_c)
    Cint_HO(j) = Omega_HO(j)
  end do

  ss = 0d0
  do j = 1,Num_HO
    ss = ss + Cint_HO(j)**2/Omega_HO(j)**2
  end do

  r = Omega_max/Omega_c
  f1 = Xi_Kondo*omega_c*(1d0-exp(-r))

  M_HO = ss/f1

  if(myrank == 0)then
    write(*,"(A,2x,e26.16e3)")"Bath GS energy:",0.5d0*sum(omega_ho)
  end if

end subroutine setting_bath_parameters
