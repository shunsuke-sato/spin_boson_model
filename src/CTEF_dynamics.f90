!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine CTEF_dynamics
  use global_variables
  use CTEF_module
  implicit none
  real(8) :: Sz_av
  integer :: it
  real(8) :: X_Cint_av
  complex(8) :: zpsi_t(2,2),zhpsi_t(2,2)


  zpsi_CTEF(1,:) = 1d0; zpsi_CTEF(2,:) = 0d0

  Szt_t = 0d0

  Sz_av = abs(zpsi(1))**2 - abs(zpsi(2))**2
  Szt_t(0) = Sz_av

!Back propagation for the initial condition
  a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
  a_HO = a_HO/M_HO
  X_HO_old = X_HO + (-dt)*V_HO + 0.5d0*(-dt)**2*a_HO

  do it = 0,Nt-1

    call set_hamiltonian_spin(zHO_CTEF,zHO_dot_CTEF)
    zpsi_t = zpsi_CTEF
    call zhpsi_CTEF(zpsi_t,zhpsi_t)
    call set_hamiltonian_bath(zpsi_t,zhpsi_t,zHO_CTEF, zHO_dot_CTEF)

  end do


end subroutine CTEF_dynamics
