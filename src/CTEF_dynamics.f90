!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine CTEF_dynamics
  use global_variables
  implicit none
  real(8) :: Sz_av
  integer :: it
  real(8) :: X_Cint_av


  zpsi_CTEF(1,:) = 1d0; zpsi_CTEF(2,:) = 0d0

  Szt_t = 0d0

  Sz_av = abs(zpsi(1))**2 - abs(zpsi(2))**2
  Szt_t(0) = Sz_av

!Back propagation for the initial condition
  a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
  a_HO = a_HO/M_HO
  X_HO_old = X_HO + (-dt)*V_HO + 0.5d0*(-dt)**2*a_HO

  do it = 0,Nt-1
! Propagate X_HO
    a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
    a_HO = a_HO/M_HO
    X_HO_new = 2d0*X_HO -X_HO_old  + dt**2*a_HO

! Propagate spin (Enforced time-reversal symmetry scheme)
    X_Cint_av = sum(X_HO*Cint_HO)
    H_spin = eps_SP*Sz + delta_SP*Sx -X_Cint_av*Sz
    call spin_propagation(zpsi,H_spin,dt*0.5d0)
    X_Cint_av = sum(X_HO_new*Cint_HO)
    H_spin = eps_SP*Sz + delta_SP*Sx -X_Cint_av*Sz
    call spin_propagation(zpsi,H_spin,dt*0.5d0)

    X_HO_old = X_HO; X_HO = X_HO_new

    Sz_av = abs(zpsi(1))**2 - abs(zpsi(2))**2
    Szt_t(it+1) = Sz_av    

  end do


end subroutine CTEF_dynamics
