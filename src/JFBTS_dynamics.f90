!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine JFBTS_dynamics(zweight0)
  use global_variables
  implicit none
  complex(8),intent(inout) :: zweight0
  real(8) :: Sz_av
  integer :: it
  real(8) :: X_Cint_av
  complex(8) :: zs
  logical :: if_jump

  Szt_t = 0d0

  z_m = x_m + zI*p_m; z_n = x_n + zI*p_n
  Sz_av = 0.5d0*(x_m(1)**2 + p_m(1)**2 + x_n(1)**2 + p_n(1)**2 &
    -x_m(2)**2 - p_m(2)**2 - x_n(2)**2 - p_n(2)**2 )

  zSzt_t(0) = conjg(z_m(1))*z_n(1) - conjg(z_m(2))*z_n(2)
  zSzt_t(0) = zSzt_t(0)*zweight0

!!Back propagation for the initial condition
!  a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
!  a_HO = a_HO/M_HO
!  X_HO_old = X_HO + (-dt)*V_HO + 0.5d0*(-dt)**2*a_HO

  do it = 0,Nt-1
! Propagate X_HO
    a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
    a_HO = a_HO/M_HO
    V_HO = V_HO + a_HO*0.5d0*dt

! Propagate spin (Enforced time-reversal symmetry scheme)
    X_Cint_av = sum(X_HO*Cint_HO)
    H_spin = eps_SP*Sz + delta_SP*Sx -X_Cint_av*Sz

    z_m = x_m + zI*p_m; z_n = x_n + zI*p_n
    call spin_propagation(z_m,H_spin,dt*0.5d0)
    call spin_propagation(z_n,H_spin,dt*0.5d0)

    X_HO = X_HO + V_HO*dt
    X_Cint_av = sum(X_HO*Cint_HO)
    H_spin = eps_SP*Sz + delta_SP*Sx -X_Cint_av*Sz

    call spin_propagation(z_m,H_spin,dt*0.5d0)
    call spin_propagation(z_n,H_spin,dt*0.5d0)
    x_m = real(z_m); p_m = aimag(z_m)
    x_n = real(z_n); p_n = aimag(z_n)

    Sz_av = 0.5d0*(x_m(1)**2 + p_m(1)**2 + x_n(1)**2 + p_n(1)**2 &
      -x_m(2)**2 - p_m(2)**2 - x_n(2)**2 - p_n(2)**2 )

    a_HO = -M_HO*Omega_HO**2*X_HO + Cint_HO*Sz_av
    a_HO = a_HO/M_HO
    V_HO = V_HO + a_HO*0.5d0*dt

    zSzt_t(it+1) = conjg(z_m(1))*z_n(1) - conjg(z_m(2))*z_n(2)
    zSzt_t(it+1) = zSzt_t(it+1) *zweight0

    if_jump = .false.
    if(it /= 0 .and. mod(it,Mjump_FBTS) ==0 ) if_jump = .true.

    if(if_jump)then

      z_m = x_m + zI*p_m; z_n = x_n + zI*p_n
      call FBTS_initial_distribution
      zs = conjg(z_m(1))*(x_m(1) + zI* p_m(1)) + conjg(z_m(2))*(x_m(2) + zI*p_m(2))
      zs = zs* ( &
        z_n(1)*conjg(x_n(1) + zI*p_n(1)) + z_n(2)*conjg(x_n(2) + zI*p_n(2)) )
      zweight0  = zweight0 * zs
    end if


  end do


end subroutine JFBTS_dynamics
