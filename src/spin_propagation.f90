!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine spin_propagation(zpsi,H_spin,dt)
  implicit none
  complex(8),parameter :: zI = (0d0,1d0)
  complex(8),intent(inout) :: zpsi(2)
  real(8),intent(in) :: H_spin(2,2)
  real(8),intent(in) :: dt
  real(8) :: vec(2,2),lambda(2),alpha,delta
  complex(8) :: zC(2)
  real(8) :: ss

  alpha = H_spin(1,1)
  delta = H_spin(2,1)
! Assuming the following relations
! H_spin(1,1) = - H_spin(2,2)
! H_spin(2,1) = H_spin(1,2)

! eigen states of H_spin
  lambda(1) = sqrt(alpha**2+delta**2)
  lambda(2) = -sqrt(alpha**2+delta**2)
  vec(1,1) = 1d0; vec(2,1) = (lambda(1)-alpha)/delta
  vec(1,2) = 1d0; vec(2,2) = (lambda(2)-alpha)/delta
  ss = sum(vec(:,1)**2); vec(:,1) = vec(:,1)/sqrt(ss)
  ss = sum(vec(:,2)**2); vec(:,2) = vec(:,2)/sqrt(ss)

! projection to eigen states of H_spin
  zC(1) = sum(vec(:,1)*zpsi(:))
  zC(2) = sum(vec(:,2)*zpsi(:))
  
! propagation of eigenstates
  zC(1) = zC(1)*exp(-zI*dt*lambda(1))
  zC(2) = zC(2)*exp(-zI*dt*lambda(2))

! reconstraction of zpsi
  zpsi(:) = zC(1)*vec(:,1) + zC(2)*vec(:,2)



end subroutine spin_propagation
