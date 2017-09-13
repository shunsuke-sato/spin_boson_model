!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module exact_smsb
  implicit none
! mathamtical constant
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)

! model parameters
  real(8),parameter :: delta_SP = 1d0, eps_SP = delta_SP
  real(8),parameter :: Sz(2,2) = reshape( (/1d0, 0d0, 0d0, -1d0/), (/2,2/) )
  real(8),parameter :: Sx(2,2) = reshape( (/0d0, 1d0, 1d0, 0d0/), (/2,2/) )
  real(8),parameter :: hs(2,2) = reshape( (/eps_SP, delta_SP, delta_SP, -eps_SP/), (/2,2/) )
  real(8),parameter :: mass_p = 1d0
  real(8),parameter :: omega_p = 1d0
  real(8),parameter :: cint_p = 1d0
  real(8),parameter :: gamma_p = 1d0 ! gamma should be computed by cint and mass

! numerical parameter
  integer,parameter :: Nmax = 10
  real(8),parameter :: dt = 0.01d0
  real(8),parameter :: Tprop = 20d0
  integer,parameter :: Nt = aint(Tprop/dt)+1


! wavefunctions
  complex(8) :: zpsi(0:Nmax,2)
  complex(8) :: ztpsi(0:Nmax,2)
  complex(8) :: zhpsi(0:Nmax,2)




  contains

!-----------------------------------------------------------------------------------------
    subroutine dynamics
      integer :: it

      call init_wf

      do it = 0,Nt
        write(*,"(A,2x,I9,A,I9)")"it=",it,"/",Nt

        call dt_prop(zpsi)

      end do

    end subroutine dynamics
!-----------------------------------------------------------------------------------------
    subroutine init_wf
      
      zpsi = 0d0
      zpsi(0,1) = 1d0

    end subroutine init_wf
!-----------------------------------------------------------------------------------------
    subroutine dt_prop(zpsi_t)
      complex(8),intent(inout) :: zpsi_t(0:Nmax,2)
      complex(8) :: zpsi_s(0:Nmax,2),zhpsi_s(0:Nmax,2)
      integer,parameter :: nTaylor = 6
      complex(8) :: zfact
      integer :: iexp

      zpsi_s = zpsi_t
      zfact = 1d0
      
      do iexp = 1,nTaylor

        zfact = zfact*(-zI*dt)/iexp
        call hamiltonian_full(zpsi_s,zhpsi_s)
        zpsi_t = zpsi_t + zfact*zhpsi_s
        zpsi_s = zhpsi_s
        
      end do
      
    end subroutine dt_prop
!-----------------------------------------------------------------------------------------
    subroutine hamiltonian_full(zpsi_t,zhpsi_t)
      complex(8),intent(in) :: zpsi_t(0:Nmax,2)
      complex(8),intent(out) :: zhpsi_t(0:Nmax,2)
      integer :: ib

! spin
      zhpsi_t(:,1) = hs(1,1)*zpsi_t(:,1) + hs(1,2)*zpsi_t(:,2)
      zhpsi_t(:,2) = hs(2,1)*zpsi_t(:,1) + hs(2,2)*zpsi_t(:,2)

! bath
      do ib = 0,Nmax
        zhpsi_t(ib,1)  = zhpsi_t(ib,1) + omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,1)
        zhpsi_t(ib,2)  = zhpsi_t(ib,2) + omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,2)
      end do

! coupling
!! S_z
      zhpsi_t(0,1)  = zhpsi_t(0,1) -gamma_p*sqrt(1d0)*zpsi_t(1,1)
      zhpsi_t(0,2)  = zhpsi_t(0,2) +gamma_p*sqrt(1d0)*zpsi_t(1,2)
      do ib = 1,Nmax-1
        zhpsi_t(ib,1)  = zhpsi_t(ib,1) -gamma_p*( &
          sqrt(dble(ib))*zpsi_t(ib-1,1) &
         +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,1) )

        zhpsi_t(ib,2)  = zhpsi_t(ib,2) +gamma_p*( &
          sqrt(dble(ib))*zpsi_t(ib-1,2) &
         +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,2) )

      end do

      zhpsi_t(Nmax,1)  = zhpsi_t(Nmax,1) -gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,1)
      zhpsi_t(Nmax,2)  = zhpsi_t(Nmax,2) +gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,2)
      

    end subroutine hamiltonian_full

!-----------------------------------------------------------------------------------------
end module exact_smsb
!-----------------------------------------------------------------------------------------
program main
  use exact_smsb
  implicit none

  call dynamics

end program main
!-----------------------------------------------------------------------------------------
