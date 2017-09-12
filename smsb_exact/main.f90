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
  real(8),parameter :: mass = 1d0
  real(8),parameter :: omega = 1d0
  real(8),parameter :: cint = 1d0
  real(8),parameter :: gamma = 1d0 ! gamma should be computed by cint and mass

! numerical parameter
  integer,parameter :: Nmax = 10
  complex(8) :: zpsi(0:Nmax,2)
  complex(8) :: ztpsi(0:Nmax,2)
  complex(8) :: zhpsi(0:Nmax,2)


  contains

end module exact_smsb
!-----------------------------------------------------------------------------------------
program main
  use exact_smsb
  implicit none

  call initialization
  call dynamics
  call 


end program main
!-----------------------------------------------------------------------------------------
