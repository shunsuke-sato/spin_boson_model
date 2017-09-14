!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module global_variables
! mathamtical constant
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)

! control parameter
!  character(32) :: calc_mode = 'JFBTS' !'FBTS_mod,' 'FBTS' 'MTEF','PBME','PBME_mod'
  character(32) :: calc_mode = 'CTEF' !'FBTS_mod,' 'FBTS' 'MTEF','PBME','PBME_mod'

! Spin
  real(8),parameter :: delta_SP = 0.5d0, eps_SP = 0d0 !delta_SP
  complex(8) :: zpsi(2)
  real(8),parameter :: Sz(2,2) = reshape( (/1d0, 0d0, 0d0, -1d0/), (/2,2/) )
  real(8),parameter :: Sx(2,2) = reshape( (/0d0, 1d0, 1d0, 0d0/), (/2,2/) )
  real(8),parameter :: hs_m(2,2) = reshape( (/eps_SP, delta_SP, delta_SP, -eps_SP/), (/2,2/) )
  real(8) :: H_spin(2,2)

! Spin-Boson model parameters
  real(8),parameter :: Xi_Kondo = 0.1d0
  real(8),parameter :: Omega_c = 2d0*delta_SP, Omega_max = 5d0*Omega_c

!! PBME and FBTS
  real(8) :: x_m(2),p_m(2),x_n(2),p_n(2)
  complex(8) :: z_m(2), z_n(2)
  integer,parameter :: Mjump_FBTS = 100

! Harmonic oscillator
  integer,parameter :: Num_HO = 1 !400
  real(8) :: X_HO(Num_HO),V_HO(Num_HO),a_HO(Num_HO)
  real(8) :: X_HO_old(Num_HO),V_HO_old(Num_HO),a_HO_old(Num_HO)
  real(8) :: X_HO_new(Num_HO),V_HO_new(Num_HO),a_HO_new(Num_HO)
  real(8) :: Omega_HO(Num_HO),Cint_HO(Num_HO)
  real(8) :: M_HO


! Whole system
  real(8),parameter :: beta_kB = 1d10 !5d0/delta_SP !5d0/delta_SP

! Parameters for time-propagation
  real(8),parameter :: Tprop = 20d0,dt = 0.01d0
  integer,parameter :: Nt = aint(Tprop/dt)+1
! 
  integer,parameter :: Ntraj= 1000 !1000
  real(8) :: Szt(0:Nt),Szt_t(0:Nt),Szt_l(0:Nt)
  complex(8) :: zSzt_t(0:Nt)

! I/O parameters
! MTEF
  character(64) :: file_MTEF_Sz="MTEF_Sz.out"
  integer :: nfile_MTEF_Sz=41
! PBME
  character(64) :: file_PBME_Sz="PBME_Sz.out"
  integer :: nfile_PBME_Sz=42
! FBTS
  character(64) :: file_FBTS_Sz="FBTS_Sz.out"
  integer :: nfile_FBTS_Sz=43
! CTEF
  character(64) :: file_CTEF_Sz="CTEF_Sz.out"
  integer :: nfile_CTEF_Sz=44

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr

end module global_variables
