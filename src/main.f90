!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
program main
  use global_variables
  implicit none

  call initialize_mpi

  select case(calc_mode)
  case('MTEF')
    call MTEF
  case('PBME','PBME_mod')
    call PBME
  case('FBTS','FBTS_mod','FBTS_approx','JFBTS')
    call FBTS
  case default
    call err_finalize('Invalid calc_mode in main')
  end select

  call MPI_finalize(ierr)

end program main
