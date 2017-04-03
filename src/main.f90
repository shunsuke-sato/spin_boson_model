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
  case default
    call err_finalize('Invalid calc_mode')
  end select

  call MPI_finalize(ierr)

end program main
