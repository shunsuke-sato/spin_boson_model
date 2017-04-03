!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
Subroutine err_finalize(err_message)
  use Global_Variables
  implicit none
  character(*),intent(in) :: err_message

  if(myrank == 0)write(*,*) err_message
  call MPI_finalize(ierr)

  stop
end Subroutine err_finalize
