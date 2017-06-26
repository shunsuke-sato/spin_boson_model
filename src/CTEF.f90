!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine CTEF
  use global_variables
  implicit none
  integer :: itraj,it


  Szt_t=0d0; Szt_l = 0d0
  call setting_bath_parameters

  do itraj = 1,Ntraj

    call thermal_bath_distribution


    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj

    call CTEF_dynamics
    Szt_l = Szt_l + Szt_t

  end do
  call MPI_ALLREDUCE(Szt_l,Szt,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Szt = Szt/Ntraj

  if(myrank == 0)then
    open(nfile_CTEF_Sz,file=file_CTEF_Sz)
    do it = 0,Nt
      write(nfile_MTEF_Sz,"(999e26.16e3)")dt*it,Szt(it)
    end do
    close(nfile_MTEF_Sz)
  end if


end subroutine CTEF
