!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine FBTS
  use global_variables
  implicit none
  integer :: itraj,it
  complex(8) :: zweight0

  Szt_t=0d0; Szt_l = 0d0
  call setting_bath_parameters

  do itraj = 1,Ntraj

    call thermal_bath_distribution
    call FBTS_initial_distribution

    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj
    zweight0 = (x_m(1) + zI * p_m(1))*(x_n(1) - zI * p_n(1))

    select case(calc_mode)
    case('FBTS','FBTS_mod','FBTS_approx')
      call FBTS_dynamics
      Szt_l = Szt_l + zSzt_t*zweight0
    case('JFBTS')
      call JFBTS_dynamics(zweight0)
      Szt_l = Szt_l + zSzt_t
    case default
      call err_finalize('Invalid calc_mode in FBTS')
    end select



  end do
  call MPI_ALLREDUCE(Szt_l,Szt,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Szt = Szt/Ntraj

  if(myrank == 0)then
    open(nfile_FBTS_Sz,file=file_FBTS_Sz)
    do it = 0,Nt
      write(nfile_FBTS_Sz,"(999e26.16e3)")dt*it,Szt(it)
    end do
    close(nfile_FBTS_Sz)
  end if


end subroutine FBTS
