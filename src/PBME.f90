!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine PBME
  use global_variables
  implicit none
  integer :: itraj,it
  real(8) :: weight0

  Szt_t=0d0; Szt_l = 0d0
  call setting_bath_parameters

  do itraj = 1,Ntraj

    call thermal_bath_distribution
    call PBME_initial_distribution

    if(mod(itraj,Nprocs) /= myrank)cycle
    if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj
    select case(calc_mode)
    case('PBME')
      weight0 = 2d0*(x_m(1)**2+p_m(1)**2-0.5d0)
    case('PBME_mod')
      weight0 = x_m(1)**2 + p_m(1)**2 + x_m(2)**2 + p_m(2)**2
      weight0 = 2d0*(x_m(1)**2+p_m(1)**2-0.5d0)*2d0**4*exp(-weight0)
    case default
      call err_finalize('Invalid calc_mode')
    end select


    call PBME_dynamics
    Szt_l = Szt_l + Szt_t*weight0

  end do
  call MPI_ALLREDUCE(Szt_l,Szt,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
  Szt = Szt/Ntraj

  if(myrank == 0)then
    open(nfile_PBME_Sz,file=file_PBME_Sz)
    do it = 0,Nt
      write(nfile_PBME_Sz,"(999e26.16e3)")dt*it,Szt(it)
    end do
    close(nfile_PBME_Sz)
  end if


end subroutine PBME
