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

  if(calc_mode == 'JFBTS')then
     call random_seeds_for_parallel(Myrank,Nprocs)
     do it = 1,10000
        call FBTS_initial_distribution
     end do
  end if

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

subroutine random_seeds_for_parallel(myrank,Nprocs)
  implicit none
  integer,intent(in) :: myrank,nprocs
  integer i, seedsize
  integer,allocatable:: seed(:)

  call random_seed(size=seedsize) 
  allocate(seed(seedsize))  
  call random_seed(get=seed) 

  do i = 1,seedsize
     seed(i) = seed(i) -1**myrank* (myrank*(seed(i)/Nprocs)+myrank)
  end do

end subroutine random_seeds_for_parallel
