!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module CTEF_module
  use math_module
  use global_variables
  implicit none

  private
! phase average
  integer,parameter :: Nphi = 2
! spin
  complex(8) :: zpsi_CTEF(2,2)
! Harmonic oscillator 
  complex(8) :: ZHO_CTEF(2,Num_HO)
  complex(8) :: ZHO_dot_CTEF(2,Num_HO)
!
  complex(8) :: zHs_CTEF(2,2,2,2)
  complex(8) :: zDs_CTEF(2,2), zSs_CTEF(2,2)
  complex(8) :: zSs_inv_CTEF(2,2)

  complex(8) :: zDb_CTEF(2,2), zSb_CTEF(2,2)
  complex(8) :: zSb_inv_CTEF(2,2)
  complex(8) :: zHb_CTEF(2,2),zFb_CTEF(2)

  public :: CTEF

  contains
!-----------------------------------------------------------------------------------------
    subroutine CTEF
      implicit none
      integer :: itraj,it,iphi
      real(8) :: norm
      real(8) :: phi0,phi


      Szt_t=0d0; Szt_l = 0d0
      call setting_bath_parameters
      
      do itraj = 1,Ntraj
        
! set_forward_backward_trajectories
        call random_number(phi0); phi0 = 2d0*pi*phi0

        if(mod(itraj,Nprocs) /= myrank)cycle
        if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj
        
        do iphi = 0,Nphi-1
          phi = phi0 + 2d0*pi*dble(iphi)/Nphi
! set_initial_condition(phi,norm)

          call CTEF_dynamics
          Szt_l = Szt_l + Szt_t*norm*exp(-zI*phi)

        end do

      end do
      call MPI_ALLREDUCE(Szt_l,Szt,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      Szt = Szt/Ntraj/Nphi

      if(myrank == 0)then
        open(nfile_CTEF_Sz,file=file_CTEF_Sz)
        do it = 0,Nt
          write(nfile_MTEF_Sz,"(999e26.16e3)")dt*it,Szt(it)
        end do
        close(nfile_MTEF_Sz)
      end if


    end subroutine CTEF
!-----------------------------------------------------------------------------------------
    subroutine CTEF_dynamics
! Assuming zpsi_CTEF and zHO_CTEF are given 
      implicit none
      real(8) :: Sz_av
      integer :: it
      real(8) :: X_Cint_av
      complex(8) :: zpsi_t(2,2),zhpsi_t(2,2)

      
      Szt_t = 0d0 ! zero clear

      call calc_Szt(Szt_t(0),zpsi_CTEF, zHO_CTEF)
      zHO_dot_CTEF = 0d0
      call refine_effective_hamiltonian(zpsi_t,zHO_CTEF, zHO_dot_CTEF)

      do it = 0,Nt-1


        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF, zHO_dot_CTEF)
        call calc_Szt(Szt_t(it+1),zpsi_CTEF, zHO_CTEF)

        
      end do


    end subroutine CTEF_dynamics
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_etrs(zpsi_t,zHO_t, zHO_dot_t)
      implicit none
      complex(8),intent(inout) :: zpsi_t(2,2),zHO_t(2,Num_HO),zHO_dot_t(2,Num_HO)
      complex(8) :: zpsi_s(2,2),zHO_s(2,Num_HO),zHO_dot_s(2,Num_HO)
      complex(8) :: zpsi_tmp(2,2),zhpsi_tmp(2,2)
      integer,parameter :: Npred_corr = 2
      integer :: ipred_corr

! propagation: t => t + dt/2
      call dt_evolve_spin_Taylor(zpsi_t,0.5d0*dt)
      call dt_evolve_bath(zHO_t,0.5d0*dt)
      call refine_effective_hamiltonian(zpsi_t,zHO_CTEF, zHO_dot_CTEF)
      zpsi_s = zpsi_t
      zHO_s = zHO_t

      do ipred_corr = 1, Npred_corr
        zpsi_t = zpsi_s
        zHO_t = zHO_s

        call dt_evolve_spin_Taylor(zpsi_t,0.5d0*dt)
        call dt_evolve_bath(zHO_t,0.5d0*dt)
        if(ipred_corr /= Npred_corr)then
          call refine_effective_hamiltonian(zpsi_t,zHO_CTEF, zHO_dot_CTEF)
        end if

      end do

    end subroutine dt_evolve_etrs
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_spin_Taylor(zpsi_t,dt_t)
      implicit none
      complex(8),intent(inout) :: zpsi_t(2,2)
      real(8),intent(in) :: dt_t
      integer,parameter :: Ntaylor_t = 4
      integer :: itaylor
      complex(8) :: zpsi_s(2,2),zhpsi_s(2,2)
      complex(8) :: zfact

      zpsi_s = zpsi_t
      zfact = 1d0
      do itaylor = 1,Ntaylor_t
        zfact = zfact*(-zI*dt_t)/itaylor
        call zhpsi_CTEF(zpsi_s,zhpsi_s)
        zpsi_t = zpsi_t + zfact*zhpsi_s
        zpsi_s = zhpsi_s
      end do

    end subroutine dt_evolve_spin_Taylor
!-----------------------------------------------------------------------------------------
    subroutine dt_evolve_bath(zHO_t,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_t(2,Num_HO)
      real(8),intent(in) :: dt_t
      integer :: iho
      complex(8) :: zhHO_t(2,Num_HO),zhhHO_t(2,Num_HO)
      complex(8) :: zHeff(2,2), zFeff(2)

      zHeff = matmul(zSb_inv_CTEF, (zHB_CTEF - zDb_CTEF))
      zFeff = matmul(zSb_inv_CTEF, zFB_CTEF)

      do iho = 1,Num_HO
        zHO_t(:,iho) = zHO_t(:,iho) +(0.5d0*dt_t/zI)*zFeff(:) &
          *(sqrt(0.5d0/(M_HO*omega_HO(iho))))
      end do

      zhHO_t = matmul(zHeff,zHO_t)
      zhhHO_t = matmul(zHeff,zhHO_t)
      
      zHO_t = zHO_t + (-zI*dt_t)*zhHO_t + (-zI*dt_t)**2/2*zhhHO_t

      do iho = 1,Num_HO
        zHO_t(:,iho) = zHO_t(:,iho) +(0.5d0*dt_t/zI)*zFeff(:) &
          *(sqrt(0.5d0/(M_HO*omega_HO(iho))))
      end do


    end subroutine dt_evolve_bath

!-----------------------------------------------------------------------------------------
    subroutine zhpsi_CTEF(zpsi_t,zhpsi_t)
      implicit none
      complex(8),intent(in) :: zpsi_t(2,2)
      complex(8),intent(out) :: zhpsi_t(2,2)
      complex(8) :: zhpsi_s(2,2)
      real(8) :: det

      zhpsi_s(:,1) = zDs_CTEF(1,1)*zpsi_t(:,1) + zDs_CTEF(1,2)*zpsi_t(:,2)
      zhpsi_s(:,2) = zDs_CTEF(2,1)*zpsi_t(:,1) + zDs_CTEF(2,2)*zpsi_t(:,2)

      zhpsi_s(1,1) = zhpsi_s(1,1) + zHs_CTEF(1,1,1,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(1,2,1,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(1,1,1,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(1,2,1,2)*zpsi_t(2,2) 

      zhpsi_s(2,1) = zhpsi_s(2,1) + zHs_CTEF(2,1,1,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(2,2,1,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(2,1,1,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(2,2,1,2)*zpsi_t(2,2) 

      zhpsi_s(1,2) = zhpsi_s(1,2) + zHs_CTEF(1,1,2,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(1,2,2,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(1,1,2,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(1,2,2,2)*zpsi_t(2,2) 

      zhpsi_s(2,2) = zhpsi_s(2,2) + zHs_CTEF(2,1,2,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(2,2,2,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(2,1,2,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(2,2,2,2)*zpsi_t(2,2) 


      zhpsi_t(:,1) = zSs_inv_CTEF(1,1)*zhpsi_s(:,1) + zSs_inv_CTEF(1,2)*zhpsi_s(:,2)
      zhpsi_t(:,2) = zSs_inv_CTEF(2,1)*zhpsi_s(:,1) + zSs_inv_CTEF(2,2)*zhpsi_s(:,2)


    end subroutine zhpsi_CTEF
!-----------------------------------------------------------------------------------------
    subroutine refine_effective_hamiltonian(zpsi_t,zHO_t,zHO_dot_t)
      implicit none
      complex(8),intent(in) :: zpsi_t(2,2), zHO_t(2,Num_HO)
      complex(8),intent(inout) :: zHO_dot_t(2,Num_HO) 
      integer,parameter :: Nscf_CTEF = 2
      complex(8) :: zs_ab(2,2) ! overlap matrix for spin
      complex(8) :: zs_cd(2,2) ! overlap matrix for bath
      complex(8) :: zd_ab(2,2) ! time-derivative-overlap matrix for bath
      complex(8) :: zd_cd(2,2) ! time-derivative-overlap matrix for bath
      complex(8) :: zHb_tot_t(2,2)
      complex(8) :: zX_Cint_av(2,2), zEb(2,2)
      complex(8) :: zSz_ab(2,2)
      complex(8) :: zs,zvec(2)
      complex(8) :: zhpsi_t(2,2)
      integer :: iho
      integer :: iter
      integer :: i,j

      zs_ab(1,1) = 1d0
      zs_ab(2,2) = 1d0
      zs_ab(1,2) = sum(conjg(zpsi_t(:,1))*zpsi_t(:,2))
      zs_ab(2,1) = conjg(zs_ab(1,2))

      zs_cd(1,1) = 1d0
      zs_cd(2,2) = 1d0
      zs_cd(1,2) = sum(exp(&
            -0.5d0*abs(zHO_t(1,:))**2 -0.5d0*abs(zHO_t(2,:))**2  &
            + conjg(zHO_t(1,:))*zHO_t(2,:) &
            ))
      zs_cd(2,1) = conjg(zs_cd(1,2))

! partially prepare spin-Hamiltonian
      zSs_CTEF(:,:) = zs_ab(:,:)
      call inverse_2x2_matrix(zSs_CTEF,zSs_inv_CTEF)

      do i = 1,2
        do j =1,2
          zX_Cint_av(i,j) = sum( &
            Cint_HO(:) &
            *sqrt(0.5d0/(M_HO*Omega_HO(:)))*(conjg(zHO_t(i,:))+zHO_t(j,:)) )
          zX_Cint_av(i,j) = zX_Cint_av(i,j) *zSs_CTEF(i,j) 
          zEb(i,j) = sum(conjg(zHO_t(i,:))*zHO_t(j,:))+dble(Num_HO)/2d0
        end do
      end do
      zX_Cint_av = zX_Cint_av * zs_cd
      zEb = zEb * zs_cd
      

      do i = 1,2
        do j =1,2
          zHs_CTEF(:,:,i,j) = (eps_SP*Sz(:,:) + delta_SP*Sx(:,:))*zSs_CTEF(i,j) &
            -zX_Cint_av(i,j)*Sz(:,:) + zEb(i,j)
        end do
      end do

! partially prepare bath-Hamiltonian
      zSb_CTEF = zs_ab*zs_cd
      call inverse_2x2_matrix(zSb_CTEF,zSb_inv_CTEF)

      zvec(:) = matmul(zHs_CTEF(:,:,1,2),zpsi_t(:,2))
      zs = sum(conjg(zpsi_t(:,1))*zvec(:))
      zHb_CTEF(1,1) = -real(zs)
      zHb_CTEF(1,2) = zs
      zvec(:) = matmul(zHs_CTEF(:,:,2,1),zpsi_t(:,1))
      zs = sum(conjg(zpsi_t(:,2))*zvec(:))
      zHb_CTEF(2,2) = -real(zs)
      zHb_CTEF(2,1) = zs

      zvec(:) = matmul(Sz(:,:),zpsi_t(:,1))
      zSz_ab(1,1) = sum(conjg(zpsi_t(:,1))*zvec(:))
      zSz_ab(2,1) = sum(conjg(zpsi_t(:,2))*zvec(:))
      zSz_ab(1,2) = conjg(zSz_ab(2,1) )
      zvec(:) = matmul(Sz(:,:),zpsi_t(:,2))
      zSz_ab(2,2) = sum(conjg(zpsi_t(:,2))*zvec(:))

      zFb_CTEF(1) = zSz_ab(1,1) + zSz_ab(1,2)*zs_cd(1,2)
      zFb_CTEF(2) = zSz_ab(2,2) + zSz_ab(2,1)*zs_cd(2,1)

      do iter = 1,Nscf_CTEF


        zd_cd(1,1) = 0.5d0*sum(conjg(zHO_t(1,:))*zHO_dot_t(1,:) &
          - conjg(zHO_dot_t(1,:))*zHO_t(1,:))

        zd_cd(2,1) = sum(-0.5d0*(conjg(zHO_t(1,:))*zHO_dot_t(1,:) &
                               + conjg(zHO_dot_t(1,:))*zHO_t(1,:) ) &
                               +conjg(zHO_t(2,:))*zHO_dot_t(1,:) )

        zd_cd(1,2) = sum(-0.5d0*(conjg(zHO_t(2,:))*zHO_dot_t(2,:) &
                               + conjg(zHO_dot_t(2,:))*zHO_t(2,:) ) &
                               +conjg(zHO_t(1,:))*zHO_dot_t(2,:) )

        zd_cd(2,2) = 0.5d0*sum(conjg(zHO_t(2,:))*zHO_dot_t(2,:) &
          - conjg(zHO_dot_t(2,:))*zHO_t(2,:))

        zd_cd = zd_cd * zs_cd

! partially prepare spin-Hamiltonian
        zDs_CTEF(1,1) = real(zI*zd_cd(1,1))
        zDs_CTEF(1,2) = zI*zd_cd(1,2)
        zDs_CTEF(2,1) = zI*zd_cd(2,1)
        zDs_CTEF(2,2) = real(zI*zd_cd(2,2))

        call zhpsi_CTEF(zpsi_t,zhpsi_t)

        zd_ab(1,1) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,1)/zI )
        zd_ab(1,2) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,2)/zI )
        zd_ab(2,1) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,1)/zI )
        zd_ab(2,2) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,2)/zI )

! partially prepare spin-Hamiltonian
        zDb_CTEF(1,1) = zi*real(zd_ab(1,1)) &
          - real(zI*zd_ab(1,2)*zs_cd(1,2) + zI*zs_ab(1,2)*zd_cd(1,2))
        zDb_CTEF(1,2) = zI*zd_ab(1,2)*zs_cd(1,2) + zI*zs_ab(1,2)*zd_cd(1,2)
        zDb_CTEF(2,1) = zI*zd_ab(2,1)*zs_cd(2,1) + zI*zs_ab(2,1)*zd_cd(2,1)
        zDb_CTEF(2,2) = zi*real(zd_ab(2,2)) &
          - real(zI*zd_ab(2,1)*zs_cd(2,1) + zI*zs_ab(2,1)*zd_cd(2,1))


        do iho = 1,Num_HO
          zHb_tot_t = zHb_CTEF-zDb_CTEF
          zHb_tot_t(1,1) = zHb_tot_t(1,1) + omega_HO(iho)
          zHb_tot_t(2,2) = zHb_tot_t(2,2) + omega_HO(iho)
          zHO_dot_t(:,iho) = matmul(zHb_tot_t,zHO_t(:,iho)) &
            -sqrt(0.5d0/(M_HO*omega_HO(iho)))*zFb_CTEF(:)
        end do
        zHO_dot_t = matmul(zSb_inv_CTEF,zHO_dot_t)/zI

      end do


! partially prepare spin-Hamiltonian
        zDs_CTEF(1,1) = real(zI*zd_cd(1,1))
        zDs_CTEF(1,2) = zI*zd_cd(1,2)
        zDs_CTEF(2,1) = zI*zd_cd(2,1)
        zDs_CTEF(2,2) = real(zI*zd_cd(2,2))

        call zhpsi_CTEF(zpsi_t,zhpsi_t)

        zd_ab(1,1) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,1)/zI )
        zd_ab(1,2) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,2)/zI )
        zd_ab(2,1) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,1)/zI )
        zd_ab(2,2) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,2)/zI )

! partially prepare spin-Hamiltonian
        zDb_CTEF(1,1) = zi*real(zd_ab(1,1)) &
          - real(zI*zd_ab(1,2)*zs_cd(1,2) + zI*zs_ab(1,2)*zd_cd(1,2))
        zDb_CTEF(1,2) = zI*zd_ab(1,2)*zs_cd(1,2) + zI*zs_ab(1,2)*zd_cd(1,2)
        zDb_CTEF(2,1) = zI*zd_ab(2,1)*zs_cd(2,1) + zI*zs_ab(2,1)*zd_cd(2,1)
        zDb_CTEF(2,2) = zi*real(zd_ab(2,2)) &
          - real(zI*zd_ab(2,1)*zs_cd(2,1) + zI*zs_ab(2,1)*zd_cd(2,1))


  end subroutine refine_effective_hamiltonian
!-----------------------------------------------------------------------------------------
  subroutine calc_Szt(tSz,zpsi_t,zHO_t)
    implicit none
    real(8),intent(out) :: tSz
    complex(8),intent(in) :: zpsi_t(2,2), zHO_t(2,Num_HO)
    complex(8) :: zs, zSz_ab(2,2), zvec(2)

    zs = sum(exp(&
            -0.5d0*abs(zHO_t(1,:))**2 -0.5d0*abs(zHO_t(2,:))**2  &
            + conjg(zHO_t(1,:))*zHO_t(2,:) &
            ))
    
    tSz = abs(zpsi_t(1,1))**2 - abs(zpsi_t(2,1))**2 &
        + abs(zpsi_t(1,2))**2 - abs(zpsi_t(2,2))**2
    zvec(:) = matmul(Sz,zpsi_t(:,2))
    tSz = tSz + 2d0*real( sum(conjg(zpsi_t(:,1))*zvec(:)) )

  end subroutine calc_Szt
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
end module CTEF_module
