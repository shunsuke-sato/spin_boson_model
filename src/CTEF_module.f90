!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module CTEF_module
  use math_module
  use global_variables
  use random_number_module
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
  complex(8) :: zSzs_CTEF(2,2), zEs_CTEF(2,2)
  complex(8) :: zDs_CTEF(2,2), zSs_CTEF(2,2)
  complex(8) :: zSs_inv_CTEF(2,2)

  complex(8) :: zDb_CTEF(2,2), zSb_CTEF(2,2)
  complex(8) :: zSb_inv_CTEF(2,2)
  complex(8) :: zHb_CTEF(2,2),zFb_CTEF(2,Num_HO)
  complex(8) :: zXb_cint_CTEF(2,2),zEb_CTEF(2,2)

  complex(8) :: zSsb_CTEF(2,2), zSsb_inv_CTEF(2,2), zEc_CTEF(2,2)

  public :: CTEF

  contains
!-----------------------------------------------------------------------------------------
    subroutine CTEF
      implicit none
      integer :: itraj,it,iphi
      real(8) :: norm
      real(8) :: phi0,phi
      complex(8) :: zpsi_stored(2,2), zHO_stored(2,num_HO)
      complex(8) :: zweight
      real(8) :: Szt_c(0:Nt),norm_c(0:Nt),Eb_c(0:Nt),Ec_c(0:Nt)
      real(8) :: Szt_cl(0:Nt),norm_cl(0:Nt),Eb_cl(0:Nt),Ec_cl(0:Nt)
      real(8) :: Szt_ct(0:Nt),norm_ct(0:Nt),Eb_ct(0:Nt),Ec_ct(0:Nt)

      Szt_cl=0d0; norm_cl=0d0;Eb_cl=0d0;Ec_cl=0d0
      call setting_bath_parameters
      
      do itraj = 1,Ntraj
        call set_forward_backward_trajectries(zpsi_stored,zHO_stored,zweight)
        call random_number(phi0); phi0 = 2d0*pi*phi0

        if(mod(itraj,Nprocs) /= myrank)cycle
        if(myrank == 0)write(*,*)"itraj=",itraj,"/",Ntraj
        do iphi = 0,Nphi-1
          phi = phi0 + 2d0*pi*dble(iphi)/Nphi
          call set_initial_condition(zpsi_stored,zHO_stored, &
                                     zpsi_CTEF,  zHO_CTEF, phi, norm)

          call CTEF_dynamics(norm_ct, Szt_ct, Eb_ct, Ec_ct)
          Szt_cl  = Szt_cl  + norm*exp(-zI*phi)*zweight*Szt_ct
          norm_cl = norm_cl + norm*exp(-zI*phi)*zweight*norm_ct
          Eb_cl   = Eb_cl   + norm*exp(-zI*phi)*zweight*Eb_ct
          Ec_cl   = Ec_cl   + norm*exp(-zI*phi)*zweight*Ec_ct


        end do

      end do

      Szt_cl  = Szt_cl/(Ntraj*Nphi)
      norm_cl = norm_cl/(Ntraj*Nphi)
      Eb_cl   = Eb_cl/(Ntraj*Nphi)
      Ec_cl   = Ec_cl/(Ntraj*Nphi)
      call MPI_ALLREDUCE(Szt_cl,Szt_c,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(norm_cl,norm_c,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Eb_cl,Eb_c,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE(Ec_cl,Ec_c,Nt+1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ierr)

      if(myrank == 0)then
        open(nfile_CTEF_Sz,file=file_CTEF_Sz)
        do it = 0,Nt
          write(nfile_CTEF_Sz,"(999e26.16e3)")dt*it,norm_c(it), &
                                                    Szt_c(it),  &
                                                    Eb_c(it),   &
                                                    Ec_c(it) 
        end do
        close(nfile_CTEF_Sz)
      end if


    end subroutine CTEF
!-----------------------------------------------------------------------------------------
    subroutine CTEF_dynamics(norm_ct, Szt_ct, Eb_ct, Ec_ct)
! Assuming zpsi_CTEF and zHO_CTEF are given 
      implicit none
      real(8),intent(out) :: Szt_ct(0:Nt),norm_ct(0:Nt),Eb_ct(0:Nt),Ec_ct(0:Nt)
      real(8) :: Sz_av
      integer :: it
      real(8) :: X_Cint_av

      
      Szt_ct   = 0d0
      norm_ct  = 0d0
      Eb_ct    = 0d0
      Ec_ct    = 0d0



      zHO_dot_CTEF = 0d0
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF, zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF, zHO_dot_CTEF)
      call refine_effective_hamiltonian(zpsi_CTEF,zHO_CTEF, zHO_dot_CTEF)
      call calc_output(norm_ct(0),Szt_ct(0),Eb_ct(0),Ec_ct(0))

      do it = 0,Nt-1

        call dt_evolve_etrs(zpsi_CTEF,zHO_CTEF, zHO_dot_CTEF)
        call calc_output(norm_ct(it+1),Szt_ct(it+1),Eb_ct(it+1),Ec_ct(it+1))

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
      call dt_evolve_bath_direct(zHO_t, zHO_dot_t, 0.5d0*dt)
      call refine_effective_hamiltonian(zpsi_t,zHO_CTEF, zHO_dot_CTEF)
      zpsi_s = zpsi_t
      zHO_s = zHO_t

      do ipred_corr = 1, Npred_corr
        zpsi_t = zpsi_s
        zHO_t = zHO_s

        call dt_evolve_spin_Taylor(zpsi_t,0.5d0*dt)
        call dt_evolve_bath_direct(zHO_t, zHO_dot_t, 0.5d0*dt)
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
    subroutine dt_evolve_bath_direct(zHO_t,zHO_dot_t,dt_t)
      implicit none
      complex(8),intent(inout) :: zHO_t(2,Num_HO)
      complex(8),intent(in) :: zHO_dot_t(2,Num_HO)
      real(8),intent(in) :: dt_t

      zHO_t = zHO_t + dt_t * zHO_dot_t

    end subroutine dt_evolve_bath_direct

!-----------------------------------------------------------------------------------------
    subroutine zhpsi_CTEF(zpsi_t,zhpsi_t)
      implicit none
      complex(8),intent(in) :: zpsi_t(2,2)
      complex(8),intent(out) :: zhpsi_t(2,2)
      complex(8) :: zhpsi_s(2,2)
      real(8) :: det


      zhpsi_s(1,1) =  zHs_CTEF(1,1,1,1)*zpsi_t(1,1) &
                    + zHs_CTEF(1,2,1,1)*zpsi_t(2,1) &
                    + zHs_CTEF(1,1,1,2)*zpsi_t(1,2) &
                    + zHs_CTEF(1,2,1,2)*zpsi_t(2,2) 

      zhpsi_s(2,1) =  zHs_CTEF(2,1,1,1)*zpsi_t(1,1) &
                    + zHs_CTEF(2,2,1,1)*zpsi_t(2,1) &
                    + zHs_CTEF(2,1,1,2)*zpsi_t(1,2) &
                    + zHs_CTEF(2,2,1,2)*zpsi_t(2,2) 

      zhpsi_s(1,2) =  zHs_CTEF(1,1,2,1)*zpsi_t(1,1) &
                    + zHs_CTEF(1,2,2,1)*zpsi_t(2,1) &
                    + zHs_CTEF(1,1,2,2)*zpsi_t(1,2) &
                    + zHs_CTEF(1,2,2,2)*zpsi_t(2,2) 

      zhpsi_s(2,2) =  zHs_CTEF(2,1,2,1)*zpsi_t(1,1) &
                    + zHs_CTEF(2,2,2,1)*zpsi_t(2,1) &
                    + zHs_CTEF(2,1,2,2)*zpsi_t(1,2) &
                    + zHs_CTEF(2,2,2,2)*zpsi_t(2,2) 


      zhpsi_s(:,1) = zhpsi_s(:,1) - zDb_CTEF(1,1)*zpsi_t(:,1) - zDb_CTEF(1,2)*zpsi_t(:,2)
      zhpsi_s(:,2) = zhpsi_s(:,2) - zDb_CTEF(2,1)*zpsi_t(:,1) - zDb_CTEF(2,2)*zpsi_t(:,2)

      zhpsi_t(:,1) = zSb_inv_CTEF(1,1)*zhpsi_s(:,1) + zSb_inv_CTEF(1,2)*zhpsi_s(:,2)
      zhpsi_t(:,2) = zSb_inv_CTEF(2,1)*zhpsi_s(:,1) + zSb_inv_CTEF(2,2)*zhpsi_s(:,2)


    end subroutine zhpsi_CTEF
!-----------------------------------------------------------------------------------------
    subroutine refine_effective_hamiltonian(zpsi_t,zHO_t,zHO_dot_t)
      implicit none
      complex(8),intent(in) :: zpsi_t(2,2), zHO_t(2,Num_HO)
      complex(8),intent(inout) :: zHO_dot_t(2,Num_HO) 
      integer,parameter :: Nscf_CTEF = 2
      complex(8) :: zHb_tot_t(2,2)
      complex(8) :: zs,zvec(2)
      complex(8) :: zhpsi_t(2,2)
      integer :: iho
      integer :: iter
      integer :: i,j

! zSs
      zSs_CTEF(1,1) = sum(abs(zpsi_t(:,1))**2)
      zSs_CTEF(2,1) = sum(conjg(zpsi_t(:,2))*zpsi_t(:,1))
      zSs_CTEF(1,2) = sum(conjg(zpsi_t(:,1))*zpsi_t(:,2))
      zSs_CTEF(2,2) = sum(abs(zpsi_t(:,2))**2)

! zSb
      zSb_CTEF(1,1) = 1d0
      zSb_CTEF(2,2) = 1d0
      zSb_CTEF(1,2) = exp(sum(&
            -0.5d0*abs(zHO_t(1,:))**2 -0.5d0*abs(zHO_t(2,:))**2  &
            + conjg(zHO_t(1,:))*zHO_t(2,:) &
            ))
      zSb_CTEF(2,1) = conjg(zSb_CTEF(1,2))

      call inverse_2x2_matrix(zSb_CTEF,zSb_inv_CTEF)

      zSsb_CTEF = zSs_CTEF*zSb_CTEF
      call inverse_2x2_matrix(zSsb_CTEF,zSsb_inv_CTEF)


      do i = 1,2
        do j =1,2
          zXb_cint_CTEF(i,j) = sum( &
            Cint_HO(:) &
            *sqrt(0.5d0/(M_HO*Omega_HO(:)))*(conjg(zHO_t(i,:))+zHO_t(j,:)) )
          zXb_cint_CTEF(i,j) = zxb_cint_CTEF(i,j) *zSs_CTEF(i,j) 
          zEb_CTEF(i,j) = sum(Omega_HO(:)*conjg(zHO_t(i,:))*zHO_t(j,:))
        end do
      end do
      zXb_cint_CTEF = zXb_cint_CTEF *zSb_CTEF
!      zEb_CTEF = (zEb_CTEF + sum(omega_HO)*0.5d0)*zSb_CTEF ! with zero-point energy
      zEb_CTEF = zEb_CTEF*zSb_CTEF ! without zero-point energy
      
! zEs
      zvec(:) = matmul(hs_m,zpsi_t(:,1))
      zEs_CTEF(1,1) = sum( conjg(zpsi_t(:,1))*zvec(:))
      zEs_CTEF(2,1) = sum( conjg(zpsi_t(:,2))*zvec(:))
      zEs_CTEF(1,2) = conjg(zEs_CTEF(2,1))
      zvec(:) = matmul(hs_m,zpsi_t(:,2))
      zEs_CTEF(2,2) = sum( conjg(zpsi_t(:,2))*zvec(:))

! zSzs
      zvec(:) = matmul(Sz,zpsi_t(:,1))
      zSzs_CTEF(1,1) = sum( conjg(zpsi_t(:,1))*zvec(:))
      zSzs_CTEF(2,1) = sum( conjg(zpsi_t(:,2))*zvec(:))
      zSzs_CTEF(1,2) = conjg(zSzs_CTEF(2,1))
      zvec(:) = matmul(Sz,zpsi_t(:,2))
      zSzs_CTEF(2,2) = sum( conjg(zpsi_t(:,2))*zvec(:))

! zEc
      zEc_CTEF = -zXb_cint_CTEF*zSzs_CTEF

! effective hamiltonian for spin
      do i = 1,2
        do j =1,2
          zHs_CTEF(:,:,i,j) = hs_m(:,:)*zSb_CTEF(i,j) &
            -zXb_cint_CTEF(i,j)*Sz(:,:) + zEb_CTEF(i,j)
        end do
      end do

      do i = 1, Num_HO
        zFb_CTEF(1,i) = - (zSzs_CTEF(1,1)*zSb_CTEF(1,1)  + zSzs_CTEF(1,2)*zSb_CTEF(1,2)) &
          *cint_HO(i)*sqrt(0.5d0/(M_HO*omega_HO(i)))
        zFb_CTEF(2,i) = - (zSzs_CTEF(2,2)*zSb_CTEF(2,2)  + zSzs_CTEF(2,1)*zSb_CTEF(2,1)) &
          *cint_HO(i)*sqrt(0.5d0/(M_HO*omega_HO(i)))
      end do


      do iter = 1,Nscf_CTEF


        zDb_CTEF(1,1) = real( 0.5d0*zI*sum(conjg(zHO_t(:,1))*zHO_dot_t(:,1) &
          -conjg(zHO_dot_t(:,1))*zHO_t(:,1) ) )
        zDb_CTEF(2,2) = real( 0.5d0*zI*sum(conjg(zHO_t(:,2))*zHO_dot_t(:,2) &
          -conjg(zHO_dot_t(:,2))*zHO_t(:,2) ) )
        zDb_CTEF(1,2) = zI*sum(-0.5d0*( &
          conjg(zHO_dot_t(:,2))*zHO_t(:,2) &
          +zHO_dot_t(:,2)*conjg(zHO_t(:,2)) ) &
          +conjg(zHO_t(:,1))*zHO_dot_t(:,2) )*zSb_CTEF(1,2)
        zDb_CTEF(2,1) = zI*sum(-0.5d0*( &
          conjg(zHO_dot_t(:,1))*zHO_t(:,1) &
          +zHO_dot_t(:,1)*conjg(zHO_t(:,1)) ) &
          +conjg(zHO_t(:,2))*zHO_dot_t(:,1) )*zSb_CTEF(2,1)


! implementing here
        call zhpsi_CTEF(zpsi_t,zhpsi_t)

        zDs_CTEF(1,1) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,1) )
        zDs_CTEF(2,2) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,2) )
        zDs_CTEF(1,2) = sum( conjg(zpsi_t(:,1))*zhpsi_t(:,2) )
        zDs_CTEF(2,1) = sum( conjg(zpsi_t(:,2))*zhpsi_t(:,1) )



        do iho = 1,Num_HO
          zHb_tot_t = omega_HO(iho)*zSsb_CTEF
          zHb_tot_t(1,1) = zHb_tot_t(1,1) - zi*real(zDs_CTEF(1,1)/zI) &
            + real(zDs_CTEF(1,2)*zSb_CTEF(1,2) + zSs_CTEF(1,2)*zDb_CTEF(1,2)) &
            - real(zEs_CTEF(1,2)*zSb_CTEF(1,2) + zEc_CTEF(1,2) + zSs_CTEF(1,2)*zEb_CTEF(1,2))

          zHb_tot_t(2,2) = zHb_tot_t(2,2) - zi*real(zDs_CTEF(2,2)/zI) &
            + real(zDs_CTEF(2,1)*zSb_CTEF(2,1) + zSs_CTEF(2,1)*zDb_CTEF(2,1)) &
            - real(zEs_CTEF(2,1)*zSb_CTEF(2,1) + zEc_CTEF(2,1) + zSs_CTEF(2,1)*zEb_CTEF(2,1))

          zHb_tot_t(1,2) = zHb_tot_t(1,2) &
            -zDs_CTEF(1,2)*zSb_CTEF(1,2) - zSs_CTEF(1,2)*zDb_CTEF(1,2) &
            + zEs_CTEF(1,2)*zSb_CTEF(1,2) + zEc_CTEF(1,2) + zSs_CTEF(1,2)*zEb_CTEF(1,2)

          zHb_tot_t(2,1) = zHb_tot_t(2,1) &
            -zDs_CTEF(2,1)*zSb_CTEF(2,1) - zSs_CTEF(2,1)*zDb_CTEF(2,1) &
            + zEs_CTEF(2,1)*zSb_CTEF(2,1) + zEc_CTEF(2,1) + zSs_CTEF(2,1)*zEb_CTEF(2,1)

          zHO_dot_t(:,iho) = matmul(zHb_tot_t,zHO_t(:,iho)) + zFb_CTEF(:,iho)
        end do
        zHO_dot_t = matmul(zSs_inv_CTEF,zHO_dot_t)/zI

      end do

  end subroutine refine_effective_hamiltonian
!-----------------------------------------------------------------------------------------
  subroutine set_forward_backward_trajectries(zpsi_t,zHO_t,zweight)
    implicit none
    complex(8),intent(out) :: zpsi_t(2,2), zHO_t(2,Num_HO), zweight
    integer :: iho
    real(8) :: q1,p1,q2,p2
    real(8) :: beta_ww, exp_beta_ww
    complex(8) :: z1,z2

    zpsi_t(1,:) = 1d0
    zpsi_t(2,:) = 0d0
    
    do iho = 1,num_HO
      call gaussian_random_number(q1,p1)
      call gaussian_random_number(q2,p2)
      z1 = q1 + zi*p1
      z2 = q2 + zi*p2
      zHO_t(1,iho) = z1 * sqrt(2d0/3d0)
      zHO_t(2,iho) = z2*sqrt(0.5d0) + 0.5d0*zHO_t(1,iho)
    end do


    zweight = 1d0
    do iho = 1, num_ho
      
      beta_ww = beta_kB * omega_ho(iho)
      exp_beta_ww = exp(-beta_ww)
      zweight = zweight *(4d0/3d0)* (1d0-exp_beta_ww)*exp(exp_beta_ww &
        *conjg(zHO_t(2,iho))*zHO_t(1,iho)) &
        *exp(0.5d0*abs(zHO_t(1,iho)-zHO_t(2,iho))**2)
    end do
!    if(myrank == 0)write(*,*)"zweight0=",zweight
    
  end subroutine set_forward_backward_trajectries
!-----------------------------------------------------------------------------------------
  subroutine set_initial_condition(zpsi_in,zHO_in, zpsi_out,zHO_out, phi, norm)
    implicit none
    complex(8),intent(in) :: zpsi_in(2,2),zHO_in(2,num_HO)
    complex(8),intent(out) :: zpsi_out(2,2),zHO_out(2,num_HO)
    real(8),intent(in) :: phi
    real(8),intent(out) :: norm
    complex(8) :: zs_spin, zs_bath


    zpsi_out(:,1) = zpsi_in(:,1)
    zpsi_out(:,2) = exp(zI*phi)*zpsi_in(:,2)
    zHO_out = zHO_in

    zs_spin = sum(conjg(zpsi_out(:,1))*zpsi_out(:,2))
    zs_bath = exp(sum(& 
            -0.5d0*abs(zHO_out(1,:))**2 -0.5d0*abs(zHO_out(2,:))**2  &
            + conjg(zHO_out(1,:))*zHO_out(2,:) &
            ))


    norm = sum( abs(zpsi_out(:,:))**2 ) + 2d0*real( zs_spin* zs_bath)
    
    zpsi_out = zpsi_out/ sqrt(norm)


  end subroutine set_initial_condition
!-----------------------------------------------------------------------------------------
  subroutine calc_Szt(tSz,zpsi_t,zHO_t)
    implicit none
    real(8),intent(out) :: tSz
    complex(8),intent(in) :: zpsi_t(2,2), zHO_t(2,Num_HO)
    complex(8) :: zs, zSz_ab(2,2), zvec(2)

    zs = exp(sum(&
            -0.5d0*abs(zHO_t(1,:))**2 -0.5d0*abs(zHO_t(2,:))**2  &
            + conjg(zHO_t(1,:))*zHO_t(2,:) &
            ))
    
    tSz = abs(zpsi_t(1,1))**2 - abs(zpsi_t(2,1))**2 &
        + abs(zpsi_t(1,2))**2 - abs(zpsi_t(2,2))**2
    zvec(:) = matmul(Sz,zpsi_t(:,2))
    tSz = tSz + 2d0*real( sum(conjg(zpsi_t(:,1))*zvec(:))*zs ) 

  end subroutine calc_Szt
!-----------------------------------------------------------------------------------------
  subroutine calc_output(norm_ct0,Szt_ct0,Eb_ct0,Ec_ct0)
    implicit none
    real(8),intent(out) :: Szt_ct0,norm_ct0,Eb_ct0,Ec_ct0

    norm_ct0 = sum(zSsb_CTEF)
    Szt_ct0  = sum(zEs_CTEF*zSb_CTEF)
    Eb_ct0   = sum(zSs_CTEF*zEb_CTEF)
    Ec_ct0 = sum(zEc_CTEF)

  end subroutine calc_output
!-----------------------------------------------------------------------------------------
end module CTEF_module
