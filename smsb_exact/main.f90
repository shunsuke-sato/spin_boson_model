!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module exact_smsb
  implicit none
! mathamtical constant
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zI = (0d0,1d0)

! model parameters
  real(8),parameter :: eps_SP = 0d0, delta_SP = 0.5d0
  real(8),parameter :: Sz(2,2) = reshape( (/1d0, 0d0, 0d0, -1d0/), (/2,2/) )
  real(8),parameter :: Sx(2,2) = reshape( (/0d0, 1d0, 1d0, 0d0/), (/2,2/) )
  real(8),parameter :: hs(2,2) = reshape( (/eps_SP, delta_SP, delta_SP, -eps_SP/), (/2,2/) )
  real(8),parameter :: mass_p = 1d0
  real(8),parameter :: omega_p = 1d0
  real(8),parameter :: cint_p = 1d0
  real(8),parameter :: gamma_p = 1d0 ! gamma should be computed by cint and mass

! numerical parameter
  integer,parameter :: Nmax = 20
  real(8),parameter :: dt = 0.01d0
  real(8),parameter :: Tprop = 20d0
  integer,parameter :: Nt = aint(Tprop/dt)+1


! wavefunctions
  complex(8) :: zpsi(0:Nmax,2)


! I/O
  character(64),parameter :: file_quantity = "td_quantity.dat"
  integer,parameter :: nfile_quantity = 21

  contains

!-----------------------------------------------------------------------------------------
    subroutine dynamics
      integer :: it
      real(8) :: norm_t(0:Nt),Sz_t(0:Nt),Es_t(0:Nt),Ec_t(0:Nt),Eb_t(0:Nt)

      call init_wf




      do it = 0,Nt
        write(*,"(A,2x,I9,A,I9)")"it=",it,"/",Nt

        call dt_prop(zpsi)
        call calc_quantities(zpsi, &
                             norm_s=norm_t(it), &
                             Sz_s=Sz_t(it), &
                             Es_s=Es_t(it), &
                             Ec_s=Ec_t(it), &
                             Eb_s=Eb_t(it) )

      end do

      
      open(nfile_quantity,file=file_quantity)
      do it = 0,Nt
        write(nfile_quantity,"(999e26.16e3)")dt*it, norm_t(it), &
                                                    Sz_t(it), &
                                                    Es_t(it), &
                                                    Ec_t(it), &
                                                    Eb_t(it), &
                                                    Es_t(it) + Ec_t(it) + Eb_t(it)
        
                                                   
      end do
      close(nfile_quantity)

    end subroutine dynamics
!-----------------------------------------------------------------------------------------
    subroutine init_wf
      
      zpsi = 0d0
      zpsi(0,1) = 1d0

    end subroutine init_wf
!-----------------------------------------------------------------------------------------
    subroutine dt_prop(zpsi_t)
      complex(8),intent(inout) :: zpsi_t(0:Nmax,2)
      complex(8) :: zpsi_s(0:Nmax,2),zhpsi_s(0:Nmax,2)
      integer,parameter :: nTaylor = 6
      complex(8) :: zfact
      integer :: iexp

      zpsi_s = zpsi_t
      zfact = 1d0
      
      do iexp = 1,nTaylor

        zfact = zfact*(-zI*dt)/iexp
        call hamiltonian_full(zpsi_s,zhpsi_s)
        zpsi_t = zpsi_t + zfact*zhpsi_s
        zpsi_s = zhpsi_s
        
      end do
      
    end subroutine dt_prop
!-----------------------------------------------------------------------------------------
    subroutine hamiltonian_full(zpsi_t,zhpsi_t)
      complex(8),intent(in) :: zpsi_t(0:Nmax,2)
      complex(8),intent(out) :: zhpsi_t(0:Nmax,2)
      integer :: ib

! spin
      zhpsi_t(:,1) = hs(1,1)*zpsi_t(:,1) + hs(1,2)*zpsi_t(:,2)
      zhpsi_t(:,2) = hs(2,1)*zpsi_t(:,1) + hs(2,2)*zpsi_t(:,2)

! bath
      do ib = 0,Nmax
        zhpsi_t(ib,1)  = zhpsi_t(ib,1) + omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,1)
        zhpsi_t(ib,2)  = zhpsi_t(ib,2) + omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,2)
      end do

! coupling
!! S_z
      zhpsi_t(0,1)  = zhpsi_t(0,1) -gamma_p*sqrt(1d0)*zpsi_t(1,1)
      zhpsi_t(0,2)  = zhpsi_t(0,2) +gamma_p*sqrt(1d0)*zpsi_t(1,2)
      do ib = 1,Nmax-1
        zhpsi_t(ib,1)  = zhpsi_t(ib,1) -gamma_p*( &
          sqrt(dble(ib))*zpsi_t(ib-1,1) &
         +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,1) )

        zhpsi_t(ib,2)  = zhpsi_t(ib,2) +gamma_p*( &
          sqrt(dble(ib))*zpsi_t(ib-1,2) &
         +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,2) )

      end do

      zhpsi_t(Nmax,1)  = zhpsi_t(Nmax,1) -gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,1)
      zhpsi_t(Nmax,2)  = zhpsi_t(Nmax,2) +gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,2)
      

    end subroutine hamiltonian_full
!-----------------------------------------------------------------------------------------
    subroutine calc_quantities(zpsi_t, norm_s, Sz_s, Es_s,Ec_s,Eb_s )
      complex(8),intent(in) :: zpsi_t(0:Nmax,2)
      real(8),intent(out),optional :: norm_s, Sz_s, Es_s,Ec_s,Eb_s 
      complex(8) :: zhpsi_s(0:Nmax,2)
      integer :: ib

      if(present(norm_s))then
        norm_s = sum(abs(zpsi_t)**2)
      end if

      if(present(Sz_s))then
        zhpsi_s(:,1) = zpsi_t(:,1) 
        zhpsi_s(:,2) = -zpsi_t(:,2)
        Sz_s = sum(conjg(zpsi_t(:,:))*zhpsi_s(:,:))
      end if

      if(present(Es_s))then
        zhpsi_s(:,1) = hs(1,1)*zpsi_t(:,1) + hs(1,2)*zpsi_t(:,2)
        zhpsi_s(:,2) = hs(2,1)*zpsi_t(:,1) + hs(2,2)*zpsi_t(:,2)
        Es_s = sum(conjg(zpsi_t(:,:))*zhpsi_s(:,:))
      end if

      if(present(Ec_s))then

        zhpsi_s(0,1)  = -gamma_p*sqrt(1d0)*zpsi_t(1,1)
        zhpsi_s(0,2)  = +gamma_p*sqrt(1d0)*zpsi_t(1,2)

        do ib = 1,Nmax-1
          zhpsi_s(ib,1)  = -gamma_p*( &
            sqrt(dble(ib))*zpsi_t(ib-1,1) &
            +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,1) )

          zhpsi_s(ib,2)  =  +gamma_p*( &
            sqrt(dble(ib))*zpsi_t(ib-1,2) &
            +sqrt(dble(ib)+1d0)*zpsi_t(ib+1,2) )

        end do

        zhpsi_s(nmax,1)  = -gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,1)
        zhpsi_s(nmax,2)  = +gamma_p*sqrt(dble(nmax))*zpsi_t(nmax-1,2)

        Ec_s = sum(conjg(zpsi_t(:,:))*zhpsi_s(:,:))

      end if

      if(present(Eb_s))then
        do ib = 0,Nmax
          zhpsi_s(ib,1)  = omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,1)
          zhpsi_s(ib,2)  = omega_p*(dble(ib)+0.5d0)*zpsi_t(ib,2)

          Eb_s = sum(conjg(zpsi_t(:,:))*zhpsi_s(:,:))
        end do
      end if


    end subroutine calc_quantities
!-----------------------------------------------------------------------------------------
end module exact_smsb
!-----------------------------------------------------------------------------------------
program main
  use exact_smsb
  implicit none

  call dynamics

end program main
!-----------------------------------------------------------------------------------------
