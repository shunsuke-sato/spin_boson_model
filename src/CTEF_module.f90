!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module CTEF_module
  use global_variables
  implicit none

! spin
  complex(8) :: zpsi_CTEF(2,2)
! Harmonic oscillator 
  complex(8) :: ZHO_CTEF(Num_HO,2)
  complex(8) :: ZHO_dot_CTEF(Num_HO,2)
!
  complex(8) :: zHs_CTEF(2,2,2,2)
  complex(8) :: zDs_CTEF(2,2), zSs_CTEF(2,2)

  complex(8) :: zDb_CTEF(2,2), zSb_CTEF(2,2)
  complex(8) :: zHb_CTEF(2,2),zFb_CTEF(2)

  contains

    subroutine zhpsi_CTEF(zpsi_t,zhpsi_t)
      implicit none
      complex(8),intent(inout) :: zpsi_t(2,2)
      complex(8),intent(out) :: zhpsi_t(2,2)
      complex(8) :: zS_inv(2,2)
      real(8) :: det

      zhpsi_t(:,1) = zDs_CTEF(1,1)*zpsi_t(:,1) + zDs_CTEF(1,2)*zpsi_t(:,2)
      zhpsi_t(:,2) = zDs_CTEF(2,1)*zpsi_t(:,1) + zDs_CTEF(2,2)*zpsi_t(:,2)

      zhpsi_t(1,1) = zhpsi_t(1,1) + zHs_CTEF(1,1,1,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(1,2,1,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(1,1,1,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(1,2,1,2)*zpsi_t(2,2) 

      zhpsi_t(2,1) = zhpsi_t(2,1) + zHs_CTEF(2,1,1,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(2,2,1,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(2,1,1,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(2,2,1,2)*zpsi_t(2,2) 

      zhpsi_t(1,2) = zhpsi_t(1,2) + zHs_CTEF(1,1,2,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(1,2,2,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(1,1,2,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(1,2,2,2)*zpsi_t(2,2) 

      zhpsi_t(2,2) = zhpsi_t(2,2) + zHs_CTEF(2,1,2,1)*zpsi_t(1,1) &
                                  + zHs_CTEF(2,2,2,1)*zpsi_t(2,1) &
                                  + zHs_CTEF(2,1,2,2)*zpsi_t(1,2) &
                                  + zHs_CTEF(2,2,2,2)*zpsi_t(2,2) 

      det = zSs_CTEF(1,1)*zSs_CTEF(2,2) - abs(zSs_CTEF(1,2))**2
      zS_inv(1,1) = zSs_CTEF(2,2)
      zS_inv(2,2) = zSs_CTEF(1,1)
      zS_inv(1,2) = -zSs_CTEF(1,2)
      zS_inv(2,1) = -zSs_CTEF(2,1)
      zS_inv = zS_inv/det
      
      zhpsi_t = matmul(zS_inv,zhpsi_t)

    end subroutine zhpsi_CTEF

    subroutine refine_zHO_dot(zpsi_t,zHO_t,zHO_dot_t)
      implicit none
      complex(8),intent(in) :: zpsi_t(2,2), zHO_t(Num_HO,2)
      complex(8),intent(inout) :: zHO_dot_t(Num_HO,2) 
      integer,parameter :: Nscf_CTEF = 2
      complex(8) :: zs_ab(2,2) ! overlap matrix for spin
      complex(8) :: zs_cd(2,2) ! overlap matrix for bath
      complex(8) :: zd_ab(2,2) ! time-derivative-overlap matrix for bath
      complex(8) :: zd_cd(2,2) ! time-derivative-overlap matrix for bath
      complex(8) :: zX_Cint_av(2,2), zEb(2,2)
      complex(8) :: zs,zvec(2)
      integer :: iter
      integer :: i,j

      zs_ab(1,1) = 1d0
      zs_ab(2,2) = 1d0
      zs_ab(1,2) = sum(conjg(zpsi_t(:,1))*zpsi_t(:,2))
      zs_ab(2,1) = conjg(zs_ab(1,2))

      zs_cd(1,1) = 1d0
      zs_cd(2,2) = 1d0
      zs_cd(1,2) = sum(exp(&
            -0.5d0*abs(zHO_t(:,1))**2 -0.5d0*abs(zHO_t(:,2))**2  &
            + conjg(zHO_t(:,1))*zHO_t(:,2) &
            ))
      zs_cd(2,1) = conjg(zs_cd(1,2))

! partially prepare spin-Hamiltonian
      zSs_CTEF(:,:) = zs_ab(:,:)

      do i = 1,2
        do j =1,2
          zX_Cint_av(i,j) = sum( &
            Cint_HO(:) &
            *sqrt(0.5d0/(M_HO*Omega_HO(:)))*(conjg(zHO_t(:,i))+zHO_t(:,j)) )
          zX_Cint_av(i,j) = zX_Cint_av(i,j) *zSs_CTEF(i,j) 
          zEb(i,j) = sum(conjg(zHO_t(:,i))*zHO_t(:,j))+dble(Num_HO)/2d0
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

      zvec(:) = matmul(zHs_CTEF(:,:,1,2),zpsi_t(:,2))
      zs = sum(conjg(zpsi_t(:,1))*zvec(:))
      zHb_CTEF(1,1) = -real(zs)
      zHb_CTEF(1,2) = zs
      zvec(:) = matmul(zHs_CTEF(:,:,2,1),zpsi_t(:,1))
      zs = sum(conjg(zpsi_t(:,2))*zvec(:))
      zHb_CTEF(2,2) = -real(zs)
      zHb_CTEF(2,1) = zs

      do iter = 1,Nscf_CTEF


        zd_cd(1,1) = 0.5d0*sum(conjg(zHO_t(:,1))*zHO_dot_t(:,1) &
          - conjg(zHO_dot_t(:,1))*zHO_t(:,1))

        zd_cd(2,1) = sum(-0.5d0*(conjg(zHO_t(:,1))*zHO_dot_t(:,1) &
                               + conjg(zHO_dot_t(:,1))*zHO_t(:,1) ) &
                               +conjg(zHO_t(:,2))*zHO_dot_t(:,1) )

        zd_cd(1,2) = sum(-0.5d0*(conjg(zHO_t(:,2))*zHO_dot_t(:,2) &
                               + conjg(zHO_dot_t(:,2))*zHO_t(:,2) ) &
                               +conjg(zHO_t(:,1))*zHO_dot_t(:,2) )

        zd_cd(2,2) = 0.5d0*sum(conjg(zHO_t(:,2))*zHO_dot_t(:,2) &
          - conjg(zHO_dot_t(:,2))*zHO_t(:,2))

        zd_cd = zd_cd * zs_cd

! partially prepare spin-Hamiltonian
        zDs_CTEF(1,1) = real(zI*zd_cd(1,1))
        zDs_CTEF(1,2) = zI*zd_cd(1,2)
        zDs_CTEF(2,1) = zI*zd_cd(2,1)
        zDs_CTEF(2,2) = real(zI*zd_cd(2,2))



      end do


    end subroutine refine_zHO_dot
end module CTEF_module
