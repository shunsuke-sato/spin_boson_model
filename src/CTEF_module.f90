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
  complex(8) :: zEb_CTEF(2,2),zFb_CTEF(2)

  contains
    subroutine set_hamiltonian_spin(zHO_tmp, zHO_dot_tmp)
      implicit none
      complex(8) :: zHO_tmp(Num_HO,2),zHO_dot_tmp(Num_HO,2)
      complex(8) :: zX_Cint_av(2,2)
      complex(8) :: zs
      integer :: i,j

      zSs_CTEF(1,1) = 1d0
      zSs_CTEF(2,2) = 1d0
      zSs_CTEF(1,2) = sum(exp(&
            -0.5d0*abs(zHO_tmp(:,1))**2 -0.5d0*abs(zHO_tmp(:,2))**2  &
            + conjg(zHO_tmp(:,2))*zHO_tmp(:,2) &
            ))
      zSs_CTEF(2,1) = conjg(zSs_CTEF(1,2))

      do i = 1,2
        do j =1,2
          zX_Cint_av(i,j) = sum( &
            Cint_HO(:) &
            *sqrt(0.5d0/(M_HO*Omega_HO(:)))*(conjg(zHO_tmp(:,i))+zHO_tmp(:,j)) )
          zX_Cint_av(i,j) = zX_Cint_av(i,j) *zSs_CTEF(i,j) 
        end do
      end do
      

      do i = 1,2
        do j =1,2
          zHs_CTEF(:,:,i,j) = (eps_SP*Sz(:,:) + delta_SP*Sx(:,:))*zSs_CTEF(i,j) &
            -zX_Cint_av(i,j)*Sz(:,:)
        end do
      end do


      zDs_CTEF(1,1) = real(zI*sum(0.5d0*(conjg(zHO_tmp(:,1))*zHO_dot_tmp(:,1) &
                          -zHO_tmp(:,1)*conjg(zHO_dot_tmp(:,1))) ))
      zDs_CTEF(2,2) = real((zI*sum(0.5d0*(conjg(zHO_tmp(:,2))*zHO_dot_tmp(:,2) &
                          -zHO_tmp(:,2)*conjg(zHO_dot_tmp(:,2))) ))
      zDs_CTEF(1,2) = zI*sum(- 0.5d0*(conjg(zHO_tmp(:,2))*zHO_dot_tmp(:,2) &
                          -zHO_tmp(:,2)*conjg(zHO_dot_tmp(:,2))) &
                          +conjg(zHO_tmp(:,1))*zHO_dot_tmp(:,2))*zSs_CTEF(1,2) 
      zDs_CTEF(2,1) = zI*sum(- 0.5d0*(conjg(zHO_tmp(:,1))*zHO_dot_tmp(:,1) &
                          -zHO_tmp(:,1)*conjg(zHO_dot_tmp(:,1))) &
                          +conjg(zHO_tmp(:,2))*zHO_dot_tmp(:,1))*zSs_CTEF(2,1) 
      
    end subroutine set_hamiltonian_spin

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

    subroutine set_hamiltonian_bath(zpsi_t,zhpsi_t,zHO_t, zHO_dot_t)
      implicit none
      complex(8) :: zHO_tmp(Num_HO,2),zHO_dot_tmp(Num_HO,2)      
      complex(8) :: zs

      zSb_CTEF(1,1) = sum(abs(zpsi_t(:,1))**2)
      zSb_CTEF(2,2) = sum(abs(zpsi_t(:,2))**2)
      zSb_CTEF(1,2) = sum(conjg(zpsi_t(:,1))*zpsi(:,2))
      zSb_CTEF(2,1) = conjg(zSb_CTEF(1,2))
      zSb_CTEF = zSb_CTEF*zSs_CTEF

      zDb_CTEF(1,1) = zI*real(sum(conjg(zpsi_t(:,1))*zhpsi_t(:,1)/zI))
      zDb_CTEF(2,2) = zI*real(sum(conjg(zpsi_t(:,2))*zhpsi_t(:,2)/zI))
      zDb_CTEF(1,2) = zI*sum(conjg(zpsi_t(:,1))*zhpsi_t(:,2)/zI)*zSs_CTEF(1,2) &
        + zI*sum(- 0.5d0*(conjg(zHO_tmp(:,2))*zHO_dot_tmp(:,2) &
                          -zHO_tmp(:,2)*conjg(zHO_dot_tmp(:,2))) &
                          *sum(conjg(zpsi_t(:,1))*zpsi_t(:,2))
      zDb_CTEF(2,1) = zI*sum(conjg(zpsi_t(:,2))*zhpsi_t(:,1)/zI)*zSs_CTEF(2,1) &
        + zI*sum(- 0.5d0*(conjg(zHO_tmp(:,1))*zHO_dot_tmp(:,1) &
                          -zHO_tmp(:,1)*conjg(zHO_dot_tmp(:,1))) &
                          *sum(conjg(zpsi_t(:,1))*zpsi_t(:,2))


    end subroutine set_hamiltonian_bath
  
end module CTEF_module
