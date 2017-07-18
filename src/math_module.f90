!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module math_module
  implicit none

  private

  public :: inverse_2x2_matrix

  interface inverse_2x2_matrix
! real (double precision)
     module procedure inverse_2x2_matrix_double
! complex (double precision)
     module procedure inverse_2x2_matrix_complex
  end interface inverse_2x2_matrix

contains

!-----------------------------------------------------------------------------------------
  subroutine inverse_2x2_matrix_double(Amat,invAmat)
    implicit none
    real(8), intent(in) :: Amat(2,2)
    real(8), intent(out) :: invAmat(2,2)
    real(8) :: det

    det = Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)
    invAmat(1,1) = Amat(2,2)/det
    invAmat(1,2) = -Amat(1,2)/det
    invAmat(2,1) = -Amat(2,1)/det
    invAmat(2,2) = Amat(2,2)/det
    
  end subroutine inverse_2x2_matrix_double
!-----------------------------------------------------------------------------------------
  subroutine inverse_2x2_matrix_complex(Amat,invAmat)
    implicit none
    complex(8), intent(in) :: Amat(2,2)
    complex(8), intent(out) :: invAmat(2,2)
    complex(8) :: det

    det = Amat(1,1)*Amat(2,2) - Amat(1,2)*Amat(2,1)
    invAmat(1,1) = Amat(2,2)/det
    invAmat(1,2) = -Amat(1,2)/det
    invAmat(2,1) = -Amat(2,1)/det
    invAmat(2,2) = Amat(2,2)/det
    
  end subroutine inverse_2x2_matrix_complex

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
end module math_module
