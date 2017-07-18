!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
module random_number_module
  implicit none
  private

  public :: gaussian_random_number, &
            correlated_gaussian_random_number
contains
  subroutine gaussian_random_number(x1,x2)
    implicit none
    real(8),parameter :: pi = 4d0*atan(1d0)
    real(8), intent(out) :: x1,x2
    real(8) :: r1,r2,tmp
    
    call random_number(r1)
    call random_number(r2)
    
    if(r1 == 0d0)then
      x1 = 0d0
      x2 = 0d0
    else 
      tmp = sqrt(-2d0*log(r1))
      x1 = tmp*cos(2d0*pi*r2)
      x2 = tmp*sin(2d0*pi*r2)
    end if
  end subroutine gaussian_random_number
!-----------------------------------------------------------------------------------------
  subroutine correlated_gaussian_random_number(x1,x2,width, normalization_factor)
    implicit none
    real(8),parameter :: pi = 4d0*atan(1d0)
    real(8), intent(out) :: x1,x2
    real(8),intent(in),optional :: width
    real(8),intent(out),optional ::normalization_factor
    real(8) :: r1,r2,tmp, sigma

    if(present(width))then
      sigma = width
    else
      sigma = 1d0
    end if

    if(sigma >= 1d0)then

      do 
        call random_number(r1)
        call random_number(r2)
      
        if(r1 == 0d0)then
          x1 = 0d0
          x2 = 0d0
        else 
          tmp = sqrt(-2d0*log(r1))
          x1 = tmp*cos(2d0*pi*r2)
          x2 = tmp*sin(2d0*pi*r2)
        end if
      
        tmp = x1 -x2
        r1 = exp(-0.5d0*tmp**2/sigma)
        call random_number(r2)
        if(r2 < r1)exit
      end do

    else

      do 
        call random_number(r1)
        call random_number(r2)
      
        if(r1 == 0d0)then
          x1 = 0d0
          x2 = 0d0
        else 
          tmp = sqrt(-2d0*log(r1))
          x1 = tmp*cos(2d0*pi*r2)
          x2 = tmp*sin(2d0*pi*r2)
        end if
      
        x2 = x1 + x2*sqrt(sigma)
        r1 = exp(-0.5d0*x2**2)
        call random_number(r2)
        if(r2 < r1)exit
      end do
    end if

    if(present(normalization_factor))then
      normalization_factor = sqrt(2d0*pi)/ ( &
        sqrt(1d0+1d0/sigma)*sqrt((sigma+2d0)/(2d0*pi*(sigma+1d0)))&
        )
    end if
    
  end subroutine correlated_gaussian_random_number
!-----------------------------------------------------------------------------------------  
end module random_number_module
