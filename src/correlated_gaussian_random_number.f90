!---------------------------------------------------!
! Copyright (c) 2017 Shunsuke A. Sato               !
! Released under the MIT license                    !
! https://opensource.org/licenses/mit-license.php   !
!---------------------------------------------------!
subroutine correlated_gaussian_random_number(x1,x2)
  implicit none
  real(8),parameter :: pi = 4d0*atan(1d0)
  real(8), intent(out) :: x1,x2
  real(8) :: r1,r2,tmp

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
    r1 = exp(-0.5d0*tmp**2)
    call random_number(r2)
    if(r2 < r1)exit
  end do


  return
end subroutine correlated_gaussian_random_number
