!======================================================
!  This subroutine computes the latitude and longitude 
!  of the middle point between two given ponits.
!
!  There are two formulae available to compute it.
!  
!  One derived from a more general m-sect formula:
!
!  xyz = sin((1-f)*theta) / sin(theta) * xyz1 +
!        sin(f*theta) /sin(theta) * xyz2 ;
!  where theta is the angle of xyz1, and xyz2.
!
!  xyz = 0.5 / sqrt[(1+dot(xyz1,xyz2))/2] * (xyz1+xyz2)
!
!  and the other one is the normalized middle point of
!  the two end points:
!
!  xyz = 0.5 * (xyz1+xyz2), xyz = xyz / sqrt(dot(xyz,xyz))
!
!  Author: Ning Wang,   March, 2006
!======================================================
SUBROUTINE middle(p1,p2,p)
     IMPLICIT NONE

     REAL*8, INTENT(IN) :: p1(2),p2(2)
     REAL*8, INTENT(OUT) :: p(2)

     REAL*8 :: xyz1(3),xyz2(3),xyz(3)
     REAL*8, PARAMETER :: pi = acos(-1.0D0)

     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     ! middle point:
     xyz = 0.5D0 * (xyz1 + xyz2)
     xyz = xyz / sqrt(DOT_PRODUCT(xyz,xyz))

     ! Convert the middle point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))) 
     p(2) = atan2(xyz(2), xyz(1)) 

END SUBROUTINE middle

SUBROUTINE pt_at_r(p1,p2,p,r)
     IMPLICIT NONE

     REAL*8, INTENT(IN) :: p1(2),p2(2)
     REAL*8, INTENT(OUT) :: p(2)
     REAL*8, INTENT(IN) :: r

     REAL*8 :: xyz1(3),xyz2(3),xyz(3)
     REAL*8 :: theta, sin_theta, sin_r_theta, sin_1mr_theta
     REAL*8, PARAMETER :: pi = acos(-1.0D0)

     ! Convert them into Cardesian coor:
     xyz1(1) = cos(p1(1)) * cos(p1(2))
     xyz1(2) = cos(p1(1)) * sin(p1(2))
     xyz1(3) = sin(p1(1))

     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     theta = acos(cos(p1(1))*cos(p2(1))*cos(p1(2)-p2(2))+sin(p1(1))*sin(p2(1)))
     sin_theta = sin(theta)
     ! point at r (0 <= r <= 1):
     sin_1mr_theta = sin((1.0-r)*theta); sin_r_theta = sin(r*theta)
     xyz = sin_1mr_theta*xyz1 + sin_r_theta*xyz2
     xyz = xyz / sin_theta 

     ! Convert the point to lat/lon coor:
     p(1) = atan2(xyz(3), sqrt(xyz(1) * xyz(1) + xyz(2) * xyz(2))) 
     p(2) = atan2(xyz(2), xyz(1)) 

END SUBROUTINE pt_at_r

SUBROUTINE sym_ext(p1,p2,p)
     IMPLICIT NONE

     REAL*8, INTENT(IN) :: p(2),p2(2)
     REAL*8, INTENT(OUT) :: p1(2)

     REAL*8 :: xyz1(3),xyz2(3),xyz(3)
     REAL*8 :: cos_a

     ! Convert p and p2 into Cardesian coor:

     xyz(1) = cos(p(1)) * cos(p(2))
     xyz(2) = cos(p(1)) * sin(p(2))
     xyz(3) = sin(p(1))
     
     xyz2(1) = cos(p2(1)) * cos(p2(2))
     xyz2(2) = cos(p2(1)) * sin(p2(2))
     xyz2(3) = sin(p2(1))

     cos_a = xyz(1)*xyz2(1)+xyz(2)*xyz2(2)+xyz(3)*xyz2(3)
     xyz = xyz*cos_a

     ! the left symmetric extension point:
     xyz1 = 2.0D0 * xyz - xyz2
     xyz1 = xyz1 / sqrt(DOT_PRODUCT(xyz1,xyz1))

     ! Convert the left symmetric extension point to lat/lon coor:
     p1(1) = atan2(xyz1(3), sqrt(xyz1(1) * xyz1(1) + xyz1(2) * xyz1(2))) 
     p1(2) = atan2(xyz1(2), xyz1(1)) 

END SUBROUTINE sym_ext

REAL*8 FUNCTION gamma_v(ll)
     IMPLICIT NONE

     REAL*8, INTENT(IN) :: ll(2)
     REAL*8 :: xyz(3)
     REAL*8 :: org(2), dist2org, theta

     org = 0.0
     theta = dist2org(ll(:))  
     xyz(1) = cos(ll(1))*cos(ll(2))-1.0*cos(theta)
     xyz(2) = cos(ll(1))*sin(ll(2)) 
     xyz(3) = sin(ll(1))

     xyz = xyz / sqrt(DOT_PRODUCT(xyz,xyz))
     gamma_v = acos(xyz(2))
     IF (ll(1) < 0.0) THEN
       gamma_v = -gamma_v
     ENDIF 

END FUNCTION gamma_v

REAL*8 FUNCTION dist2org(p2)
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: p2(2)

    REAL*8 dlatov2, dlonov2, a

    dlatov2 = p2(1)/2.0
    dlonov2 = p2(2)/2.0
    a = sin(dlatov2) * sin(dlatov2) + cos(p2(1))*sin(dlonov2)*sin(dlonov2)
    dist2org = 2.0 * atan2(sqrt(a), sqrt(1.0-a))

END FUNCTION dist2org

#ifdef PROG
PROGRAM test_pt_r
  REAL*8 p1(2), p2(2), p(2), p_new(2), r(4)

  REAL*8 d2r
  INTEGER :: i 

  d2r = acos(-1.0)/180.0;
  p1(1) = 10.0*d2r; p1(2) = 20.0*d2r 
  p2(1) = 50.0*d2r; p2(2) = 20.0*d2r 
  
 ! r(1) = 0.5 - 1.0 / sqrt(3.0) / 2.0
 ! r(2) = 0.5 + 1.0 / sqrt(3.0) / 2.0
 ! r(3) = 1.0 - 1.0 / sqrt(3.0)
 ! r(4) = 1.0 / sqrt(3.0)
   r(1) = 0.2; r(2) = 0.4; r(3) = 0.6; r(4) = 0.8

  DO i = 1, 4
    CALL pt_at_r(p1,p2,p,r(i))
    PRINT*,  i, r(i), p/d2r
  ENDDO

  p(1) = 10.0*d2r; p(2) = 20.0*d2r 
  p2(1) = 10.0*d2r; p2(2) = 30.0*d2r 
  CALL sym_ext(p1, p2, p)
  PRINT*, "p1 ext =", p1/d2r

  CALL middle(p1, p2, p_new)
  PRINT*, "middle pt =", p_new/d2r

END PROGRAM test_pt_r
#endif
