#ifdef testprog
program tstrh
  real p(2), p1(2), p2(2)
  real d2r
  logical righthand

  d2r = acos(-1.0)/180.0

  p(1) = 60.0*d2r; p(2) = 40.0*d2r
  p1(1) = 60.0*d2r; p1(2) = 30.0*d2r
  p2(1) = 50.0*d2r; p2(2) = 40.0*d2r

  p(1) = -60.0*d2r; p(2) = -40.0*d2r
  p2(1) = -60.0*d2r; p2(2) = -30.0*d2r
  p1(1) = -50.0*d2r; p1(2) = -40.0*d2r
  if (righthand(p,p1,p2)) then
     print*, 'right hand'
  else
     print*, 'not right hand'
  endif 

end program tstrh
#endif

logical function righthand(p, p1, p2)
  real p(2), p1(2), p2(2)
  real p1_xy(2), p2_xy(2)
  real cos_d2c, cp_z
  
  cos_d2c = sin(p(1))*sin(p1(1)) + cos(p(1))*cos(p1(1))*cos(p1(2)-p(2)) 
  p1_xy(1) = (cos(p1(1))*sin(p1(2) - p(2))) / cos_d2c
  p1_xy(2) = (cos(p(1))*sin(p1(1)) - sin(p(1))*cos(p1(1))*cos(p1(2) - p(2))) / cos_d2c
   
  cos_d2c = sin(p(1))*sin(p2(1)) + cos(p(1))*cos(p2(1))*cos(p2(2)-p(2)) 
  p2_xy(1) = (cos(p2(1))*sin(p2(2) - p(2))) / cos_d2c
  p2_xy(2) = (cos(p(1))*sin(p2(1)) - sin(p(1))*cos(p2(1))*cos(p2(2) - p(2))) / cos_d2c

  cp_z = p1_xy(1)*p2_xy(2) - p1_xy(2)*p2_xy(1)
  righthand = cp_z > 0

end function righthand


SUBROUTINE cross_product(p,v,w)

IMPLICIT NONE

REAL*8, INTENT(IN) :: v(3),w(3)
REAL*8, INTENT(OUT) :: p(3)

p(1) = v(2)*w(3)-v(3)*w(2)
p(2) = v(3)*w(1)-v(1)*w(3)
p(3) = v(1)*w(2)-v(2)*w(1)

END SUBROUTINE cross_product
