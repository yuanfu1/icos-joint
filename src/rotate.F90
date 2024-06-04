  PROGRAM testER
    IMPLICIT NONE
  
    REAL*8 :: alpha, beta, gamma
    REAL*8 :: edge(2,2), grids(2, 10)
    REAL*8, PARAMETER :: d2r = acos(-1.0) / 180.0
  
    edge(1,1) = 49.0*d2r; edge(2,1) = 21.0*d2r
    edge(1,2) = 51.0*d2r; edge(2,2) = 19.0*d2r

    grids(1,1) = 49.0*d2r; grids(2,1) = 19.0*d2r
    grids(1,1) = 51.0*d2r; grids(2,1) = 21.0*d2r

    CALL EulerRotate(edge, grids, 2)
    PRINT*, grids(:,1:2)/d2r
  
  END PROGRAM testER

!==================================================================
! The two subroutines in this file rotate the input grid points 
! given Euler angles (alpha, beta, and gamma) - positive means
! clock wise rotation.
!
! Ning Wang, June, 2019
!
!===================================================================
  SUBROUTINE EulerRotate(edge,gps,n)
    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: edge(2,2)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(INOUT) :: gps(2,n)

    REAL(KIND=8) :: latlon(2), latlon1(2), latlon2(2), mp(2)
    REAL(KIND=8) :: xyz(3), xyz1(3), xyz2(3)
    REAL(KIND=8) :: alpha, beta, gamma, gamma_v
    
    INTEGER :: i

! calculate middle point of the edge
    latlon1(1) = edge(1,1); latlon1(2) = edge(2,1)
    latlon2(1) = edge(1,2); latlon2(2) = edge(2,2)
    CALL middle(latlon1,latlon2,mp)
print*, 'after middle', mp*180.0/acos(-1.0)
! rotate the edge so that center is at the (lat,lon) = (0.0, 0.0)
    alpha = mp(2); beta = mp(1); gamma = 0.0  
    CALL rotate_alpha_beta_gamma(latlon1,alpha,beta,gamma)
print*, 'after rotate_alpha... latlon1=', latlon1*180.0/acos(-1.0)
!    CALL rotate_alpha_beta_gamma(latlon2,alpha,beta,gamma)
!print*, 'after rotate_alpha... latlon2=', latlon2*180.0/acos(-1.0)

! calculate 'gamma' using the edge end one 
    gamma = atan2(latlon1(1), latlon1(2))
    gamma = gamma_v(latlon1(:))
print*, 'gamma=', gamma*180.0/acos(-1.0)
    DO i = 1, n
      latlon(1) = gps(1,i); latlon(2) = gps(2,i)
      CALL rotate_alpha_beta_gamma(latlon, alpha, beta, gamma)
print*, 'after rotate_alpha...', i
      gps(1,i) = latlon(1); gps(2,i) = latlon(2)
    ENDDO
  END SUBROUTINE eulerRotate

  SUBROUTINE rotate_alpha_beta_gamma(latlon, alpha, beta, gamma)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(INOUT) :: latlon(2)
    REAL(KIND=8), INTENT(IN) :: alpha, beta, gamma
    
    REAL(KIND=8) :: xyz(3), xyz1(3), xyz2(3), xyz3(3)
    REAL(KIND=8) :: sin_alpha, cos_alpha, sin_beta, cos_beta
    REAL(KIND=8) :: sin_gamma, cos_gamma
    REAL*8, PARAMETER :: pi = acos(-1.0)

    sin_alpha = sin(alpha); cos_alpha = cos(alpha)
    sin_beta = sin(beta); cos_beta = cos(beta)
!    sin_gamma = sin(gamma); cos_gamma = cos(gamma)
    sin_gamma = 0.0; cos_gamma = 1.0 

    ! Convert them into Cardesian coor:
    xyz(1) = cos(latlon(1)) * cos(latlon(2))
    xyz(2) = cos(latlon(1)) * sin(latlon(2))
    xyz(3) = sin(latlon(1))
    ! rotate about Z axis, alpha degrees
    xyz1(1) = cos_alpha * xyz(1) + sin_alpha * xyz(2)
    xyz1(2) = -sin_alpha * xyz(1) + cos_alpha * xyz(2)
    xyz1(3) = xyz(3)
    ! ratate about Y' axis, beta degrees 
    xyz2(1) = cos_beta * xyz1(1) + sin_beta * xyz1(3)
    xyz2(2) = xyz1(2)
    xyz2(3) = -sin_beta * xyz1(1) + cos_beta * xyz1(3)
    ! ratate about Z' axis, gamma degrees 
    xyz3(1) = cos_gamma * xyz2(1) + sin_gamma * xyz2(2)
    xyz3(2) = -sin_gamma * xyz2(1) + cos_gamma * xyz2(2)
    xyz3(3) = xyz2(3)

    ! Convert the grid point back to lat/lon coor:
    latlon(1) = atan2(xyz3(3), sqrt(xyz3(1) * xyz3(1) + xyz3(2) * xyz3(2))) 
    latlon(2) = atan2(xyz3(2), xyz3(1)) 

    IF (latlon(2) < 0.0) THEN 
      latlon(2) = latlon(2) + 2 * pi 
    END IF

  END SUBROUTINE rotate_alpha_beta_gamma
    

