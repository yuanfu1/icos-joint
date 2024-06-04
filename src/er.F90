!==================================================================
! The two subroutines in this file rotate the input grid points 
! given Euler angles (alpha, beta, and gamma) - positive means
! clock wise rotation.
!
! Ning Wang, June, 2019
!
!===================================================================
!#define testProg
#ifdef testProg
  PROGRAM testER
    IMPLICIT NONE
  
    REAL*8 :: alpha, beta, gamma
    REAL*8 :: edge(2,2), grids(2, 10)
    REAL*8, PARAMETER :: d2r = acos(-1.0) / 180.0
  
    edge(1,1) = 31.0*d2r; edge(2,1) = 21.0*d2r
    edge(1,2) = 29.0*d2r; edge(2,2) = 19.0*d2r

    grids(1,1) = 25.0*d2r; grids(2,1) = 15.0*d2r
    grids(1,2) = 35.0*d2r; grids(2,2) = 25.0*d2r

    PRINT*,'before ER', grids(:,1:2)/d2r
    CALL EulerRotate1(edge, grids, 2)
    PRINT*,'after ER', grids(:,1:2)/d2r
  
  END PROGRAM testER
#endif

  SUBROUTINE EulerRotate1(edge,gps,gps_rot,n, dbg_prt)
    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: edge(2,2)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN) :: gps(2,n)
    REAL(KIND=8), INTENT(OUT) :: gps_rot(2,n)
    LOGICAL, INTENT(IN) :: dbg_prt

    REAL(KIND=8) :: latlon(2), latlon1(2), latlon2(2), mp(2)
    REAL(KIND=8) :: alpha, beta, gamma
    REAL(KIND=8) :: gamma_v, temp
    REAL(KIND=8), PARAMETER :: r2d = 180.0/acos(-1.0)
    
    INTEGER :: i, j 
    INTEGER, PARAMETER :: nsbk = 0

! calculate middle point of the edge
    latlon1(1) = edge(1,1); latlon1(2) = edge(2,1)
    latlon2(1) = edge(1,2); latlon2(2) = edge(2,2)
    CALL middle(latlon1,latlon2,mp)
! rotate the edge so that center is at (lat,lon) = (0.0, 0.0)
    alpha = mp(2); beta = mp(1); gamma = 0.0  
    CALL rotate_alpha_beta_gamma(latlon1,alpha,beta,gamma)
    IF (dbg_prt) THEN
      PRINT*, 'Origin (mp) = ', mp * r2d
    ENDIF

! calculate 'gamma' using the one edge end  
!    gamma = atan2(latlon1(1), latlon1(2))
    gamma = gamma_v(latlon1(:))
    IF (dbg_prt) THEN
       PRINT*, 'gamma = ', gamma
    ENDIF
    DO i = 1, n
      latlon(:) = gps(:,i)
      CALL rotate_alpha_beta_gamma(latlon, alpha, beta, gamma)
      gps_rot(:,i) = latlon(:)
!      print*, 'gps_rot:', i, gps_rot(:,i)*r2d
    ENDDO
! scaling back 2**(-nsbk)
    mp(:) = 0.0
    DO i = 1, n
      DO j  = 1, nsbk
         latlon(:) = gps_rot(:,i)
         CALL middle(mp,latlon,gps_rot(:, i)) 
      ENDDO
    ENDDO
    DO i = 1, n
      temp = gps_rot(1, i)
      gps_rot(1,i) = gps_rot(2,i)
      gps_rot(2,i) = temp
    ENDDO
  END SUBROUTINE eulerRotate1

  SUBROUTINE EulerRotate2(cnt,gps,gps_rot,n)
    IMPLICIT NONE

    REAL(KIND=8), INTENT(IN) :: cnt(2)
    INTEGER, INTENT(IN) :: n
    REAL(KIND=8), INTENT(IN) :: gps(2,n)
    REAL(KIND=8), INTENT(OUT) :: gps_rot(2,n)

    REAL(KIND=8) :: latlon(2), mp(2)
    REAL(KIND=8) :: alpha, beta, gamma
    REAL(KIND=8) :: temp
    
    INTEGER :: i, j 
    INTEGER, PARAMETER :: nsbk = 0

! rotate the stencil so that center is at the (lat,lon) = (0.0, 0.0)
    alpha = cnt(2); beta = cnt(1); gamma = 0.0  

    DO i = 1, n
      latlon(:) = gps(:,i)
      CALL rotate_alpha_beta_gamma(latlon, alpha, beta, gamma)
      gps_rot(:,i) = latlon(:)
    ENDDO
 ! scaling back 2**(-nsbk)
    mp(:) = 0.0
    DO i = 1, n
      DO j  = 1, nsbk
         latlon(:) = gps_rot(:,i)
         CALL middle(mp,latlon,gps_rot(:, i)) 
      ENDDO
    ENDDO
    DO i = 1, n
      temp = gps_rot(1, i)
      gps_rot(1,i) = gps_rot(2,i)
      gps_rot(2,i) = temp
    ENDDO
  END SUBROUTINE eulerRotate2

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
    sin_gamma = sin(gamma); cos_gamma = cos(gamma)
!    sin_gamma = 0.0; cos_gamma = 1.0 

    ! Convert them into Cardesian coor:
    xyz(1) = cos(latlon(1)) * cos(latlon(2))
    xyz(2) = cos(latlon(1)) * sin(latlon(2))
    xyz(3) = sin(latlon(1))
    ! rotate about Z axis, alpha degrees
    xyz1(1) = cos_alpha*xyz(1)+sin_alpha*xyz(2)
    xyz1(2) = -sin_alpha*xyz(1)+cos_alpha*xyz(2)
    xyz1(3) = xyz(3)
    ! ratate about Y' axis, beta degrees 
    xyz2(1) = cos_beta*xyz1(1)+sin_beta*xyz1(3)
    xyz2(2) = xyz1(2)
    xyz2(3) = -sin_beta*xyz1(1)+cos_beta*xyz1(3)
    ! ratate about X' axis, gamma degrees 
    xyz3(1) = xyz2(1)
    xyz3(2) = cos_gamma*xyz2(2)+sin_gamma*xyz2(3)
    xyz3(3) = -sin_gamma*xyz2(2)+cos_gamma*xyz2(3)

    ! Convert the grid point back to lat/lon coor:
    latlon(1) = atan2(xyz3(3),sqrt(xyz3(1)*xyz3(1)+xyz3(2)*xyz3(2))) 
    latlon(2) = atan2(xyz3(2),xyz3(1)) 

  END SUBROUTINE rotate_alpha_beta_gamma
    

