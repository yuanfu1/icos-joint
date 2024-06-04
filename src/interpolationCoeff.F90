!!--------------------------------------------------------------------------------------------------
! PROJECT           : Model debugging tools and operators
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : Beta 0.0
! HISTORY           :
!   Created by Yuanfu Xie (yuanfu_xie@yahoo.com), 2021/03/17, @GBA-MWF, Shenzhen
!!--------------------------------------------------------------------------------------------------

!!===================================================================
!> @brief
!! # Calculating the interpolation coefficients of a set of gridpoints
!!
!!  *This module provides a routine for the calculation*
!!
!!  \author Yuanfu Xie
!!  \b history: Yuanfu Xie created 2021-03
!!     
!!
!!  Reference:
!!    
!!    under
!!    on my laptop*
!!
!!  \b history: Yuanfu Xie modified 2021-03-17
!!
!! @author Yuanfu Xie
!! @copyright (C) 2021 GBA-MWF, All rights reserved.
!! @test This is a beta testing
!!
!!===================================================================
MODULE interpolationCoeff_m
  USE kinds_m, ONLY : i_kind,r_kind

  IMPLICIT NONE

  PUBLIC comp_cf2021,leastNormIntp

  CONTAINS

    !>
    !!===============================================================
    !!
    !!  This routine calculates interpolation coefficients based on
    !!  the intplType for a given set of gridpoints, locPoints.
    !!
    !!  Inputs:
    !!    numPoints:    number of gridpoints, or stencils
    !!    locPoints:    gridpoint locations
    !!    intplType:    interpolation type, e.g., 1 for a second
    !!                  order interpolation of function value
    !!    printInfo:    print options
    !!
    !!  Output:
    !!    outStatus:    interpolation outStatus
    !!    outCoeffs:    interpolation coefficients
    !!    
    !!===============================================================
    !
    SUBROUTINE comp_cf2021(numPoints,locPoints,outCoeffs,outResidu, &
                          intplType,printInfo,outStatus)
      INTEGER(i_kind), INTENT(IN) :: numPoints,intplType,printInfo
      REAL(r_kind),    INTENT(IN) :: locPoints(2,numPoints)

      INTEGER(i_kind), INTENT(OUT) :: outStatus
      REAL(r_kind),    INTENT(OUT) :: outResidu,outCoeffs(numPoints)

      ! LINPACK QR variables:
      INTEGER :: ipvt(numPoints),job,info,row
      DOUBLE PRECISION :: qy(numPoints),qty(numPoints),rsd(numPoints), &
                          xkb(numPoints),amx(numPoints,numPoints), &
                          aux(numPoints),rhs(numPoints),wrk(numPoints), &
                          sol(numPoints)

      ! Local variables:
      INTEGER(i_kind) :: i,krank
      REAL(r_kind), PARAMETER :: threshold = 1.2D-16,funcBound = 5.0D2
      REAL(r_kind) :: dervBound,closestDs

      closestDs = MINVAL(locPoints(1,:)**2 + locPoints(2,:)**2)
      IF (closestDs .LT. threshold) THEN
        WRITE(*,1) closestDs
1       FORMAT('comp_cf/interpolationCoeff_m: Too close to interpolate: ', &
        D16.8)
        STOP
      END IF

      dervBound = 1.0D1/closestDs

      ! Default status:
      outStatus = 0

      ! Interpolation setting:
      amx = 0.0D0
      rhs = 0.0D0
      amx(1,:) = 1.0D0      ! First term in a Taylor expansion always 1 do i=1,np

      ! Setting number of rows:
      IF (intplType .LT. 4 .AND. intplType .GE. 1) THEN
        IF (numPoints .LT. 6) THEN
          WRITE(*,2) intplType,numPoints
2         FORMAT('comp_cf2021/interpolationCoeff_m: too few points ', &
                 ' interpolation type: ',I2,' number of points: ',I3)
          STOP
        END IF
        row = 6
        DO i=1,numPoints
          amx(2,i) = locPoints(1,i)
          amx(3,i) = locPoints(2,i)
          amx(4,i) = locPoints(1,i)**2
          amx(5,i) = locPoints(1,i)*locPoints(2,i)
          amx(6,i) = locPoints(2,i)**2
        END DO
      ELSE IF (intplType .GE. 4 .AND. intplType .LE. 6 .AND. intplType .GE. 1) THEN
        IF (numPoints .LT. 10) THEN
          WRITE(*,2) intplType,numPoints
          STOP
        END IF
        row = 10
        DO i=1,numPoints
          amx( 2,i) = locPoints(1,i)
          amx( 3,i) = locPoints(2,i)
          amx( 4,i) = locPoints(1,i)**2
          amx( 5,i) = locPoints(1,i)*locPoints(2,i)
          amx( 6,i) = locPoints(2,i)**2
          amx( 7,i) = locPoints(1,i)**3
          amx( 8,i) = locPoints(1,i)**2*locPoints(2,i)
          amx( 9,i) = locPoints(1,i)*locPoints(2,i)**2
          amx(10,i) = locPoints(2,i)**3
        END DO
      ELSE ! Undefined type
        WRITE(3,*) intplType
3       FORMAT('comp_cf/interpolationCoeff_m: undefined type: ',I2)
        STOP
      END IF

      ! Setting interpolation matrix and right hand side:
      rhs(MOD(intplType-1,3)+1) = 1.0D0

      ! QR factorization:
      ! parameters:
      job = 1   ! Pivoting 1 yes, 0 no
      ipvt = 0  ! 0: free column; 1: initial column; -1: final column
      ipvt(1:6) = 1
      CALL dqrdc(amx,numPoints,row,numPoints,aux,ipvt,wrk,job)

      ! Default rank:
      krank = row
      DO i=1,row
        IF (printInfo .GT. 0) PRINT*,'Diagonal elements: ',amx(i,i),i,row

        IF (ABS(amx(i,i)) .LT. threshold) THEN
          krank = i-1
          EXIT
        END IF
      END DO

      IF (printInfo .GT. 0) WRITE(*,4) krank,row,ipvt(1:numPoints)
4     FORMAT('Rank: ',I3,' Number of rows: ',I2' IPVT: ',14I8)

      IF (krank .LT. row) THEN
        ! Exit with report of rank deficiency:
        outStatus = 1
        ! RETURN
      END IF

      ! QR solves for the interpolation coefficients:
      job = 110
      sol = 0.0D0
      rsd = 0.0D0
      CALL dqrsl(amx,numPoints,row,krank,aux,rhs,qy,qty,sol,rsd,xkb,job,info)

      ! Check if the residue is sufficiently small as the matrix may be singular:
      outResidu = MAXVAL(ABS(rsd))

      ! Check solution:
      IF (intplType .EQ. 1 .OR. intplType .EQ. 4) THEN
        IF (MAXVAL(ABS(sol)) .GT. funcBound) outStatus = 2
      ELSE
        IF (MAXVAL(ABS(sol)) .GT. dervBound) outStatus = 2
      END IF

7     CONTINUE
      ! Pivoting the solution:
      outCoeffs = 0.0D0
      DO i=1,row
        outCoeffs(ipvt(i)) = sol(i)
      END DO

      ! Check the interpolation:
#ifdef DIAG
      ! Actually, dqrsl returns residue. There is no need to check in the following call:
      ! The following call can be used to check each of the interpolation equations:
      !CALL check_intp(outCoeffs,numPoints,locPoints)
      !stop
#endif
      
    END SUBROUTINE comp_cf2021

    SUBROUTINE check_intp(c,n,p)
      IMPLICIT NONE

      INTEGER(i_kind), INTENT(IN) :: n
      REAL(r_kind), INTENT(IN) :: c(n),p(2,n)

      ! Local variables:
      INTEGER(i_kind) :: i
      REAL(r_kind) :: amx(n,n),residues(n)

      amx = 0.0D0
      amx(1,:) = 1.0D0
      DO i=1,n
        amx( 2,i) = p(1,i)
        amx( 3,i) = p(2,i)
        amx( 4,i) = p(1,i)**2
        amx( 5,i) = p(1,i)*p(2,i)
        amx( 6,i) = p(2,i)**2
        amx( 7,i) = p(1,i)**3
        amx( 8,i) = p(1,i)**2*p(2,i)
        amx( 9,i) = p(1,i)*p(2,i)**2
        amx(10,i) = p(2,i)**3
      END DO

      ! Check interpolation:
      WRITE(*,2)
      DO i=1,10
        residues(i) = DOT_PRODUCT(amx(i,:),c)
        WRITE(*,1) i,residues(i)
1       FORMAT('|  Equation ',I2: ' Residue: ',D20.12,'  |')
      END DO
      WRITE(*,3)
2     FORMAT('+---------- Check the interpolation ----------+')
3     FORMAT('+-------- End heck the interpolation ---------+')

    END SUBROUTINE check_intp

    !>
    !!===============================================================
    !!
    !!  The leastNormIntp routine calculates a set of interpolation 
    !!  coefficients with a minimial L-2 norm, based on the intplType 
    !!  for a given set of gridpoints, locPoints.
    !!
    !!  Inputs:
    !!    numPoints:    number of gridpoints, or stencils
    !!    locPoints:    gridpoint locations
    !!    intplType:    interpolation type, e.g., 1 for a second
    !!                  order interpolation of function value
    !!    infoPrint:    print options
    !!
    !!  Output:
    !!    outStatus:    interpolation outStatus
    !!    outCoeffs:    interpolation coefficients with a least norm
    !!
    !!  Documentation:
    !!    
    !!    
    !!===============================================================
    !
    SUBROUTINE leastNormIntp(numPoints,locPoints,outCoeffs, &
                              intplType,infoPrint,outStatus)

      INTEGER(i_kind), INTENT(IN) :: numPoints,intplType,infoPrint
      REAL(r_kind),    INTENT(IN) :: locPoints(2,numPoints)
                        
      INTEGER(i_kind), INTENT(OUT) :: outStatus
      REAL(r_kind),    INTENT(OUT) :: outCoeffs(numPoints)
                        
      ! LINPACK QR variables:
      INTEGER :: ipvt(numPoints),jpvt(numPoints),job,info,row
      DOUBLE PRECISION :: qy(numPoints),qty(numPoints),rsd(numPoints), &
                          xkb(numPoints),amx(numPoints,numPoints), &
                          aux(numPoints),rhs(numPoints),wrk(numPoints), &
                          sol(numPoints)
                        
      ! Local variables:
      INTEGER(i_kind) :: i,j,k,krank
      REAL(r_kind), PARAMETER :: threshold = 1.2D-16,funcBound = 5.0D2
      REAL(r_kind) :: dervBound,closestDs

      ! Least norm variables:
      REAL(r_kind) :: anx(numPoints,numPoints),cnn(numPoints,numPoints), &
                      cnn_origin(numPoints,numPoints),rhs_tem(numPoints), &
                      xn(numPoints),gradient(numPoints)

      closestDs = MINVAL(locPoints(1,:)**2 + locPoints(2,:)**2)
      IF (intplType .NE. 1 .AND. closestDs .LT. threshold) THEN
        WRITE(*,1) closestDs
1       FORMAT('comp_cf/interpolationCoeff_m: Too close to interpolate: ', &
        D16.8)
        STOP
      END IF

      dervBound = 1.0D1
      IF (intplType .NE. 1) dervBound = dervBound/closestDs

      ! Default status:
      outStatus = 0

      ! Interpolation setting:
      amx = 0.0D0
      rhs = 0.0D0
      amx(1,:) = 1.0D0      ! First term in a Taylor expansion always 1 do i=1,np

      ! Setting number of rows:
      IF (intplType .LT. 4 .AND. intplType .GE. 1) THEN
        IF (numPoints .LT. 6) THEN
          WRITE(*,2) intplType,numPoints
2         FORMAT('leastNormIntp/interpolationCoeff_m: too few points ', &
                 ' interpolation type: ',I2,' number of points: ',I3)
          STOP
        END IF
        row = 6
        DO i=1,numPoints
          amx(2,i) = locPoints(1,i)
          amx(3,i) = locPoints(2,i)
          amx(4,i) = locPoints(1,i)**2
          amx(5,i) = locPoints(1,i)*locPoints(2,i)
          amx(6,i) = locPoints(2,i)**2
        END DO
      ELSE IF (intplType .GE. 4 .AND. intplType .LE. 6 .AND. intplType .GE. 1) THEN
        IF (numPoints .LT. 10) THEN
          WRITE(*,2) intplType,numPoints
          STOP
        END IF
        row = 10
        DO i=1,numPoints
          amx( 2,i) = locPoints(1,i)
          amx( 3,i) = locPoints(2,i)
          amx( 4,i) = locPoints(1,i)**2
          amx( 5,i) = locPoints(1,i)*locPoints(2,i)
          amx( 6,i) = locPoints(2,i)**2
          amx( 7,i) = locPoints(1,i)**3
          amx( 8,i) = locPoints(1,i)**2*locPoints(2,i)
          amx( 9,i) = locPoints(1,i)*locPoints(2,i)**2
          amx(10,i) = locPoints(2,i)**3
        END DO
      ELSE ! Undefined type
        WRITE(3,*) intplType
3       FORMAT('comp_cf/interpolationCoeff_m: undefined type: ',I2)
        STOP
      END IF

      ! Setting interpolation matrix and right hand side:
      ! PRINT*,'Right index: ',MOD(intplType-1,3)+1
      rhs = 0.0D0
      rhs(MOD(intplType-1,3)+1) = 1.0D0

      ! QR factorization:
      ! parameters:
      job = 1   ! Pivoting 1 yes, 0 no
      ipvt = 0  ! 0: free column; 1: initial column; -1: final column
      ipvt(1:6) = 1 !; ipvt(7) = 1; ipvt(10) = 1
      !ipvt = 0
      CALL dqrdc(amx,numPoints,row,numPoints,aux,ipvt,wrk,job)

      ! Default rank:
      krank = row
      DO i=1,row
        IF (infoPrint .GT. 0) WRITE(*,8) amx(i,i),i,row
 8      FORMAT('Diagonal elements: ',D20.12,' row: ',I2,' of: ',I2)

        IF (ABS(amx(i,i)) .LT. threshold) THEN
          krank = i-1
          EXIT
        END IF
      END DO

      IF (infoPrint .GT. 0) THEN
        WRITE(*,4) krank,row
4       FORMAT('Rank: ',I3,' Number of rows: ',I2)
        WRITE(*,5) ipvt(1:numPoints)
5       FORMAT('IPVT: ',6I5)
      END IF

      IF (krank .LT. row) THEN
        ! Exit with report of rank deficiency:
        outStatus = 1
        ! RETURN
      END IF

      ! Solve the least norm solution of Ax=b:

      ! 1. B^(-1) Q^T b: => sol
      job = 1100   ! sol only
      sol = 0.0D0
      qty = 0.0D0
      CALL dqrsl(amx,numPoints,row,krank,aux,rhs,qy,qty,sol,rsd,xkb,job,info)
      IF (krank .LT. row .AND. MAXVAL(ABS(qty(krank+1:row))) .GT. threshold) THEN
          WRITE(*,7) MAXVAL(ABS(qty(krank+1:row)))
7         FORMAT('The QTY has non-zero component in the null of A: ',D20.12)
          outStatus = 2
      END IF

      ! 2. B^(-1)N: => anx
      anx = 0.0D0
      DO i=krank+1,numPoints
        job = 10000 ! qy
        qy = 0.0D0
        CALL dqrsl(amx,numPoints,row,krank,aux,amx(:,i),qy,qty,anx(:,i-krank),rsd,xkb,job,info)
        job = 100 ! solution only: B^(-1) Q^T Q amx(:,i) = B^(-1)N
        CALL dqrsl(amx,numPoints,row,krank,aux,qy,qy,qty,anx(:,i-krank),rsd,xkb,job,info)
      END DO

      ! 3. N^T B^{-T) B^(-1) Q^T b:
      rhs = 0.0D0
      DO i=1,numPoints-krank
        DO j=1,krank
          rhs(i) = rhs(i) + anx(j,i)*sol(j)
        END DO
      END DO
      rhs_tem = rhs

      ! 4. (I+N^T B^(-T) B^(-1) N):
      cnn = 0.0D0
      DO i=1,numPoints-krank
        cnn(i,i) = 1.0D0  ! cnn = I
      END DO
      DO i=1,numPoints-krank
        DO j=1,numPoints-krank
          DO k=1,krank
            cnn(i,j) = cnn(i,j) + anx(k,i)*anx(k,j)
          END DO
        END DO
      END DO
      ! Save the original matrix I+N^TB^(-T) B^(-1) N:
      cnn_origin = cnn

      ! 5. QR factorization of cnn, a positive definite matrix:
      job = 1   ! Pivoting 1 yes, 0 no
      jpvt = 1  ! This matrix is positive definite, no need to pivot
      CALL dqrdc(cnn,numPoints,numPoints-krank,numPoints-krank,aux,jpvt,wrk,job)

      ! 6. Solve for x_N:
      job = 100 ! Solution only
      CALL dqrsl(cnn,numPoints,numPoints-krank,numPoints-krank,aux,rhs_tem,qy,qty,xn,rsd,xkb,job,info)

      ! 7. Solve for x_B: B^(-1) QT b - B^(-1) N xn = sol - B^(-1) N xn
      xkb = sol
      DO i=1,krank
        DO j=1,numPoints-krank
          xkb(i) = xkb(i)-anx(i,j)*xn(jpvt(j))
        END DO
      END DO

      ! Check the gradient:
      IF (infoPrint .GT. 0) THEN
        gradient = rhs
        DO i=1,numPoints-krank
          DO j=1,numPoints-krank
            gradient(i) = gradient(i)-cnn_origin(i,j)*xn(jpvt(j))
          END DO
        END DO
        WRITE(*,6) SQRT(DOT_PRODUCT(gradient(1:krank),gradient(1:krank)))
6       FORMAT('Norm of the gradient: ',D20.12)
      END IF

      outCoeffs = 0.0D0
      DO i=1,krank
        outCoeffs(ipvt(i)) = xkb(i)
      END DO
      DO i=1,numPoints-krank
        outCoeffs(ipvt(i+krank)) = xn(jpvt(i))
      END DO

      ! Checking the norm of function interpolation coefficients:
      IF (intplType .EQ. 4) THEN
        IF (DOT_PRODUCT(outCoeffs,outCoeffs) .GT. 1.0D0) THEN
          PRINT*,'Least normal coefficients of func intp greater than 1: ', &
		DSQRT(DOT_PRODUCT(outCoeffs,outCoeffs))
          STOP
        END IF
      END IF

      ! Check the interpolation:
      CALL check_intp(outCoeffs,numPoints,locPoints)
    END SUBROUTINE leastNormIntp
END MODULE interpolationCoeff_m
