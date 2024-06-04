!>
!!=============================================================================
!!  Modified from comp_cf.F90 to calculate the interpolation coefficients of 
!!  function, normal and tangential derivative according to the cf_tp option.
!!
!!  \author Yuanfu Xie
!!  Created in Mar. 2019.
!!  Modified in May 2019 by Yuanfu Xie adding more options to calculate 3rd
!!                                        accurate function value, normal and
!!                                        tangential derivatives
!!  Modified in Aug. 2020 by Yuanfu Xie adding a iflag to show the status
!!  Modified in Oct. 2021 by Yuanfu Xie removing minimization and add a search
!!                                      a set of interpolation coefficients 
!!                                      with minimal L-2 norm.
!!      Note that: the norm of interpolation coefficients is not differentiable
!!                  function of rotation angles.
!!
!!  Routine specification:
!!  This routine computes the interpolation coefficients according
!!    to cf_tp for function value or normal derivative. Tangential
!!    derivative is assumed to a finite difference of two 3rd order
!!    function values.
!!
!!  Facts: 1. directly interpolating tangential derivatives of 10
!!            cells yields a second order accurate scheme but it
!!            cannot guarantee tangent line integral of a cell to
!!            be zero, conservation?
!!         2. To meet the requirement of tangent line integral is
!!            zero, we have to use the same stencial surrounding a
!!            vertex no matter which cell it belongs to.
!!         3. To have the same stencil for a vertex, we have to use
!!            5 or 6 cells attached to interpolate the function
!!            value at the vertex. Then the tangential derivative
!!            can be obtained by a center finite difference scheme
!!            at the edge center. However, this scheme is a first
!!            order accurate at most.
!!         4. To have the same stencil and guarantee a second order
!!            accurate scheme, we have to use 10 or 12 cell 
!!            adjacent to the 5 or 6 cells sharing one edge.
!!        Thus, we have three options interpolating tangential
!!        derivative. 1. Directly interpolation; 2. first order
!!        accurate scheme but conserve using 5/6 stencil; 3. second
!!        order accurate scheme conserving and using 10/12 stencil.
!!        These options is read in from a interpolation namelist in
!!        the calling routine but not here this routine!
!!
!!  Inputs:
!!    np:     Number of interpolation points
!!    p:      np 2D grid points assumed on a Cartesian coordinate
!!            tangent to the sphere at origin
!!
!!    NOTE!!!: For invalid points when actual points is less than
!!            np, p needs to set sufficiently large!!!
!!
!!    cf_tp:  interpolation type:
!!            1: a second order accurate function value
!!            2: a second order accurate normal derivative
!!            3: a second order accurate tangential derivative
!!            4: a third order accurate function value
!!            5: a third order accurate normal derivative
!!            6: a third order accurate tangential derivative
!!  Output:
!!    cf:     the np interpolation coefficients
!!    iflag:  status of the coefficient calculation
!!=============================================================================
!
SUBROUTINE comp_cf4(np, p, cf, cf_tp, iprint,iflag)
  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: np, &     ! Number of interpolation points
                         cf_tp     ! Interpolation type
  INTEGER, INTENT(IN) :: iprint    ! Printing option for debugging
  REAL*8,  INTENT(IN) :: p(2,np)   ! Interpolation points
  INTEGER, INTENT(OUT) :: iflag    ! 0 successful; 1 projection error;
                                   ! 2 coefficients too large error
  REAL*8, INTENT(OUT) :: cf(np)    ! Interpolation coefficients

  ! Local variables:
  DOUBLE PRECISION :: threshold,amx(np,np),rhs(np)
  INTEGER :: i,job,ipvt(np),nrow,krank,info

  ! Coefficients of derivative bound:
  DOUBLE PRECISION :: dervBound,closestDs,distances(np)

  ! QR factorization parameters:
  DOUBLE PRECISION :: aux(np),sol(np),wrk(np),rsd(np),qy(np),qty(np),xkb(np)

  ! Default returning flag:
  iflag = 0

  ! QR factorization smallest diagonal element threshold:
  job = 1  ! Pivoting 1 yes, 0 no
  ipvt = 0  ! 0: free column; 1: initial column; -1: final column
  threshold = 1.2D-16

  closestDs = MINVAL(p(1,:)**2 + p(2,:)**2)
  IF (cf_tp .NE. 1 .AND. closestDs .LT. threshold) THEN
    WRITE(*,1) closestDs,cf_tp,np
1       FORMAT('comp_cf4/interpolationCoeff_m: Too close to interpolate: ', &
    D16.8,' interpolation type:',2I2)
    STOP
  END IF

  dervBound = 1.0D1
  IF (cf_tp .NE. 1) dervBound = dervBound/closestDs

 ! Interpolation setting:
  amx = 0.0D0
  rhs = 0.0D0
  amx(1,:) = 1.0D0      ! First term in a Taylor expansion always 1

  ! Setting number of rows:
  IF (cf_tp .EQ. 1) THEN
    IF (np .LT. 3) THEN
      WRITE(*,2) cf_tp,np
2         FORMAT('comp_cf4/interpolationCoeff_m: too few points ', &
             ' interpolation type: ',I2,' number of points: ',I3)
      STOP
    END IF
    nrow = 3
    DO i=1,np
      amx(2,i) = p(1,i)
      amx(3,i) = p(2,i)
    END DO
  ELSE IF (cf_tp .LT. 4 .AND. cf_tp .GE. 2) THEN
    IF (np .LT. 6) THEN
      WRITE(*,2) cf_tp,np
      STOP
    END IF
    nrow = 6
    DO i=1,np
      amx(2,i) = p(1,i)
      amx(3,i) = p(2,i)
      amx(4,i) = p(1,i)**2
      amx(5,i) = p(1,i)*p(2,i)
      amx(6,i) = p(2,i)**2
    END DO
  ELSE IF (cf_tp .GE. 4 .AND. cf_tp .LE. 6 .AND. cf_tp .GE. 1) THEN
    IF (np .LT. 10) THEN
      WRITE(*,2) cf_tp,np
      STOP
    END IF
    nrow = 10
    DO i=1,np
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
    ipvt = 0 ! 0: free column; 1: initial column; -1: final column
    ipvt(1:6) = 1
  ELSE ! Undefined type
    WRITE(3,*) cf_tp
3       FORMAT('comp_cf/interpolationCoeff_m: undefined type: ',I2)
    STOP
  END IF

  ! Setting the right hand side for interpolation:
  rhs(MOD(cf_tp-1,3)+1) = 1.0D0
        
  ! QR factorization:
  CALL dqrdc(amx,np,nrow,np,aux,ipvt,wrk,job)
        
  ! Check the rank:
  krank = nrow
  DO i=1,nrow
    IF (iprint .GT. 0) PRINT*,'Diagonal elements: ',amx(i,i),i,nrow

    IF (ABS(amx(i,i)) .LT. threshold) THEN
      krank = i-1
      EXIT
    END IF
  END DO

  IF (iprint .GT. 0) WRITE(*,4) krank,nrow,ipvt(1:np)
4     FORMAT('Rank: ',I3,' Number of rows: ',I2' IPVT: ',14I8)

  IF (krank .LT. nrow) THEN
    ! Exit with report of rank deficiency:
    iflag = 1
    RETURN
  END IF

  ! QR solves for the interpolation coefficients:
  job = 110
  sol = 0.0D0
  CALL dqrsl(amx,np,nrow,krank,aux,rhs,qy,qty,sol,rsd,xkb,job,info)

  ! Check solution:
  IF (cf_tp .EQ. 1 .OR. cf_tp .EQ. 4) THEN
    IF (MAXVAL(ABS(sol)) .GT. 1.0D1) iflag = 2
  ELSE
    IF (MAXVAL(ABS(sol)) .GT. dervBound) iflag = 2
  END IF

  ! Pivoting the solution:
  cf = 0.0D0
  DO i=1,nrow
    cf(ipvt(i)) = sol(i)
  END DO

END SUBROUTINE comp_cf4