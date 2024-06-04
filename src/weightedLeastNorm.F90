!!---------------------------------------------------------------------------------------
! PROJECT           : Utility.weightedLeastNorm and Utility.weightedNormIntp
! AFFILIATION       : Guangdong-HongKong-Macao Greater Bay Area Weather Research Center 
!                     for Monitoring Warning and Forecasting (GBA-MWF)
! AUTOHR(S)         : Yuanfu Xie
! VERSION           : V 0.1
! HISTORY           :
!   Created  by Yuanfu Xie, 2023/05 @GBA-MWF, Shenzhen
!       It is under $(HOME)/developments/MOTOR-Z/Utility
!       other copies of this file are not supposed to modify but simply copy it over.
!!---------------------------------------------------------------------------------------

!> @brief
!! A direct solution of a weighted least norm problem.
!! It is mainly used for wider interpolation stencil problem that is underdetermined.
!! This package allows users to solve
!!      min 1/2 x^T D x
!!      s.t. Ax = b
!! where A is a full rank matrix.
!! @author Yuanfu Xie
!! @copyright (C) 2023 GBA-MWF, All rights reserved.
!! @warning
!! @attention
!
! Input:
!   m,n:    are integers for the dimension of A,b and D
!   D:      a symmetric matrix in R(nxn)
!   A:      a matrix in R(mxn)
!   b:       vector in R(m)
! Output:
!   x:      the weighted least norm solution in R(n)
SUBROUTINE weightedLeastNorm(m,n,D,A,b,x)
    USE kinds_m
    USE parameters_m, ONLY : machineEps

    IMPLICIT NONE

    INTEGER(i_kind) :: m,n
    REAL(r_kind), INTENT(IN) :: D(n,n),A(m,n),b(m)
    REAL(r_kind), INTENT(OUT) :: x(n)

    ! QR decomposition:
    INTEGER :: job,info,jpvt(n)
    DOUBLE PRECISION, ALLOCATABLE :: cmatrix(:,:),tmatrix(:,:),tvector(:,:)
    DOUBLE PRECISION :: aux(n),work(n),qy(n),qty(n),sol(n),xb(n),rsd(n)

    ! Local variables:
    INTEGER(i_kind) :: i,j,k

    x = 0.0D0
    ALLOCATE(cmatrix(n,n),tmatrix(n,m),tvector(n,2)) ! working matrix

    ! See the note "A QR Factorization Based Interpolation Scheme" under 
    ! $(HOME)/developments/models/doc
    ! There is a section describing the method, see A new approach ...
    ! Solve for lambda first: lamda = -(A D^-1 A^T)^-1 b
    ! D^-1 A^T:
    cmatrix = D
    job = 1 ! Pivoting
    jpvt = 0
    CALL dqrdc(cmatrix,n,n,n,aux,jpvt,work,job)
    k=1
    DO i=1,n
        IF (ABS(cmatrix(i,i)) .LT. machineEps) THEN
            EXIT
        END IF
    END DO
    IF (i .LT. n) THEN
        WRITE(*,1) i,cmatrix(i,i),machineEps
1       FORMAT('The D matrix is not full rank as required by this software, exit!', &
            I4,' R element/eps: ',2D12.4)
        STOP
    END IF

    job = 100 ! Just compute the xbb only
    rsd = 0.0D0
    DO j=1,m
        ! tmatrix holds D^-1 A^T
        tvector(:,1) = A(j,1:n)
        CALL dqrsl(cmatrix,n,n,n,aux,tvector(:,1),qy,qty,tvector(:,2),rsd,xb,job,info)
        DO i=1,n
            tmatrix(jpvt(i),j) = tvector(i,2)
        END DO
    END DO

    ! A D^-1 A^T:
    cmatrix = 0.0D0
    DO j=1,m
        DO i=1,m
            DO k=1,n
                cmatrix(i,j) = cmatrix(i,j)+A(i,k)*tmatrix(k,j)
            END DO
        END DO
    END DO

    ! (A D^-1 A^T)^-1 b:
    job = 11 ! Pivoting
    jpvt = 0
    CALL dqrdc(cmatrix,n,m,m,aux,jpvt,work,job)
    job = 100 ! Just compute the xbb only
    rsd = 0.0D0
    tvector(1:m,1) = b
    CALL dqrsl(cmatrix,n,m,m,aux,tvector(:,1),qy,qty,tvector(:,2),rsd,xb,job,info)
    tvector(:,1) = 0.0D0
    DO i=1,m
        tvector(jpvt(i),1) = tvector(i,2)
    END DO

    ! x = D^-1 A^T (A D^-1 A^T)^-1 b
    x = 0.0D0
    DO i=1,n
        DO j=1,m
            x(i) = x(i)+tmatrix(i,j)*tvector(j,1)
        END DO
    END DO

    ! Check the residual:
    ! IF (DOT_PRODUCT(rsd,rsd) .GT. 1.0D-10) THEN
    !     PRINT*,'The solution is not exact, needed to check!'
    !     STOP
    ! END IF
    ! IF (ABS(SUM(x)-1.0D0) .GT. 1.0D-8) THEN
    !     PRINT*,'Function value interpolated incorrect: ',SUM(x)
    !     STOP
    ! END IF

    DEALLOCATE(cmatrix)
END SUBROUTINE weightedLeastNorm

!> @brief
!! Preprocess of a weighted least norm problem.
!! It is mainly used for wider interpolation stencil problem that is underdetermined.
!! This package allows users to solve
!!      min 1/2 x^T D x
!!      s.t. Ax = b
!! where A is a full rank matrix.
!! @author Yuanfu Xie
!! @copyright (C) 2023 GBA-MWF, All rights reserved.
!! @warning
!! @attention
!
! Input:
!   n:      are integer of stencil points
!   x:      n stencil points
!   t:      interpolation type: 1 - 2nd order function value; 2 - 2nd normal; 3 - 2nd tangent;
!                               5 - 3rd order function value; 6 - 3rd normal; 7 - 3rd tangent;
!   p:      printing option
!   f:      flag showing the status of the calculation
! Output:
!   c:      coefficients of the interpolation requested
SUBROUTINE weightedNormIntp(n,x,c,t,sigma,iprint,f)
    USE kinds_m

    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN) :: n,t,iprint
    INTEGER(i_kind), INTENT(OUT) :: f
    REAL(r_kind), INTENT(IN) :: x(2,n)
    REAL(r_kind), INTENT(IN) :: sigma
    REAL(r_kind), INTENT(OUT) :: c(n)

    ! Local variables:
    INTEGER(i_kind), PARAMETER :: me = 10   ! Maximum number of equations
    INTEGER(i_kind) :: m,i,j,k,nx   ! nx is the number of points to interpolate. 
    REAL(r_kind) :: D(n,n),A(me,n),b(me),dis2,sigmaaa

    ! In the test_weightedLeastNorm.F90, sigma is set to 0.15*rd for triangle grids
    ! assume the grid is triangle, we can reversely calculate rd and set sigma
    ! in the functionBasedNormal_figs.pptx, x = rd/2, where rd is the length of an edge. 
    ! A Gauss-Legendre point is 1/SQRT(3)/2*rd
    ! The minimal distance d of these stencil to a Gauss-Legendre point is rd*(1-1/SQRT(3)/)/2
    ! rd = 2d/(1-1/SQRT(3)), thus, after rotation of the xy in that PPTX, the minimal distance
    ! is the y direction, x(2,:):
    ! sigma = 5.0D-1*2.0D0/(1.0D0-1.0D0/DSQRT(3.0D0))*MINVAL(ABS(x(2,:)))
    ! sigma = 1.6D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2)) ! This choice causes G5 big errors
    ! sigma = 1.57D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))   ! current best choice for G6
    !! sigma = 1.58D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))   ! current best choice for G6 *** these two lines being testing
    !! sigma = 1.d0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))
    !sigma = 1.58D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))   ! ???
    ! sigma = 1.58D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))   ! This is unstable!
    ! sigma = 1.5D0*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))   ! Try 1.5 for G7 failed at 1.57 and 1.55
    ! PRINT*,'Sigma = ',sigma

    ! On 2024-05-02, Yuanfu adds a nx to control the number of points to be use.
    ! for function values, he intends to use 10 or 12 points for a third order accurate scheme.
    nx = n
    sigmaaa = sigma*DSQRT(MINVAL(x(1,:)**2+x(2,:)**2))

    ! Calculate the weights:
    D = 0.0D0
    DO i=1,nx
        dis2 = x(1,i)**2+x(2,i)**2
        D(i,i) = EXP((dis2)/sigmaaa/sigmaaa)
        IF (iprint .EQ. 1) WRITE(*,3) i,D(i,i),1.0D0/D(i,i),dis2,sigma,t
3       FORMAT('D at: ',I3,2D12.4,' distance to itpl: ',D12.4,' sigma: ',D12.4,' t:',I2)

        IF (dis2 .GT. 1.0D0) THEN
            print*,'Found dis2 larger than 1: ',dis2
            STOP
        END IF
        IF (MINVAL(x(1,:)**2+x(2,:)**2) .GT. 1.0D0) THEN
            PRINT*,'Found min-x^2 larger than 1: ',MINVAL(x(1,:)**2+x(2,:)**2)
            STOP
        END IF
        IF (dis2 .LT. 1.0D0) THEN
        END IF
    END DO

    ! Load the interpolation coefficient matrix:
    A = 0.0D0
    b = 0.0D0
    c = 0.0D0
    A(1,1:nx) = 1.0D0
    SELECT CASE (t)
        CASE(1)
            m = 3; b(1) = 1.0D0
            DO i=1,n
                A(2,i) = x(1,i)
                A(3,i) = x(2,i)
            END DO
        CASE(2,3)
            m = 6; 
            IF (t .LE. 3) THEN
                b(t) = 1.0D0
            ELSE
                b(1) = 1.0D0
            END IF
            DO i=1,nx
                A(2,i) = x(1,i)
                A(3,i) = x(2,i)
                A(4,i) = x(1,i)*x(1,i)
                A(5,i) = x(1,i)*x(2,i)
                A(6,i) = x(2,i)*x(2,i)
            END DO
        CASE(4,5,6)
            m = 10;
            IF (t .EQ. 4) THEN
                b(1) = 1.0D0
                ! m = 6   ! Third order accurate would be sufficient for function interpolations Yuanfu 2024-05-02
                ! nx = 12
            ELSE IF (t .EQ. 5) THEN
                b(2) = 1.0D0
            ELSE
                b(3) = 1.0D0
            END IF
            DO i=1,nx
                A(2,i) = x(1,i)
                A(3,i) = x(2,i)
                A(4,i) = x(1,i)*x(1,i)
                A(5,i) = x(1,i)*x(2,i)
                A(6,i) = x(2,i)*x(2,i)
                A(7,i) = x(1,i)*x(1,i)*x(1,i)
                A(8,i) = x(1,i)*x(1,i)*x(2,i)
                A(9,i) = x(1,i)*x(2,i)*x(2,i)
                A(10,i) = x(2,i)*x(2,i)*x(2,i)
            END DO
        CASE DEFAULT
            WRITE(*,4) t
4           FORMAT('weightedNormIntp: no such interpolation type: ',I2,' Check and rereun!')
            STOP
    END SELECT

    ! Find the interpolation coefficients:
    CALL weightedLeastNorm(m,nx,D(1:nx,1:nx),A(1:m,1:nx),b,c)

    f = 0 ! Always sucessful call?
    CALL solutionCheck(m,nx,A(1:m,1:nx),b,c,f)

END SUBROUTINE weightedNormIntp

SUBROUTINE solutionCheck(m,n,A,b,c,f)
    USE kinds_m

    IMPLICIT NONE

    INTEGER(i_kind), INTENT(IN) :: m,n
    INTEGER(i_kind), INTENT(OUT) :: f
    REAL(r_kind), INTENT(IN) :: A(m,n),b(m),c(n)

    ! Local variables:
    INTEGER(i_kind) :: i,j
    REAL(r_kind) :: y(m)

    ! y = Ac:
    y = 0.0D0
    DO i=1,m
      DO j=1,n
        y(i) = y(i) + A(i,j)*c(j)
      END DO
    END DO

    ! Check:
    f = 0
    DO i=1,m
        IF (ABS(y(i)-b(i)) .GT. 1.0D-10) THEN
            PRINT*,'This is not a solution: residual - ',y(i)-b(i)
            f = 1
        END IF
    END DO
END SUBROUTINE solutionCheck
