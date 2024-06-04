  !>
  !!=================================================================
  !! This testing code is to test interpolation coefficients by three
  !! different approaches, least square, least norm and Gaussian 
  !! elimination to see which is varying continuously with locations
  !! of non-singular coefficient matrix to locations of singular.
  !!
  !! +----+----+----+
  !! |
  !! |
  !! +----+----+----+
  !!
  !! We set the ultimate interpolation point at the origin, (0,0).
  !! The 6 interpolation stencil are uniformly distributed around the
  !! interpolation point for a normal directional derivative pointing
  !! to the right (east).
  !!
  !!  \author Yuanfu Xie
  !!  \history 
  !!   2018-06-10 
  !!     Orignal version
  !!   2019-01-26 
  !!     Changed to a test program with a public subroutine
  !1     comp_cf(). N.W.
  !!=================================================================
  !
#ifdef PROG
PROGRAM coeff_optm

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: m=10
  DOUBLE PRECISION :: p(2,m),b(m),dx,dy,sigma,cof(m)
  INTEGER :: cof_type

  dx = 1.0D0
  dy = 1.0D0
  sigma = 1.0D0
  
  p(1,1) = -0.6D0*dx
  p(2,1) = -1.0D0*dy*sigma
  p(1,2) = -0.5D0*dx
  p(2,2) =  0.0D0*dy
  p(1,3) = -0.6D0*dx
  p(2,3) =  1.0D0*dy
  p(1,6) =  0.6D0*dx
  p(2,6) = -1.0D0*dy
  p(1,5) =  0.5D0*dx
  p(2,5) =  0.0D0*dy
  p(1,4) =  0.6D0*dx
  p(2,4) =  1.0D0*dy
  
  p(1,7) = -0.5D0*dx
  p(2,7) = -2.0D0*dy*sigma
  p(1,8) = -0.5D0*dx
  p(2,8) =  2.0D0*dy
  p(1,9) =  0.5D0*dx
  p(2,9) =  2.0D0*dy
  p(1,10) =  0.5D0*dx
  p(2,10) = -2.0D0*dy

  CALL comp_cf(p, cof, cof_type)
  PRINT*, 'cof='
  PRINT*, cof
 
END PROGRAM coeff_optm
  
#endif

SUBROUTINE comp_cf(p, cf, cf_tp)
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: p(2,10)
  REAL*8, INTENT(OUT) :: cf(10)
  INTEGER :: cf_tp

  INTEGER, PARAMETER :: m=10
  INTEGER :: i,j,ipvt(m),lda,n,job,isolve,info,k,l
  DOUBLE PRECISION :: qy(m),qty(m),rsd(m),xb(m),threshold
  
  DOUBLE PRECISION :: a(m,m),aa(m,m),b(m),dx,dy,eps(2),rcond, &
                      x(m),z(m),sigma,w(m),dis(m)

  ! Optimization variables:
  DOUBLE PRECISION :: f,g(m),gb(m),gn(m),rtb(m),xbb(m),xnb(m),rhs(m)
  
  ! LBFGS variables:
  INTEGER ,PARAMETER :: MM=10
  CHARACTER(LEN=60) :: TA,CS
  LOGICAL           :: LS(4)
  INTEGER    :: I0,IC,IP,IT,ISBMN,O,S,T,NO,ER
  INTEGER    :: IS(44)
  INTEGER, ALLOCATABLE :: NB(:),IW(:)
  DOUBLE PRECISION, ALLOCATABLE :: LB(:), UB(:)
  DOUBLE PRECISION :: FA, PG, DS(29)
  DOUBLE PRECISION,ALLOCATABLE :: WA(:)

  eps(1) =  0.0D0*dx
  eps(2) =  0.0D0*dy
  threshold = 1.0D-10
  
  n = 10
  lda = 10
  job = 1

  a = 0.0D0
  a(1,1:10) = 1.0D0
  DO i=1,n
    ipvt = 0
    a(2,i) = p(1,i)-eps(1)
    a(3,i) = p(2,i)-eps(2)
    a(4,i) = (p(1,i)-eps(1))**2
    a(5,i) = (p(1,i)-eps(1))*(p(2,i)-eps(2))
    a(6,i) = (p(2,i)-eps(2))**2
  END DO

  ! Save an undecomposed coefficient matrix:
  aa = a
  
  ! Compute distances between interpolation points to the point of interpolation
  DO i=1,n
    dis(i) = (p(1,i)-eps(1))**2+(p(2,i)-eps(2))**2
    dis(i) = 1.0D0/dis(i)**8 !*dis(i)*dis(i)
  END DO
  dis(2) = 0.0D0
  dis(5) = 0.0D0
  
  b = 0.0D0
  b(2) = 1.0D0 ! 1.0D0 ! function value
#ifdef DIAG
  DO i=1,n
    write(*,11) a(i,:)
  END DO
#endif
11 FORMAT('A: ',10e10.2)
  
  ! QR factorization:
  CALL dqrdc(a,lda,n,n,z,ipvt,w,job)
  k = n
  DO i=1,n
#ifdef DIAG
    PRINT*,'A: ',a(i,i),i
#endif
    IF (ABS(a(i,i)) .LT. threshold) THEN
      k =i-1
      EXIT
    END IF
  END DO
#ifdef DIAG
  print*,'IPVT: ',ipvt(1:10),k
#endif

  job = 1111
#ifdef DIAG
  print*,'before K = ',k,n
#endif
    
  ! Start minimization:
  
  IP=-1
  FA=1.0d-15
  PG=1.0d-20
  ISBMN = 1     ! SUBSPACE MINIMIZATION
  ALLOCATE(WA(2*MM*(n-k)+5*(n-k)+12*MM*MM+12*MM), &
           LB(n-k),UB(n-k),NB(n-k),IW(3*(n-k)),STAT=ER)
  DO i=1,n-k
    NB(i)= 0
  ENDDO
  ub(1) =  0.2D0
  lb(1) =  0.0D0
  ub(2) =  0.0D0
  lb(2) = -8.0D0
  ub(3) =  0.0D0
  lb(3) = -1.0D0
  ub(4) =  0.2D0
  lb(4) =  0.0D0

  TA='START'
  I0=0
  IC=0
    
  xnb = 0.0D0
  xnb(2) = -1.0D0
  
  ! Q^T b: saved in qty:
  CALL dqrsl(a,lda,n,n,z,b,qy,qty,xbb,rsd,xb,job,info)
  print*,'QTB: ',qty(k+1:n)
        
  PRINT*,'Minimizing......'
  DO WHILE(ta(1:2) .EQ. 'FG' .OR. ta .EQ. 'NEW_X' .OR. &
           ta .EQ. 'START')
    
    ! Minimization:
    CALL SETULB(n-k,MM,xnb,LB,UB,NB,f,g,FA,PG,WA,IW, &
                TA,IP,CS,LS,IS,DS)
    !print*,'Task: ',ta
    IF (ta(1:2) .EQ. 'FG') THEN
      ! Calculate xbb = b-C_n xnb:
      rhs = qty
      DO i=1,k
        DO j=k+1,n
          rhs(i) = rhs(i)-a(i,j)*xnb(j-k)
        END DO
      END DO
   
      ! Solve xbb by inverting R:
      DO i=k,1,-1
        xbb(i) = rhs(i)
        DO j=i+1,k
          xbb(i) = xbb(i)-a(i,j)*xbb(j)
        END DO
        xbb(i) = xbb(i)/a(i,i)
      END DO
      
      w = 0.0D0
      DO i=1,k
        w(ipvt(i)) = xbb(i)
      END DO
      DO i=k+1,n
        w(ipvt(i)) = xnb(i-k)
      END DO
#ifdef DIAG
      WRITE(*,22) w
#endif
22    FORMAT('X: ',10E12.4)
  
#ifdef DIAG
      DO i=1,n
        PRINT*,'Solution: ',w(i),i,ipvt(i)
      END DO
#endif

      DO i=1,k
        x(ipvt(i)) = xbb(i)
      END DO
      DO i=k+1,n
        x(ipvt(i)) = xnb(i-k)
      END DO
    
      ! Get function value and gradient:
      CALL funcgrad(n,k,ipvt,x,dis,f,gb,gn)
#ifdef DIAG
      write(*,13) gb(1:k)
      write(*,14) gn(1:n-k)
#endif
13    format('Gradient base: ',3e12.4)
14    format('Gradient nbas: ',3e12.4)
    
      ! Solving R^(-T)gb by forward replacing:
      DO i=1,k
        rtb(i) = gb(i)
        DO j=1,i-1
          rtb(i) = rtb(i)-a(j,i)*rtb(j)
        END DO
        rtb(i) = rtb(i)/a(i,i)
      END DO
    
      ! Multiplying by N^T from the left and plus gn:
      DO i=1,n-k
        g(i) = gn(i)
        DO j=1,k
          g(i) = g(i)-a(j,i+k)*rtb(j)
        END DO
      END DO

#ifdef DIAG
      WRITE(*,12) xnb(1:3),f,g(1:1)
#endif
12    FORMAT('X/F/G: ',10E12.4)

      ! stop
        
    END IF
      
    IF (TA(1:5).EQ.'NEW_X') THEN
      IT = IT+1
    END IF
  END DO ! End of LBFGS loop
    
#ifdef DIAG
  ! Check solution:
  rsd(1) = 0.0D0
  DO i=1,n
    rsd(1) = rsd(1)+w(i)
  END DO
  PRINT*,'Equation 1: ',rsd(1)
  
  rsd(2) = 0.0D0
  DO i=1,n
    rsd(2) = rsd(2)+w(i)*(p(1,i)-eps(1))
  END DO
  PRINT*,'Equation 2: ',rsd(2)
  
  rsd(3) = 0.0D0
  DO i=1,n
    rsd(3) = rsd(3)+w(i)*(p(2,i)-eps(2))
  END DO
  PRINT*,'Equation 3: ',rsd(3)
  
  rsd(4) = 0.0D0
  DO i=1,n
    rsd(4) = rsd(4)+w(i)*(p(1,i)-eps(1))**2
  END DO
  PRINT*,'Equation 4: ',rsd(4)
  
  rsd(5) = 0.0D0
  DO i=1,n
    rsd(5) = rsd(5)+w(i)*(p(1,i)-eps(1))*(p(2,i)-eps(2))
  END DO
  PRINT*,'Equation 5: ',rsd(5)
  
  rsd(6) = 0.0D0
  DO i=1,n
    rsd(6) = rsd(6)+w(i)*(p(2,i)-eps(2))**2
  END DO
  PRINT*,'Equation 6: ',rsd(6)
#endif
  
END SUBROUTINE comp_cf

SUBROUTINE funcgrad(n,k,ipvt,x,dis,f,gb,gn)
  INTEGER,          INTENT(IN) :: n,k,ipvt(n)
  DOUBLE PRECISION, INTENT(IN) :: x(n),dis(n)
  DOUBLE PRECISION, INTENT(OUT) :: f,gb(k),gn(n-k)
  
  f = 0.0D0
  gb = 0.0D0
  gn = 0.0D0
  
  ! New distance weighted function:
  f = 0.0D0
  DO i=1,n
    f = f+dis(i)*x(i)*x(i)
  END DO
  
  f = 0.5D0*f
  
  ! New distance weighted gradients:
  gb = 0.0D0
  gn = 0.0D0
  DO i=1,k
    gb(i) = gb(i)+dis(ipvt(i))*x(ipvt(i))
  END DO
  DO i=k+1,n
    gn(i-k) = gn(i-k)+dis(ipvt(i))*x(ipvt(i))
  END DO
END SUBROUTINE funcgrad

SUBROUTINE funcgrad1(n,k,ipvt,x,f,gb,gn)
  INTEGER,          INTENT(IN) :: n,k,ipvt(n)
  DOUBLE PRECISION, INTENT(IN) :: x(n)
  DOUBLE PRECISION, INTENT(OUT) :: f,gb(k),gn(n-k)
  
  f = 0.0D0
  gb = 0.0D0
  gn = 0.0D0
  
  f = f+(x(2)+1.0D0)**2
  f = f+(x(5)-1.0D0)**2
  f = f+(x(2)-x(5))**2
  
  DO i=1,n
    IF (i .NE. 2 .AND. i .NE. 5) f = f+x(i)*x(i)*1.0D0
  END DO
  
  f = 0.5D0*f
  
  DO i=1,k
    IF (ipvt(i) .EQ. 2) gb(i) = x(ipvt(i))+1.0D0+x(2)-x(5)
    IF (ipvt(i) .EQ. 5) gb(i) = x(ipvt(i))-1.0D0+x(5)-x(2)
    IF (ipvt(i) .NE. 2 .AND. ipvt(i) .NE. 5) &
        gb(i) = x(ipvt(i))*1.0D0
  END DO
  DO i=k+1,n
    IF (ipvt(i) .EQ. 2) gn(i-k) = x(ipvt(i))+1.0D0+x(2)-x(5)
    IF (ipvt(i) .EQ. 5) gn(i-k) = x(ipvt(i))-1.0D0+x(5)-x(2)
    IF (ipvt(i) .NE. 2 .AND. ipvt(i) .NE. 5) &
        gn(i-k) = x(ipvt(i))*1.0D0
  END DO
END SUBROUTINE funcgrad1
