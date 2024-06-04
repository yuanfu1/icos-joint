!>
!!=============================================================================
!!  Calculate the interpolation coefficients of function, normal and tangential
!!  derivatives.
!!
!!  This is an include file in precal_itmesh.F90.
!!
!!  \author Yuanfu Xie
!!  Created in Mar. 2019.
!!  Modified in May 2019 by Yuanfu Xie and Ning Wang converting the original
!!                          file into a routine.
!!                       by Yuanfu Xie making all interpolation into 3rd order
!!                          for Z-grid finite volume scheme
!!  Modified in Dec 2019 by Yuanfu Xie cleaning up unneeded changes, like the
!!                          rotation. But the debugging information is still
!!                          kept inside the DIAG blocks.
!!  Modified in Jul 2020 by Yuanfu Xie moving the rotation of coordinate from
!!                          pre_itmesh.F90 to here as many calls of this routine
!!                          so we unify the rotation in this routine.
!!                          Change xy1_ssz/stcl as IN variable from INOUT.
!!                          So far, stcl seems not being used.
!!
!!  References:
!!    "Finite Volume Scheme Based on Function Interpolation" and
!!    "A QR Factorization Based Interpolation Scheme"
!!    under developments/models/doc
!!=============================================================================
!
#define EULER_ROT

#ifdef MID_EDGE
SUBROUTINE comp_intp_coef_mid_edge (itmesh, i, j, stcl, latlon1, latlon2, &
                           mp, xy1_ssz, dbg_prt, glvl)
  USE precal_data
  IMPLICIT NONE

  TYPE(ITmeshGeoQnty), INTENT(INOUT) :: itmesh
  INTEGER, INTENT(IN) :: i, j, glvl
  LOGICAL, INTENT(IN) :: dbg_prt
  INTEGER, INTENT(INOUT) :: stcl(edge_stcl_sz)
  REAL*8, INTENT(IN) :: latlon1(2), latlon2(2), mp(2)
  REAL*8, INTENT(INOUT) :: xy1_ssz(2,edge_stcl_sz)
  
  ! Local variables:
  REAL*8, PARAMETER :: earth_radius = 6371.220D3
          
  REAL*8 :: cf(edge_stcl_sz), ll(2,edge_stcl_sz),xy(2,edge_stcl_sz)
  INTEGER :: k, ncell_vrtx1, ncell_vrtx2
  
  ! Debugging variables:
  INTEGER, PARAMETER :: i_debug = 1, j_debug = 2, mlvl = 2
  INTEGER :: iprint,iflag
  REAL*8 :: cf4(edge_stcl_sz,2)

  ! Yuanfu changes the coordinates from Ning's xy1_10
  ! as his x is set opposite to the tangential and y the normal.
  ! so x' = y and y' = -x
  DO k=1,edge_stcl_sz
    xy(1,k) =  xy1_ssz(2,k)
    xy(2,k) = -xy1_ssz(1,k)
  END DO
  
  iprint = 0
  iflag = 0
#ifdef DIAG
  IF (i .EQ. i_debug .AND. j .EQ. j_debug .AND. glvl .EQ. mlvl)  iprint = 1
#endif

  IF (dbg_prt) iprint = 1
  print*, 'iprint=', iprint
! Function value:
  CALL comp_cf4(edge_stcl_sz,xy,cf,1,iprint,iflag)

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    PRINT*,'When interpolating a function...'
    PRINT*,'Null space is not empty at cell:',i,' edge:',j
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    PRINT*,'Out of interpolation bound at: ',i,j
    STOP
  END IF

  itmesh%itpl_cf_fv(:,j,i) = cf
#ifdef DIAG
  IF (iprint .EQ. 1) THEN
    DO k=1,edge_stcl_sz
      PRINT*,'DEBUG: ', k, 'function cf= ', cf(k)
    END DO
  END IF
#endif
  
! Normal derivative:
  ! See the document: A QR Factorization Based Interpolation Scheme
  ! under developments/models/doc

  CALL comp_cf4(edge_stcl_sz,xy,cf,2,iprint,iflag)
  IF (iprint .EQ. 1) then
    ! For plotting: YX
    IF (i .EQ. i_debug .AND. j .EQ. j_debug) write(11,*) 10
    DO k = 1, edge_stcl_sz
      PRINT*, k, 'cf2= ', cf(k)
      IF (j .EQ. j_debug) write(11,*) xy(1,k),xy(2,k)
    END DO
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    PRINT*,'When interpolation normal derivative...'
    PRINT*,'Null space is not empty at cell:',i,' edge:',j
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    PRINT*,'Out of interpolation bound at: ',i,j
    STOP
  END IF

  itmesh%itpl_cf_nd(:,j,i) = cf/earth_radius
  
! For tangential derivatives, there are 3 options now:
  CALL comp_cf4(edge_stcl_sz,xy,cf,3,iprint,iflag)

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    PRINT*,'When interpolation tangential derivative...'
    PRINT*,'Null space is not empty at cell:',i,' edge:',j
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    PRINT*,'Out of interpolation bound at: ',i,j
    STOP
  END IF

  itmesh%itpl_cf_td(:,j,i) = cf/earth_radius

END SUBROUTINE comp_intp_coef_mid_edge

#else

SUBROUTINE comp_intp_coef_gqp (itmesh, i, j, stcl, &
           xy1_1_ssz, xy2_1_ssz, dbg_prt,glvl)
  USE precal_data
  USE interpolationCoeff_m

  USE weights_m, ONLY: sigma
  IMPLICIT NONE

  TYPE(ITmeshGeoQnty), INTENT(INOUT) :: itmesh
  INTEGER, INTENT(IN) :: i, j, glvl
  LOGICAL, INTENT(IN) :: dbg_prt
  INTEGER, INTENT(INOUT) :: stcl(edge_stcl_sz)
  REAL*8, INTENT(INOUT) :: xy1_1_ssz(2,edge_stcl_sz), xy2_1_ssz(2,edge_stcl_sz)
  
  ! Local variables:
  REAL*8, PARAMETER :: earth_radius = 6371.220D3
          
  REAL*8 :: cf(edge_stcl_sz), ll(2,edge_stcl_sz)
  REAL*8 :: xy1(2,edge_stcl_sz),xy2(2,edge_stcl_sz)
  INTEGER :: k, ncell_vrtx1, ncell_vrtx2,iflag
  
  ! Debugging variables:
  INTEGER, PARAMETER :: i_debug = 1, j_debug = 1, f_debug = 2, mlvl = 4
  INTEGER :: iprint
  REAL*8 :: cf4(edge_stcl_sz,2)
 
  ! Yuanfu changes the coordinates from Ning's xy1_10
  ! as his x is set opposite to the tangential and y the normal.
  ! so x' = y and y' = -x
  DO k=1,edge_stcl_sz
    xy1(1,k) =  xy1_1_ssz(2,k)
    xy1(2,k) = -xy1_1_ssz(1,k)
    xy2(1,k) =  xy2_1_ssz(2,k)
    xy2(2,k) = -xy2_1_ssz(1,k)
  END DO
  
  iprint = 0
#ifdef DIAG
  IF (i .EQ. i_debug .AND. j_debug .EQ. j_debug .AND. glvl .EQ. mlvl) THEN
    iprint = 1
  END IF
#endif
  ! Function value 1:
  CALL comp_cf4(edge_stcl_sz,xy1,cf,4,iprint,iflag)
  IF (iprint .EQ. 1 .AND. f_debug .EQ. 1) THEN ! .AND. l_debug .EQ. 1) THEN
    write(11,*) edge_stcl_sz
    DO k=1,edge_stcl_sz
      PRINT*,k, 'function cf-1= ', cf(k)
      ! Write xy for plotting:
      WRITE(11,*) xy1(1,k),xy1(2,k)
    END DO
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,210) i,j,glvl
210 FORMAT('1-2 Function: null space non-zero at ',I10,2I3)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,215) i,j,glvl,MAXVAL(ABS(cf))
215 FORMAT('1-2 Function: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
           print*,'Stencil size: ',edge_stcl_sz,dbg_prt
    STOP
  END IF
  itmesh%itpl_cf_fv(1:edge_stcl_sz,1,j,i) = cf

  ! Function value 2:
  CALL comp_cf4(edge_stcl_sz,xy2,cf,4,iprint,iflag)
  IF (iprint .EQ. 1 .AND. f_debug .EQ. 1) THEN !  .AND. l_debug .EQ. 2) THEN
    write(11,*) edge_stcl_sz
    DO k=1,edge_stcl_sz
      PRINT*,k, 'function cf-2= ', cf(k)
      ! Write xy for plotting:
      WRITE(11,*) xy2(1,k),xy2(2,k)
    END DO
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,220) i,j,glvl
220 FORMAT('2-2 Function: null space non-zero at ',I10,2I2)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,225) i,j,glvl,MAXVAL(ABS(cf))
225 FORMAT('2-2 Function: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
    STOP
  END IF
  itmesh%itpl_cf_fv(1:edge_stcl_sz,2,j,i) = cf

! Normal derivative:
  ! See the document: A QR Factorization Based Interpolation Scheme
  ! under developments/models/doc
  cf = 0.0D0
  ! CALL comp_cf4(edge_stcl_sz,xy1,cf,5,iprint,iflag)
  ! CALL leastNormIntp(edge_stcl_sz,xy1,cf,5,iprint,iflag)
  CALL weightedNormIntp(edge_stcl_sz,xy1,cf,5,sigma(glvl),iprint,iflag)

  ! Debugging
  IF (ABS(SUM(cf)-0.0D0) .GT. 1.0D-8) THEN
    PRINT*,'Normal derivaitve 1 - Interpolation violated: ',i,glvl,SUM(cf)
    STOP
  END IF

  DO k=1,edge_stcl_sz
  IF (isNaN(cf(k))) THEN
    PRINT*,'Found a NaN: ',k
    STOP
  END IF
  END DO

  IF (iprint .EQ. 1) then
    ! For plotting: YX
    IF (iprint .EQ. 1 .AND. f_debug .EQ. 2) THEN !  .AND. l_debug .EQ. 1) THEN
      write(11,*) edge_stcl_sz,' Norm 1 ',glvl
      DO k = 1, edge_stcl_sz
        PRINT*, k, 'Normal cf2-1= ', cf(k),j_debug
        ! Write xy for plotting:
        WRITE(11,*) xy1(1,k),xy1(2,k)
      END DO
    END IF
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,230) i,j,glvl
230 FORMAT('1-2 normal derivative: null space non-zero at ',I10,2I2)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,235) i,j,glvl,MAXVAL(ABS(cf))
235 FORMAT('1-2 Normal: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
    STOP
  END IF
  itmesh%itpl_cf_nd(1:edge_stcl_sz,1,j,i) = cf/earth_radius

  cf = 0.0D0
  ! CALL comp_cf4(edge_stcl_sz,xy2,cf,5,iprint,iflag)
  ! CALL leastNormIntp(edge_stcl_sz,xy2,cf,5,iprint,iflag)
  CALL weightedNormIntp(edge_stcl_sz,xy2,cf,5,sigma(glvl),iprint,iflag)

  ! Debugging
  IF (ABS(SUM(cf)-0.0D0) .GT. 1.0D-8) THEN
    PRINT*,'Normal derivaitve 2 - Interpolation violated: ',i,glvl,SUM(cf)
    STOP
  END IF

  DO k=1,edge_stcl_sz
  IF (isNaN(cf(k))) THEN
    PRINT*,'Found a nan: ',k
    STOP
  END IF
  END DO

  IF (iprint .EQ. 1) then
    ! For plotting: YX
    IF (iprint .EQ. 1 .AND. f_debug .EQ. 2) THEN !  .AND. l_debug .EQ. 2) THEN
      write(11,*) edge_stcl_sz
      DO k = 1, edge_stcl_sz
        PRINT*, k, 'Normal cf2-2= ', cf(k)
        ! Write xy for plotting:
        WRITE(11,*) xy2(1,k),xy2(2,k)
      END DO
    END IF
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,231) i,j,glvl
231 FORMAT('2-2 normal derivative: null space non-zero at ',I10,2I2)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,232) i,j,glvl,MAXVAL(ABS(cf))
232 FORMAT('2-2 Normal: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
    STOP
  END IF
  itmesh%itpl_cf_nd(1:edge_stcl_sz,2,j,i) = cf/earth_radius

! For tangential derivatives, there are 2 needded to interpolate:
  ! First point:
  CALL comp_cf4(edge_stcl_sz,xy1,cf,6,iprint,iflag)
  IF (iprint .EQ. 1) then
    ! For plotting: YX
    IF (iprint .EQ. 1 .AND. f_debug .EQ. 3) THEN !  .AND. l_debug .EQ. 1) THEN
      write(11,*) edge_stcl_sz
      DO k = 1, edge_stcl_sz
        PRINT*, k, 'Tangential cf2-1= ', cf(k)
        WRITE(11,*) xy1(1,k),xy1(2,k)
      END DO
    END IF
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,240) i,j,glvl
240 FORMAT('1-2 tangential derivative: null space non-zero at ',I10,2I2)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,241) i,j,glvl,MAXVAL(ABS(cf))
241 FORMAT('1-2 Tangential: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
    STOP
  END IF

  itmesh%itpl_cf_td(1:edge_stcl_sz,1,j,i) = cf/earth_radius

  ! Second point:
  CALL comp_cf4(edge_stcl_sz,xy2,cf,6,iprint,iflag)
  IF (iprint .EQ. 1) then
    ! For plotting: YX
    IF (iprint .EQ. 1 .AND. f_debug .EQ. 3) THEN !  .AND. l_debug .EQ. 2) THEN
      write(11,*) edge_stcl_sz
      DO k = 1, edge_stcl_sz
        PRINT*, k, 'Tangential cf2-2= ', cf(k)
        WRITE(11,*) xy2(1,k),xy2(2,k)
        !write(*,350) maxval(ABS(cf(:))),iflag, &
        !1.0D0/MINVAL(ABS(xy2(1,:))),1.0D0/MINVAL(ABS(xy2(2,:)))
350     format('Max cf: ',D20.12,' iflag: ',I2,' bound: ',2D20.12)
      END DO
    END IF
  END IF

  ! Check if the null space is empty:
  IF (iflag .EQ. 1) THEN
    WRITE(*,250) i,j,glvl
250 FORMAT('2-2 tangential derivative: null space non-zero at ',I10,2I2)
    STOP
  ELSE IF (iflag .EQ. 2) THEN
    WRITE(*,251) i,j,glvl,MAXVAL(ABS(cf))
251 FORMAT('2-2 Tangential: interpolation coefficients too large: ', &
           I10,2I3,' cf val: ',D20.12)
    STOP
  END IF

  itmesh%itpl_cf_td(1:edge_stcl_sz,2,j,i) = cf/earth_radius

END SUBROUTINE comp_intp_coef_gqp


SUBROUTINE comp_intp_coef (itmesh, ipn, j, stcl, stcl_sz, xy1_1_ssz, dbg_prt, glvl)
  USE precal_data
  USE interpolationCoeff_m

  USE weights_m, ONLY: sigma

  IMPLICIT NONE

  TYPE(ITmeshGeoQnty), INTENT(INOUT) :: itmesh
  INTEGER, INTENT(IN) :: ipn, j, stcl_sz, glvl
  LOGICAL, INTENT(IN) :: dbg_prt
  INTEGER, INTENT(INOUT) :: stcl(stcl_sz)
  REAL*8, INTENT(INOUT) :: xy1_1_ssz(2,stcl_sz)
  
  ! Local variables:
  ! Yuanfu adds num_pts to control interpolation points to be use:
  INTEGER :: num_pts
  REAL*8, PARAMETER :: earth_radius = 6371.220D3
          
  REAL*8 :: cf(stcl_sz), residue
  
  ! Debugging variables:
  INTEGER, PARAMETER :: ipn_debug = 1081, & ! Vertex to be debugged
                          j_debug = 1,    & ! j-th G-L along its emitting edge
                        lvl_debug = 4       ! G level to debug
  INTEGER :: iprint, k, iflag
  
  iprint = 0
!#ifdef DIAG
  IF (ipn .EQ. ipn_debug .AND. j .EQ. j_debug .AND. glvl .EQ. lvl_debug) iprint = 1
  if (iprint .eq. 1) then
    write(*,1) ipn,j,glvl,dbg_prt
1   format('Diagnose vertex: ',I6,' edge: ',I2,' glvl: ',I2,' dbg_prt: ',L)
    write(*,2) stcl
2   format('Stencil indices: ',20I6)
  end if
  
!#endif

! Function value:
  ! Switch to the new least norm interpolation scheme:
  ! CALL comp_cf4(stcl_sz,xy1_1_ssz, cf, 4, iprint, iflag)
  ! CALL leastNormIntp(stcl_sz,xy1_1_ssz,cf,4,iprint,iflag)
  ! IF (ipn .EQ. 8721 .AND. j .EQ. 7 .AND. glvl .EQ. 5) THEN
  ! IF (ipn .EQ. 2328 .AND. j .EQ. 5 .AND. glvl .EQ. 5) THEN
  IF (ipn .EQ. 8976 .AND. glvl .EQ. 5) THEN
    iprint = 1 ! Test the least square weights
    print*,'IPRINT 1: ',ipn,j,glvl,sigma(glvl)
  END IF

  ! Use 10 or 12 points instead of 20 or 24 2024-05-02 by Yuanfu
  ! num_pts = 2*itmesh%num_tri(ipn) 
  num_pts = stcl_sz

  ! CALL weightedNormIntp(stcl_sz,xy1_1_ssz,cf,4,sigma(glvl),iprint,iflag)
  cf = 0.0D0  ! We may pass the first num_pts only so we have to initialize it
  CALL weightedNormIntp(num_pts,xy1_1_ssz,cf,4,sigma(glvl),iprint,iflag)
  IF (iflag .NE. 0) THEN
    PRINT*,'Cannot find a weighted norm solution: ',ipn,j,glvl
    STOP
  END IF

  IF (ABS(SUM(cf(1:num_pts))-1.0D0) .GT. 1.0D-8) THEN
    PRINT*,'Interpolation violated: ',ipn,glvl,SUM(cf(1:num_pts))
    STOP
  END IF

  itmesh%itpl_v_cf_fv(1:stcl_sz,j,ipn) = cf

  ! There is a possible bug at 1815 or 4311 as the compiled code gave different answers,
  ! some time, it went through fine but some did not for the same compiled executive file.
  ! IF (ipn .EQ. 1815 .AND. glvl .EQ. 5) THEN
  !   PRINT*,'Sum of cf: ',SUM(cf)
  !   DO k=1,24
  !     PRINT*,'Coef: ',k,cf(k)
  !   END DO
  !   !STOP
  ! END IF

  ! IF (ipn .EQ. 8721 .AND. j .EQ. 7 .AND. glvl .EQ. 5) THEN
  IF (ipn .EQ. 2328 .AND. j .EQ. 5 .AND. glvl .EQ. 5) THEN
    PRINT*,'||X||: ',DOT_PRODUCT(cf,cf),SUM(cf)
    DO k=1,4
      WRITE(*,4321) cf(1+6*(k-1):6+6*(k-1))
4321  FORMAT('IntpCoef: ',6E12.4)
    END DO
    DO k=1,24
      WRITE(*,4322) k,cf(k)
4322  FORMAT('In sequence: ',I3,E12.4)
    END DO
  END IF

#ifdef DIAG
  IF (iprint .EQ. 1) THEN
print*, "New scheme DEBUG, ipn = ", ipn
    WRITE(31,*) stcl_sz
    DO k=1,stcl_sz
      WRITE(33,33) k,xy1_1_ssz(:,k),cf(k)
33    FORMAT('Coef at: ',I2,' XY: ',2D20.12,' F: ',D20.12)
      write(31,*) xy1_1_ssz(:,k)
    END DO
  END IF
#endif
! print*,'XXYYFF: ',stcl_sz,ipn,j,ipn_debug,j_debug,SQRT(DOT_PRODUCT(cf(1:stcl_sz),cf(1:stcl_sz)))

! Check if the null space is empty:
IF (iflag .EQ. 1) THEN
  WRITE(*,250) ipn,j,glvl
250 FORMAT('Function: null space non-zero at ',I10,2I2)
  STOP
ELSE IF (iflag .EQ. 2) THEN
  WRITE(*,251) ipn,j,glvl,MAXVAL(ABS(cf))
251 FORMAT('Function: interpolation coefficients too large: ', &
         I10,2I3,' cf val: ',D20.12)
print*,'IPRINT: ',iprint
WRITE(11,*) stcl_sz
DO k=1,stcl_sz
  write(11,*) xy1_1_ssz(:,k)
END DO
  STOP
END IF

IF (ipn .EQ. 6 .AND. j .EQ. 1 .AND. glvl .EQ. 5) PRINT*,'NNOORRMM: ',DOT_PRODUCT(cf,cf)

  ! Check if the interpolation coefficient norm is greater than 1:
  IF (DOT_PRODUCT(cf,cf) .GT. 1.5D0) THEN
    WRITE(*,31) DOT_PRODUCT(cf,cf),ipn,j,glvl
31  FORMAT('Least norm interpolation coefficients too large: ',D20.12,' at: ',I8,2I2)
    WRITE(*,36) SUM(cf)
36  FORMAT('Sum of the coef: ',D12.4)
    WRITE(111,35) stcl_sz,xy1_1_ssz(1,:)
    WRITE(111,35) stcl_sz,xy1_1_ssz(2,:)
    WRITE(111,35) stcl_sz,cf(:)
35  FORMAT(I2,24D12.4)
    STOP
  END IF
  
END SUBROUTINE comp_intp_coef
#endif
