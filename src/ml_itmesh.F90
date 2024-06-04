!=================================================================
!  This program uses module precal_itmesh to calculate multiple
!  resolutions (levels) of icosahedral-triangular mesh for multi-
!  grid solver.
! 
!  Ning Wang, March 2019, original version
!  
!==================================================================
PROGRAM ml_itmesh
  USE module_readnl,only: ReturnGLVL, ReturnGridType
  USE precal_data
  USE precal
  USE ml_precal_util
  IMPLICIT NONE

  TYPE(ITmeshGeoQnty), DIMENSION(8) :: itmeshes
  TYPE(ITmeshGeoQnty), DIMENSION(8) :: itmeshes2

  ! variables for namelist
  CHARACTER(len=256) :: icos_tri_grid_pc_file
  CHARACTER(len=2) :: gls
  INTEGER :: glvl,ierr

  INTEGER :: i, j, k, l, jp, s_idx, stcl_cndt
  INTEGER :: jp1, jm1, jp_p1, jp_m1, nb_sidx, nb_sidx_p1, nb_sidx_m1

  INTEGER, PARAMETER :: i_debug = 2

  CALL ReturnGLVL(glvl)
  PRINT*, 'Highest level of glvl=', glvl
  
  CALL alloc_ml_itmesh(glvl, itmeshes)
  DO i = 2, glvl
  PRINT*, 'Current level of glvl=', i
    CALL precal_itmesh(i,itmeshes(i))
  ENDDO

  PRINT*, 'Create multi-level links ...'
  DO i = 2, glvl-1
    PRINT*, 'Between level',i,'and',i+1
    CALL create_ml_link(i, itmeshes(i), itmeshes(i+1))
  END DO

  PRINT*, 'Create multi-level interpolation stencil ...'
  DO i = 2, glvl-1
    CALL create_ml_c2f_intp_stcl(i, itmeshes(i), itmeshes(i+1))
  ENDDO

  PRINT*, 'Create multi-level interpolation coeffs ...'
  DO i = 2, glvl-1
    CALL comp_c2f_intp_coff(i, itmeshes(i), itmeshes(i+1))
  ENDDO

  PRINT*, 'Write multi-level data file ...'
  CALL write_ml_data(glvl, itmeshes)
!  CALL alloc_ml_itmesh(glvl, itmeshes2)
!  CALL read_ml_data(glvl, itmeshes2)

END PROGRAM ml_itmesh

SUBROUTINE create_ml_link(glvl, itmesh_i, itmesh_ip1)
  USE precal_data
  USE kdt, ONLY: init_kd_tree
  IMPLICIT NONE
  TYPE(ITmeshGeoQnty) :: itmesh_i, itmesh_ip1
  INTEGER :: glvl

  INTEGER :: nip_i, ntmp_i, nip_ip1, ntmp_ip1, side_len
  INTEGER, PARAMETER :: subdiv(20) = 2

  REAL*8 :: v(2,3), v_m(2,3)
  REAL, ALLOCATABLE :: llpts(:,:)
  INTEGER :: i, j, k, jp1, idx, tri_idx(4)
  INTEGER, ALLOCATABLE :: i2t(:,:), i2t_ct(:)
  INTEGER :: nn(6)

  side_len = 1
  DO i = 1, glvl
   side_len = side_len*subdiv(i) 
  ENDDO
  nip_i = 10*(side_len*side_len)+2
  ntmp_i = 2*(nip_i-2)
  side_len = side_len*subdiv(glvl+1) 
  nip_ip1 = 10*(side_len*side_len)+2
  ntmp_ip1 = 2*(nip_ip1-2)
  
  ALLOCATE(llpts(2,nip_ip1))
  llpts(:,:) = itmesh_ip1%icos_grid(:,:)
  ALLOCATE(i2t(6,nip_ip1), i2t_ct(nip_ip1))

  i2t_ct = 0; i2t = 0
  DO i = 1, ntmp_ip1
    DO j = 1, 3
      idx = itmesh_ip1%cv_vrtx(j, i) 
      i2t_ct(idx) = i2t_ct(idx) + 1
      i2t(i2t_ct(idx),idx) = i
    ENDDO
  ENDDO 

!         1
!        /\
!       /  \
!    4 /____\6
!     /\    /\
!    /  \  /  \
!   /____\/____\
!  2     5      3
!
  CALL init_kd_tree(llpts,nip_ip1,1)
  DO i = 1, ntmp_i
    DO j = 1, 3
      v(:,j) = itmesh_i%icos_grid(:,itmesh_i%cv_vrtx(j,i)) 
    ENDDO
    DO j = 1, 3
      jp1 = mod(j,3)+1
      CALL middle(v(:,j),v(:,jp1), v_m(:,j))
!print*, 'v_m', v(:,j), v(:,jp1), v_m(:,j) 
    ENDDO
    nn = 0
    DO j = 1, 3
      CALL nn_search(v(:,j), nn(j)) 
!    PRINT*, i, nn 
!print*, v(:,j), itmesh_ip1%icos_grid(:,nn(j))
    ENDDO
    DO j = 1, 3
      CALL nn_search(v_m(:,j), nn(j+3)) 
!    PRINT*, i, nn 
!print*, v_m(:,j), itmesh_ip1%icos_grid(:,nn(j+3))
    ENDDO

    CALL find_common_triangle(nn(1), nn(4), nn(6), i2t, nip_ip1, tri_idx(1))
    CALL find_common_triangle(nn(4), nn(2), nn(5), i2t, nip_ip1, tri_idx(2))
    CALL find_common_triangle(nn(6), nn(5), nn(3), i2t, nip_ip1, tri_idx(3))
    CALL find_common_triangle(nn(4), nn(6), nn(5), i2t, nip_ip1, tri_idx(4))
!    print*, 'tri_idx=', tri_idx
    itmesh_i%c_tc_idx(1:4,i) = tri_idx(1:4)
    itmesh_ip1%p_tc_idx(tri_idx(1)) = i 
    itmesh_ip1%p_tc_idx(tri_idx(2)) = i 
    itmesh_ip1%p_tc_idx(tri_idx(3)) = i 
    itmesh_ip1%p_tc_idx(tri_idx(4)) = i 
  ENDDO
END SUBROUTINE create_ml_link  

SUBROUTINE find_common_triangle(ipn1, ipn2, ipn3, i2t, nip, tmesh_idx)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ipn1, ipn2, ipn3, nip, i2t(6,nip)
  INTEGER, INTENT(OUT) :: tmesh_idx
 
  INTEGER :: i, j, k, comm(2)
  
  k = 0
  DO i = 1, 6
    DO j = 1, 6
      IF (i2t(i, ipn1) == i2t(j, ipn2)) THEN
        k = k + 1
        comm(k) = i2t(i, ipn1)
        EXIT
      END IF  
    ENDDO
  ENDDO

  DO i = 1, 2
    DO j = 1, 6
      IF (comm(i) == i2t(j, ipn3)) THEN
        tmesh_idx = comm(i)
        RETURN   
      ENDIF
    ENDDO
  ENDDO
  
END SUBROUTINE find_common_triangle

SUBROUTINE create_ml_c2f_intp_stcl(glvl, itmesh_i, itmesh_ip1)
  USE precal_data
  IMPLICIT NONE
  TYPE(ITmeshGeoQnty) :: itmesh_i, itmesh_ip1
  INTEGER :: glvl
 
  INTEGER :: i, j, side_len, nip_i, ntmp_i
  INTEGER :: stcl(4), f_idx(4)
  INTEGER, PARAMETER :: subdiv(20) = 2

  side_len = 1
  DO i = 1, glvl
   side_len = side_len*subdiv(i) 
  ENDDO
  nip_i = 10*(side_len*side_len)+2
  ntmp_i = 2*(nip_i-2)
  
  DO i = 1, ntmp_i
    stcl(1) = i; stcl(2:4) = itmesh_i%adjnb(1:3,i)
    f_idx(1:4) = itmesh_i%c_tc_idx(1:4, i)
    itmesh_ip1%c2f_intp_stcl(:,f_idx(1)) = stcl(1:4)
    itmesh_ip1%c2f_intp_stcl(:,f_idx(2)) = stcl(1:4)
    itmesh_ip1%c2f_intp_stcl(:,f_idx(3)) = stcl(1:4)
    itmesh_ip1%c2f_intp_stcl(:,f_idx(4)) = stcl(1:4)
  ENDDO

END SUBROUTINE create_ml_c2f_intp_stcl  

SUBROUTINE comp_c2f_intp_coff(glvl, itmesh_i, itmesh_ip1)
  USE precal_data
  IMPLICIT NONE
  TYPE(ITmeshGeoQnty) :: itmesh_i, itmesh_ip1
  INTEGER :: glvl

  INTEGER :: i, j, side_len, nip_i, ntmp_i, nip_ip1, ntmp_ip1
  INTEGER :: stcl(4) 
  INTEGER, PARAMETER :: subdiv(20) = 2
  REAL*8 :: mp(2), ll(2,4), xy1_4(2,4), cf(4)
  INTEGER :: istat
 
  side_len = 1
  DO i = 1, glvl
   side_len = side_len*subdiv(i) 
  ENDDO
  nip_i = 10*(side_len*side_len)+2
  ntmp_i = 2*(nip_i-2)
  side_len = side_len*subdiv(glvl+1) 
  nip_ip1 = 10*(side_len*side_len)+2
  ntmp_ip1 = 2*(nip_ip1-2)

  DO i = 1, ntmp_ip1 
    stcl(1:4) = itmesh_ip1%c2f_intp_stcl(:,i) 
  ! compute the xy for the coarse grid cell centers
    mp(:) = itmesh_ip1%center(:,i)
    DO j = 1, 4
      ll(:,j) = itmesh_i%center(:, stcl(j))
    ENDDO
  ! translate to the required location 
    CALL EulerRotate2(mp,ll,xy1_4,4)
    CALL comp_cf4(4, xy1_4, cf, 1, 0, istat)
    itmesh_ip1%c2f_intp_cf(:,i) = cf(:) 
  ENDDO

  ! call coefficient computation routine

END SUBROUTINE comp_c2f_intp_coff

SUBROUTINE nn_search(q_ll, nn)
  USE kdt, ONLY:knn_search_ts
  IMPLICIT NONE
  REAL*8, INTENT(IN) :: q_ll(2)
  INTEGER, INTENT(OUT) :: nn(1) 

  REAL :: q_ll4(2), nn_sr, min_dists(1)
  REAL :: hp1(3), hp2(3)
  INTEGER :: num_nn, k

  hp1 = 0.0; hp2 = 0.0; nn_sr = 1.0
  k = 1; num_nn = 0
  q_ll4(:) = q_ll(:)
  CALL knn_search_ts(q_ll4, nn, min_dists, hp1, hp2, nn_sr, k, num_nn)
!  PRINT*, '--------------------------------------------'
!  PRINT*, 'min_dists, num_nn', min_dists, num_nn

END SUBROUTINE nn_search









