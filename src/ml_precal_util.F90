MODULE ml_precal_util
  USE precal_data
  IMPLICIT NONE
CONTAINS

SUBROUTINE write_ml_data(glvl, itmeshes)
  INTEGER :: glvl
  TYPE(ITmeshGeoQnty) :: itmeshes(glvl)

  CHARACTER(len=256) :: ml_icos_tri_grid_pc_file
  CHARACTER(len=2) gls
  INTEGER :: side_len, i, j
  INTEGER, PARAMETER :: subdiv(20) = 2

  WRITE (gls, "(I2.2)") glvl
  ml_icos_tri_grid_pc_file = "ml_icos_tri_grid_pc_g" // gls //  ".dat"
  PRINT*, 'write ', trim(ml_icos_tri_grid_pc_file)
  OPEN(10,file=ml_icos_tri_grid_pc_file,form='unformatted')
  print*,'VRTX_STCL: ',itmeshes(glvl)%vrtx_stcl(:,1)
  WRITE(*,1) glvl,MAXVAL(itmeshes(glvl)%lnth),MINVAL(itmeshes(glvl)%lnth)
1 FORMAT('G',I2,' Edge max-min lengths: ',2D12.4)
  DO i = 2, glvl
    WRITE(10)itmeshes(i)%icos_grid
  ! save memory for all control-volumes
    WRITE(10)itmeshes(i)%area,itmeshes(i)%center
    WRITE(10)itmeshes(i)%cv_stcl,itmeshes(i)%cv_vrtx
    WRITE(10)itmeshes(i)%adjnb,itmeshes(i)%adjnb_si
  ! save memory for all sides of control-volumes
    WRITE(10)itmeshes(i)%lnth,itmeshes(i)%nb
    WRITE(10)itmeshes(i)%edge_stcl,itmeshes(i)%edge_vrtx
    WRITE(10)itmeshes(i)%vrtx_stcl
    WRITE(10)itmeshes(i)%num_tri
    WRITE(10)itmeshes(i)%t_vec,itmeshes(i)%n_vec
    WRITE(10)itmeshes(i)%t_vec2d,itmeshes(i)%n_vec2d
    WRITE(10)itmeshes(i)%itpl_cf_nd,itmeshes(i)%itpl_cf_td
    WRITE(10)itmeshes(i)%itpl_cf_fv
    WRITE(10)itmeshes(i)%edge_itpl_idx
    WRITE(10)itmeshes(i)%itpl_v_cf_fv
  ! save memory for parent-child triangle cell indexes
    WRITE(10)itmeshes(i)%p_tc_idx
    WRITE(10)itmeshes(i)%c_tc_idx
    WRITE(10)itmeshes(i)%c2f_intp_stcl
    WRITE(10)itmeshes(i)%c2f_intp_cf
  ENDDO
  CLOSE(10)
  PRINT*, 'Write ml precal data done.'

END SUBROUTINE write_ml_data

SUBROUTINE read_ml_data(glvl, itmeshes)
  INTEGER :: glvl
  TYPE(ITmeshGeoQnty) :: itmeshes(glvl)

  CHARACTER(len=256) :: ml_icos_tri_grid_pc_file
  CHARACTER(len=2) gls
  INTEGER :: side_len, nip, ntmp, i, j
  INTEGER, PARAMETER :: subdiv(20) = 2

  WRITE (gls, "(I2.2)") glvl
  ml_icos_tri_grid_pc_file = "ml_icos_tri_grid_pc_g" // gls //  ".dat"
  PRINT*, 'read ', trim(ml_icos_tri_grid_pc_file)
  OPEN(10,file=ml_icos_tri_grid_pc_file,status='old', form='unformatted')
  DO i = 2, glvl
    READ(10)itmeshes(i)%icos_grid
  ! save memory for all control-volumes
    READ(10)itmeshes(i)%area,itmeshes(i)%center
    READ(10)itmeshes(i)%cv_stcl,itmeshes(i)%cv_vrtx
    READ(10)itmeshes(i)%adjnb,itmeshes(i)%adjnb_si
  ! save memory for all sides of control-volumes
    READ(10)itmeshes(i)%lnth,itmeshes(i)%nb
    READ(10)itmeshes(i)%edge_stcl,itmeshes(i)%edge_vrtx
    READ(10)itmeshes(i)%vrtx_stcl
    READ(10)itmeshes(i)%num_tri
    READ(10)itmeshes(i)%t_vec,itmeshes(i)%n_vec
    READ(10)itmeshes(i)%t_vec2d,itmeshes(i)%n_vec2d
    READ(10)itmeshes(i)%itpl_cf_nd,itmeshes(i)%itpl_cf_td
    READ(10)itmeshes(i)%itpl_cf_fv
    READ(10)itmeshes(i)%edge_itpl_idx
    READ(10)itmeshes(i)%itpl_v_cf_fv
  ! save memory for parent-child triangle cell indexes
    READ(10)itmeshes(i)%p_tc_idx
    READ(10)itmeshes(i)%c_tc_idx
    READ(10)itmeshes(i)%c2f_intp_stcl
    READ(10)itmeshes(i)%c2f_intp_cf
  ENDDO
  CLOSE(10)
  PRINT*, 'Read ml precal data done.'

END SUBROUTINE read_ml_data

SUBROUTINE alloc_ml_itmesh(glvl, itmeshes)
  USE precal_data
  IMPLICIT NONE  
  INTEGER, INTENT(IN) :: glvl 
  TYPE(ITmeshGeoQnty) :: itmeshes(glvl)
   
  INTEGER :: side_len, i, j, nip, ntmp
  INTEGER, PARAMETER :: subdiv(20) = 2
  
  DO i = 2, glvl
  ! compute nip and ntmp for each level
    side_len = 1
    DO j = 1, i
      side_len = side_len*subdiv(j) 
    ENDDO
    nip = 10*(side_len*side_len)+2
    ntmp = 2*(nip-2)
    ALLOCATE(itmeshes(i)%icos_grid(2,nip))
    ALLOCATE(itmeshes(i)%vrtx_stcl(24,nip))
    ALLOCATE(itmeshes(i)%num_tri(nip))
    ALLOCATE(itmeshes(i)%itpl_v_cf_fv(24,7,nip))
  ! allocate memory for all control-volumes
    ALLOCATE(itmeshes(i)%area(ntmp),itmeshes(i)%center(2,ntmp))
    ALLOCATE(itmeshes(i)%cv_stcl(24,ntmp),itmeshes(i)%cv_vrtx(3,ntmp))  ! Yuanfu temporarily sets it to 24
    ALLOCATE(itmeshes(i)%adjnb(3,ntmp),itmeshes(i)%adjnb_si(3,ntmp))
  ! allocate memory for all sides of control-volumes
    ALLOCATE(itmeshes(i)%lnth(3,ntmp),itmeshes(i)%nb(3,ntmp))
    ALLOCATE(itmeshes(i)%edge_vrtx(2,3,ntmp))
    ALLOCATE(itmeshes(i)%edge_stcl(edge_stcl_sz,3,ntmp))
    ALLOCATE(itmeshes(i)%t_vec(3,3,ntmp),itmeshes(i)%n_vec(3,3,ntmp))
    ALLOCATE(itmeshes(i)%t_vec2d(2,3,ntmp),itmeshes(i)%n_vec2d(2,3,ntmp))
!    ALLOCATE(itmeshes(i)%edge_stcl_idx(2,2,3,ntmp)   ! vertex and edge index to stcls of right and left parts of edge  
#ifdef MID_EDGE
    ALLOCATE(itmeshes(i)%itpl_cf_nd(edge_stcl_sz,3,ntmp))
    ALLOCATE(itmeshes(i)%itpl_cf_td(edge_stcl_sz,3,ntmp))
    ALLOCATE(itmeshes(i)%itpl_cf_fv(edge_stcl_sz,3,ntmp))
#else
    ALLOCATE(itmeshes(i)%itpl_cf_nd(edge_stcl_sz,2,3,ntmp))
    ALLOCATE(itmeshes(i)%itpl_cf_td(edge_stcl_sz,2,3,ntmp))
    ALLOCATE(itmeshes(i)%itpl_cf_fv(edge_stcl_sz,2,3,ntmp))
    ALLOCATE(itmeshes(i)%edge_itpl_idx(2,4,3,ntmp)) ! 10 interp coffs, function val
#endif
  ! allocate memory for parent-child triangle cell indexes
    ALLOCATE(itmeshes(i)%p_tc_idx(ntmp))
    ALLOCATE(itmeshes(i)%c_tc_idx(4,ntmp))
    ALLOCATE(itmeshes(i)%c2f_intp_stcl(4,ntmp))
    ALLOCATE(itmeshes(i)%c2f_intp_cf(4,ntmp))
  ENDDO

END SUBROUTINE alloc_ml_itmesh

SUBROUTINE dealloc_ml_itmesh(glvl, itmeshes)
  USE precal_data
  IMPLICIT NONE  
  INTEGER, INTENT(IN) :: glvl
  TYPE(ITmeshGeoQnty) :: itmeshes(glvl)

  INTEGER :: i

  DO i = 2, glvl
    DEALLOCATE(itmeshes(i)%icos_grid)
    DEALLOCATE(itmeshes(i)%num_tri)
    DEALLOCATE(itmeshes(i)%area,itmeshes(i)%center)
    DEALLOCATE(itmeshes(i)%cv_stcl,itmeshes(i)%cv_vrtx)
    DEALLOCATE(itmeshes(i)%adjnb,itmeshes(i)%adjnb_si)
    DEALLOCATE(itmeshes(i)%lnth,itmeshes(i)%nb)
    DEALLOCATE(itmeshes(i)%edge_stcl,itmeshes(i)%edge_vrtx)
    DEALLOCATE(itmeshes(i)%vrtx_stcl)
    DEALLOCATE(itmeshes(i)%t_vec,itmeshes(i)%n_vec)
    DEALLOCATE(itmeshes(i)%t_vec2d,itmeshes(i)%n_vec2d)
    DEALLOCATE(itmeshes(i)%itpl_cf_nd,itmeshes(i)%itpl_cf_td)
    DEALLOCATE(itmeshes(i)%itpl_cf_fv)
    DEALLOCATE(itmeshes(i)%edge_itpl_idx)
    DEALLOCATE(itmeshes(i)%itpl_v_cf_fv)
    DEALLOCATE(itmeshes(i)%p_tc_idx,itmeshes(i)%c_tc_idx)
    DEALLOCATE(itmeshes(i)%c2f_intp_stcl)
    DEALLOCATE(itmeshes(i)%c2f_intp_cf)
  ENDDO

END SUBROUTINE dealloc_ml_itmesh

END MODULE ml_precal_util

