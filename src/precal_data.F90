MODULE precal_data 
  INTEGER, PARAMETER :: edge_stcl_sz = 14
  TYPE ITmeshGeoQnty
    REAL*8, ALLOCATABLE ::icos_grid(:,:) ! of the size (2, nip)
    INTEGER, ALLOCATABLE :: vrtx_stcl(:,:) ! of the size (12, nip)
    REAL, ALLOCATABLE :: itpl_v_cf_fv(:,:,:) ! 12 interp coffs for vertex, and neighboring GL pts
    INTEGER, ALLOCATABLE ::num_tri(:) ! of the size (nip)
    ! The following arrays are for tri-CV, of the size (...,ntmp)
    REAL,ALLOCATABLE :: area(:)  ! area
    REAL,ALLOCATABLE :: center(:,:) ! center, (lat/lon)
    INTEGER,ALLOCATABLE :: cv_stcl(:,:) ! 13 stcl indices (tri-CV seqnum)
    INTEGER,ALLOCATABLE :: cv_vrtx(:,:) ! 3 vertices (icos grid point seqnum) 
    INTEGER,ALLOCATABLE :: adjnb(:,:) ! 3 adjacent neighbors (tri-CV seqnum)
    INTEGER,ALLOCATABLE :: adjnb_si(:,:) ! 3 adjacent neighbors side index
    ! The following arrays are for side of tri-CV, of the size (...,3,ntmp)
    REAL, ALLOCATABLE :: lnth(:,:)        ! the length of each CV side 
    INTEGER, ALLOCATABLE :: nb(:,:)       ! the neighbor on the opposite side 
    INTEGER, ALLOCATABLE :: edge_stcl(:,:,:)   ! 10 stcl indices for each edge
    INTEGER, ALLOCATABLE :: edge_stcl_idx(:,:,:,:)   ! vertex index for right and left edge stcls 
    INTEGER, ALLOCATABLE :: edge_vrtx(:,:,:) ! 2 vertices of each CV side
    REAL, ALLOCATABLE :: t_vec(:,:,:)     ! unit tngnt vector of each side
    REAL, ALLOCATABLE :: n_vec(:,:,:)     ! unit norml vector of each
    REAL, ALLOCATABLE :: t_vec2d(:,:,:)     ! unit tngnt vector of each side
    REAL, ALLOCATABLE :: n_vec2d(:,:,:)     ! unit norml vector of each
#ifdef MID_EDGE
    REAL, ALLOCATABLE :: itpl_cf_nd(:,:,:) ! 10 interp coffs, normal drvtv  
    REAL, ALLOCATABLE :: itpl_cf_td(:,:,:) ! 10 interp coffs, tengtl drvtv
    REAL, ALLOCATABLE :: itpl_cf_fv(:,:,:) ! 10 interp coffs, function val
#else
    REAL, ALLOCATABLE :: itpl_cf_nd(:,:,:,:) ! 10 interp coffs, normal drvtv  
    REAL, ALLOCATABLE :: itpl_cf_td(:,:,:,:) ! 10 interp coffs, tengtl drvtv
    REAL, ALLOCATABLE :: itpl_cf_fv(:,:,:,:) ! 10 interp coffs, function val
    INTEGER, ALLOCATABLE :: edge_itpl_idx(:,:,:,:) ! 10 interp coffs, function val
#endif
    INTEGER, ALLOCATABLE :: p_tc_idx(:)    ! parent triangle cell index 
    INTEGER, ALLOCATABLE :: c_tc_idx(:,:)  ! child triangle cell indexes
    INTEGER, ALLOCATABLE :: c2f_intp_stcl(:,:)  ! coarse to fine resolution interpolation stencil 
    REAL, ALLOCATABLE :: c2f_intp_cf(:,:)  ! coarse to fine resolution interpolation coefficient 
    LOGICAL, ALLOCATABLE :: comp_vrtx(:,:) ! A flag to indicate the vertex needs to conduct computation 
    REAL, ALLOCATABLE :: vrtx_intp_cf(:,:) ! 24 interp coeffs for vertex, function val
  END TYPE ITmeshGeoQnty
END MODULE precal_data 

