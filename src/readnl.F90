module module_readnl
 implicit none

 character(8)   :: ComputeTasks = '2'        ! Number of compute tasks; 'S' means Serial
 integer        :: glvl                      ! The grid level defined in the Makefile
 integer        :: gtype                     ! The grid type for the icosaheral grid: 0 standare recursive,
                                             ! 2 -- modified recursive, 3 -- modified nonrecursive.
 integer        :: curve
 integer        :: subDivNum(20)
 contains

subroutine Readnamelist
 namelist /Gridspec/ ComputeTasks,glvl,gtype,curve,subDivNum
 OPEN (10,file="gridspec.nl")
 READ (10,NML=Gridspec)
 close(10)
end subroutine Readnamelist

subroutine GetNprocs  (nprocs)
 integer,intent(OUT) :: nprocs
 call Readnamelist
 if(ComputeTasks=='S'.or.ComputeTasks=='s') then
   nprocs=1
 else
   read(ComputeTasks,*) nprocs
 endif
end subroutine GetNprocs

subroutine ReturnGLVL (glvlout)
 integer,intent(OUT) :: glvlout
 call Readnamelist
 glvlout = glvl
end subroutine ReturnGLVL

subroutine ReturnGridType (gtout)
 integer,intent(OUT) :: gtout
 call Readnamelist
 gtout = gtype
end subroutine ReturnGridType

subroutine ReturnSubDivNums (sdvn)
 integer,intent(OUT) :: sdvn(20)
 call Readnamelist
 sdvn = subDivNum
end subroutine ReturnSubDivNums

end module module_readnl


