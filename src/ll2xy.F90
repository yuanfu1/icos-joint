SUBROUTINE ll2xy(lonc,latc,lon,lat,x,y)
    IMPLICIT NONE

    REAL*8, INTENT(IN) :: lonc,latc,lon,lat
    REAL*8, INTENT(OUT) :: x, y

    REAL*8 mf

    mf = 2.0/(1.0+sin(lat)*sin(latc)+cos(lat)*cos(latc)*cos(lon-lonc))

    x = mf*(cos(lat)*sin(lon-lonc))
    y = mf*((sin(lat)*cos(latc)-cos(lat)*sin(latc)*cos(lon-lonc)))

END SUBROUTINE ll2xy

