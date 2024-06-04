SUBROUTINE xyz2ll(e, p)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: e(3)
    REAL*8, INTENT(OUT) :: p(2)

    p(1) = atan2(e(3), SQRT(e(1) * e(1) + e(2) * e(2)))
    p(2) = atan2(e(2), e(1))

END SUBROUTINE xyz2ll

SUBROUTINE ll2xyz(p, e)
    IMPLICIT NONE
    REAL*8, INTENT(IN) :: p(2)
    REAL*8, INTENT(OUT) :: e(3)

    e(1) = cos(p(1)) * cos(p(2))
    e(2) = cos(p(1)) * sin(p(2))
    e(3) = sin(p(1))

END SUBROUTINE ll2xyz

SUBROUTINE normalize(e)
    IMPLICIT NONE
    REAL*8, INTENT(INOUT) :: e(3)

    REAL*8 mag

    mag = sqrt(e(1)*e(1)+e(2)*e(2)+e(3)*e(3))
    e = e / mag

END SUBROUTINE normalize
    
