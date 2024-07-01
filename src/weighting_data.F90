!>
!!=============================================================================
!!  This file contains a module defining weights of the weighted norm interpolation:
!!  weightedNormIntp.
!!  Created by Yuanfu Xie on April 27, 2024
!!=============================================================================

MODULE weights_m
  DOUBLE PRECISION, PARAMETER :: &
    ! for G5, the value 1.4D0 results in larger streamfunction trunction errors: 7620.0995279026420
    ! but smaller Jacobian value: 4.8E-21
    ! sigma(12) = (/1.8D0, 1.7D0, 1.6D0, 1.5D0, 1.4D0, 1.3D0, 1.2D0, 1.1D0, 1.0D0, 0.9D0, 0.8D0, 0.7D0/)

    ! for G5, the value 0.18D-1 results in smaller streamfunction trunction errors: 1946.1462206929832
    ! but larger Jacobian value: 1.16E-20
    !
    ! See my note: "A QR Factorization Based Interpolation Scheme" under ${HOME}/developments/models/doc
    !
    ! But this weighting changes the other operators, causing streamfunction with NaN.
    ! sigma(12) = (/1.8D0, 1.7D0, 1.6D0, 1.5D0, 0.18D-1, 1.3D0, 1.2D0, 1.1D0, 1.0D0, 0.9D0, 0.8D0, 0.7D0/)
    ! sigma(12) = (/1.8D0, 1.7D0, 1.6D0, 1.5D0, 1.6D0, 1.3D0, 1.2D0, 1.1D0, 1.0D0, 0.9D0, 0.8D0, 0.7D0/)
    sigma(12) = (/1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0, 1.57D0/)
    ! sigma(12) = (/1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0, 1.2D0/)
END MODULE weights_m