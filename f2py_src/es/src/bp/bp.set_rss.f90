SUBROUTINE set_rss_alat(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(1) = rss_in 
END SUBROUTINE set_rss_alat

SUBROUTINE set_rss_e0(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(2) = rss_in 
END SUBROUTINE set_rss_e0

SUBROUTINE set_rss_b0(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(3) = rss_in 
END SUBROUTINE set_rss_b0

SUBROUTINE set_rss_ec(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(4) = rss_in 
END SUBROUTINE set_rss_ec

SUBROUTINE set_rss_g(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(5) = rss_in 
END SUBROUTINE set_rss_g

SUBROUTINE set_rss_e(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(6) = rss_in 
END SUBROUTINE set_rss_e

SUBROUTINE set_rss_v(rss_in)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: rss_in
!###########################################################
rss_weights(7) = rss_in 
END SUBROUTINE set_rss_v












