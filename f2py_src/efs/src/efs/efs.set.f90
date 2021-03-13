SUBROUTINE set_weights(config_weight, energy_weight, force_weight, stress_weight)
!###########################################################
REAL(kind=DoubleReal), INTENT(IN) :: config_weight
REAL(kind=DoubleReal), INTENT(IN) :: energy_weight
REAL(kind=DoubleReal), INTENT(IN) :: force_weight
REAL(kind=DoubleReal), INTENT(IN) :: stress_weight
!###########################################################
rss_weights(1) = config_weight
rss_weights(2) = energy_weight
rss_weights(3) = force_weight
rss_weights(4) = stress_weight
END SUBROUTINE set_weights










