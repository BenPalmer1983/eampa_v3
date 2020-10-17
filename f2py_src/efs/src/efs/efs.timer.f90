SUBROUTINE start_t()
CALL CPU_TIME(timer_temp)
END SUBROUTINE start_t

SUBROUTINE end_t(time_out)
REAL(kind=DoubleReal), INTENT(OUT) :: time_out
REAL(kind=DoubleReal) :: now
CALL CPU_TIME(now)
time_out = now - timer_temp
END SUBROUTINE end_t