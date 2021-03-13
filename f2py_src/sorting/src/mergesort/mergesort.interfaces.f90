INTERFACE ms_set
  MODULE PROCEDURE set_1d_int, set_2d_int, set_1d_dp, set_2d_dp
END INTERFACE ms_set

INTERFACE ms_match_key_table
  MODULE PROCEDURE match_key_table_1d_int, match_key_table_2d_int, match_key_table_1d_dp, match_key_table_2d_dp
END INTERFACE ms_match_key_table

INTERFACE ms_return
  MODULE PROCEDURE ms_return_1d_dp, ms_return_2d_dp, ms_return_1d_int, ms_return_2d_int
END INTERFACE ms_return


