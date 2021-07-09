! USE CACHE
LOGICAL ::																							function_cache = .TRUE.

INTEGER(kind=StandardInteger), PARAMETER ::             mpf = 100   ! Maximum number of potential functions
INTEGER(kind=StandardInteger), PARAMETER ::             mps = 100  
INTEGER(kind=StandardInteger), PARAMETER ::             mpfs = 100  



INTEGER(kind=StandardInteger) ::                        pot_count = 0
INTEGER(kind=StandardInteger) ::                        label_max = 0
INTEGER(kind=StandardInteger) ::                        pair_max = 0
INTEGER(kind=StandardInteger) ::                        fgroup_max = 0
INTEGER(kind=StandardInteger) ::                        dens_key_offset = 0
INTEGER(kind=StandardInteger) ::                        embe_key_offset = 0
INTEGER(kind=StandardInteger) ::                        key_max = 0




INTEGER(kind=StandardInteger) ::                        pot_f_type(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                        pot_option_a(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                        pot_option_b(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                        pot_f_group(1:mpf) = 0
CHARACTER*64 ::                                         pot_f_name(1:mpf)
INTEGER(kind=StandardInteger) ::                        p_count(1:mpf) = 0
INTEGER(kind=StandardInteger) ::                        pf_count(1:mpf) = 0
REAL(kind=DoubleReal) ::                                params(1:mpf, 1:mps) = 0.0D0
REAL(kind=DoubleReal) ::                                paramsfixed(1:mpf, 1:mpfs) = 0.0D0
INTEGER(kind=StandardInteger) ::                        f_groups(100,100) = 0













