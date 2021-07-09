SUBROUTINE make_nl()
!###########################################################
INTEGER(kind=StandardInteger) :: cn, an, gn, nn
REAL(kind=DoubleReal) :: r(1:3), rsq, rsq_cut, rmag
!###########################################################

nn = 1
DO cn = 1, cc
  nl_key(cn,1) = nn
  rsq_cut = rcut(cc)**2
  DO an = coords_key(cn, 1), coords_key(cn, 2)
    DO gn = halo_key(cn, 1), halo_key(cn, 2)
      IF(atom_id(an) .LT. atom_halo_id(gn))THEN
        r(:) = coords_halo(gn,:) - coords_real(an,:)
        rsq = SUM(r(:)**2)
        IF(rsq .LE. rsq_cut)THEN
          rmag = sqrt(rsq) 
          nl_atoms(nn, 1) = atom_id(an)
          nl_atoms(nn, 2) = atom_halo_id(gn)
          nl_atoms(nn, 3) = atom_type(an)
          nl_atoms(nn, 4) = atom_halo_type(gn)
          nl_r(nn) = rmag
          nl_rvec(nn,1:3) = r(1:3) / rmag
          nl_in_halo(nn) = in_halo(gn)
	        nn = nn + 1
        END IF
      END IF
    END DO
  END DO
	nl_key(cn,2) = nn - 1
	nl_key(cn,3) = nn - nl_key(cn,1)
END DO


END SUBROUTINE make_nl