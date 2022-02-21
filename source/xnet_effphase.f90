  Real(8) Function xeffc (xene)
        Use xnet_integ_phase, Only: qcap, t9me, xchem
        Use xnet_types, Only: dp
        Implicit none
        Real (dp) :: xene
        Real (dp) :: p2
        Real (dp) :: zzz, fdist
        if (xene<1.0d0) then
                xeffc=1d-99
                return
        endif

        zzz=(xene-xchem)/t9me
        if ( zzz < 36.0d0 ) then
        fdist = 1.0d0/(1.0d0+exp(zzz))
        else if ( zzz < 708.0d0 ) then
        fdist = exp (-zzz)
        else
        fdist = 0.0d0
        end if
        xeffc=xene**2*(qcap+xene)**2*fdist
 !       p2 = xene**2-1.d0
 !       xeffc=xene*dsqrt(p2)*(qcap+xene)**2*fdist
  End Function xeffc
