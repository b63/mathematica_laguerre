
:Evaluate: LG::error = "`1`"


void laguerre_guassian_mag P((int, int, double, double,
        int, int, double, double, double, int));
:Begin:
:Function:       laguerre_guassian_mag
:Pattern:        LaguerreGuassianMag[w_Integer, h_Integer, wscale_Real, hscale_Real, l_Integer, p_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, wscale, hscale, l, p, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Real, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianMag::usage = "get mag"


void laguerre_guassian_phase P((int, int, double, double,
        int, int, double, double, double, int));
:Begin:
:Function:       laguerre_guassian_phase
:Pattern:        LaguerreGuassianPhase[w_Integer, h_Integer, wscale_Real, hscale_Real, l_Integer, p_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, wscale, hscale, l, p, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Real, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianPhase::usage = "get phase"




void two_laguerre_guassian_phase P((int, int, double, double,
        int, int, int, int, double, double, double, int));
:Begin:
:Function:       two_laguerre_guassian_phase
:Pattern:        TwoLaguerreGuassianPhase[w_Integer, h_Integer, wscale_Real, hscale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, wscale, hscale, l0, p0, l1, p1, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianPhase::usage = "get phase for sum of two LG modes"


void two_laguerre_guassian_mag P((int, int, double, double,
        int, int, int, int, double, double, double, int));
:Begin:
:Function:       two_laguerre_guassian_mag
:Pattern:        TwoLaguerreGuassianMag[w_Integer, h_Integer, wscale_Real, hscale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, wscale, hscale, l0, p0, l1, p1, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianMag::usage = "get mag for sum of two LG modes"
