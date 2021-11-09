
:Evaluate: LG::error = "`1`"


void laguerre_guassian_mag P((int, int, double,
        int, int, double, double, double, int));
:Begin:
:Function:       laguerre_guassian_mag
:Pattern:        LaguerreGuassianMag[w_Integer, h_Integer, scale_Real, l_Integer, p_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, scale, l, p, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianMag::usage = "get mag"


void laguerre_guassian_phase P((int, int, double,
        int, int, double, double, double, int));
:Begin:
:Function:       laguerre_guassian_phase
:Pattern:        LaguerreGuassianPhase[w_Integer, h_Integer, scale_Real, l_Integer, p_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, scale, l, p, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianPhase::usage = "get phase"




void two_laguerre_guassian_phase P((int, int, double,
        int, int, int, int, double, double, double, int));
:Begin:
:Function:       two_laguerre_guassian_phase
:Pattern:        TwoLaguerreGuassianPhase[w_Integer, h_Integer, scale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, scale, l0, p0, l1, p1, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianPhase::usage = "get phase for sum of two LG modes"


void two_laguerre_guassian_mag P((int, int, double,
        int, int, int, int, double, double, double, int));
:Begin:
:Function:       two_laguerre_guassian_mag
:Pattern:        TwoLaguerreGuassianMag[w_Integer, h_Integer, scale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, z_Real, clamp_Integer]
:Arguments:      { w, h, scale, l0, p0, l1, p1, k, w0, z, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianMag::usage = "get mag for sum of two LG modes"



void WS_blaze P((double, double, double, double, double, int, int, int, int, int, int));

:Begin:
:Function:       WS_blaze
:Pattern:        Blaze[angle_Real, scale_Real, offset_Real, amplitude_Real, ampPhase_Real, x_Integer, y_Integer, w_Integer, h_Integer, fullCycle_, clamped_, data_]
:Arguments:      { angle, scale, offset, amplitude, ampPhase, x, y, w, h, fullCycle, clamped, data }
:ArgumentTypes:  { Real, Real, Real, Real, Real, Integer, Integer, Integer, Integer, Integer, Integer, Manual  }
:ReturnType:     Manual
:End:

:Evaluate: Blaze::usage = "Blaze[angle, scale, amp, ampPhase, x, y, w, h, fullCycle, clamped, data]"
