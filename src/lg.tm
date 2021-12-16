
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



void laguerre_guassian_mag_zx P((int, int, double,
        int, int, double, double, double, int, int));

:Begin:
:Function:       laguerre_guassian_mag_zx
:Pattern:        LaguerreGuassianMagZX[w_Integer, h_Integer, scale_Real, l_Integer, p_Integer, k_Real, w0_Real, x_Real, usex_Integer, clamp_Integer]
:Arguments:      { w, h, scale, l, p, k, w0, x, usex, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Real, Real, Real, Integer, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianMagZX::usage = "get mag along ZX or ZY plane"


void laguerre_guassian_phase_zx P((int, int, double,
        int, int, double, double, double, int, int));

:Begin:
:Function:       laguerre_guassian_phase_zx
:Pattern:        LaguerreGuassianPhaseZX[w_Integer, h_Integer, scale_Real, l_Integer, p_Integer, k_Real, w0_Real, x_Real, usex_Integer, clamp_Integer]
:Arguments:      { w, h, scale, l, p, k, w0, x, usex, clamp}
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Real, Real, Real, Integer, Integer}
:ReturnType:     Manual
:End:

:Evaluate: LaguerreGuassianPhaseZX::usage = "get phase along ZX or ZY plane"




void two_laguerre_guassian_phase P((int, int, double,
        int, int, int, int, double, double, double,int));
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


void two_laguerre_guassian_mag_zx P((int, int, double,
        int, int, int, int, double, double, double, int, int));
:Begin:
:Function:       two_laguerre_guassian_mag_zx
:Pattern:        TwoLaguerreGuassianMagZX[w_Integer, h_Integer, scale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, x_Real, usex_Integer, clamp_Integer]
:Arguments:      { w, h, scale, l0, p0, l1, p1, k, w0, x, usex, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianMagZX::usage = "get mag for sum of two LG modes - along the ZX or ZY plane"


void two_laguerre_guassian_phase_zx P((int, int, double,
        int, int, int, int, double, double, double, int, int));
:Begin:
:Function:       two_laguerre_guassian_phase_zx
:Pattern:        TwoLaguerreGuassianPhaseZX[w_Integer, h_Integer, scale_Real, l0_Integer, p0_Integer, l1_Integer, p1_Integer, k_Real, w0_Real, x_Real, usex_Integer, clamp_Integer]
:Arguments:      { w, h, scale, l0, p0, l1, p1, k, w0, x, usex, clamp }
:ArgumentTypes:  { Integer, Integer, Real, Integer, Integer, Integer, Integer, Real, Real, Real, Integer, Integer}
:ReturnType:     Manual
:End:

:Evaluate: TwoLaguerreGuassianPhaseZX::usage = "get phase for sum of two LG modes - along the ZX or ZY plane"


void WS_blaze P((double, double, double, int, int, int, int, int));

:Begin:
:Function:       WS_blaze
:Pattern:        Blaze[angle_Real, scale_Real, offset_Real, x_Integer, y_Integer, w_Integer, h_Integer, clamped_, data_]
:Arguments:      { angle, scale, offset, x, y, w, h, clamped, data }
:ArgumentTypes:  { Real, Real, Real, Integer, Integer, Integer, Integer, Integer, Manual  }
:ReturnType:     Manual
:End:

:Evaluate: Blaze::usage = "Blaze[angle, scale, offset, x, y, w, h, clamped, data]"

void WS_spherical P((double, double, int, int));

:Begin:
:Function:       WS_spherical
:Pattern:        SphericalPhase[scale_Real, offset_Real, w_Integer, h_Integer]
:Arguments:      { scale, offset, w, h }
:ArgumentTypes:  { Real, Real, Integer, Integer }
:ReturnType:     Manual
:End:

:Evaluate: SphericalPhase::usage = "SphericalPhase[scale, offset, w, h]"
