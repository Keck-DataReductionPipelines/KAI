import math

def refco(HM, TDK, PMB, RH, WL, PHI, TLR, EPS):
    """
    Determine the constants A and B in the atmospheric refraction
    model dZ = A tan Z + B tan**3 Z.
 
    Z is the "observed" zenith distance (i.e. affected by refraction)
    and dZ is what to add to Z to give the "topocentric" (i.e. in vacuo)
    zenith distance.
 
    Given:
      HM      d     height of the observer above sea level (metre)
      TDK     d     ambient temperature at the observer (deg K)
      PMB     d     pressure at the observer (millibar)
      RH      d     relative humidity at the observer (range 0-1)
      WL      d     effective wavelength of the source (micrometre)
      PHI     d     latitude of the observer (radian, astronomical)
      TLR     d     temperature lapse rate in the troposphere (degK/metre)
      EPS     d     precision required to terminate iteration (radian)
      
    Returned:
      REFA    d     tan Z coefficient (radian)
      REFB    d     tan**3 Z coefficient (radian)

    Called:  refco

    Notes:

    1  Typical values for the TLR and EPS arguments might be 0.0065D0 and
       1D-10 respectively.

    2  The radio refraction is chosen by specifying WL > 100 micrometres.

    3  The routine is a slower but more accurate alternative to the
       slRFCQ routine.  The constants it produces give perfect
       agreement with slRFRO at zenith distances arctan(1) (45 deg)
       and arctan(4) (about 76 deg).  It achieves 0.5 arcsec accuracy
       for ZD < 80 deg, 0.01 arcsec accuracy for ZD < 60 deg, and
       0.001 arcsec accuracy for ZD < 45 deg.

    P.T.Wallace   Starlink   3 June 1997

    Copyright (C) 1997 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    
    Migrated from SLALIB (fortran) to python:
    J. R. Lu -- 2015-05-13
    """
    # Sample zenith distances: arctan(1) and arctan(4)
    atn1 = 0.7853981633974483
    atn4 = 1.325817663668033

    # Determine refraction for the two sample zenith distances
    R1 = refro(atn1, HM, TDK, PMB, RH, WL, PHI, TLR, EPS)
    R2 = refro(atn4, HM, TDK, PMB, RH, WL, PHI, TLR, EPS)

    # Solve for refraction constants
    refa = (64.0 * R1 - R2) / 60.0
    refb = (R2 - 4.0 * R1) / 60.0

    return refa, refb

def refro(ZOBS, HM, TDK, PMB, RH, WL, PHI, TLR, EPS):
    """
    Atmospheric refraction for radio and optical/IR wavelengths.
  
    Given:
      ZOBS    d  observed zenith distance of the source (radian)
      HM      d  height of the observer above sea level (metre)
      TDK     d  ambient temperature at the observer (deg K)
      PMB     d  pressure at the observer (millibar)
      RH      d  relative humidity at the observer (range 0-1)
      WL      d  effective wavelength of the source (micrometre)
      PHI     d  latitude of the observer (radian, astronomical)
      TLR     d  temperature lapse rate in the troposphere (degK/metre)
      EPS     d  precision required to terminate iteration (radian)
  
    Returned:
      REF     d  refraction: in vacuo ZD minus observed ZD (radian)
  
    Notes:
  
    1  A suggested value for the TLR argument is 0.0065.  The
       refraction is significantly affected by TLR, and if studies
       of the local atmosphere have been carried out a better TLR
       value may be available.
  
    2  A suggested value for the EPS argument is 1D-8.  The result is
       usually at least two orders of magnitude more computationally
       precise than the supplied EPS value.
  
    3  The routine computes the refraction for zenith distances up
       to and a little beyond 90 deg using the method of Hohenkerk
       and Sinclair (NAO Technical Notes 59 and 63, subsequently adopted
       in the Explanatory Supplement, 1992 edition - see section 3.281).
  
    4  The code is a development of the optical/IR refraction subroutine
       AREF of C.Hohenkerk (HMNAO, September 1984), with extensions to
       support the radio case.  Apart from merely cosmetic changes, the
       following modifications to the original HMNAO optical/IR refraction
       code have been made:
  
       .  The angle arguments have been changed to radians.
  
       .  Any value of ZOBS is allowed (see note 6, below).
  
       .  Other argument values have been limited to safe values.
  
       .  Murray's values for the gas constants have been used
          (Vectorial Astrometry, Adam Hilger, 1983).
  
       .  The numerical integration phase has been rearranged for
          extra clarity.
  
       .  A better model for Ps(T) has been adopted (taken from
          Gill, Atmosphere-Ocean Dynamics, Academic Press, 1982).
  
       .  More accurate expressions for Pwo have been adopted
          (again from Gill 1982).
  
       .  Provision for radio wavelengths has been added using
          expressions devised by A.T.Sinclair, RGO (private
          communication 1989), based on the Essen & Froome
          refractivity formula adopted in Resolution 1 of the
          13th International Geodesy Association General Assembly
          (Bulletin Geodesique 70 p390, 1963).
  
       .  Various small changes have been made to gain speed.
  
       None of the changes significantly affects the optical/IR results
       with respect to the algorithm given in the 1992 Explanatory
       Supplement.  For example, at 70 deg zenith distance the present
       routine agrees with the ES algorithm to better than 0.05 arcsec
       for any reasonable combination of parameters.  However, the
       improved water-vapour expressions do make a significant difference
       in the radio band, at 70 deg zenith distance reaching almost
       4 arcsec for a hot, humid, low-altitude site during a period of
       low pressure.
  
    5  The radio refraction is chosen by specifying WL > 100 micrometres.
       Because the algorithm takes no account of the ionosphere, the
       accuracy deteriorates at low frequencies, below about 30 MHz.
  
    6  Before use, the value of ZOBS is expressed in the range +/- pi.
       If this ranged ZOBS is -ve, the result REF is computed from its
       absolute value before being made -ve to match.  In addition, if
       it has an absolute value greater than 93 deg, a fixed REF value
       equal to the result for ZOBS = 93 deg is returned, appropriately
       signed.
  
    7  As in the original Hohenkerk and Sinclair algorithm, fixed values
       of the water vapour polytrope exponent, the height of the
       tropopause, and the height at which refraction is negligible are
       used.
  
    8  The radio refraction has been tested against work done by
       Iain Coulson, JACH, (private communication 1995) for the
       James Clerk Maxwell Telescope, Mauna Kea.  For typical conditions,
       agreement at the 0.1 arcsec level is achieved for moderate ZD,
       worsening to perhaps 0.5-1.0 arcsec at ZD 80 deg.  At hot and
       humid sea-level sites the accuracy will not be as good.
  
    9  It should be noted that the relative humidity RH is formally
       defined in terms of "mixing ratio" rather than pressures or
       densities as is often stated.  It is the mass of water per unit
       mass of dry air divided by that for saturated air at the same
       temperature and pressure (see Gill 1982).
  
    Called:  slDA1P, slATMT, slATMS
  
    P.T.Wallace   Starlink   3 June 1997
  
    Copyright (C) 1997 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    # Fixed parameters
    
    # 93 degrees in radians
    D93 = 1.623156204
    # Universal gas constant
    GCR = 8314.32
    # Molecular weight of dry air
    DMD = 28.9644
    # Molecular weight of water vapour
    DMW = 18.0152
    # Mean Earth radius (metre)
    S = 6378120.0
    # Exponent of temperature dependence of water vapour pressure
    DELTA = 18.36
    # Height of tropopause (metre)
    HT = 11000.0
    # Upper limit for refractive effects (metre)
    HS = 80000.0
    
    # The refraction integrand
    def refi(r, dn, rdndr):
        return rdndr / (dn + rdndr)

    # Transform ZOBS into the normal range.
    ZOBS1 = da1p(ZOBS)
    ZOBS2 = min(abs(ZOBS1), 1.0e93)

    # Keep other arguments within safe bounds.
    HMOK = min(max(HM, -1.0e3), 10.0e3)
    TDKOK = min(max(TDK, 100.0), 500.0)
    PMBOK = min(max(PMB, 0.0), 10000.0)
    RHOK = min(max(RH, 0.0), 1.0)
    WLOK = max(WL, 0.1)
    ALPHA = min(max(abs(TLR), 0.001), 0.01)

    # Tolerance for iteration.
    TOL = min(max(abs(EPS), 1.0e-12), 0.1) / 2.0

    # Decide whether optical/IR or radio case - switch at 100 microns.
    OPTIC = WLOK <= 100.0

    # Set up model atmosphere parameters defined at the observer.
    WLSQ = WLOK * WLOK
    GB = 9.784 * (1.0 - 0.0026 * math.cos(PHI + PHI) - 0.00000028 * HMOK)
    if OPTIC:
        A = (287.604 + (1.6288 + 0.0136 / WLSQ) / WLSQ) * 273.15e-6 / 1013.25
    else:
        A = 77.624e-6
    
    GAMAL = (GB * DMD) / GCR
    GAMMA = GAMAL / ALPHA
    GAMM2 = GAMMA - 2.0
    DELM2 = DELTA - 2.0
    TDC = TDKOK - 273.15
    PSAT = 10.0**((0.7859 + 0.03477 * TDC) / (1.0 + 0.00412 * TDC))
    PSAT *= (1.0 + PMBOK * (4.5e-6 + 6e-10 * TDC * TDC))
    if (PMBOK > 0.0):
        PWO = RHOK * PSAT / (1.0 - (1.0 - RHOK) * PSAT / PMBOK)
    else:
        PWO = 0.0

    W = PWO * (1.0 - DMW / DMD) * GAMMA / (DELTA - GAMMA)
    C1 = A * (PMBOK + W) / TDKOK
    if OPTIC:
        C2 = (A * W + 11.2684e-6 * PWO) / TDKOK
    else:
        C2 = (A * W + 12.92e-6 * PWO) / TDKOK

    C3 = (GAMMA - 1.0) * ALPHA * C1 / TDKOK
    C4 = (DELTA - 1.0) * ALPHA * C2 / TDKOK
    if OPTIC:
        C5 = 0.0
        C6 = 0.0
    else:
        C5 = 371897e-6 * PWO / TDKOK
        C6 = C5 * DELM2 * ALPHA / (TDKOK * TDKOK)

    # Conditions at the observer.
    R0 = S + HMOK
    TEMPO, DN0, RDNDR0 = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                              C1, C2, C3, C4, C5, C6, R0)
                              
    SK0 = DN0 * R0 * math.sin(ZOBS2)
    F0 = refi(R0, DN0, RDNDR0)

    # Conditions in the troposphere at the tropopause.
    RT = S + HT
    TT, DNT, RDNDRT = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                           C1, C2, C3, C4, C5, C6, RT)
    SINE = SK0 / (RT * DNT)
    ZT = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE, 0.0)))
    FT = refi(RT, DNT, RDNDRT)

    # Conditions in the stratosphere at the tropopause.
    DNTS, RDNDRP = atms(RT, TT, DNT, GAMAL, RT)
    SINE = SK0 / (RT * DNTS)
    ZTS = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE,0.0)))
    FTS = refi(RT, DNTS, RDNDRP)

    # Conditions at the stratosphere limit.
    RS = S + HS
    DNS, RDNDRS = atms(RT, TT, DNT, GAMAL, RS)
    SINE = SK0 / (RS * DNS)
    ZS = math.atan2(SINE, math.sqrt(max(1.0 - SINE * SINE, 0.0)))
    FS = refi(RS, DNS, RDNDRS)

    # Integrate the refraction integral in two parts;  first in the
    # troposphere (K=1), then in the stratosphere (K=2).

    # Initialize previous refraction to ensure at least two iterations.
    REFOLD = 1.0e6

    # Start off with 8 strips for the troposphere integration, and then
    # use the final troposphere value for the stratosphere integration,
    # which tends to need more strips.
    IS = 8

    # Troposphere then stratosphere.
    for K in [1,2]:
        # Start Z, Z range, and start and end values.
        if K == 1:
            Z0 = ZOBS2
            ZRANGE = ZT - Z0
            FB = F0
            FF = FT
        else:
            Z0 = ZTS
            ZRANGE = ZS - Z0
            FB = FTS
            FF = FS

        # Sums of odd and even values.
        FO = 0.0
        FE = 0.0

        # First time through the loop we have to do every point.
        N = 1
       
        # Start of iteration loop (terminates at specified precision).
        LOOP = True
        while LOOP:
            # Strip width.
            H = ZRANGE / float(IS)

            # Initialize distance from Earth centre for quadrature pass.
            if K == 1:
                R = R0
            else:
                R = RT

            # One pass (no need to compute evens after first time).
            for I in range(1, IS, N):
                # Sine of observed zenith distance.
                SZ = math.sin(Z0 + H * float(I))

                # Find R (to the nearest metre, maximum four iterations).
                if SZ > 1e-20:
                    W = SK0 / SZ
                    RG = R
                    DR = 1e6
                    J = 0
                    
                    while (abs(DR) > 1.0) and (J < 4):
                        J = J + 1
                        if K == 1:
                            TG, DN, RDNDR = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                                                 C1, C2, C3, C4, C5, C6, RG)
                        else:
                            DN, RDNDR = atms(RT, TT, DNT, GAMAL, RG)

                        DR = (RG * DN - W) / (DN + RDNDR)
                        RG = RG - DR

                    R = RG

                # Find the refractive index and integrand at R.
                if K == 1:
                    T, DN, RDNDR = atmt(R0, TDKOK, ALPHA, GAMM2, DELM2,
                                        C1, C2, C3, C4, C5, C6, R)
                else:
                    DN,RDNDR = atms(RT, TT, DNT, GAMAL, R)

                F = refi(R, DN, RDNDR)

                # Accumulate odd and (first time only) even values.
                if (N == 1) and ((I % 2) == 0):
                    FE = FE + F
                else:
                    FO = FO + F

            # Evaluate the integrand using Simpson's Rule.
            REFP = H * (FB + 4.0 * FO + 2.0 * FE + FF) / 3.0

            # Has the required precision been achieved?
            if (abs(REFP - REFOLD) > TOL):

                # No: prepare for next iteration.

                # Save current value for convergence test.
                REFOLD = REFP

                # Double the number of strips.
                IS = IS + IS

                # Sum of all current values = sum of next pass's even values.
                FE = FE + FO

                # Prepare for new odd values.
                FO = 0.0

                # Skip even values next time.
                N = 2
                
            else:
                # Yes: save troposphere component and terminate the loop.
                if (K == 1):
                    REFT = REFP

                LOOP = False
            # END IF
        # END FOR
    # END WHILE

    # Result.
    REF = REFT + REFP
    if (ZOBS1 < 0.0):
        REF = -REF

    return REF


def atmt(R0, T0, ALPHA, GAMM2, DELM2, C1, C2, C3, C4, C5, C6, R):
    """
    Internal routine used by REFRO
  
    Refractive index and derivative with respect to height for the
    troposphere.
  
    Given:
      R0      d    height of observer from centre of the Earth (metre)
      T0      d    temperature at the observer (deg K)
      ALPHA   d    alpha          )
      GAMM2   d    gamma minus 2  ) see HMNAO paper
      DELM2   d    delta minus 2  )
      C1      d    useful term  )
      C2      d    useful term  )
      C3      d    useful term  ) see source
      C4      d    useful term  ) of slRFRO
      C5      d    useful term  )
      C6      d    useful term  )
      R       d    current distance from the centre of the Earth (metre)
  
    Returned:
      T       d    temperature at R (deg K)
      DN      d    refractive index at R
      RDNDR   d    R    rate the refractive index is changing at R
  
    Note that in the optical case C5 and C6 are zero.
  
    P.T.Wallace   Starlink   30 May 1997
  
    Copyright (C) 1997 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    T = max(min(T0 - ALPHA * (R - R0), 320.0), 100.0)
    TT0 = T / T0
    TT0GM2 = TT0**GAMM2
    TT0DM2 = TT0**DELM2
    DN = 1.0 + (C1 * TT0GM2 - (C2 - C5 / T) * TT0DM2) * TT0
    RDNDR = R * (-C3 * TT0GM2 + (C4 - C6 / TT0) * TT0DM2)

    return T, DN, RDNDR


def atms(RT, TT, DNT, GAMAL, R):
    """
  
    Internal routine used by REFRO
  
    Refractive index and derivative with respect to height for the
    stratosphere.
  
    Given:
      RT      d    height of tropopause from centre of the Earth (metre)
      TT      d    temperature at the tropopause (deg K)
      DNT     d    refractive index at the tropopause
      GAMAL   d    constant of the atmospheric model = G  MD/R
      R       d    current distance from the centre of the Earth (metre)
  
    Returned:
      DN      d    refractive index at R
      RDNDR   d    R    rate the refractive index is changing at R
  
    P.T.Wallace   Starlink   14 July 1995
  
    Copyright (C) 1995 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    B = GAMAL / TT
    W = (DNT - 1.0) * math.exp(-B * (R - RT))
    DN = 1.0 + W
    RDNDR = -R * B * W

    return DN, RDNDR
        

def da1p(ANGLE):
    """
    Normalize angle into range +/- pi  (double precision)
  
    Given:
       ANGLE     dp      the angle in radians
  
    The result (double precision) is ANGLE expressed in the range +/- pi.
  
    P.T.Wallace   Starlink   23 November 1995
  
    Copyright (C) 1995 Rutherford Appleton Laboratory
    Copyright (C) 1995 Association of Universities for Research in Astronomy Inc.
    """
    DPI = 3.141592653589793238462643
    D2PI = 6.283185307179586476925287

    slDA1P = ANGLE % D2PI
    if (abs(slDA1P) >= DPI):
        slDA1P = slDA1P - math.copysign(D2PI, ANGLE)

    return slDA1P
    
