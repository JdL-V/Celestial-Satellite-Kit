function astroConstants(in::Union{Int64,Vector{Int64}})

# astroConstants.jl - Returns astrodynamic-related physical constants.
#
# PROTOTYPE:
#   out = astro_constants(in)
#
# DESCRIPTION:
#   Returns a row vector of constants, in which there is the corresponding
#   constant for each element of the input vector.
#
#   List of identifiers:
#       Generic astronomical constants:
#           1   Universal gravity constant (G) (from DITAN and Horizon) [km^3/(kg*s^2)]
#           2   Astronomical Unit (AU) (from DE405) [km]
#               Note:  The value for 1 au is from the IAU 2012 Resolution B1.
#       Sun related:
#           3   Sun mean radius (from DITAN) [km]
#           4   Sun planetary constant (mu = mass * G) (from DE405)
#               [km^3/s^2]
#           31  Energy flux density of the Sun (from Wertz,SMAD)
#               [W/m2 at 1 AU]
#       Other:
#           5   Speed of light in the vacuum (definition in the SI and Horizon) [km/s]
#           6   Standard free fall (the acceleration due to gravity on the
#               Earth's surface at sea level) (from Wertz,SMAD) [m/s^2]
#           7   Mean distance Earth-Moon (from Wertz,SMAD) [km]
#           8   Obliquity (angle) of the ecliptic at Epoch 2000 (from
#               Horizon) [rad]
#           9   Gravitatonal field constant of the Earth (from Wertz,SMAD,
#               taken from JGM-2). This should be used in conjunction to
#               Earth radius = 6378.1363 km
#           10  Gravitatonal field constant of the Earth (taken from JGM-3).
#           32  Days in a Julian year y = 365.25 d  (from Horizon)
#       Planetary constants of the planets (mu = mass * G) [km^3/s^2]:
#           11  Me      (from DE405)
#           12  V       (from DE405)
#           13  E       (from DE405)
#           14  Ma      (from DE405)
#           15  J       (from DE405)
#           16  S       (from DE405)
#           17  U       (from DE405)
#           18  N       (from DE405)
#           19  P       (from DE405)
#           20  Moon    (from DE405)
#       Mean radius of the planets [km]:
#           21  Me      (from Horizon)
#           22  V       (from Horizon)
#           23  E       (from Horizon)
#           24  Ma      (from Horizon)
#           25  J       (from Horizon)
#           26  S       (from Horizon)
#           27  U       (from Horizon)
#           28  N       (from Horizon)
#           29  P       (from Horizon)
#           30  Moon    (from Horizon)
#
#   Notes for upgrading this function:
#       It is possible to add new constants.
#       - DO NOT change the structure of the function, as well as its
#           prototype.
#       - DO NOT change the identifiers of the constants that have already
#           been defined in this function. If you want to add a new
#           constant, use an unused identifier.
#       - DO NOT add constants that can be easily computed starting form
#           other ones (avoid redundancy).
#       Contact the author for modifications.
#
# INPUT:
#   in      Vector of identifiers of required constants.
#
# OUTPUT:
#   out     Vector of constants.
#
# EXAMPLE:
#   astroConstants([2, 4, 26])
#      Returns a row vector in which there is the value of the AU, the Sun
#      planetary constant and the mean radius of Saturn.
#
#   astroConstants(10 + [1:9])
#      Returns a row vector with the planetary constant of each planet.
#
# REFERENCES:
#   - DITAN (Direct Interplanetary Trajectory Analysis), Massimiliano
#       Vasile, 2006.
#	- Wertz J. R., Larson W. J., "Space Mission Analysis and Design", Third
#       Edition, Space Technology Library 2003.
#   [A]   DE405 - http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf
#   [B]   Explanatory Supplement to the Astronomical Almanac. 1992. K. P.
#         Seidelmann, Ed., p.706 (Table 15.8) and p.316 (Table 5.8.1),
#         University Science Books, Mill Valley, California. 
#   [C]   Tholen, D.J. and Buie, M.W. 1990. "Further Analysis of
#         Pluto-Charon Mutual Event Observations" BAAS 22(3):1129.
#   [D]   Seidelmann, P.K. et al. 2007. "Report of the IAU/IAG Working
#         Group on cartographic coordinates and rotational elements: 2006"
#         Celestial Mech. Dyn. Astr. 98:155-180. 
#   [F]   Anderson, J.D., et al. 1987. "The mass, gravity field, and
#         ephemeris of Mercury" Icarus 71:337-349.
#   [G]   Konopliv, A.S., et al. 1999. "Venus gravity: 180th degree and
#         order model" Icarus 139:3-18.
#   [H]   Folkner, W.M. and Williams, J.G. 2008. "Mass parameters and
#         uncertainties in planetary ephemeris DE421." Interoffice Memo.
#         343R-08-004 (internal document), Jet Propulsion Laboratory,
#         Pasadena, CA. 
#   [I]   Jacobson, R.A. 2008. "Ephemerides of the Martian Satellites -
#         MAR080" Interoffice Memo. 343R-08-006 (internal document),
#         Jet Propulsion Laboratory, Pasadena, CA. 
#   [J]   Jacobson, R.A. 2005. "Jovian Satellite ephemeris - JUP230"
#         private communication. 
#   [K]   Jacobson, R.A., et al. 2006. "The gravity field of the Saturnian
#         system from satellite observations and spacecraft tracking data"
#         AJ 132(6):2520-2526. 
#   [L]   Jacobson, R.A. 2007. "The gravity field of the Uranian system and
#         the orbits of the Uranian satellites and rings" BAAS 39(3):453. 
#   [M]   Jacobson, R.A. 2008. "The orbits of the Neptunian satellites and
#         the orientation of the pole of Neptune" BAAS 40(2):296. 
#   [N]   Jacobson, R.A. 2007. "The orbits of the satellites of Pluto -
#         Ephemeris PLU017" private communication.
#   [W1]  http://ssd.jpl.nasa.gov/?planet_phys_par Last retrieved
#         20/03/2013
#   [W2]  http://ssd.jpl.nasa.gov/?sat_phys_par Last retrieved
#         20/03/2013
#   [W3]  http://ssd.jpl.nasa.gov/horizons.cgi Last retrieved
#         20/03/2013
#   [M1]  Bills, B.G. and Ferrari, A.J. 1977. ``A Harmonic Analysis of
#         Lunar Topography'', Icarus 31, 244-259.
#   [M2]  Standish, E. M. 1998. JPL Planetary and Lunar Ephemerides,
#         DE405/LE405.
#   [M3]  Lunar Constants and Models Document, Ralph B. Roncoli, 23 Sept 2005,
#         JPL Technical Document D-32296 
#
#
# CALLED FUNCTIONS:
#   (none)
#
# AUTHOR:
#   Matteo Ceriotti, 2006, MATLAB, astroConstants.jl
#
# PREVIOUS VERSION:
#   Matteo Ceriotti, 2006, MATLAB, astro_constants.jl, Ver. 1.2
#       - Header and function name in accordance with guidlines.
#
# CHANGELOG:
#   26/10/2006, Camilla Colombo: Updated.
#   22/10/2007, Camilla Colombo: astroConstants(8) added (Obliquity (angle)
#       of the ecliptic at Epoch 2000).
#   02/10/2009, Camilla Colombo: Header and function name in accordance
#       with guidlines.
#   12/11/2010, Camilla Colombo: astroConstants(9) added (J2) Note: the
#       present value of J2 is not consistent with the value of the Earth
#       radius. This value of J2 should be used in conjunction to Earth
#       radius = 6378.1363 km
#   19/03/2013, Camilla Colombo: constants updated to NASA JPL website.
#       References added.
#   20/03/2013, REVISION, Francesca Letizia.
#   22/03/2013, Francesca Letizia: all GM from DE405.
#
# -------------------------------------------------------------------------

# 9: J2
# 32: 365.25
    out = zeros(length(in));
    for i = 1:length(in)
        if in[i] == 1
            out[i]=6.67259e-20; # From DITAN and Horizon
        elseif in[i] == 2
            out[i]=149597870.691; # From DE405
        elseif in[i] == 3
            # out[i]=700000; # From DITAN
            out[i]=6.955*10^5; # From Horizon [W3]
        elseif in[i] == 4
            # out[i]=0.19891000000000E+31*6.67259e-20; # From DITAN
            out[i]=1.32712440017987E+11; # From DE405 [A]
        elseif in[i] == 5
            out[i]=299792.458; # Definition in the SI, Horizon, DE405
        elseif in[i] == 6
            out[i]=9.80665; # Definition in Wertz, SMAD
        elseif in[i] == 7
            # out[i]=384401; # Definition in Wertz, SMAD
            out[i]=384400; # From Horizon [W3]
        elseif in[i] == 8
            # out[i]=23.43928111*pi/180; # Definition in Wertz, SMAD
            out[i]=84381.412/3600*pi/180; # Definition in Horizon
            # obliquity of ecliptic (J2000)    epsilon = 84381.412 (� 0.005) arcsec 
        elseif in[i] == 9
            out[i]=0.1082626925638815e-2; # Definition in Wertz, SMAD
        elseif in[i] == 10
            out[i]=-0.25324353457544e-05; # Definition taken from JGM-3
        elseif in[i] == 11
            # out[i]=0.33020000000000E+24*6.67259e-20; # From DITAN
            #out[i]=0.330104E+24*6.67259e-20;    # From Horizon [F]
            out[i]=2.203208E+4;    # From DE405
        elseif in[i] == 12
            # out[i]=0.48685000000000E+25*6.67259e-20; # From DITAN
            #out[i]=4.86732E+24*6.67259e-20;     # From Horizon [G]
            out[i]=3.24858599E+5; # From DE405
        elseif in[i] == 13
            # out[i]=0.59736990612667E+25*6.67259e-20; # From DITAN
            # out[i]=5.97219E+24*6.67259e-20;     # From Horizon [H]
            out[i] = 3.98600433e+5; # From DE405
        elseif in[i] ==  14
            # out[i]=0.64184999247389E+24*6.67259e-20; # From DITAN
            #out[i]=0.641693E+24*6.67259e-20; 	# From Horizon [I]
            out[i] = 4.2828314E+4; #Frome DE405
        elseif in[i] ==  15
            # out[i]=0.18986000000000E+28*6.67259e-20; # From DITAN
            #out[i]=1898.13E+24*6.67259e-20; 	# From Horizon [J]
            out[i] = 1.26712767863E+08; # From DE405
        elseif in[i] ==  16
            # out[i]=0.56846000000000E+27*6.67259e-20; # From DITAN
            # out[i]=568.319E+24*6.67259e-20;     # From Horizon [k]
            out[i] = 3.79406260630E+07; # From DE405
        elseif in[i] ==  17
            # out[i]=0.86832000000000E+26*6.67259e-20; # From DITAN
            # out[i]=86.8103E+24*6.67259e-20;     # From Horizon [L]
            out[i]= 5.79454900700E+06; # From DE405
        elseif in[i] ==  18
            # out[i]=0.10243000000000E+27*6.67259e-20; # From DITAN
            # out[i]=102.410E+24*6.67259e-20;     # From Horizon [M]
            out[i] = 6.83653406400E+06; # From DE405
        elseif in[i] ==  19
            # out[i]=0.14120000000000E+23*6.67259e-20; # From DITAN
            #out[i]=.01309E+24*6.67259e-20;     # From Horizon [N]
            out[i] = 9.81601000000E+02; #From DE405
        elseif in[i] ==  20
            # out[i]=0.73476418263373E+23*6.67259e-20; # From DITAN
             out[i]=4902.801;                 # From Horizon  [M2]
            #out[i]=4902.801076;                # From Horizon  [M3]
        elseif in[i] ==  21
            # out[i]=0.24400000000000E+04; # From DITAN
            out[i]=2439.7; # From Horizon [D]
        elseif in[i] ==  22
            # out[i]=0.60518000000000E+04; # From DITAN
            out[i]=6051.8; # From Horizon [D]
        elseif in[i] ==  23
            # out[i]=0.63781600000000E+04; # From DITAN
            # out[i]=6371.00; # From Horizon [B]
            out[i]=6371.01; # From Horizon [W3]
        elseif in[i] ==  24
            # out[i]=0.33899200000000E+04; # From DITAN
            # out[i]=3389.50; # From Horizon [D]
            out[i]=3389.9; # From Horizon [W3]            
        elseif in[i] ==  25
            # out[i]=0.69911000000000E+05; # From DITAN
            out[i]=69911;   # From Horizon [D]
        elseif in[i] ==  26
            # out[i]=0.58232000000000E+05; # From DITAN
            out[i]=58232;   # From Horizon [D]
        elseif in[i] ==  27
            # out[i]=0.25362000000000E+05; # From DITAN
            out[i]=25362;   # From Horizon [D]
        elseif in[i] ==  28
            # out[i]=0.24624000000000E+05; # From DITAN
            # out[i]=24622;   # From Horizon [D]
            out[i]= 24624; # From Horizon [W3]            
        elseif in[i] ==  29
            # out[i]=0.11510000000000E+04; # From DITAN
            out[i]=1151; 	# From Horizon [C]
        elseif in[i] ==  30
            # out[i]=0.17380000000000E+04; # From DITAN
            # out[i]=1737.5;  # From Horizon [M1]
            out[i]=1738.0;    # From Horizon  [M3]
        elseif in[i] ==  31
            out[i]=1367; # From Wertz, SMAD
            # out[i]=1367.6;  # From Horizon  [W3]
        elseif in[i] ==  32
            out[i]=365.25; # From Horizon
        # Add an identifier and constant here. Prototype:
        # elseif in[i] ==  $identifier$
        #     out[i]=$constant_value$;
        else
            println("Constant identifier ", in[i], " is not defined!");
            out[i]=0;
        end
    end
    return out
end


# http://www.csr.utexas.edu/publications/statod/TabD.3.new.txt

# l	m	            C                                S

# 2	0	  -0.10826360229840D-02	           0.0
# 3	0	   0.25324353457544D-05	           0.0
# 4	0	   0.16193312050719D-05	           0.0
# 5	0	   0.22771610163688D-06	           0.0
# 6	0	  -0.53964849049834D-06	           0.0
# 7	0	   0.35136844210318D-06	           0.0
# 8	0	   0.20251871520885D-06 	       0.0
# 2	1	  -0.24140000522221D-09	           0.15430999737844D-08
# 3	1	   0.21927988018965D-05	           0.26801189379726D-06
# 4	1	  -0.50872530365024D-06	          -0.44945993508117D-06
# 5	1	  -0.53716510187662D-07	          -0.80663463828530D-07
# 6	1	  -0.59877976856303D-07	           0.21164664354382D-07
# 7	1	   0.20514872797672D-06	   	       0.69369893525908D-07
# 8	1	   0.16034587141379D-07	   	       0.40199781599510D-07 
# 2	2	   0.15745360427672D-05	  	      -0.90386807301869D-06
# 3	2	   0.30901604455583D-06	 	      -0.21140239785975D-06
# 4	2	   0.78412230752366D-07	    	   0.14815545694714D-06
# 5	2	   0.10559053538674D-06	   	      -0.52326723987632D-07
# 6	2	   0.60120988437373D-08		      -0.46503948132217D-07
# 7	2	   0.32844904836492D-07	    	   0.92823143885084D-08
# 8	2	   0.65765423316743D-08		       0.53813164055056D-08
# 3	3 	   0.10055885741455D-06	    	   0.19720132389889D-06
# 4	3	   0.59215743214072D-07	  	      -0.12011291831397D-07
# 5	3	  -0.14926153867389D-07	  	      -0.71008771406986D-08
# 6	3	   0.11822664115915D-08	      	   0.18431336880625D-09
# 7	3	   0.35285405191512D-08	  	      -0.30611502382788D-08
# 8	3	  -0.19463581555399D-09	  	      -0.87235195047605D-09
# 4	4	  -0.39823957404129D-08	   	       0.65256058113396D-08
# 5	4	  -0.22979123502681D-08	   	       0.38730050770804D-09
# 6	4	  -0.32641389117891D-09	   	      -0.17844913348882D-08
# 7	4	  -0.58511949148624D-09	   	      -0.26361822157867D-09
# 8	4	  -0.31893580211856D-09	    	   0.91177355887255D-10
# 5	5	   0.43047675045029D-09	    	  -0.16482039468636D-08
# 6	5	  -0.21557711513900D-09	    	  -0.43291816989540D-09
# 7	5	   0.58184856030873D-12	      	   0.63972526639235D-11
# 8	5	  -0.46151734306628D-11	     	   0.16125208346784D-10
# 6	6	   0.22136925556741D-11	     	  -0.55277122205966D-10
# 7	6	  -0.24907176820596D-10	       	   0.10534878629266D-10
# 8	6	  -0.18393642697634D-11	       	   0.86277431674150D-11
# 7	7	   0.25590780149873D-13 	       0.44759834144751D-12 
# 8 7	   0.34297618184624D-12 	       0.38147656686685D-12 
# 8	8	  -0.15803322891725D-12	           0.15353381397148D-12