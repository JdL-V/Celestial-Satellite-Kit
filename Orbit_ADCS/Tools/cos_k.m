function kelly = cos_k(angle)
        b = 5.3;
        c = 31;
        d = 86;
        e = 1.000016658987724;

        x = wrapTo360(angle);
        kelly = (e*cosd(x + b*exp(-(d - x)/c)).*(x >= 0).*(x < 90)...
           - e*cosd(180 - x + b*exp(-(d - 180 + x)/c)).*(x >= 90).*(x < 180)...
           - e*cosd(x - 180 + b*exp(-(d - x + 180)/c)).*(x >= 180).*(x < 270)...
           + e*cosd(360 - x + b*exp(-(d - 360 + x)/c)).*(x >= 270).*(x <= 360)).*(abs(cosd(x)) > cosd(85));
    end