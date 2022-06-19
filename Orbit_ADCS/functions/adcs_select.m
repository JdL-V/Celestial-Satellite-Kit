function pointing = adcs_select(case_adcs)
    if     case_adcs == 1   %% Magnetic orientation
        pointing = @ parallel2mag;

    elseif case_adcs == 2   %% Earth pointing case. 
        pointing =  @ point2earth;

    elseif case_adcs == 3   %% Rotation around the perpendicular axis to the orbital plane case.
        pointing = @ perp2orbit;

    elseif case_adcs == 4   %% Normal to sun
        pointing = @ norm2sun;
    end
end