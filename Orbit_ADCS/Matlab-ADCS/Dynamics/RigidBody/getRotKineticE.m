function var = getRotKineticE(Io, om)
    om = om(:);
    var = 0.5*om'*Io*om;
end