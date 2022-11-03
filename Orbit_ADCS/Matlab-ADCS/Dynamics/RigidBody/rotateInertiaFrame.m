function var = rotateInertiaFrame(Io, dcm)
    var = dcm*Io*dcm';
end