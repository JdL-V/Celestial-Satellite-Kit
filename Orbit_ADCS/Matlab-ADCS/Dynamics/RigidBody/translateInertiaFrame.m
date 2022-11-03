function var = translateInertiaFrame(Ig, M, Rp)
    var = Ig + sum(M)*skewsym(Rp)*skewsym(Rp)';
end