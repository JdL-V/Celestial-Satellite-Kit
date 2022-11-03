function var = PRV2om(prv)
    tau = skewsym(prv.x*prv.Angle);
    var = (eye(3,3) + 0.5.*tau + 1/prv.Angle^2*(1 - prv.Angle/2*cot(prv.Angle/2)).*tau^2);
end