function var = dcm2EP(dcm, eta)
    if isstruct(dcm)
        Mat = dcm.Mat;
        outFrame = dcm.outFrame;
        inFrame = dcm.inFrame;
    else
        checkMat(dcm, [3,3]);
        Mat = dcm;
        outFrame = "unk" + num2str(randi(5000));
        inFrame = "unk" + num2str(randi(5000));
    end

    if nargin == 1
        eta = 0;
    end

    Trace = Mat(1,1) + Mat(2,2) + Mat(3,3);
    check = [Trace,
             2*Mat(1,1) - Trace,
             2*Mat(2,2) - Trace,
             2*Mat(3,3) - Trace];

    q = quat(zeros(4,1), outFrame, inFrame);
    for i = 1:4
        k0 = 1:3;
        k = k0(1:end ~= (i-1));
        sgn = sign((Mat(k(1),k(2)) - Mat(k(2),k(1)))*(-1).^i)*(i ~= 1) + (i == 1);

        if check(i) > eta
            q.x(i) = 0.5*sqrt(1 + check(i))*sgn;
        else
            q.x(i) = 0.5*sqrt(( (Mat(3,2) + Mat(2,3)*(-1)^(i == 2 || i == 1))^2 ... 
                              + (Mat(1,2) + Mat(2,1)*(-1)^(i == 4 || i == 1))^2 ... 
                              + (Mat(3,1) + Mat(1,3)*(-1)^(i == 3 || i == 1))^2) / (3 - check(i)))*sgn;
        end
    end
    var = q;
end