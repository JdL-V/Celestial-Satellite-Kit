function var = sheppard(dcm)
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
    
    Trace = Mat(1,1) + Mat(2,2) + Mat(3,3);
    q = quat(sqrt([1 + Trace
                            1 - Trace + 2*Mat(1,1)
                            1 - Trace + 2*Mat(2,2)
                            1 - Trace + 2*Mat(3,3)]./4), ...
                       outFrame, inFrame);

    index = [1 2 3 4];
    index = index(min(find(q.x == max(q.x))));
    vec = 1:3;
    if index == 1
        for j = vec
            i = vec(1:end ~= j);
            q.x(j+1)  = (Mat(i(1),i(2)) - Mat(i(2),i(1)))/(4*q.x(index))*(-1).^(j + 1);
        end
    else
        i = vec(1:end ~= (index-1));
        q.x(1) = (Mat(i(1),i(2)) - Mat(i(2),i(1)))/(4*q.x(index))*(-1).^index;
        for j = vec(1:end ~= index-1)
            q.x(j+1) = (Mat(index-1,j) + Mat(j,index-1))/(4*q.x(index));
        end
        q.x = (q.x).*sign(q.x(1));
    end
    var = q;
end