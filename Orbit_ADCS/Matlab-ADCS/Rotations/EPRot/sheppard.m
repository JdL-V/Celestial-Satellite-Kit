function var = sheppard(dcm)
    Trace = dcm.Mat(1,1) + dcm.Mat(2,2) + dcm.Mat(3,3);
    q = quat(sqrt([1 + Trace
                            1 - Trace + 2*dcm.Mat(1,1)
                            1 - Trace + 2*dcm.Mat(2,2)
                            1 - Trace + 2*dcm.Mat(3,3)]./4), ...
                       dcm.outFrame, dcm.inFrame);

    index = [1 2 3 4];
    index = index(min(find(q.x == max(q.x))));
    vec = 1:3;
    if index == 1
        for j = vec
            i = vec(1:end ~= j);
            q.x(j+1)  = (dcm.Mat(i(1),i(2)) - dcm.Mat(i(2),i(1)))/(4*q.x(index))*(-1).^(j + 1);
        end
    else
        i = vec(1:end ~= (index-1));
        q.x(1) = (dcm.Mat(i(1),i(2)) - dcm.Mat(i(2),i(1)))/(4*q.x(index))*(-1).^index;
        for j = vec(1:end ~= index-1)
            q.x(j+1) = (dcm.Mat(index-1,j) + dcm.Mat(j,index-1))/(4*q.x(index));
        end
        q.x = (q.x).*sign(q.x(1));
    end
    var = q;
end