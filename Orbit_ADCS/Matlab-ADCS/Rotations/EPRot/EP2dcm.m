function var = EP2dcm(q)
    dcm = DCM([q.x(1)^2 + q.x(2)^2 - q.x(3)^2 - q.x(4)^2   2*(q.x(2)*q.x(3) + q.x(1)*q.x(4))           2*(q.x(2)*q.x(4) - q.x(1)*q.x(3))
            2*(q.x(2)*q.x(3) - q.x(1)*q.x(4))           q.x(1)^2 - q.x(2)^2 + q.x(3)^2 - q.x(4)^2   2*(q.x(3)*q.x(4) + q.x(1)*q.x(2)) 
            2*(q.x(2)*q.x(4) + q.x(1)*q.x(3))           2*(q.x(3)*q.x(4) - q.x(1)*q.x(2))           q.x(1)^2 - q.x(2)^2 - q.x(3)^2 + q.x(4)^2],...
            q.outFrame, q.inFrame);
    var = dcm;
end