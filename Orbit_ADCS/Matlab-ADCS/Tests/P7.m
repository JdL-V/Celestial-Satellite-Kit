format compact

prvBN = PRV(rad2deg(45), [0.,0.,1.], "B", "N");
dcmBN = PRV2dcm(prvBN);
neN = [-1.,0.,0.]';
nsN = [0.,1.,0.]';
neB = dcmBN.Mat*neN;
nsB = dcmBN.Mat*nsN;
dcmBN2 = TRIAD(neB, nsB, neN, nsN);
dcmBN.Mat
dcmBN2.Mat
dcmBN3 = qMethod([neB nsB], [neN nsN], [1.,1.]);
dcmBN4 = QUEST([neB nsB], [neN nsN], [1.,1.]);
dcmBN5 = OLAE([neB nsB], [neN nsN], [1.,1.]);
dcmBN3.Mat
dcmBN4.Mat
dcmBN5.Mat

dcmBN6 = TRIAD([0.8273, 0.5541, -0.0920], [-0.8285, 0.5522, -0.0955], [-0.1517, -0.9669, 0.2050], [-0.8393, 0.4494, -0.3044]);
dcmBN6.Mat
dcmBN7 = OLAE([[0.8273, 0.5541, -0.0920]', [-0.8285, 0.5522, -0.0955]'], [[-0.1517, -0.9669, 0.2050]', [-0.8393, 0.4494, -0.3044]'],  [1.,1.]);
dcmBN7.Mat

u01 = deg2rad([80, 30, 40]');
optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);
warning off
tic
    [tout1, uout1] = dop853(@feuler, linspace(0,60,1e3), u01, optiondop);
    figure; plot(tout1,uout1)
toc

u0 = [1, 0, 0, 0]';
tic
    [tout, uout] = EPdiff(@fq2, linspace(0,60,1e3), u0, optiondop);
    figure; plot(tout,uout)
toc
warning on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function du = feuler(t,u)
    global ROT
    global TP
    % 321
    w = @(t) deg2rad(5).*[sin(0.1*t), 0.01, cos(0.1*t)]';
    M = om2euler(EulerAng([u(1), u(2), u(3)], [3 2 1], 'A', 'B'));
    du = M*w(t);
end

function du = fq2(t,u)
    global ROT
    w = @(t) deg2rad(50).*[sin(0.1*t), 0.01, cos(0.1*t)]';
    M = om2EP(u);
    du = M*w(t);
end