Init
warning off
optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);


%% Initilize variables
TimeVec  = 0 : 0.1 : 60;
SunVecsI = rand(3, length(TimeVec));
MagVecsI = rand(3, length(TimeVec));

SunVecsB = NaN (3, length(TimeVec));
MagVecsB = NaN (3, length(TimeVec));

%% Solve differential equation to get quaternion from simulated attitude
tic
[tout, uout] = EPdiff(@qutDyn, TimeVec, [1 0 0 0], optiondop);
figure;plot(TimeVec, uout)
toc

%% Rotate inertial sun and magnetic vectors to body
tic
for i = 1 : length(TimeVec)
    qutBI = uout(i,:);
    DCMBI = EP2dcm(qutBI);
    SunVecsB(:,i) = normalize(ChangeFrame(SunVecsI(:,i), DCMBI) .* (1 + rand(3,1)*1e-20), 'norm');
    MagVecsB(:,i) = normalize(ChangeFrame(MagVecsI(:,i), DCMBI) .* (1 + rand(3,1)*1e-20), 'norm');
end
toc

%% Simulate attitude determination
tic
qut = simMagpie(TimeVec, SunVecsI, SunVecsB, MagVecsI, MagVecsB);
qut = EPsmooth(qut);
toc
figure;plot(TimeVec, qut)

global qutpp1 qutpp2 qutpp3 qutpp4
qutpp1 = makima(TimeVec, qut(1,:));
qutpp2 = makima(TimeVec, qut(2,:));
qutpp3 = makima(TimeVec, qut(3,:));
qutpp4 = makima(TimeVec, qut(4,:));

tic
[tout2, uout2] = dop853(@omDyn, TimeVec, [0 0 0]', optiondop);
toc

figure;plot(TimeVec, uout2)

function var = qutDyn(t, u)
    w   = @(t) deg2rad(6).*[0 0 1]';
    M   = om2EP(u);
    var = M * w(t);
end

function om = omDyn(t, u)
    global qutpp1 qutpp2 qutpp3 qutpp4

    q = @(t) [(ppval(t, qutpp1));
              (ppval(t, qutpp2));
              (ppval(t, qutpp3));
              (ppval(t, qutpp4))];

    qdot = @(t) [(ppval(t+1e-8, qutpp1) - ppval(t, qutpp1));
                 (ppval(t+1e-8, qutpp2) - ppval(t, qutpp2));
                 (ppval(t+1e-8, qutpp3) - ppval(t, qutpp3));
                 (ppval(t+1e-8, qutpp4) - ppval(t, qutpp4))]/1e-8 ;
                 
    M    = EP2om(q(t));
    om0  = M * qdot(t);
    om   = om0(2:4);
end

function var = simMagpie(TimeVec, SunVecsI, SunVecsB, MagVecsI, MagVecsB)

    weight = [1, 1];
    qut = NaN(4, length(TimeVec));

    for i = 1 : length(TimeVec)
        vkb = [SunVecsB(:,i) MagVecsB(:,i)];
        vki = [SunVecsI(:,i) MagVecsI(:,i)];

        var(:,i) = sheppard(qMethod(vkb, vki, weight)).x;
    end

end

function vecB = ChangeFrame(vecI, DCMBI)
    vecB = DCMBI.Mat * vecI;
end