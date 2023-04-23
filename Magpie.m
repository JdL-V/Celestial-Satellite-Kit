Init
seed(-1)
warning off
optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);


%% Initilize variables
TimeVec  = 0 : 1 : 60;
SunVecsI = normalize(rand(3, length(TimeVec)), 'norm');
MagVecsI = normalize(rand(3, length(TimeVec)), 'norm');

SunVecsB = NaN (3, length(TimeVec));
MagVecsB = NaN (3, length(TimeVec));

LightFlag = [zeros(1,20) zeros(1,21) zeros(1,20)];

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
    SunVecsB(:,i) = normalize(ChangeFrame(SunVecsI(:,i), DCMBI) .* (1 + (rand(3,1) - 0.5)*1e-100), 'norm');
    MagVecsB(:,i) = normalize(ChangeFrame(MagVecsI(:,i), DCMBI) .* (1 + (rand(3,1) - 0.5)*1e-100), 'norm');
end
toc

%% Simulate attitude determination
tic
qut = simMagpie(TimeVec, SunVecsI, SunVecsB, MagVecsI, MagVecsB, LightFlag, 8);
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
    w   = @(t) deg2rad(6).*normalize([0.1 0.1 1], 'norm')';
    M   = om2EP(u);
    var = M * w(t);
end

function om = omDyn(t, u)
    global qutpp1 qutpp2 qutpp3 qutpp4

    q = @(t) [(ppval(t, qutpp1));
              (ppval(t, qutpp2));
              (ppval(t, qutpp3));
              (ppval(t, qutpp4))];

    qdot = @(t) [(ppval(t+1e-0, qutpp1) - ppval(t-1e-0, qutpp1));
                 (ppval(t+1e-0, qutpp2) - ppval(t-1e-0, qutpp2));
                 (ppval(t+1e-0, qutpp3) - ppval(t-1e-0, qutpp3));
                 (ppval(t+1e-0, qutpp4) - ppval(t-1e-0, qutpp4))]/2e-0 ;
                 
    M    = EP2om(q(t));
    om0  = M * qdot(t);
    om   = om0(2:4);
end

function var = simMagpie(TimeVec, SunVecsI, SunVecsB, MagVecsI, MagVecsB, LightFlag, tol)

    weight = [1, 1];
    var = NaN(4, length(TimeVec));

    for i = 1 : length(TimeVec)
        switch LightFlag(i)
        case 1
            vkb = [SunVecsB(:,i) MagVecsB(:,i)];
            vki = [SunVecsI(:,i) MagVecsI(:,i)];
    
            DCM = qMethod(vkb, vki, weight);
            var(:,i) = sheppard(DCM).x;
        case 0
            vkb = [[0 0 1]' MagVecsB(:,i)];
            vki = [[0 0 1]' MagVecsI(:,i)];

            DCM = qMethod(vkb, vki, weight);
            var(:,i) = sheppard(DCM).x;

            shift = 1 + abs(var(1,i));
%             err = norm(DCM.Mat * MagVecsI(:,i) - MagVecsB(:,i)) * (abs(var(1,i)));
%             err = acosd(dot(DCM.Mat * MagVecsI(:,i), MagVecsB(:,i)));
            err = acosd(dot([0 0 1], dcm2PRV(DCM).x)) * (abs(var(1,i)) < 1 - 1d-5);
            err = min(abs(err - 0), abs(err - 180)) * shift * 0.8;
            if err >= tol
                var(:,i) = zeros(4,1);
            end
        end
    end

end

function vecB = ChangeFrame(vecI, DCMBI)
    vecB = DCMBI.Mat * vecI;
end