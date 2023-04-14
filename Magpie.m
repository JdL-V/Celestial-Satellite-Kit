clc
clear
close all
Init

TimeVec  = 0 : 1 : 60;
SunVecsI = rand(3, length(TimeVec));
MagVecsI = rand(3, length(TimeVec));

SunVecsB = NaN (3, length(TimeVec));
MagVecsB = NaN (3, length(TimeVec));

warning off
optiondop = rdpset('RelTol',1e-7,'AbsTol',1e-7,'Refine',10);
[tout, uout] = EPdiff(@qutDyn, TimeVec, [1 0 0 0], optiondop);

for i = 1 : length(TimeVec)
    DCMBI = EP2dcm(quat(uout(i,:), "B", "I"));
    SunVecsB(:,i) = ChangeFrame(SunVecsI(:,i), DCMBI);
    MagVecsB(:,i) = ChangeFrame(MagVecsI(:,i), DCMBI);
end

qut = simMagpie(TimeVec, SunVecsI, SunVecsB, MagVecsI, MagVecsB);
figure; plot(qut')
qut = EPsmooth(qut);
figure; plot(qut')
% for i = 1:4; figure; plot(qut(i,:)); figure; plot(uout(:,i)); end

% qutpp1 = polyfit(TimeVec, qut(1,:), 4)
% plot(0:0.01:25, polyval(qutpp1, 0:0.01:25))

global qutpp1 qutpp2 qutpp3 qutpp4
qutpp1 = spline(TimeVec, qut(1,:));
qutpp2 = spline(TimeVec, qut(2,:));
qutpp3 = spline(TimeVec, qut(3,:));
qutpp4 = spline(TimeVec, qut(4,:));

[tout2, uout2] = dop853(@omDyn, TimeVec, [0 0 0]', optiondop);

figure;plot(uout2)

function var = qutDyn(t, u)
    w   = @(t) deg2rad(6).*[0 0 1]';
    M   = om2EP(u);
    var = M * w(t);
end

function var = omDyn(t, u)
    global qutpp1 qutpp2 qutpp3 qutpp4
    q = @(t) [(ppval(t, qutpp1)); (ppval(t, qutpp2)); (ppval(t, qutpp3)); (ppval(t, qutpp4))];
    qdot = @(t) [(ppval(t+1e-8, qutpp1) - ppval(t, qutpp1)); (ppval(t+1e-8, qutpp2) - ppval(t, qutpp2)); (ppval(t+1e-8, qutpp3) - ppval(t, qutpp3)); (ppval(t+1e-8, qutpp4) - ppval(t, qutpp4))]/1e-8 ;
    M    = EP2om(q(t));
    var  = M * qdot(t);
    var = var(2:4);
end

function var = simMagpie(tTime, SunVecsI, SunVecsB, MagVecsI, MagVecsB)

    weight = [1, 1];
    qut = NaN(4, length(tTime));

    for i = 1 : length(tTime)
        vkb = [SunVecsB(:,i) MagVecsB(:,i)];
        vki = [SunVecsI(:,i) MagVecsI(:,i)];

        var(:,i) = sheppard(qMethod(vkb, vki, weight)).x;
    end

end

function var = ChangeFrame(vecI, DCMBI)
    var = DCMBI.Mat * vecI;
end