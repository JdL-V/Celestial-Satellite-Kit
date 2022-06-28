function plot_rotation(i, j, k, npoints, tspan, yrot, plot_speed)

O = [0 0 0]';
figure
quiver3(O,O,O,i,j,k,'k')
xlabel('x')
ylabel('y')
zlabel('z')
axis equal
xlim([-1 1]);
ylim([-1 1]);
zlim([-1 1]);
[az, el] = view;
view([az-90 el])
text(i,j,k, {'i','j','k'})
hold on

filename    = 'POCKET_CUBE.STL';
[p,t,~] = import_stl_fast2(filename,1);
p = p/1e2; % scaling
p = [p(:,1) - max(p(:,1))/2,  p(:,2) - max(p(:,2))/2,  p(:,3) - max(p(:,3))/2]; % centering
p2 = p;

for kk = [1:plot_speed:length(tspan) npoints]

I = yrot(kk, 1:3)';
J = yrot(kk, 4:6)';
K = yrot(kk, 7:9)';

for kn = 1:size(p,1)
    p2(kn,:) = ([I J K]'*p(kn,:)')';
end

try
    delete(ppp)
    delete(hQuiver)
catch
end

ppp = patch('vertices',p2,'faces',t,'FaceColor','c','edgecolor','k');
hQuiver = quiver3(O,O,O,I,J,K,'b');
drawnow
title(['t = ' num2str(round(tspan(kk))) ' s'])
end
end

