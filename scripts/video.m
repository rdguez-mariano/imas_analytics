addpath('extrafuns');
addpath('plotball');


figure;
%set (gcf, 'Units', 'normalized', 'Position', [0,0,1,1]);
optionsdraw1.drawtitle = false;
optionsdraw1.tilt = 8;
optionsdraw1.draw_circle_eachtilt= false;
optionsdraw1.drawtitle=false;
optionsdraw1.disk_random_color=false;

r=log(sqrt(2));

loops = 300;

v = VideoWriter('newfile.avi');
open(v)
t = exp([linspace(log(1),log(sqrt(2)),loops/2) linspace(log(sqrt(2)),log(2*sqrt(2)),loops/2)])- 1;
t = unique(t);
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    x = t(j).*cos(t(j));
    y = t(j).*sin(t(j));
    
    clf
    subplot(1,2,1);
    draw_ball_pol((t(j)+1),{10*t(j)},r,optionsdraw1);
    subplot(1,2,2);
    draw_ball3d((t(j)+1),10*t(j),r);
    drawnow
    
    F(j) = getframe(gcf);
    
    writeVideo(v,F(j));
end

%movie(F,2)
close(v)