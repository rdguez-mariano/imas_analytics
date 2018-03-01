function [t_vect,psi_cell] = draw_ball_pol(tvec,psicell,r,opts)

global balls beta NUMPOINTSINTERACTIVE cont_simu

%colors
 green = [51/255 150/255 0 ];
 red = [204/255 0 51/255];
 DRAW_CENTER_MarkerSize = 20;
 DRAW_LineSize = 1.5;
 DRAW_DashedLineSize = 2;

%Default Parameters
DRAWOPT = false;
PLOTTILT = [4*sqrt(2) 8];
PLOTTILT_NAME = {'\Gamma' ' '};%{'\Gamma' '\Gamma^{\prime}'};
NUMPOINTS = 100000;
DRAWTITLE = true;
DRAW_CIRCLE_TILT = true;
INTERACTIVE = false;
NUMPOINTSINTERACTIVE = 10000;
DRAW_REGIONS = false;
DRAW_REGIONS_TILT = sqrt(2)*4;
DRAW_CENTER_DISKS_OPTS = '.';
DRAW_CENTER_DISKS_COLOR = green;
DISK_RANDOM_COLOR = false;




% Input Parameters
switch nargin
    case 4
        if (~isstruct(opts))
            error('Optional parameter is not of type struct!');
        end
        names = fieldnames(opts);
        for i=1:length(names)
            name = names{i};
            switch name
                case 'numpoints'
                    NUMPOINTS = opts.numpoints;
                case 'numpointsinteractive'
                    NUMPOINTSINTERACTIVE = opts.numpointsinteractive;                    
                case 'tilt'
                    DRAWOPT = true;
                    PLOTTILT = opts.tilt;
                case 'drawtitle'
                    DRAWTITLE = opts.drawtitle;
                case 'draw_circle_eachtilt'
                    DRAW_CIRCLE_TILT = opts.draw_circle_eachtilt;
                case 'interactive'
                    INTERACTIVE = opts.interactive;
                case 'draw_regions'
                    DRAW_REGIONS = opts.draw_regions;                    
                case 'draw_regions_tilt'
                    DRAW_REGIONS_TILT = opts.draw_regions_tilt;  
                case 'draw_center_disks_opts'
                    DRAW_CENTER_DISKS_OPTS = opts.draw_center_disks_opts;
                case 'disk_random_color'
                    DISK_RANDOM_COLOR = opts.disk_random_color;   
                case 'draw_center_markersize'
                    DRAW_CENTER_MarkerSize = opts.draw_center_markersize;
                case 'draw_linesize'
                    DRAW_LineSize = opts.draw_linesize;
                case 'draw_dashedlinesize'
                    DRAW_DashedLineSize = opts.draw_dashedlinesize;
            end
        end
    otherwise
        if (nargin~=3)
            error('Wrong number of input parameters!');
        end
end





gamma = exp(r);
beta = ( (gamma^2)+1 )/(2*gamma);
balls = {};
cont_simu = 0;
phivect = linspace(0,2*pi,NUMPOINTS);
hold on;
for i=1:length(tvec)
    t=tvec(i);
    psivec=psicell{i};
    cont_simu = cont_simu + length(psivec);
    for j=1:length(psivec)
        psi = psivec(j);


        [pointinf,pointsup] = edges_one_ball(phivect,t,psi,beta);

        
        if (~INTERACTIVE)            
            if (DISK_RANDOM_COLOR)
                colorvec = .9*rand(1,3);
            else
                colorvec =red;
            end
            ballind = length(balls)+1;
            tempball.ballcenter = [t psi];
            plot(log(t).*cos(psi),log(t).*sin(psi),DRAW_CENTER_DISKS_OPTS,'MarkerSize',DRAW_CENTER_MarkerSize,'Color',DRAW_CENTER_DISKS_COLOR);
            tempball.curvsup = plot(log(pointsup).*cos(phivect),log(pointsup).*sin(phivect),'-','LineWidth',DRAW_LineSize,'Color',colorvec);
            tempball.curvinf = plot(log(pointinf).*cos(phivect),log(pointinf).*sin(phivect),'-','LineWidth',DRAW_LineSize,'Color',colorvec);
            tempball.colorvec = colorvec;
            balls{ballind} = tempball;
        else
            if (DISK_RANDOM_COLOR)
                colorvec = .9*rand(1,3);
            else
                colorvec = red;
            end
            h = impoint(gca,log(t).*cos(psi),log(t).*sin(psi));
            setColor(h,'r');
            ballind = length(balls)+1;
            tempball.ballcenter = [t psi];
            
            addNewPositionCallback(h,@(h) MyPositionCallback(h,ballind));
            
            tempball.curvsup = plot(log(pointsup).*cos(phivect),log(pointsup).*sin(phivect),'-','Color',colorvec);
            tempball.curvinf = plot(log(pointinf).*cos(phivect),log(pointinf).*sin(phivect),'-','Color',colorvec);
            tempball.colorvec = colorvec;
            balls{ballind} = tempball;
        end
        
    end
    
    if (DRAW_CIRCLE_TILT)
        plot(log(t)*cos(phivect),log(t)*sin(phivect),':','Color',.3*rand(1,3));
    end
    
end
if (DRAWOPT)
    for i=1:length(PLOTTILT)
        tilt2plot = PLOTTILT(i);
        plot(log(tilt2plot)*cos(phivect),log(tilt2plot)*sin(phivect),'--k','LineWidth',DRAW_DashedLineSize);
        %text((log(tilt2plot)+0.01)*cos(pi/4),(log(tilt2plot)+0.01)*sin(pi/4),['\leftarrow Tilt = ' num2str(tilt2plot) ]);
        textpointer = text((log(tilt2plot)+0.05)*cos(pi/4),(log(tilt2plot)+0.05)*sin(pi/4),[PLOTTILT_NAME{i}]);
        textpointer.FontWeight='bold';
        textpointer.FontSize = DRAW_CENTER_MarkerSize;
    end
end


if (DRAWTITLE)
if ((length(tvec)==1) && (length(psicell)==1))
    title(['d( T_{' num2str(t) '}R_{' num2str(psi) '} , T_sR_{\phi-' num2str(psi) '} ) \leq ' num2str(r)]);
else
    textpointer = title([num2str(cont_simu) ' affine simulations']);
    %textpointer.FontSize = 14;
end
end
    
xlabel('log(\tau) cos(\phi)','FontSize',DRAW_CENTER_MarkerSize);
ylabel('log(\tau) sin(\phi)','FontSize',DRAW_CENTER_MarkerSize);
axis equal




if (INTERACTIVE)
    set(gcf,'WindowbuttonDownFcn',@doubleclickcallback)
    waitfor(gcf);
else
    if (DRAW_REGIONS)
        discrete_rho = 500;
        discrete_phi = 1000;
        [s,phi] = meshgrid(linspace(1,DRAW_REGIONS_TILT,discrete_rho),linspace(0,pi,discrete_phi));
        M = zeros(size(s));
        for i=1:length(balls)
            t = balls{i}.ballcenter(1);
            psi = balls{i}.ballcenter(2);
            G = ( (t ./s + s ./t).*(cos(phi-psi).^2) + ((1 ./(s.*t)) +s.*t).*(sin(phi-psi).^2) )/2;
            M((G<beta)) = M((G<beta)) + 1;
        end
        Z = NaN(size(s));
        Z(M==0)=-.1;
        h = surf(log(s).*cos(phi),log(s).*sin(phi),Z,'EdgeColor','none','LineStyle','none');
        set(h,'FaceColor',[.5 .5 .5]);
        h = surf(log(s).*cos(phi+pi),log(s).*sin(phi+pi),Z,'EdgeColor','none','LineStyle','none');
        set(h,'FaceColor',[.5 .5 .5]);
        
        Z = NaN(size(s));
        Z(M>1)=-.1;
        h = surf(log(s).*cos(phi),log(s).*sin(phi),Z,'EdgeColor','none','LineStyle','none');
        set(h,'FaceColor',[.871 .922 .98]);
        h = surf(log(s).*cos(phi+pi),log(s).*sin(phi+pi),Z,'EdgeColor','none','LineStyle','none');
        set(h,'FaceColor',[.871 .922 .98]);
    end
end

t_vect = [];
psi_cell = {};
for i=1:length(balls)
    t_vect = [t_vect balls{i}.ballcenter(1)];
    psi_cell = [psi_cell [balls{i}.ballcenter(2)]];
end

end

function doubleclickcallback(obj,evt)
global balls cont_simu
switch get(obj,'SelectionType')
    case 'normal'
    case 'open'
        CurrentPoint = get(gca,'CurrentPoint');
        point = CurrentPoint(1,1:2);
        colorvec = .9*rand(1,3);
        ballind = length(balls)+1;
        tempball.ballcenter = [exp(sqrt(point(1)^2+point(2)^2)),atan(point(2)/point(1))];
        tempball.curvsup = plot(point(1),point(2),'.r');
        tempball.curvinf = plot(point(1),point(2),'.r');
        tempball.colorvec = colorvec;
        h = impoint(gca,point(1),point(2));
        setColor(h,'r');

        balls{ballind} = tempball;
        cont_simu = cont_simu + 1;
        title([num2str(cont_simu) ' Simulations']);
        
        addNewPositionCallback(h,@(h) MyPositionCallback(h,ballind));
        
end
end



function MyPositionCallback(h,ballind)
%hfun = get(gcf,'WindowButtonUpFcn');
%set(gcf,'WindowButtonUpFcn',@(src,event) mybuttonup(src,event,hfun,h));

global balls beta NUMPOINTSINTERACTIVE
phivect = linspace(0,2*pi,NUMPOINTSINTERACTIVE);
        [pointinf,pointsup] = edges_one_ball(phivect,exp(sqrt(h(1)^2+h(2)^2)),atan(h(2)/h(1)),beta);
tempball = balls{ballind}; 
        tempball.ballcenter = [exp(sqrt(h(1)^2+h(2)^2)),atan(h(2)/h(1))];
            delete(tempball.curvsup);
            delete(tempball.curvinf);
            
            colorvec = tempball.colorvec;
            tempball.curvsup = plot(log(pointsup).*cos(phivect),log(pointsup).*sin(phivect),'-','Color',colorvec);
            tempball.curvinf = plot(log(pointinf).*cos(phivect),log(pointinf).*sin(phivect),'-','Color',colorvec);
            tempball.colorvec = colorvec;
            balls{ballind} = tempball;
       

%title(sprintf('indexed at %i -> (%1.1f,%1.1f)',ind,h(1),h(2)));
end
% les A tels que d( T_2, A) <= log(sqrt(2)) 
% draw_ball(2,0,log(sqrt(2)))



function [pointinf,pointsup] = edges_one_ball(phivect,t,psi,beta)
        pointsup = [];
        pointinf = [];
        for phi = phivect
            Gphi = cos(psi-phi)^2;
            dis = beta^2 - ( Gphi/t + t*(1-Gphi) )*( (1-Gphi)/t + t*Gphi );
            switch (logical(true))
                case dis>0
                    tupper = ( beta + sqrt(dis) )/( Gphi/t + (1-Gphi)*t );
                    tlower = ( beta - sqrt(dis) )/( Gphi/t + (1-Gphi)*t );
                    if (tupper<1)
                        pointsup = [pointsup NaN];
                        pointinf = [pointinf NaN];
                    else
                        pointsup = [pointsup tupper];
                        if (tlower<1)
                            pointinf = [pointinf 1];
                        else
                            pointinf = [pointinf tlower];
                        end
                    end
                case dis==0
                    tupper = beta/( Gphi/t + (1-Gphi)*t );
                    if (tupper<1)
                        pointsup = [pointsup NaN];
                        pointinf = [pointinf NaN];
                    else
                        pointsup = [pointsup tupper];
                        pointinf = [pointinf tupper];
                    end
                case dis<0
                    pointsup = [pointsup NaN];
                    pointinf = [pointinf NaN];
            end
        end
end