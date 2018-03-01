function draw_ball3d(tvec,psicell,r,opts)

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
                case 'tilt'
                    DRAWOPT = true;
                    PLOTTILT = opts.tilt;
                case 'drawtitle'
                    DRAWTITLE = opts.drawtitle;   
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


discretise = 1600;%1000;
phivect = linspace(0,2*pi,NUMPOINTS);
[X,Y] = meshgrid(linspace(-1,1,discretise),linspace(-1,1,discretise));
gamma = exp(r);
beta = ( (gamma^2)+1 )/(2*gamma);
phi = atan2(Y, X);
rho2 = X.^2 + Y.^2;
s = 1 ./sqrt(1 - rho2);

%figure;
colormap parula%parula
Z = sqrt(max(0,1-rho2));
h=surf(X,Y,Z,'EdgeColor','none','LineStyle','none');
shading interp
%set(h,'FaceColor',[0 .7 .7])
camlight left
% lighting phong, light
axis equal
hold on
incrpatch = 1.01;
incrpoint = 1.012;
for i=1:length(tvec)
    t=tvec(i);
    psivec=psicell{i};
    for j=1:length(psivec)
        psi = psivec(j);
        G = ( (t ./s + s ./t).*(cos(phi-psi).^2) + ((1 ./(s.*t)) +s.*t).*(sin(phi-psi).^2) )/2;
        
        %imshow((G<beta)&(rho2<1));
        

        Z = sqrt(max(0,1-rho2));
        Z((G>beta)) = NaN;
        Z(rho2>1) = NaN;
        %incrpatch = incrpatch+0.01;
        %h = surf(X*incrpatch,Y*incrpatch,Z*incrpatch,'EdgeColor','none','LineStyle','none');'FaceColor','interp','FaceLighting','gouraud'
        h = surf(X*incrpatch,Y*incrpatch,Z*incrpatch,'EdgeColor','none');
        set(h,'FaceColor',[.0 .0 .0]);
        %set(h,'FaceColor',[.58 .39 .39]);
        
        %draw borders in red
        
        [tinf,tsup] = edges_one_ball(phivect,t,psi,beta);
        thetainf = (acos(1 ./tinf));
        thetasup = (acos(1 ./(tsup)));
        thetainf(thetainf == Inf) = NaN;
        thetasup(thetasup == Inf) = NaN;
        %thetasup = mod(pi/2-thetasup,pi/2);
        tempball.curvsup = plot3(sin(thetainf).*cos(phivect)*incrpatch,sin(thetainf).*sin(phivect)*incrpatch,cos(thetainf)*incrpatch,'-','LineWidth',DRAW_LineSize,'Color',red);
        tempball.curvinf = plot3(sin(thetasup).*cos(phivect)*incrpatch,sin(thetasup).*sin(phivect)*incrpatch,cos(thetasup)*incrpatch,'-','LineWidth',DRAW_LineSize,'Color',red);
        
        %draw centers
        rho = sqrt(1-1/(t^2));
        psi= pi-psi;
        A = rho*[cos(psi) sin(psi);-sin(psi) cos(psi)]*[1;0];
        x=A(1);
        y=A(2);
        z=1/t;
        
        plot3(x*incrpoint,y*incrpoint,z*incrpoint,'.','MarkerSize',DRAW_CENTER_MarkerSize,'Color',green);
    end
end
%incrpatch = incrpatch + 0.01;
incrpatch1 = incrpatch +0.2; 
if (DRAWOPT)
    for i=1:length(PLOTTILT)
        theta = (acos(1 ./PLOTTILT(i)));
        plot3(sin(theta)*cos(phivect)*incrpatch,sin(theta)*sin(phivect)*incrpatch,cos(theta).*ones(1,length(phivect))*incrpatch,'--w','LineWidth',DRAW_DashedLineSize);
        
        phivect = -pi/2-pi/4;
        theta1 = theta;
        textpointer = text(sin(theta1)*cos(phivect)*incrpatch1,sin(theta1)*sin(phivect)*incrpatch1,cos(theta1).*ones(1,length(phivect))*incrpatch1,[PLOTTILT_NAME{i}]);
        textpointer.FontWeight='bold';
        textpointer.FontSize = DRAW_CENTER_MarkerSize;
        textpointer.Color = [.8 .8 .8];
    end
end
%view([90 0])
%title(['d( T_{' num2str(t) '}R_{' num2str(psi) '} , T_sR_{\phi+' num2str(psi) '} ) \leq ' num2str(r)]);
end

% les A tels que d( T_2, A) <= log(sqrt(2))
% draw_ball(2,0,(sqrt(2)))

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