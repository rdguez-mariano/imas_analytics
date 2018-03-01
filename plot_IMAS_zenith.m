function [ output_args ] = plot_IMAS_zenith( data_matches, opts )
%   PLOT_IMAS_ZENITH Summary of this function goes here
%   Detailed explanation goes here

global UNIT_BALL_RADIUS im1_center im2_center

threshold = 15;
lw = 2;
UNIT_BALL_RADIUS = 8;
markersize = 15;

im1_center = [-3 0];
im2_center = [3 0];
GRAY = false;

tilts_num = 7;
tilts_vect = 1;
for i=2:tilts_num
    tilts_vect =[tilts_vect tilts_vect(end)*sqrt(2)];
end


switch nargin
    case 2
        if (~isstruct(opts))
            error('Optional parameter is not of type struct!');
        end
        names = fieldnames(opts);
        for i=1:length(names)
            name = names{i};
            switch name
                case 'threshold'
                    threshold = opts.threshold;
                case 'linewidth'
                    lw = opts.linewidth;
                case 'markersize'
                    markersize = opts.markersize;                    
                case 'tilts_num'
                    tilts_num = opts.tilts_num;                    
                case 'im1_center'
                    im1_center = opts.im1_center;  
                case 'im2_center'
                    im2_center = opts.im2_center;  
                case 'tilts_vect'
                    tilts_vect = opts.tilts_vect; 
                case 'gray'
                    GRAY = opts.gray;                     
                case 'tilts_num'
                    tilts_num = opts.tilts_num;  
                    tilts_vect = 1;
                    for i=2:tilts_num
                        tilts_vect =[tilts_vect tilts_vect(end)*sqrt(2)];
                    end
            end
        end
    otherwise
        if (nargin~=1)
            error('Wrong number of input parameters!');
        end
end  



figure; hold on;
axis equal tight;
title({'IMAS matches per simulation' 'Zenith view'});
plot_structure(tilts_vect);


    data = data_matches; % data_matches
    tilts1_data = data(6,:);
    tilts1 = unique(tilts1_data);
 plotdata.phi1=[];   
 plotdata.phi2=[];   
 plotdata.t1=[];
 plotdata.t2=[];
 plotdata.score=[];
 plotdata.phi2=[];
 mscore = 0;
  cont = 0;
    for i=1:length(tilts1)
        global_active_data = data(repmat(tilts1(i)==tilts1_data,[14 1]));
        global_active_data = reshape( global_active_data, [14 size(global_active_data,1)/14]);
        tilts2_data = global_active_data(13,:);
        tilts2 = unique(tilts2_data);

        for j=1:length(tilts2)
            active_data = global_active_data(repmat(tilts2(j)==tilts2_data,[14 1]));
            active_data = reshape( active_data, [14 size(active_data,1)/14]);
            
                    %prendre le maximum de matches pour un pair (phi1,phi2)
                    data_phi1 = active_data(7,:);
                    data_phi2 = active_data(14,:);
                    for phi1=unique(data_phi1)
                        for phi2=unique(data_phi2)
                            
                            data_sim = active_data(repmat( (data_phi1==phi1).*(data_phi2==phi2)==1 ,[14 1]));
                            data_sim = reshape( data_sim, [14 size(data_sim,1)/14]);
                            
                            if (( length(unique(data_sim(6,:)))<=1)&&(length(unique(data_sim(13,:)))<=1)&&...
                                    (length(unique(data_sim(7,:)))<=1)&&(length(unique(data_sim(14,:)))<=1) )
                                if ((( length(unique(data_sim(6,:)))==1)&&(length(unique(data_sim(13,:)))==1)&&...
                                        (length(unique(data_sim(7,:)))==1)&&(length(unique(data_sim(14,:)))==1) ))
                                    if (~( (unique(data_sim(6,:))==tilts1(i))&&(unique(data_sim(13,:))==tilts2(j))&&...
                                            (unique(data_sim(7,:))==phi1)&&(unique(data_sim(14,:))==phi2) ) )
                                        error('Assertion Error: Data is not comming from expected fixed t1, t2, phi1, phi2')
                                    end
                                end
                            else                                
                                error('Assertion Error: Handling mixed data!')
                            end
                            cont = cont + size(data_sim,2);
                            
                            score =  size(data_sim,2);
                            
                            if (mscore<score)
                                mscore = score;
                            end
                            
                            plotdata.t1 = [plotdata.t1 log(tilts1(i))];
                            plotdata.t2 = [plotdata.t2 log(tilts2(j))];
                            plotdata.phi1 = [plotdata.phi1 phi1];
                            plotdata.phi2 = [plotdata.phi2 phi2];
                            plotdata.score = [plotdata.score score];
                                                        
                        end
                    end
            

        end
    end
if (GRAY)    
    cvalues = colormap(flipud(gray(mscore)));
else
    cvalues = colormap(jet(mscore));
end
caxis([0 mscore])
h = colorbar;
set(get(h,'title'),'string','Number of matches');
    for i=1:length(plotdata.t1)

        if (plotdata.score(i)>threshold)
            [x1,y1] = coord(plotdata.t1(i),plotdata.phi1(i),im1_center,'right');
            [x2,y2] = coord(plotdata.t2(i),plotdata.phi2(i),im2_center,'left');
            plot3(x1,y1,plotdata.score(i),'.k','MarkerSize',markersize);
            plot3(x2,y2,plotdata.score(i),'.k','MarkerSize',markersize);
            %contour3([x1,x2],[y1,y2],[score,score]);
            line([x1,x2],[y1,y2],[plotdata.score(i),plotdata.score(i)],'LineWidth',lw,'Color',cvalues(plotdata.score(i),:));
        end
    end
    
if (cont~=size(data,2))
    error('Assertion Error: All data not being covered!')
end
end

function [x,y] = coord(t,phi,center,direction)
phi = phi*pi/180;
switch direction
    case 'right'
        ang = phi-pi/2; 
    case 'left'
        ang = (pi-phi)+pi/2; 
end
x=t*cos(ang)+center(1);
y=t*sin(ang)+center(2);
end


function semicircle(vec,r,direction)
%vec = [x,y] where x and y are the coordinates of the center of the circle
%r is the radius of the circle
%you might notice imperfections (not very smooth)
x=vec(1);y=vec(2);
switch direction
    case 'right'
        ang = linspace(-pi/2,pi/2,100); 
    case 'left'
        ang = linspace(pi/2,3*pi/2,100); 
end
xp=r*cos(ang);
yp=r*sin(ang);
plot(x+xp,y+yp,'-k');
end

function drawangles(vec,r,direction)
%vec = [x,y] where x and y are the coordinates of the center of the circle
%r is the radius of the circle
%you might notice imperfections (not very smooth)
switch direction
    case 'left'
        da = 'right';
    case 'right'
        da = 'left';
end
x=vec(1);y=vec(2);
for phi= linspace(0,180,10)
    [xp,yp] = coord(r,phi,vec,direction);
    plot([x xp],[y yp],':','Color',[.1 .1 .1]);
    text(xp,yp,[num2str(phi)],'HorizontalAlignment',da);
end

end

function plot_structure(tilts_vect)
global im1_center im2_center
plot(im1_center(1),im1_center(2),'.k');
plot(im2_center(1),im2_center(2),'.k');
ytick = []; yticklabel = [];
for t=log(tilts_vect)
    semicircle(im1_center,t,'right');
    semicircle(im2_center,t,'left');
    if (t~=0)
        ytick = [ytick t -t];
        yticklabel = [yticklabel {[ num2str(exp(t))] [ num2str(exp(t))]}];
    else
        ytick = [ytick t];
        yticklabel = [yticklabel {[ num2str(exp(t))]}];
    end
end
    t = log(tilts_vect(length(tilts_vect))+0.5);
    drawangles(im1_center,t,'right');
    drawangles(im2_center,t,'left');
    ax = gca;
    [~,ind] = sort(ytick);
    
    ax.YTick = ytick(ind);
    ax.YTickLabel = {yticklabel{ind}};
    ylabel('Applied Tilts');
end