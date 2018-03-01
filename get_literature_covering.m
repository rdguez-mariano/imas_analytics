function [ tvec, psicell, radius, region ] = get_literature_covering( name )
%GET_LITERATURE_COVERING Summary of this function goes here
%   In here you will find the 



switch name
    case 'ASIFT'
        radius = 1.8;
        region = 5.5;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2)];
        fixed_angle = 72;
    case 'ASIFT8'
        radius = 1.8;
        region = 7.7;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2) 8];
        fixed_angle = 72;
    case 'FAIR-SURF fixed tilts covering'
        radius = 1.5;
        region = 1.7;
        fixed_tilts = [1 2*sqrt(3)/3 sqrt(2) 2 4];
        fixed_angle = 72;
    case 'FAIR-SURF simulated tilts'
        radius = 1.5;
        region = 1.65;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2)];
        fixed_angle = 72;
    case 'MODS DOG-SIFT HARD'
        radius = 1.8;
        region = 9.6;
        fixed_tilts = [1 2 4 6 8];
        fixed_angle = 60;
    case 'MODS DOG-SIFT MEDIUM'
        radius = 1.8;
        region = 1.8;
        fixed_tilts = [1 2 3 4 5 6 7 8 9];
        fixed_angle = 180;
    case 'MODS DOG-SIFT EASY'
        radius = 1.8;
        region = 1.8;
        fixed_tilts = [1 5 9];
        fixed_angle = 360;
    case 'MODS SURF-SURF EASY'
        radius = 1.5;
        region = 1.5;
        fixed_tilts = [1 5 9];
        fixed_angle = 360;
    case 'MODS SURF-SURF MEDIUM'
        radius = 1.5;
        region = 1.5;
        fixed_tilts = [1 2 3 4 5 6 7 8 9];
        fixed_angle = 360;
    case 'MODS SURF-SURF HARD'
        radius = 1.5;
        region = 1.5;
        fixed_tilts = [1 2 3 4 5 6 7 8 9];
        fixed_angle = 72;
    case 'MODS SURF-FREAK EASY'
        radius = 1.8;
        region = 9.6;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2) 8];
        fixed_angle = 60; % must be 360
    case 'MODS SURF-FREAK MEDIUM'
        radius = 1.8;
        region = 2.95;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2) 8];
        fixed_angle = 90;
    case 'MODS SURF-FREAK HARD'
        radius = 1.8;
        region = 7.4;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2) 8];
        fixed_angle = 72;
    case 'MODS ORB EASY'
        radius = 1.5;
        region = 1.5;
        fixed_tilts = [1 5 9];
        fixed_angle = 360;
    case 'MODS ORB MEDIUM'
        radius = 1.5;
        region = 1.5;
        fixed_tilts = [1 2 3 4 5 6 7 8 9];
        fixed_angle = 90;
    case 'MODS ORB HARD'
        radius = 1.5;
        region = 1.6;
        fixed_tilts = [1 sqrt(2) 2 2*sqrt(2) 4 4*sqrt(2) 8];
        fixed_angle = 72;
        
end

tvec = [];
psicell = {};

angle_stop = pi;

for n=1:length(fixed_tilts)
    t=fixed_tilts(n);
    tvec = [tvec t];
    theta = [];
    if (t~=1)
        delta_theta = (pi*fixed_angle/180)/t;
        theta0 = 0:delta_theta:(angle_stop-5*pi/180);
        for rr = 1:length(theta0)
            theta = [theta theta0(rr)];
        end
    else
        theta = 0;
    end
    psicell = [psicell theta];
end

end

