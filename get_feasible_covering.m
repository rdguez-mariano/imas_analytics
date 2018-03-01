function [ TVEC1, PSICELL1, REGION ] = get_feasible_covering( radius )
%GET_FEASIBLE_COVERING Summary of this function goes here
%   Detailed explanation goes here

switch radius
    case 1.4
        [TVEC1,PSICELL1,~,~] = covering_asift_next2optimal2(1.4,2); REGION = 5;
    case 1.55
        t1=2.28479; phi1=0.318792; t2=4.63662; phi2=0.16015; REGION=5.5;
    case 1.6
        t1=2.37907; phi1=0.35292; t2=4.94169; phi2=0.176161;REGION=5.6;
    case 1.7
        t1=2.61211; phi1=0.3961; t2=5.33812; phi2=0.196909; REGION=5.8;
        t1=2.61212; phi1=0.392893; t2=5.37923; phi2=0.196769; REGION=5.8; % update 15/08/2017
    case 1.8
        t1=2.96171; phi1=0.35041; t2=6.91615; phi2=0.175084; REGION=7;
        t1=2.88447; phi1=0.394085; t2=6.2197; phi2=0.196389; REGION=6; % update 15/08/2017
    case 1.9
        t1=3.14247; phi1=0.393272; t2=6.73059; phi2=0.174594;REGION=8;
        t1=3.17743; phi1=0.392951; t2=6.74304; phi2=0.174692;REGION=8; % update 15/08/2017
    case 2
        t1=2.85821; phi1=0.525828; t2=5.76899; phi2=0.243276; REGION=6;
        t1=3.54845; phi1=0.357496; t2=8.76615; phi2=0.142867; REGION=10; % update 15/08/2017
end



if (radius~=1.4)
    TVEC1 = [1];
    PSICELL1 = {[0]};
    
    TVEC1 = [TVEC1 t1];
    theta = 0:phi1:(pi);
    PSICELL1 = [PSICELL1 theta];
    
    TVEC1 = [TVEC1 t2];
    theta = 0:phi2:(pi);
    PSICELL1 = [PSICELL1 theta]; 
end
end
