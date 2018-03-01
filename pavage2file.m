function [ output_args ] = pavage2file( tvec1,tvec2, psicell1,psicell2 )
%PAVAGE2FILE Summary of this function goes here
%   Detailed explanation goes here
    
if (exist('/tmp/pavage.txt')==2)
    [status, ~] = system('rm /tmp/pavage.txt');
end

tvec = tvec1;
psicell = psicell1;
pavage = [];
for j=1:length(tvec)
    t=tvec(j);
    for i=1:length( psicell{j} )
        psi = psicell{j}(i);
        pavage = [pavage [t;psi]];
    end
end
dlmwrite('/tmp/pavage.txt',pavage,'delimiter',',','precision',10);


tvec = tvec2;
psicell = psicell2;
pavage = [];
for j=1:length(tvec)
    t=tvec(j);
    for i=1:length( psicell{j} )
        psi = psicell{j}(i);
        pavage = [pavage [t;psi]];
    end
end
%csvwrite('pavage.txt',pavage)
dlmwrite('/tmp/pavage.txt',pavage,'-append','delimiter',',','precision',10);

end

