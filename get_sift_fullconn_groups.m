function [ Y, newW, newVec, newArr ] = get_sift_fullconn_groups(W,vec_representants )
%SIFT_GROUPS Summary of this function goes here
%   Detailed explanation goes here
visual = true;
verbose = false;

tic;
Y = fullconnections(W);

classes = unique(Y);
newW = zeros(size(W));

groupmax = 0;
groupmin = Inf;

% New arrangement positions for the W matrix -> will be output in newW
pos = 1;
newArr = {};
for i = 1:length(classes)
    finalpos = sum(Y == classes(i))+pos-1;
    newArr = [newArr pos:finalpos];
    pos = finalpos+1;
end
newVec = -1*ones(size(vec_representants));
% Create newW
for i = 1:length(classes)
    newVec(:,newArr{i}) = vec_representants(:,Y == classes(i));
    
    groupmax = max(size(newArr{i},2),groupmax);
    groupmin = min(size(newArr{i},2),groupmin);
    newW(newArr{i},newArr{i}) = W(Y == classes(i),Y == classes(i));
    
    for j=(i+1):length(classes)
        newW(newArr{i},newArr{j}) = W(Y == classes(i),Y == classes(j));
        newW(newArr{j},newArr{i}) = W(Y == classes(j),Y == classes(i));
    end
end

if(visual)
    figure;
    subplot(1,2,1);imshow(newW>0);title('New arrangement')
    subplot(1,2,2);imshow(W>0);title('Old arrangement')
    suptitle('Similarity matrix');
end

end


function [Y] = fullconnections(W)
tic;
Y = -1*ones(1,size(W,2));

current_group_id = 1;
cind = find(Y==-1,1);
current_group = [cind];
while (~isempty(cind))
    Y(cind) = current_group_id;
    avail_conn = find(W(cind,:)>0);
    
    for a=avail_conn
        if verify_fullyconn(W,current_group,a)
            current_group = [current_group a];
            Y(a) = current_group_id;
        end
    end
    
    current_group_id = current_group_id + 1;
    cind = find(Y==-1,1);
    current_group = [cind];
end

end

function [ok] = verify_fullyconn(W,current_group,a)
inds = [current_group a];
n = length(inds);
maxval = n*(n-1);

% Condition to have a fully connected graph
if (sum(sum(W(inds,inds)>0))==maxval)
    ok = true;
else
    ok = false;
end
end

