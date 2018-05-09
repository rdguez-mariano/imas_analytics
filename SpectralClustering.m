function SpectralClustering(W,vec,type)
thisW = W;
thisvec = vec;
[ thisY,newthisW, newthisvec, newArr] = get_sift_groups(thisW,thisvec,type);
subclasses = unique(thisY);

% by Clusters
mapnewW = newthisW;
for subi = 1:length(newArr)
    mapnewW(newArr{subi},newArr{subi}) = 0;
end

colorgroups = 0.8*rand(length(newArr),3);
countgroups = zeros(length(newArr)+1,1);
pause(0.5)
h=figure;
hold on;
for subi=1:length(newArr)
    thisdata = newthisvec(1:2,newArr{subi});
    subnewW = newthisW(newArr{subi},newArr{subi});
    if(length(thisdata(1,:))==1)
        continue;
    end
    
    plot(thisdata(1,:),thisdata(2,:),'o','MarkerSize',10,'MarkerEdgeColor','black','MarkerFaceColor',colorgroups(subi,:))
    for x=1:length(thisdata(1,:))
        for y=x:length(thisdata(1,:))
            if (subnewW(x,y)>0)
                line(thisdata(1,[x y]),thisdata(2,[x y]),'LineWidth',2*subnewW(x,y),'Color',colorgroups(subi,:));
                countgroups(subi) = countgroups(subi) + 1;
            end
        end
    end
end

thisdata = newthisvec(1:2,:);

for x=1:length(mapnewW)
    for y=x:length(mapnewW)
        if (mapnewW(x,y)>0)
            line(thisdata(1,[x y]),thisdata(2,[x y]),'LineStyle',':','LineWidth',2*mapnewW(x,y),'Color','black');
            countgroups(length(countgroups)) = countgroups(length(countgroups)) + 1;
        end
    end
end
leg = [''];
for j=1:(length(countgroups)-1)
    thiscolor = mat2str(colorgroups(j,:)); thiscolor = thiscolor(2:(length(thiscolor)-1));
    leg = [leg ' {\color[rgb]{' num2str(thiscolor) '}' [num2str(countgroups(j)) ' links}']];
end
leg = [leg ' {\color{black}' [num2str(countgroups(length(countgroups))) ' links}']];
title(leg);
end