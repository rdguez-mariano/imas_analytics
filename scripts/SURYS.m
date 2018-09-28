im_ac = double(imread('./image_BD/a-contrario-image/im3_sub.png'));
if (size(im_ac, 3) ~= 1)
    im_ac = sum(im_ac,3)/3;
end

listing = dir('textures/coca.png');
for f=1:1%length(listing)
    
    im1 = double(imread(['./textures/' listing(f).name]));
    if (size(im1, 3) ~= 1)
        im1 = sum(im1,3)/3;
    end
    [m,n]= size(im1);
    
    [data,simi,W,vec] = get_sift_stats( im1,im_ac,0.3, 10);
    
    Y = vec(4,:);
    classes = unique(Y);
    for i=1:length(classes)
        thisW = W(Y == classes(i),Y == classes(i));
        thisvec = vec(:,Y == classes(i));
        N = sum(Y == classes(i)); %number of elements in the class
        N_1 = N/2; %number of elements in the class
        Nconn_0 = N*(N+1)/2 - N; %max number of connections between N elements
        Nconn_1 = N_1*(N_1+1)/2 - N_1; %%max number of connections between N-1 elements
        if (N>4)
            % Spectral Clustering groups
            %[ thisY,newthisW, newthisvec, newArr] = get_sift_groups(thisW,thisvec,'sym');
            
            % Fully connected groups
            [ thisY,newthisW, newthisvec, newArr] = get_sift_fullconn_groups(thisW,thisvec);
            
            subclasses = unique(thisY);
            
            if (length(newArr)>1)
                % All links
                thisdata = thisvec;
                pause(0.5)
                h=figure;
                count = 0;
                count2 = 0;
                imshow(im1/255);hold on;
                plot(thisdata(1,:),thisdata(2,:),'o','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
                for x=1:length(thisdata(1,:))
                    for y=x:length(thisdata(1,:))
                        if (thisW(x,y)>0)
                            line(thisdata(1,[x y]),thisdata(2,[x y]),'LineWidth',2*thisW(x,y),'Color','b');
                            count = count + 1;
                        end
                        if(newthisW(x,y)>0)
                            count2 = count2 + 1;
                        end
                    end
                end
                title([num2str(N) ' nodes with ' num2str(count) ' connections'])
                
                % by Clusters
                mapnewW = newthisW;
                for subi = 1:length(newArr)
                    mapnewW(newArr{subi},newArr{subi}) = 0;
                end
                
                colorgroups = 0.8*rand(length(newArr),3);
                countgroups = zeros(length(newArr)+1,1);
                pause(0.5)
                h=figure;
                imshow(im1/255);hold on;
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
                pause;
                %close(h);
            end
            close all
        end
        
        ext_links = W(Y ~= classes(i),Y == classes(i));
        ext_links = sum(sum(ext_links>0));
        if (ext_links>0)
            disp('problem');
        end
        
    end
    
    
    
    
    datax = data(1,:);
    datay = data(2,:);
    dataangle = data(3,:);
    dataoct = data(4,:);
    
    octmin = -1;
    octmax = 4;
    
    
    x_div = 3;
    y_div = 3;
    
    x_block = n/x_div;
    y_block = m/y_div;
    h = figure;
    hold on;
    for i= 0:(x_div-1)
        for j=0:(y_div-1)
            indfilter = true(1,size(data,2));
            xmin = round(i*x_block + 1);
            xmax = round((i+1)*x_block);
            ymin = round(j*y_block + 1);
            ymax = round((j+1)*y_block);
            
            % [x_min,x_max]
            indfilter = indfilter .* (datax>=xmin) .* (datax<=xmax);
            %[y_min,y_max]
            indfilter = indfilter .* (datay>=ymin) .* (datay<=ymax);
            
            otcdensity = zeros(1,octmax-octmin+1);
            for oct = unique(dataoct(indfilter==true))
                otcdensity(round(oct+2)) = (4^(oct+1))*1000*(sum(dataoct(indfilter==true)==oct))/(m*n);
            end
            subplot(y_div,x_div, 1+i+j*y_div );
            if (length(otcdensity)>length(octmin:octmax))
                otcdensity = otcdensity(1:length(octmin:octmax));
            end
            bar(octmin:octmax,otcdensity);
            
            
            ylabel('density');
            xlabel('octave');
            
            
        end
    end
end
