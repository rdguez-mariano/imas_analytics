im_ac = double(imread('./image_BD/a-contrario-image/im3_sub.png'));
    if (size(im_ac, 3) ~= 1)
        im_ac = sum(im_ac,3)/3;
    end

listing = dir('textures/*.png');
for f=1:1%length(listing)

    im1 = double(imread(['./textures/' listing(f).name]));
    if (size(im1, 3) ~= 1)
        im1 = sum(im1,3)/3;
    end
    [m,n]= size(im1);

    [data,simi,W,vec] = get_sift_stats( im1,im_ac,0.6, 4);

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
    
[ groups, Y,newW ] = get_sift_groups(W,vec,'unn');

pause(1)
classes = unique(Y);
for i=1:length(groups)
    thisdata = groups{i};
    if(length(thisdata(1,:))==1)
        continue;
    end
    thisW = W(Y == classes(i),Y == classes(i))
    ext_links = W(Y ~= classes(i),Y == classes(i));
    ext_links = sum(sum(ext_links>0));
    
    if (ext_links>0)
        h=figure;
        imshow(im1/255);hold on;
        plot(thisdata(1,:),thisdata(2,:),'o','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor',[1 .6 .6])
        title([num2str(ext_links) ' external links'])
        for x=1:length(thisdata(1,:))
            for y=1:length(thisdata(1,:))
                if (thisW(x,y)>0)
                    line(thisdata(1,[x y]),thisdata(2,[x y]),'LineWidth',2*thisW(x,y),'Color','b');
                end
            end
        end
        
        pause
        close(h);
    end
end

%     drawnow;
%     set(get(handle(gcf),'JavaFrame'),'Maximized',1);
%     drawnow;
% 
%     print(['./textures/data/' listing(f).name],'-dpng','-r100')
%     pause(1)
%     close(h);
%     
%     h = figure;
%     %histogram(g_cont(g_cont~=2));
%     gg = unique(g_cont);
%     gc = zeros(1,size(gg,2));
%     for g=1:length(gg)
%         gc(g) = sum(g_cont==gg(g));
%     end
%     plot(gg,gc,'*:');
%     title(['mean = ' num2str(mean(g_cont)) ' / median = ' num2str(median(g_cont))])
%     xlabel('cardinality of groups')
%     ylabel('number of groups of similar size')
%     
%     drawnow;
%     set(get(handle(gcf),'JavaFrame'),'Maximized',1);
%     drawnow;
% 
%     print(['./textures/data2/grouphist_' listing(f).name],'-dpng','-r100')
%     pause(1)
%     close(h)
end
