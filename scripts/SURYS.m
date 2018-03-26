listing = dir('textures/*.png');
for f=1:length(listing)
    
    im1 = double(imread(['./textures/' listing(f).name]));
    if (size(im1, 3) ~= 1)
        continue
    end
    [m,n]= size(im1);
    
    data = get_sift_stats( im1);
    
    datax = data(1,:);
    datay = data(2,:);
    dataangle = data(3,:);
    dataoct = data(4,:);
    dataratio = data(8,:);
    datagroups = data(9,:);
    
    
    g_vec=unique(datagroups);
    gind = size(g_vec,2);
    g_cont = zeros(1,gind);
    for g = 1:gind
        g_cont(g) = sum(datagroups==g_vec(g));
    end
    g_cont = g_cont(g_vec~=0);
    
    octmin = -1;
    octmax = 8;
    
    
    x_div = 3;
    y_div = 3;
    
    x_block = n/x_div;
    y_block = m/y_div;
    figure; 
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
                otcdensity(round(oct+2)) = 1000*(sum(dataoct(indfilter==true)==oct))/(m*n);
            end
            subplot(y_div,x_div, 1+i+j*y_div );
            bar(octmin:octmax,otcdensity);

            
            ylabel('density');
            xlabel('octave');
            
            
        end
    end
    suptitle(['groups_{min} = ' num2str(min(g_cont)) ' / groups_{mean} ' num2str(mean(g_cont)) ' / groups_{max} = ' num2str(max(g_cont)) ]);
 
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print(['./textures/data/' listing(f).name],'-dpng','-r300')
end
