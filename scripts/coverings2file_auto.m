% ---- Setup ------

savepath = './'; % enter here where to save image coverings

mkdir([savepath '2D']);
mkdir([savepath '3D']);

milipause = 1; % avoid matlab from crashing cauze of too many auto handling

% use this for disabling unwanted printings
covering = { };
filename = { };

% ------------ LITERATURE COVERINGS ---------------
covering = {'MODS SURF-SURF HARD', 'MODS DOG-SIFT HARD', 'MODS DOG-SIFT MEDIUM', 'ASIFT', 'FAIR-SURF simulated tilts', 'FAIR-SURF fixed tilts covering' };
filename = {'MODS_SURF_SURF_HARD', 'MODS_DOG_SIFT_HARD', 'MODS_DOG_SIFT_MEDIUM', 'ASIFT', 'FAIR_SURF_simulated_tilts', 'FAIR_SURF_fixed_tilts_covering' };

covering = {'MODS SURF-SURF HARD' };
filename = {'MODS_SURF_SURF_HARD'};

for icov=1:length(covering)
    
    [ tvec, psicell, radius, region ] = get_literature_covering(covering{icov});
    val = 0; count =0;
    for i=1:length(tvec)
        t=tvec(i);
        numphi=length(psicell{i});
        count = count + numphi;
        val = val + numphi/t;
    end
    
    
    % 2D drawing
    optionsdraw.tilt = [region];
    optionsdraw.draw_regions = true;
    optionsdraw.draw_regions_tilt = max(tvec(:))*radius;
    
    h = figure;
    draw_ball_pol(tvec,psicell,log(((radius))),optionsdraw);
    htitle = get(gca,'Title');
    title( [get(htitle,'String') ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '2D/' filename{icov}],'-dpng','-r300')
    print([savepath '2D/' filename{icov}],'-depsc','-r300')
    close(h)
    pause(milipause);
    
    % 3D drawing
    optionsdraw.tilt = [region];
    
    h = figure;
    draw_ball3d(tvec,psicell,log(radius),optionsdraw);
    title( ['simulations = ' num2str(count) ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '3D/' filename{icov}],'-dpng','-r300')
    print([savepath '3D/' filename{icov}],'-depsc','-r300')
    close(h)
    pause(milipause);
end




% ------------ NEAR OPTIMAL COVERINGS ---------------
covering = { 1.6 1.7 1.8 1.9 2 };
filename = { 'near_optimal_1_6' 'near_optimal_1_7' 'near_optimal_1_8' 'near_optimal_1_9' 'near_optimal_2' };

covering = {  2};
filename = { 'near_optimal_2'};

for icov=1:length(covering)
    radius = covering{icov};
    [ tvec, psicell, region ] = get_feasible_covering(radius);
    val = 0; count =0;
    for i=1:length(tvec)
        t=tvec(i);
        numphi=length(psicell{i});
        count = count + numphi;
        val = val + numphi/t;
    end
    
    
    % 2D drawing
    optionsdraw.tilt = [region];
    optionsdraw.draw_regions = true;
    optionsdraw.draw_regions_tilt = max(tvec(:))*radius;
    
    h = figure;
    draw_ball_pol(tvec,psicell,log(((radius))),optionsdraw);
    htitle = get(gca,'Title');
    title( [get(htitle,'String') ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '2D/' filename{icov}],'-dpng','-r300')
    print([savepath '2D/' filename{icov}],'-depsc','-r300')
    close(h)
    pause(milipause);
    
    % 3D drawing
    optionsdraw.tilt = [region];
    
    h = figure;
    draw_ball3d(tvec,psicell,log(radius),optionsdraw);
    title( ['simulations = ' num2str(count) ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '3D/' filename{icov}],'-dpng','-r300')
    print([savepath '3D/' filename{icov}],'-depsc','-r300')
    close(h)
    pause(milipause);
end

% ------------ ISOLATED DISKS ---------------

covering = { [1,0] [sqrt(2),0] [2,0] [4,0] };
filename = { 'B_1_0' 'B_sqrt2_0' 'B_2_0' 'B_4_0' };

covering = { };
filename = { };

for icov=1:length(covering)
    radius = sqrt(2);
    tvec = covering{icov}(1);
    psicell = {covering{icov}(2)};
    region = 4*sqrt(2);
    
    
    val = 0; count =0;
    for i=1:length(tvec)
        t=tvec(i);
        numphi=length(psicell{i});
        count = count + numphi;
        val = val + numphi/t;
    end
    
    
    % 2D drawing
    optionsdraw.tilt = [region];
    optionsdraw.draw_regions = true;
    optionsdraw.draw_regions_tilt = region;
    
    h = figure;
    draw_ball_pol(tvec,psicell,log(((radius))),optionsdraw);
    htitle = get(gca,'Title');
    title( [get(htitle,'String') ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '2D/' filename{icov}],'-dpng','-r300')
    print([savepath '2D/' filename{icov}],'-depsc','-r300')
    close(h)
    
    pause(milipause);
    % 3D drawing
    optionsdraw.tilt = [region];
    
    h = figure;
    draw_ball3d(tvec,psicell,log(radius),optionsdraw);
    title( ['simulations = ' num2str(count) ' / radius = ' num2str(radius) ' / Area ratio = ' num2str(val) ' / Tilt \leq ' num2str(region)]);
    
    drawnow;
    set(get(handle(gcf),'JavaFrame'),'Maximized',1);
    drawnow;
    
    print([savepath '3D/' filename{icov}],'-dpng','-r300')
    print([savepath '3D/' filename{icov}],'-depsc','-r300')
    close(h)
    pause(milipause);
end
