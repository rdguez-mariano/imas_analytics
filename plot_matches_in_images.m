function plot_matches_in_images(im1,im2,data_matches,opts)

% Default Parameters
HORI = false;
add_title = 0;
BCOLOR = 255;%150;
name_algorithm = '';


% Input Parameters
switch nargin
    case 4
        if (~isstruct(opts))
            error('Optional parameter is not of type struct!');
        end
        names = fieldnames(opts);
        for i=1:length(names)
            name = names{i};
            switch name
                case 'subtitle'
                    add_title = opts.subtitle;
                case 'name_algorithm'
                    name_algorithm = opts.name_algorithm;
                case 'hori'
                    HORI = opts.hori;
                case 'background'
                    BCOLOR = opts.background;
            end
        end
    otherwise
        if (nargin~=3)
            error('Wrong number of input parameters!');
        end
end    
    
    
    [w1,h1] = size(im1);
    [w2,h2] = size(im2);
    band_w = 20;
    
    if HORI
        wo =  max(w1,w2);
        ho = h1+h2+band_w;
    else 
        ho =  max(h1,h2);
        wo = w1+w2+band_w;        
    end
    
    % images in horizontal
    imv = BCOLOR*ones(wo,ho)/255;
    imv(1:w1,1:h1) = im1/max([im1(:); im2(:)]);
    
    if HORI
        imv( 1:w2 ,(h1+band_w+1):end)=im2/max([im1(:); im2(:)]);
    else
        imv( (w1+band_w+1):end ,1:h2)=im2/max([im1(:); im2(:)]);
    end
    
    if ((size(data_matches,1)==14)) % not 8 is then 1 -> ASIFT
        hold on;
        imshow(imv);
        if isstr(add_title)
            title({[num2str(size(data_matches,2)) ' matches ' add_title]});
        else
            title([name_algorithm num2str(size(data_matches,2)) ' matches']);
        end
       
        if isempty(data_matches)
            return;
        end
        %matches are of the type (x1,y1)->(x2,y2)
        if HORI
            Y = [ data_matches(1,:); data_matches(8,:)+h1+band_w];%[x1;x2]
            X = [ data_matches(2,:); data_matches(9,:)];%[y1;y2]
        else
            Y = [ data_matches(1,:); data_matches(8,:)];%[x1;x2]
            X = [ data_matches(2,:); data_matches(9,:)+w1+band_w];%[y1;y2]
        end
        
        line(Y,X);%,'Color',[.0 .0 .9]);
        pause(0.2)
    end
end
