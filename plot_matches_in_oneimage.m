function plot_matches_in_oneimage(im1,data_matches,opts)

% Default Parameters
HORI = false;
add_title = 0;
BCOLOR = 255;%150;
name_algorithm = '';


% Input Parameters
switch nargin
    case 3
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
        if (nargin~=2)
            error('Wrong number of input parameters!');
        end
end    
    
    
    [w1,h1] = size(im1);
    
        wo =  w1;
        ho = h1;
    
    % images in horizontal
    imv = im1/max(im1(:));
    
    if ((size(data_matches,1)==14))
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
            Y = [ data_matches(1,:); data_matches(8,:)];%[x1;x2]
            X = [ data_matches(2,:); data_matches(9,:)];%[y1;y2]
        
        line(Y,X);%,'Color',[.0 .0 .9]);
        pause(0.2)
    end
end
