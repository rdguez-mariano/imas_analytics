function [Mbool, Mval, RADIUS ] = tolerance_tests_SIIMS(tvec,psivec,opts )

RADIUS = 1;
REGION = 4*sqrt(2);
DRAW_ON = true;
LIBRARYPATH = 'LD_LIBRARY_PATH= ';
MATCHER_PATH = './siim-matcher';
PERFORM_IMAS_MEX = true;

switch nargin
    case 3
        if (~isstruct(opts))
            error('Optional parameter is not of type struct!');
        end
        names = fieldnames(opts);
        for i=1:length(names)
            name = names{i};
            switch name
                case 'draw_on'
                    DRAW_ON = opts.draw_on;
                case 'radius'
                    RADIUS = opts.radius;
                case 'region'
                    REGION = opts.region;
            end
        end
    otherwise
        if (nargin~=2)
            error('Wrong number of input parameters!');
        end
end




optionsdraw.tilt = [REGION];
optionsdraw.draw_regions = true;
optionsdraw.tilt = [REGION];
optionsdraw.draw_regions_tilt = REGION;
optionsdraw.interactive = false;
optionsdraw.draw_center_disks_opts = '.k';

listing = dir('image_BD/*.png');
Mbool = ones(length(tvec),length(psivec))*true;
Mval = zeros(length(tvec),length(psivec));
for i=1:length(listing)
    
    im = double(imread(['./image_BD/' listing(i).name]));
    if (size(im,3)==3)
        im = sum(im,3)/3;
    end
    im = 255*im/max(im(:));
    
    disp('**************************************');
    disp(['********* Selected image is ' listing(i).name]);
    disp('**************************************');
    
    if (~PERFORM_IMAS_MEX)
        imwrite(im/max(im(:)),[MATCHER_PATH '/im2.png']);
    end
    
    for ti=1:length(tvec)
        for psi1i=1:length(psivec)
            t = tvec(ti);
            psi1 = psivec(psi1i);
            
            im1 = perform_tilt_on_image(im,t,0);
            im1 = perform_tilt_on_image(im1,1,psi1);
           
            if (~PERFORM_IMAS_MEX)
                imwrite(im1/max(im1(:)),[MATCHER_PATH '/im1.png']);
            end
            
            num_matches = 0;
            if (~PERFORM_IMAS_MEX)
                % EXECUTION OF THE SIIM-MATCHER
                % ----------- START -------------
                delete([MATCHER_PATH '/data_matches.csv']);
                % weights system([LIBRARYPATH ' cd ' MATCHER_PATH ' && ' LIBRARYPATH ' ./main -im1 im1.png -im2 im2.png -applyfilter 2 -match_type 1 -covering 1 -AC_type 1' ]);
                % quantised system([LIBRARYPATH ' cd ' MATCHER_PATH ' && ' LIBRARYPATH ' ./main -im1 im1.png -im2 im2.png -applyfilter 2 -match_type 1 -covering 1 -AC_type 2 -new_ori_size 18 -step_sigma 1.6 -quant_prec 0.28' ]);
                
                %system([LIBRARYPATH ' cd ' MATCHER_PATH ' && ' LIBRARYPATH ' ./main -im1 im1.png -im2 im2.png -applyfilter 2 -match_type 1 -covering 1 -AC_type 1' ]);
                system([LIBRARYPATH ' cd ' MATCHER_PATH ' && ' LIBRARYPATH ' ./main -im1 im1.png -im2 im2.png -applyfilter 2 -match_type 1 -covering 1 -AC_type 2 -new_ori_size 22 -step_sigma 1.5 -quant_prec 0.32' ]);
                [~,cmdout] = system([LIBRARYPATH 'cat ' MATCHER_PATH '/data_matches.csv | wc -l' ]);
                
                if (exist([MATCHER_PATH '/data_matches.csv'], 'file') == 2)
                    num_matches = str2num(cmdout)-1;
                end
            else    % INTERNAL MEX FUNCTION
                opts_IMAS.covering = 1;
                opts_IMAS.desc_type = 11;
                opts_IMAS.applyfilter = 2;
                opts_IMAS.match_ratio = 0.8;
                opts_IMAS.showon = false;
                
                
                data_matches = perform_IMAS(im1, im,opts_IMAS);
                num_matches = size(data_matches,2);
            end
            
            
            if (num_matches>0)
                Mbool(ti,psi1i) = true*Mbool(ti,psi1i);
            else
                Mbool(ti,psi1i) = false*Mbool(ti,psi1i);
            end
            Mval(ti,psi1i) = num_matches+Mval(ti,psi1i);
        end
        
    end
end

maxti = 1;
for ti=1:length(tvec)
    if (maxti<tvec(ti))
        if ( Mbool(ti,:)*Mbool(ti,:)'==length(Mbool(ti,:)) )
            maxti = tvec(ti);
        else
            break;
        end
    end
end





%save(FILENAME);
if (maxti>0)
    RADIUS = maxti;
end

if (DRAW_ON)
    hold on;
    TVEC1 = 1;
    PSICELL1 = {0};
    draw_ball_pol(TVEC1,PSICELL1,log(RADIUS),optionsdraw);
    for ti=1:length(tvec)
        for psi1i=1:length(psivec)
            t = tvec(ti);
            psi1 = psivec(psi1i);
            psi1 = mod(psi1 ,180.1);
            if (Mbool(ti,psi1i)==true)
                plot(log(t).*cos(pi*psi1/180),log(t).*sin(pi*psi1/180),'.b','MarkerSize',12);
            else
                plot(log(t).*cos(pi*psi1/180),log(t).*sin(pi*psi1/180),'xr','MarkerSize',8);
            end
        end
    end
end

end
