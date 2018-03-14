function [ goodm, badm ] = similarity_test( im1, im3, patch_dim, N_rep, t, phi )
%GENERATE_TEST Summary of this function goes here
%   Detailed explanation goes here


DESCVEC = 4;
descname = {'SIFT' 'RootSIFT' 'Rafa_weights' 'Rafa_quantised' 'SURF' 'SIFT_0.6' 'RootSIFT_0.6'};
goodm = zeros(1,max(DESCVEC));
badm = zeros(1,max(DESCVEC));

LIBRARYPATH = 'LD_LIBRARY_PATH= ';
MATCHER_PATH_global = './siim-matcher';
VARIANCE = 0;%0.005;

SHOWON = true;
PERFORM_IMAS_MEX = true;
DO_MOD_CHECK = true;

im1 = sum(im1,3)/3; im1 = im1*( 255/max(im1(:)) );
im3 = sum(im3,3)/3; im3 = im3*( 255/max(im3(:)) );

[x_im1,y_im1] = size(im1);
[x_im3,y_im3] = size(im3);

x = randi([1 x_im3-patch_dim]);
y = randi([1 y_im3-patch_dim]);

patch = im3(x:x+patch_dim-1, y:y+patch_dim-1);
patch = repmat(patch, N_rep);
%imshow(patch/max(patch(:)));


% image 1
x1 = randi([1 x_im1-N_rep*patch_dim]);
y1 = randi([1 y_im1-N_rep*patch_dim]);
x1_top = x1+N_rep*patch_dim-1;
y1_top = y1+N_rep*patch_dim-1;

%figure;
im1(x1:x1_top,y1:y1_top) = patch;
im1 = perform_tilt_on_image(im1,t,phi);

im2 = im1;

if (VARIANCE>0)
    im1 = imnoise(im1/max(im1(:)),'gaussian',0,VARIANCE);
    im1 = im1*( 255/max(im1(:)) );
    
    im2 = imnoise(im2/max(im2(:)),'gaussian',0,VARIANCE);
    im2 = im2*( 255/max(im2(:)) );
end

%imshow(im1/max(im1(:)));

for DESC=DESCVEC
    
    opts_IMAS.uno = 1;
    data_matches = [];
    
    MATCHER_PATH = [MATCHER_PATH_global '/' descname{DESC}];
    mkdir(MATCHER_PATH);
    
    num_matches = 0;
    
    opts_IMAS.applyfilter = 0;
    opts_IMAS.showon = false;
    switch DESC
        case 1 % SIFT L1 0.8
            opts_IMAS.covering = 1.7;
            opts_IMAS.desc_type = 1;
            opts_IMAS.match_ratio = 0.8;
        case 2 % RootSIFT 0.8
            opts_IMAS.covering = 1.9;
            opts_IMAS.desc_type = 11;
            opts_IMAS.match_ratio = 0.8;
        case 6 % SIFT L1 0.6
            opts_IMAS.covering = 1.7;
            opts_IMAS.desc_type = 1;
            opts_IMAS.match_ratio = 0.6;
        case 7 % RootSIFT 0.6
            opts_IMAS.covering = 1.9;
            opts_IMAS.desc_type = 11;
            opts_IMAS.match_ratio = 0.6;
        case 3 % Rafa weights
            opts_IMAS.covering = 1.6;
            opts_IMAS.desc_type = 31;
        case 4 % Rafa quantised
            opts_IMAS.covering = 1.7;
            opts_IMAS.desc_type = 32;
        case 5 % SURF
            opts_IMAS.covering = 1.7;
            opts_IMAS.desc_type = 2;
    end
    opts_IMAS.covering = 1;
    
    if (opts_IMAS.desc_type > 0)
        data_matches = perform_IMAS(im1, im2,opts_IMAS);
        num_matches = size(data_matches,2);
    else
        num_matches = 0;
    end
    
    
    
    
    
    
    
    if (num_matches>0)
        data_x1 = data_matches(2,:);
        data_y1 = data_matches(1,:);
        
        data_x2 = data_matches(9,:);
        data_y2 = data_matches(8,:);
        
        x1p=0; y1p=0; x2p=0 ;y2p =0;
        good = true(1,length(data_x1));
        
        
        % change of coordinates x <-> y
        ux1=[];ux2=[];uy1=[];uy2=[];
        for i=1:length(data_x1)
            
            if (norm([data_x1(i)-data_x2(i); data_y1(i)-data_y2(i)])<5)
                good(i) = false*good(i);
            end

%             tx1 = data_x1(i);
%             ty1 = data_y1(i);
%             tx2 = data_x2(i);
%             ty2 = data_y2(i);
%             
%             [tx1,ty1,tx2,ty2] = myorder(tx1, ty1, tx2, ty2);
%             
%             for j=1:length(ux1)
%                 if ( norm([tx1-ux1;ty1-uy1])<3 && norm([tx2-ux2;ty2-uy2])<3 )
%                     good(i) = false*good(i);
%                     display('aquí')
%                 else
%                     ux1 = [ux1; tx1];
%                     ux2 = [ux2; tx2];
%                     uy1 = [uy1; ty1];
%                     uy2 = [uy2; ty2];
%                 end
%             end
%             if (length(ux1)==0)
%                 ux1 = [ux1; tx1];
%                 ux2 = [ux2; tx2];
%                 uy1 = [uy1; ty1];
%                 uy2 = [uy2; ty2];
%             end
            
            [y0,x0,inside] = get_tilted_coor( data_y1(i), data_x1(i), y_im1, x_im1, t, phi*pi/180, -1);
            if (inside)
                if ( (x1<=x0)&&(x0<=x1_top)&&(y1<=y0)&&(y0<=y1_top) )
                    x1p = mod(x0-x1,patch_dim);
                    y1p = mod(y0-y1,patch_dim);
                    good(i) = true*good(i);
                else
                    good(i) = false*good(i);
                end
            else
                good(i) = false*good(i);
            end
            
            [y0,x0,inside] = get_tilted_coor( data_y2(i), data_x2(i), y_im1, x_im1, t, phi*pi/180, -1);
            if (inside)
                if ( (x1<=x0)&&(x0<=x1_top)&&(y1<=y0)&&(y0<=y1_top) )
                    x1p = mod(x0-x1,patch_dim);
                    y1p = mod(y0-y1,patch_dim);
                    good(i) = true*good(i);
                else
                    good(i) = false*good(i);
                end
            else
                good(i) = false*good(i);
            end
            
            
        end
        
        goodm(DESC) = sum(good);
        badm(DESC) = sum(1-good);
        if (SHOWON)
            pause(0.2);
            figure;
            plot_matches_in_oneimage(im1,zeros(14,1));
            title('image');
            
            pause(0.2);
            figure;
            plot_matches_in_oneimage(im1,data_matches(:,good==true));
            %title([num2str(sum(good)) ' matches ']);
            title('inside');
            
            
            pause(0.2);
            figure;
            plot_matches_in_oneimage(im1,data_matches);
            title('all');

        end
        
        
        
        
        
    end
end


%
% %im1
% imseg = zeros(size(im1));
% for x=1:10:size(im1,1)
%     for y=1:10:size(im1,2)
%         % change of coordinates x <-> y
%         x0 = x;
%         y0 = y;
%
%             if ( (x1<=x0)&&(x0<=x1_top)&&(y1<=y0)&&(y0<=y1_top) )
%                 imseg(x,y) = 1; % good matches / target region
%             else
%                 imseg(x,y) = 0.5;
%             end
%
%     end
% end
% figure;
% imshow(imseg);
%
%
%
%
%
%
%
% %im2
% imseg = zeros(size(im2));
% for x=1:5:size(im2,1)
%     for y=1:5:size(im2,2)
%         % change of coordinates x <-> y
%         [y0,x0,inside] = get_tilted_coor( y, x, y_im2, x_im2, t, phi*pi/180, -1);
%         if (inside)
%             if ( (x2<=x0)&&(x0<=x2_top)&&(y2<=y0)&&(y0<=y2_top) )
%                 imseg(x,y) = 1; % good matches / target region
%             else
%                 imseg(x,y) = 0.5;
%             end
%         end
%
%     end
% end
% figure;
% imshow(imseg);







end

function [x1_out,x2_out]=myswap(x1,x2)
x1_out = x2;
x2_out = x1;
end

function [x1,y1,x2,y2]=myorder(x1, y1, x2, y2)

    %//(x1<x2) i do nothing
    if (x1>x2)
        [x1,x2]=myswap(x1,x2);
        [y1,y2]=myswap(y1,y2);
    end
    if ((x1==x2)&&(y1>y2))
        [x1,x2]=myswap(x1,x2);
        [y1,y2]=myswap(y1,y2);
    end

end
