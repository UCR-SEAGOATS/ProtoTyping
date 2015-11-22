function [ J ] = shapeSegmentation( I )
%UNTITLED2 Summary of this function goes here
%  I = rgb color image
%  J = Binary output
%
%   TO DO: experiment with improvements namely, more accurate selection of
%   seed points for the region growing.

%impliment color constancy
% I(:,:,1) = histeq(I(:,:,1));
% I(:,:,2) = histeq(I(:,:,2));
% I(:,:,3) = histeq(I(:,:,3));

hsv = rgb2hsv(I);

% hsv(:,:,1) = histeq(hsv(:,:,1));

% hsv(:,:,1) = medfilt2(hsv(:,:,1) , [7 7], 'symmetric');

norm_hsv = uint8(round(hsv*255));

sz = size(hsv(:,:,1));

K = zeros(sz);
J = zeros(sz);
for i = 1:sz(1)
    for j = 1:sz(2)
        K(i,j) = blackWhite(norm_hsv(i,j,:));
    end
end

se = strel('diamond', 3);
K = imerode(K, se);

[x,y] =  findseed(K);

if isempty(x) || isempty(y)
    return;
end

for i=1:size(x)
    tmp_img = regiongrowing(double(hsv(:,:,1)),...
        y(i), x(i), .01 );
    J = tmp_img | J;
end


se = strel('diamond', 7);
J = imerode(J, se);
J = bwconvhull(J, 'objects');

%in order to segement letter follow these steps



    
function b = blackWhite( hsv )
%   h is a vector of the h, s, and v values at point i,j
%   b is the binary value

h = hsv(1);
s = hsv(2);
v = hsv(3);

if h == 15
    if s < 40
        b = 0;
    else
        if v < 50 || v > 230
            b = 0;
        else
            b = 1;
        end
    end
else
    b = 0;
end

function [x , y] = findseed(J)
% J is a binary image
% x & y are arrays of coordinates

w = size(J,1);
h = size(J,2);

num_w = floor(w/16);
num_h = floor(h/16);
x = [];
y = [];

for i=1:num_w
    for j=1:num_h
        if i == 1 && j == 1
            small_J = J(1:16, 1:16);
            
            num_pixels = sum(sum(small_J));
            if num_pixels >= 60
                x = 8;
                y = 8;
            end
        elseif j == 1
            small_J = J(((i-1)*16):(i*16) , 1:16);
            
            num_pixels = sum(sum(small_J));
            if num_pixels >= 60
                index_x = size(x)+1;
                index_y = size(y)+1;

                y(index_x) = (i-i)*16+8;
                x(index_y) = 8;
            end
        elseif i == 1
            small_J = J(1:16,((j-1)*16):(j*16));
            
            num_pixels = sum(sum(small_J));
            if num_pixels >= 60
                index_x = size(x)+1;
                index_y = size(y)+1;

                y(index_x) = 8;
                x(index_y) = 8+16*(j-1);

            end
        else
            small_J = J(((i-1)*16):(i*16),((j-1)*16):(j*16));
            
            num_pixels = sum(sum(small_J));
            if num_pixels >= 60
                index_x = size(x)+1;
                index_y = size(y)+1;

                y(index_x) = 8+16*(i-1);
                x(index_y) = 8+16*(j-1);
            end
        end
    end
end





function J=regiongrowing(I,x,y,reg_maxdist)
% This function performs "region growing" in an image from a specified
% seedpoint (x,y)
%
% J = regiongrowing(I,x,y,t) 
% 
% I : input image 
% J : logical output image of region
% x,y : the position of the seedpoint (if not given uses function getpts)
% t : maximum intensity distance (defaults to 0.2)
%
% The region is iteratively grown by comparing all unallocated neighbouring pixels to the region. 
% The difference between a pixel's intensity value and the region's mean, 
% is used as a measure of similarity. The pixel with the smallest difference 
% measured this way is allocated to the respective region. 
% This process stops when the intensity difference between region mean and
% new pixel become larger than a certain treshold (t)
%
% Example:
%
% I = im2double(imread('medtest.png'));
% x=198; y=359;
% J = regiongrowing(I,x,y,0.2); 
% figure, imshow(I+J);
%
% Author: D. Kroon, University of Twente

if(exist('reg_maxdist','var')==0), reg_maxdist=0.2; end
if(exist('y','var')==0), figure, imshow(I,[]); [y,x]=getpts; y=round(y(1)); x=round(x(1)); end

J = zeros(size(I)); % Output 
Isizes = size(I); % Dimensions of input image

reg_mean = double(I(x,y)); % The mean of the segmented region
reg_size = 1; % Number of pixels in region

% Free memory to store neighbours of the (segmented) region
neg_free = 10000; neg_pos=0;
neg_list = zeros(neg_free,3); 

pixdist=0; % Distance of the region newest pixel to the regio mean

% Neighbor locations (footprint)
neigb=[-1 0; 1 0; 0 -1;0 1];

% Start regiogrowing until distance between regio and posible new pixels become
% higher than a certain treshold
while(pixdist<reg_maxdist&&reg_size<numel(I))

    % Add new neighbors pixels
    for j=1:4,
        % Calculate the neighbour coordinate
        xn = x +neigb(j,1); yn = y +neigb(j,2);
        
        % Check if neighbour is inside or outside the image
        ins=(xn>=1)&&(yn>=1)&&(xn<=Isizes(1))&&(yn<=Isizes(2));
        
        % Add neighbor if inside and not already part of the segmented area
        if(ins&&(J(xn,yn)==0)) 
                neg_pos = neg_pos+1;
                neg_list(neg_pos,:) = [xn yn I(xn,yn)]; J(xn,yn)=1;
        end
    end

    % Add a new block of free memory
    if(neg_pos+10>neg_free), neg_free=neg_free+10000; neg_list((neg_pos+1):neg_free,:)=0; end
    
    % Add pixel with intensity nearest to the mean of the region, to the region
    dist = abs(neg_list(1:neg_pos,3)-reg_mean);
    [pixdist, index] = min(dist);
    J(x,y)=2; reg_size=reg_size+1;
    
    % Calculate the new mean of the region
    reg_mean= (reg_mean*reg_size + neg_list(index,3))/(reg_size+1);
    
    % Save the x and y coordinates of the pixel (for the neighbour add proccess)
    x = neg_list(index,1); y = neg_list(index,2);
    
    % Remove the pixel from the neighbour (check) list
    neg_list(index,:)=neg_list(neg_pos,:); neg_pos=neg_pos-1;
end

% Return the segmented area as logical matrix
J=J>1;










