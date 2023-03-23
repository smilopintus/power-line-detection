
%%% Reference
%%% 1.  Song, B., & Li, X. (2014). Power line detection from optical
%%%     images. Neurocomputing, 129, 350-361.
%%%
%%% 2.  Oyibo, P. O., Ajayi, O., & Abubakar, I. A. M. O. A. (2019). 
%%%     DEVELOPMENT OF POWER LINE DETECTION ALGORITHM USING FRANGI 
%%%     FILTER AND FIRST-ORDER DERIVATIVE OF GAUSSIAN. DEVELOPMENT, 8(1).

close all;clear;clc;
addpath(genpath(pwd))


%% input image
imageFile = 'img4.png';
thresholding = {'mf', 'ff'};
method = 1;  %1: matched filterthresholding , 2: frangi filter thresholding

%%
I1=imread(imageFile);

% use first order derivative of gaussian filter on the image
fim=mat2gray(I1);
[imx,imy]=gaussgradient(fim,1);
G=rgb2gray(imy);

% reference threshold based on the image response to FDOG
c = 2;
uG = mean(mean(G));
TG = c*uG;

if strcmp(thresholding(method), 'mf')
    % matched filter and first order derivative of gaussian power line detection
    % thresholding
    %get the matched filter image response
    [imx,imy]=matchedfilter(fim,1);
    M=rgb2gray(imadd(imx, imy));
    
    % adjust threshold based on the response to matched filter
    r= 10;
    R=(1/r).*ones(r);
    M_bar=imfilter(M,R);
    IM = M_bar;
    IM = IM - min(IM(:));
    IM = IM / max(IM(:));
    Mn_bar = IM;
    T = (1 + Mn_bar).*TG;
    
elseif strcmp(thresholding(method), 'ff')
    
    %use frangi vesselness filter on the image
    I2=rgb2gray(I1);
    Ivessel=FrangiFilter2D(double(I2));
    
    % adjust threshold based on the response to frangi filter
    IM = Ivessel;
    IM = IM - min(IM(:));
    IM = IM / max(IM(:));
    Ivn = IM;
    T = (1 + Ivn).*TG;
end

%apply threshold on the response to FDOG response
Q=G>T;

%apply morphological filtering to binary image
Q=bwmorph(Q, 'skel', inf);
bp=find(bwmorph(Q, 'branchpoints'));
for i=1:numel(bp)
    [i,j] = ind2sub(size(Q), bp(i));
    if any(~[i-1:i+1,j-1:j+1])==0
        Q(i-1:i+1,j-1:j+1)=0;
    end
end

% remove edge pixels with less than 30 connectedness
Q=bwareaopen(Q,30);

% filtering by thresholding the edge smoothness measurement.
CC = bwconncomp(Q);
for idx=1:CC.NumObjects
    ind = CC.PixelIdxList{idx};
    [x,y]=ind2sub(CC.ImageSize, ind);
    p = polyfit(x,y,1);
    x0 = min(x):max(x);
    y0 = polyval(p,x0, 'r');
    curvexy = [x0', y0'];
    mapxy = [x,y];
    [xy,distance,t] = distance2curve(curvexy,mapxy,'linear');
    S=var(distance)*10;
    if (S > 110)
        Q(CC.PixelIdxList{idx})=0;
    end
end

%%
W = clusterEdges(Q);

%% display detected result
dI = I1;
rInd = find(W==1);
rChan = dI(:,:,1);
rChan(rInd) = 256;
dI(:,:,1) = rChan;
figure;
h1 = subplot(1,2,1);
imshow(I1); title('original image')

h1 = subplot(1,2,2);
imshow(dI); title('detected result')

%imwrite(dI, strcat('result/', imageFile));
