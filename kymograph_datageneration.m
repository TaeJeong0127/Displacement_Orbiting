%% Jiwon Kim, Ph.D. Brown University; Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids


clear; close all; clc;
%%
timestep = 161;
thickness = 20; % thickness of the ROI (um) 
px = 1/0.65;

bar = cell(timestep,1);
post = cell(timestep,1);
savename = 'Kymograph_data';
maskpre = 'Mask';
filepre = 'Nucleus';


for count = 1:timestep
maskname = [maskpre, sprintf('%04d',count-1),'.tif'];
filename = [filepre, sprintf('%04d',count-1),'.tif'];

%% Image & info loading
M = imread(maskname);
I = imread(filename);

unitx = 60;
unity = 6;

[h, w] = size(M);
xc = floor(w/2);
yc = floor(h/2);

mask = imbinarize(M);
s = regionprops(mask, 'centroid');
centroid = cat(1,s.Centroid);

dx = xc - centroid(1);
dy = yc - centroid(2);

Ic = imtranslate(I, [dx dy]);
Mc = imtranslate(M, [dx dy]);
mask = imbinarize(Mc);

dangle = 2;     % Angle interval
dtheta = dangle/180*pi();
npieces = 360/dangle;
r = 500;        % Longer than the spheroid radius (px)

thicknesspx = floor(thickness*px);
se1 = strel('disk', thicknesspx);

%% Outermost layer
mask_in = imerode(mask, se1);
mask_layer = mask - mask_in;

%% Second outermost layer
% mask_in = imerode(mask, se1);
% mask_in_in = imerode(mask_in, se1);
% mask_layer = mask_in-mask_in_in;

%% 
mask_layer = uint8(mask_layer);

B = labeloverlay(Ic, mask_layer);
B = Ic.*mask_layer;

%% ROI check
% hold off
% imshow(B)
% pause(0.1)
% hold on

%% Slicing ROI

hp = zeros(npieces,1);
wp = zeros(npieces,1);
patch = cell(npieces, 1);
for i = 1:npieces
    theta1 = dtheta*(i-1);
    theta2 = theta1 + dtheta;

    x1 = xc + r*cos(theta1);
    y1 = yc + r*sin(theta1);
    x2 = xc + r*cos(theta2);
    y2 = yc + r*sin(theta2);

    x = [xc x1 x2 xc];
    y = [yc y1 y2 yc];

    sector = poly2mask(x,y,h,w);
    sector = uint8(sector);
    window = mask_layer.*sector;
    roi = Ic.*window;
    XX = labeloverlay(Ic, window);

    cal = (theta1+0.5*dtheta)*180/pi();
    roi_rot = imrotate(roi, cal);

    [row, col, v] = find(roi_rot);
    yy1 = min(row);
    yy2 = max(row);
    xx1 = min(col);
    xx2 = max(col);

    patch{i} = roi_rot(yy1:yy2, xx1:xx2);
    [hp(i), wp(i)] = size(patch{i});
    
end

%% Stitching ROI
patch_Cal = cell(npieces,1);
slope_Cal = cell(npieces,1);
for i = 1:npieces
    patch_temp = patch{i};
    patch_temp(patch_temp==0) = NaN;
    meancol = mean(patch_temp, 'omitNaN');
    meanrow = mean(patch_temp, 2, 'omitNaN');
    slope_temp = max(meanrow);
    slope_Cal{i} = slope_temp;
   
    a=1;
    patch_temp=zeros(unity, length(meancol));
    while a<unity+1
        patch_temp(a,:) = meancol;
        a = a+1;
    end
        
    xratio = unitx / wp(i);
    A = [xratio 0 0; 0 1 0; 0 0 1];
    tform = affinetform2d(A);
    patch_temp = imwarp(patch_temp,tform);
    [hpp, wpp] = size(patch_temp);
    if hpp > unity
        patch_temp = patch_temp(2:end, :);
    end
    if wpp > unitx
        patch_temp = patch_temp(:, 2:end);
    end

    patch_Cal{i} = patch_temp;

end
bar{count} = cat(1, patch_Cal{:});
post{count} = cat(1, slope_Cal{:});

end



%% Saving data
slope = cat(2, post{:});
slope = uint8(slope);

kymograph = cat(2,bar{:});
kymograph_output = uint8(kymograph);
kymograph_output = imadjust(kymograph_output);

savekymodata = [savename, '.mat'];
save(savekymodata, 'kymograph')

saveslopedata = [savename, '_slope.mat'];
save(saveslopedata, 'slope')


