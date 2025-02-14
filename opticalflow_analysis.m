%% Jiwon Kim, Ph.D. Brown University; Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids

clear all; close all;
%%

whichdirection = ['rad'];
if whichdirection == 'tan'
load('green4purple.mat')
elseif whichdirection == 'rad'
load('blue4red.mat')
colmap = flipud(colmap);
end

%%

modefac = 10;
maxc = 10;
dataname = 'OpticalFlow_Farneback.mat';
dx = 1;        % size of the grid of the heatmap
dx2 = 25;
df = 1;         % pick 1 data per df
dfq = 15;
um = 0.65;

dt = 0.25;
vlength = 161;
savenamepre = 'OpticalFlow_Farneback_result';
maskfilepre = 'Mask';
%%
load(dataname)
time_real = time;
mintime = 2;
maxtime = vlength;
%%
Mnew = zeros(1,9);
spheroid_displ = [0,0];
displacement2 = [0,0];
COLOR = [0 0 0];
for i=mintime:maxtime
    maskname = [maskfilepre, sprintf('%04d',i-1),'.tif'];
    maskname_before = [maskfilepre, sprintf('%04d',i-2),'.tif'];

    M = imread(maskname);

    [h, w] = size(M);
    mask = imbinarize(M);
    prop = regionprops(mask, 'centroid');   % px
    centroid = cat(1,prop.Centroid);        % px

    M_before= imread(maskname_before);
    mask_before = imbinarize(M_before);
    prop_before = regionprops(mask_before, 'centroid');   % px
    centroid2 = cat(1,prop_before.Centroid);        % px
    displacement = centroid-centroid2;
    displacement2 = displacement2+displacement;
    
    %% Optical flow data downsapling
    vxnow = Vx_cell{i}*um; % every pixel
    vynow = Vy_cell{i}*um; % every pixel

    maskscale = imresize(mask,1/df);

    rn = ceil(h/df);
    cn = ceil(w/df);
    
    unowd = zeros(rn,cn);   % every df-th pixel
    vnowd = zeros(rn,cn);   % every df-th pixel
    snowd = zeros(rn,cn);

    for p = 1:rn
        unowtemp = vxnow((p-1)*df+1,:);
        vnowtemp = vynow((p-1)*df+1,:);        
        for q = 1:cn
            unowd(p,q) = unowtemp(1,(q-1)*df+1);
            vnowd(p,q) = vnowtemp(1,(q-1)*df+1);
            snowd(p,q) = norm([unowd(p,q), vnowd(p,q)]);
        end
    end                             % every df-th pixel
    unowd = unowd.*maskscale/dt;        % px/hour    
    vnowd = vnowd.*maskscale/dt;
    snowd = snowd.*maskscale/dt;

    ind = abs(snowd) > 0.000001;           
    unowd = unowd.*ind;
    vnowd = vnowd.*ind;

    [ynow2, xnow2, unow] = find(unowd);     % matrix2vector 
    [ynow3, xnow3, vnow] = find(vnowd);     
    

    xnow = df*(xnow2-1)+1;                  % back to original scale
    ynow = df*(ynow2-1)+1;                  

    liven = length(unow);
    vrnow = zeros(liven,1);
    vrnowvec = zeros(liven,2);
    vtnow = zeros(liven,1);
    vtnowvec = zeros(liven,2);
    snow = zeros(liven,1);
    rnow = zeros(liven,1);

    for j = 1:liven
       unow(j) = unow(j)-displacement(1);
        vnow(j) = vnow(j)-displacement(2);
    end

    unowmode = mode(round(unow*modefac))/modefac;
    vnowmode = mode(round(vnow*modefac))/modefac;
    unow = unow - unowmode;
    vnow = vnow - vnowmode;     % Calibration

    for k = 1:liven
        locx = xnow(k) - centroid(1);
        locy = ynow(k) - centroid(2);
        nloc = norm([locx, locy]);

        unit = [locx locy]/nloc;
        unit_orth = [-locy/nloc locx/nloc];
        vrnow(k) = unit*[unow(k); vnow(k)];
        vrnowvec(k,:) = vrnow(k)*unit;
        vtnow(k) = unit_orth*[unow(k); vnow(k)];
        vtnowvec(k,:) = vtnow(k)*unit_orth;
        snow(k) = norm([unow(k), vnow(k)]);
        rnow(k) = nloc;

    end


     %% Plot the boundary mask and velocity profile heatmap
b = bwboundaries(mask);
bm = b{1};
bn = length(bm);
plot3(bm(:,2)-displacement2(1),bm(:,1)-displacement2(2), ones(bn,1),'Color',[0 0 0], 'LineWidth',1);
hold on


if dx == 1
    [xq,yq] = meshgrid(1:w, 1:h);
    vs = NaN(h,w);

    if whichdirection == 'tan'
    vs(sub2ind([h,w], ynow, xnow)) = vtnow;
    elseif whichdirection == 'rad'
    vs(sub2ind([h,w], ynow, xnow)) = vrnow;
    end
    surfnow = surf(xq-displacement2(1), yq-displacement2(2), vs-100, 'EdgeColor','none');


else
    [xq, yq] = meshgrid(0:dx:w, 0:dx:h);

    if whichdirection == 'tan'
        vs = griddata(xnow, ynow, vtnow, xq, yq, 'linear');     % vrnow or vtnow?
    elseif whichdirection == 'rad'
        vs = griddata(xnow, ynow, vrnow, xq, yq, 'linear'); 
    end
    surfnow = surf(xq-displacement2(1), yq-displacement2(2), vs-100, 'EdgeColor','none');
end
% 
view(2)
clim([-maxc-100 maxc-100]);
colormap(colmap)

axis equal
axis off
box on

set(gcf,'color', [1 1 1]);

hold on
    %% Plot velocity profile quiver


    [xq2, yq2] = meshgrid(0:dx2:w, 0:dx2:h);
 
    uq = griddata(xnow, ynow, unow./snow, xq2, yq2, 'natural');     % vrnow or vtnow?
    vq = griddata(xnow, ynow, vnow./snow, xq2, yq2, 'natural'); 

    vrq = griddata(xnow, ynow, vrnow, xq2, yq2, 'natural'); 
    vtq = griddata(xnow, ynow, vtnow, xq2, yq2, 'natural'); 

    xq3 = xq2(:);
    yq3 = yq2(:);
    uq3 = uq(:);
    vq3 = vq(:);
    vrq3 = vrq(:);
    vtq3 = vtq(:);
    
    ind2 = inpolygon(xq3, yq3, bm(:,2), bm(:,1));
    
    xq3 = xq3(ind2);
    yq3 = yq3(ind2);
    uq3 = uq3(ind2);
    vq3 = vq3(ind2);
    vrq3 = vrq3(ind2);
    vtq3 = vtq3(ind2);
    
    if whichdirection == 'tan'
        vqq = vtq3;
    else
        vqq = vrq3;
    end

medvqq = median(abs(vqq));

qpi = quiver(xq3-displacement2(1), yq3-displacement2(2), uq3, vq3, 0.4,'LineWidth',1,'Color','none','MaxHeadSize',0.8,'Alignment','center', 'AutoScale','off');
ind3 = find(abs(vqq)>medvqq);
xq4 = xq3(ind3);
yq4 = yq3(ind3);
uq4 = uq3(ind3);
vq4 = vq3(ind3);

qpi2 = quiver(xq4-displacement2(1), yq4-displacement2(2), uq4, vq4, 0.4,'LineWidth',1,'Color','k','MaxHeadSize',0.8,'Alignment','center', 'AutoScale','off');

X = qpi.XData;
Y = qpi.YData;
U = qpi.UData;
V = qpi.VData;

headLength = 6.5;
headWidth = 4.5;

for ii = 1:length(xq3)

    if abs(vqq(ii)) >= medvqq
        colvqq = [0, 0, 0];
        colqui = [0, 0, 0];
    
    else
        colvqq = [0.7, 0.7, 0.7];
        colqui = 'none';
    end
    
    ah2 = annotation('arrow', 'headStyle', 'cback1', 'HeadLength', headLength, 'HeadWidth', headWidth, 'Color',colvqq);
    set(ah2,'parent',gca);
    set(ah2, 'position', [X(ii)+15*U(ii) Y(ii)+15*V(ii) U(ii) V(ii)])

end

fill3([0.3 0.8 0.8 0.3]*w, [0.1 0.1 0.9 0.9]*h, [-500 -500 -500 -500], [0.9 0.9 0.88]);     % Background color
xlim([0.2*w 0.9*w])
ylim([0.1*h 0.9*h])

axis equal
view(2)
pause(0.1)

hold off
hAnnotations = findall(gcf, 'Type', 'annotation');

savefig = [savenamepre,'_',whichdirection,'_',sprintf('%04d',i-1),'.tif'];
exportgraphics(gcf,savefig,'Resolution',300);


i   % Displaying slide #


end
