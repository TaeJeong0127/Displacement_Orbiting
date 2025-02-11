clear; clc; close all
%% This code implements the methods described in:
% Jiwon Kim, Hyuntae Jeong, Carles Falc´o, Alex M. Hruska, W. Duncan 
% Martinson, Alejandro Marzoratti, Mauricio Araiza, Haiqian Yang,
% Christian Franck, Jos´e A. Carrillo, Ming Guo, and Ian Y. Wong,
% "Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids"
% The original code is written by Jiwon Kim, Brown Univ., Ian Wong Lab, 2024

tic
xratio = 4;
dangle = 0.5;
angledeg = (0:dangle:360)';
anglerad = angledeg*pi/180;

nangle = length(anglerad);
probe = 500;

n=24;
glv = 0.5;
glvc = 0.8;

redwhite = [
    linspace(glvc,1,n)', linspace(glvc,glv,n)',linspace(glvc,glv,n)';
    ones(n,1), linspace(glv,0,n)',linspace(glv,0,n)';
    linspace(1,glv,n)', zeros(n,2);];
bluewhite = rot90(redwhite,2);
redwhiteblue = [bluewhite;redwhite];


load('blue2red.mat')


meandata = zeros(1,6);
ax1 = axes;
setn = 'M';
%spheroidID = {'D7B2','E4A2','E6A2','E7C2'};
%spheroidID = {'C4A', 'C4B','C5A','C5B','C5C','C7C','D4E','D6B','D6C','D7B','D8A','E4A','E6A','E7A','E7B','E7C'};
dt = 0.5;
asp = [0.8 1 1];
tmax = 72;
nslice = tmax/dt+1;

% for count = 1:12
% for count = [1,2,3,4,5,6,9,11,13,15,16]
% for count = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]]
% for count = [2,6,9,11,15]
for count = 10
    %spheroid = spheroidID{count};
    spheroid = 'F042';
    % for G, M samples
    % maskfilepre = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
    %     ,spheroid, ' Mask/' ,spheroid, ' Mask'];
    %% Before data
    % maskfilepre = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
    % ,spheroid, ' Mask/' ,spheroid, ' Mask'];
    %% After data
    maskfilepre = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241225 peg media change/tiff2/'...
    ,spheroid, ' Mask/' ,spheroid, ' Mask'];

    % % For PEG samples
    % maskfilepre = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
    %     ,spheroid, ' Mask/' ,spheroid, ' Mask'];

    % spheroid = spheroidID{count};
    %
    % pinkshadename = [spheroid, 'vimentin_vs_angle.mat'];
    % pinkshade = load(pinkshadename);
    % pinkshade = pinkshade.Msave;
    pinkshade = zeros(36,193);
    dartfilename = [spheroid, 'DART_0109.txt'];
    dart = readtable(dartfilename);
    darttime = dart.slide;
    seq = unique(darttime);
    dartcont = dart.contraction;
    dartprot = dart.protrusion;


    dartang = (0:10:350)';
    dartangrad = deg2rad(dartang);
    dartang2 = (10:10:360)';
    dartangrad2 = deg2rad(dartang2);
    dartangnowx = cos(dartang);
    dartangnowy = -sin(dartang);

    jetn = jet(nslice);
    um = 0.65;              % 0.65um/px
    px = 1/0.65;            % 1px/0.65um

    meandistangle = zeros(nslice,1);
    maxdistangle = zeros(nslice,1);
    mindistangle = zeros(nslice,1);
    stddistangle = zeros(nslice,1);

    t = (0:dt:(nslice-1)*dt)';

    rtdata = [0 0 0 0 0 0 0 0 0];

    rel_rt = zeros(nangle,nslice);

    darts = [];
    dartsp = [];

    dartangs = [];

    for i=5:nslice
        % for i=4:nslice
        now = max(seq(seq <= i));

        dartangnow = dartang+360*(i-1);
        % dartnow = flipud(dartcont(darttime == now));
        dartnow = dartcont(darttime == now);
        darts = [darts;dartnow];

        dartnowp = dartprot(darttime == now);
        dartsp = [dartsp;dartnowp];

        dartangs = [dartangs;dartangnow];


        maskfilename = [maskfilepre, sprintf('%04d', i-1),'.tif'];
        % maskfilename = [maskfilepre, sprintf('%04d', i),'.tif'];
        M = imread(maskfilename);

        Mum = imresize(M, um);
        mask = imbinarize(Mum);       % px
        % mask = imcomplement(mask);
        [h, w] = size(mask);
        prop = regionprops(mask, 'centroid');
        centroid = cat(1,prop.Centroid);
        xc = centroid(1);
        yc = centroid(2);

        boundary = bwboundaries(mask);      % initial point = final point, loop
        boundaryvec = boundary{1};
        nboundary = length(boundaryvec);

        pboundary = polyshape(boundaryvec(:,2), boundaryvec(:,1));

        rt = zeros(nangle, 1);
        xx = zeros(nangle, 1);
        yx = zeros(nangle, 1);



        for j = 1:nangle
            lin = [xc yc; xc+probe*cos(anglerad(j)) yc-probe*sin(anglerad(j))];
            [in,out] = intersect(pboundary, lin);
            xx(j) = out(end-1,1);
            yx(j) = out(end-1,2);

            %
            % plot(pboundary)
            % hold on
            % plot([xc xc+probe*cos(anglerad(j))], [yc yc-probe*sin(anglerad(j))],'b-');
            % plot(xx(j),yx(j),'bx')

            rt(j) = sqrt((xx(j)-xc)^2 + (yx(j)-yc)^2);

            % pause(0.01)

        end


        datatemp = [ones(nangle,1)*[count, i, xc, yc], anglerad, angledeg, xx, yx, rt];
        rtdata = cat(1,rtdata, datatemp);
        meandistangle(i) = mean(rt);
        maxdistangle(i) = max(rt);
        mindistangle(i) = min(rt);
        stddistangle(i) = std(rt);

        rel_rttemp = rt-meandistangle(i)*ones(nangle,1);
        rel_rt(:,i) = rel_rttemp;

        xxx = xx-(xc-floor(w/2));
        yxx = yx-(yc-floor(h/2));

        % xnow3 = xnow-(xc-floor(w/2));
        % ynow3 = ynow-(yc-floor(h/2));
        % xnow4 = xnow2-(xc-floor(w/2));
        % ynow4 = ynow2-(yc-floor(h/2));

        % figure;

        surf([xxx(:) xxx(:)], [yxx(:) yxx(:)], [rel_rttemp(:) rel_rttemp(:)], 'FaceColor','none','EdgeColor','interp', 'LineWidth',3)
        colormap(redwhiteblue)
        % colormap(colmap)
        caxis([-25 25]);
        % set(gca, 'color', [0.8 0.8 0.8]);
        grid off
        axis equal
        pbaspect([1 1 50])


        view(2)
        % plot(angledeg, rt-meandistangle, 'Color', colortemp)

        hold on
        rectangle('Position', [1 1 w h], 'LineWidth',0.01);
        xlim([-10 w+20])
        ylim([-10 h+20])
        axis off
        % viscircles([floor(w/2),floor(h/2)],meandistangle(i), Color=[glvc glvc glvc], Linewidth= 1, EnhanceVisibility=0, LineStyle='--');

        % quiver(xnow3, ynow3, unow, vnow, 2, 'Color',[glv glv glv])
        % quiver(xnow4, ynow4, unow2, vnow2, 2, 'Color',[glv glv glv])
        % pause(0.1)
        % hold off

        xcos = cos(dartangrad);
        ysin = -sin(dartangrad);
        xcos2 = cos(dartangrad2);
        ysin2 = -sin(dartangrad2);
        % plus = meandistangle(i);
        plus = rt(6:20:720)+3;
        % plus = rt(1:20:720);
        scalar = 5;
        % scalar = 100;
        dartnow = flipud([dartnow(19:36);dartnow(1:18)]);
        dartxnow = (dartnow*scalar+plus) .* xcos;
        dartynow = (dartnow*scalar+plus) .* ysin;
        dartxnow2 = (dartnow*scalar+plus) .* xcos2;
        dartynow2 = (dartnow*scalar+plus) .* ysin2;

        dartxnow3 = probe.*xcos2;
        dartynow3 = probe.*ysin2;
        dartxnow4 = probe.*xcos;
        dartynow4 = probe.*ysin;


        dartnowp = flipud([dartnowp(19:36);dartnowp(1:18)]);
        dartxnowp = (dartnowp*scalar+plus) .* xcos;
        dartynowp = (dartnowp*scalar+plus) .* ysin;
        dartxnowp2 = (dartnowp*scalar+plus) .* xcos2;
        dartynowp2 = (dartnowp*scalar+plus) .* ysin2;



        npie = length(dartxnow);

        pinkarrayn = [pinkshade(:,nslice) zeros(36,1) pinkshade(:,nslice)];


        x3 = [dartxnow' + floor(w / 2);
            dartxnow2' + floor(w / 2);
            % dartxnow3' + floor(w / 2);
            % dartxnow4' + floor(w / 2)];
            floor(w / 2)*ones(1,36)];


        y3 = [dartynow' + floor(h / 2);
            dartynow2' + floor(h / 2);
            % dartynow3' + floor(h / 2);
            % dartynow4' + floor(h / 2)];
            floor(h / 2)*ones(1,36)];

        x4 = [dartxnowp' + floor(w / 2);
            dartxnowp2' + floor(w / 2);
            floor(w / 2)*ones(1,36)];

        y4 = [dartynowp' + floor(h / 2);
            dartynowp2' + floor(h / 2);
            floor(h / 2)*ones(1,36)];


        sector_pre = poly2mask(xxx,yxx,h,w);
        se = strel('sphere',4);
        sector = imerode(sector_pre,se);
        sectorp = imdilate(sector_pre,se);


        for ai = 1:36
            xi = x3(:, ai);
            yi = y3(:, ai);
            ci = pinkarrayn(ai,:);
            % patch(xi, yi, ci,'EdgeColor','none', 'FaceAlpha',0.3);

            bite = poly2mask(xi,yi,h,w);
            crust = bite & ~sectorp;
            boundary_c = bwboundaries(crust);      % initial point = final point, loop
            if ~isempty(boundary_c)
                boundaryvec_c = boundary_c{1};

                patch(boundaryvec_c(:,2), boundaryvec_c(:,1), 'b','EdgeColor','none', 'FaceAlpha',0.4);
                hold on
            end

            xj = x4(:, ai);
            yj = y4(:, ai);
            ci = pinkarrayn(ai,:);
            % patch(xi, yi, ci,'EdgeColor','none', 'FaceAlpha',0.3);

            piece = poly2mask(xj,yj,h,w);
            crust = piece & ~sectorp;
            boundary_p = bwboundaries(crust);      % initial point = final point, loop
            if ~isempty(boundary_p)
                boundaryvec_p = boundary_p{1};

                patch(boundaryvec_p(:,2), boundaryvec_p(:,1), 'r','EdgeColor','none', 'FaceAlpha',0.4);
                hold on
            end

            % xj = x4(:, ai);
            % yj = y4(:, ai);
            %
            % patch(xj, yj, 'r','EdgeColor','none', 'FaceAlpha',0.3);

        end

        % plot([dartxnow;dartxnow(1)]*scalar+floor(w/2), -[dartynow;dartynow(1)]*scalar+floor(h/2), 'Color',[0.3 0.3 0.3], 'LineWidth',1)
        % fill([dartxnow;dartxnow(1)]*scalar+floor(w/2), -[dartynow;dartynow(1)]*scalar+floor(h/2), [glv glv glv], 'LineStyle','none')

        % titlestring = [spheroid, ' DART at Slide # = ', num2str(now)];
        % title(titlestring);

        pause(0.1)

        savepic = [spheroid,'_overlap_vim_out_s5_p3_',num2str(i,'%.0f'),'.tif'];
        exportgraphics(gcf,savepic,'Resolution',300)
        hold off

    end

    % rtdata = rtdata(2:end,:);
    % %%
    % datatable = array2table(rtdata,'VariableNames',{'SpheroidID','time','xc','yc','angle_rad','angle_deg','inter_x','inter_y','r(t)'});

    % savename = [spheroid, '_rtdata.xlsx'];
    % writetable(datatable, savename);

end








