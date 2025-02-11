clear; close all; clc;
%% This code implements the methods described in:
% Jiwon Kim, Hyuntae Jeong, Carles Falc´o, Alex M. Hruska, W. Duncan 
% Martinson, Alejandro Marzoratti, Mauricio Araiza, Haiqian Yang,
% Christian Franck, Jos´e A. Carrillo, Ming Guo, and Ian Y. Wong,
% "Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids"
% The original code is written by Hyuntae Jeong, Brown Univ., Ian Wong Lab,
% 2024 and Mauricio Araiza, Christian Frank's Lab
% The quiver plots for the displacement fields are generated with this code.

%% Name of index for samples
% Sample_name = {'G11'; 'G13'; 'G14'; 'G15'; 'G16'; 'M1'; 'M2'; 'M3'; ...
%    'M4'; 'M5'; 'M6'; 'M7'; 'M8'; 'M9'; 'M10'; 'M11'; 'M12'};
% Sample_name = {'C4A'; 'C4B'; 'C5A'; 'C5B'; 'C5C'; 'C7C'; 'D4E'; 'D6B'; ...
%     'D6C'; 'D7B'; 'D8A'; 'E4A2'; 'E6A2'; 'E7C2'};
% Sample_name = {'N2'; 'N6'; 'N9'; 'N11'; 'N15'; 'I1'; 'I2'; 'I3'; ...
%     'I4'; 'I5'; 'I6';'I7';'I8';'I9';'I10';'I12';'I13'};

%% color map generation
% color map for the bead displacement(for the temporal)
colors = [161,215,106;
    20,20,20;
    233,163,201];

% color map for the bead deformation(for cumulative)
% colors = [255,255,255;
%     92,158,196;
%     18,49,97];

colors = colors/256;
n = 200;
x = linspace(1, size(colors, 1), n);
gradient_colormap = interp1(1:size(colors, 1), colors, x, 'linear');

%% Other type of colormaps
Sample_num = length(Sample_name);
data1 = load('blue2red.mat');
Cmap_R = data1.colmap;
data2 = load('green2purple.mat');
Cmap_A = data2.colmap;

% targeting interval for plotting
o = 8;

new_size = 256;
new_indices = linspace(1, 200, new_size);
Cmap_R_new = interp1(1:200, Cmap_R, new_indices);



%% for colormapping the spheroid max projection(here cell images)
tt = zeros(128, 3);
ttt = [139, 115, 85];
for ii = 1:128
    tt(ii,:) = [0, ii/128, 0];
    tt(ii,:) = ttt/255*ii/128/3 + ttt/255/3*2;
end
tt(1,:) = [0, 0, 0];
for ii = 2:40
    tt(ii,:) = tt(41,:)*ii/41;
end

%% main loop for the plotting
for i = 1:Sample_num
    savenameheader = 'Leader_test';
    foldername = ['Displacement_1127_withoutcellmigration',char(Sample_name(i))];
    mkdir(foldername)
    fol2 = [foldername,'/'];
    file_name = ['Results_compliation_',char(Sample_name(i)),'_new.mat'];
    load(file_name);

    time_length = size(Ucell_optic,3);
    spacing = 8; % spacing for vector maps

    for t = 1:time_length/o-1

        ti = o*(t-1)+2;

        Ucell_traj = X_pos(~isnan(Ucell_optic(:,:,ti))); Vcell_traj = Y_pos(~isnan(Vcell_optic(:,:,ti)));
        Ubead_traj = X_pos(~isnan(Udis_do(:,:,ti))); Vbead_traj = Y_pos(~isnan(Vdis_do(:,:,ti)));
        U_initial = Ucell_traj; V_initial = Vcell_traj;
        U_bead_initial = Ubead_traj; V_bead_initial = Vbead_traj;

        for t1 = 1:o
            Ucell_add = griddata(X_pos, Y_pos, Ucell_optic(:,:,ti+(t1-1)), Ucell_traj, Vcell_traj);
            Vcell_add = griddata(X_pos, Y_pos, Vcell_optic(:,:,ti+(t1-1)), Ucell_traj, Vcell_traj);
            Ubead_add = griddata(X_pos, Y_pos, Udis_do(:,:,ti+(t1-1)), Ubead_traj, Vbead_traj);
            Vbead_add = griddata(X_pos, Y_pos, Vdis_do(:,:,ti+(t1-1)), Ubead_traj, Vbead_traj);

            Ucell_add(isnan(Ucell_add)) = 0; Vcell_add(isnan(Vcell_add)) = 0;

            Ucell_traj = Ucell_traj + Ucell_add/(o/2);
            Vcell_traj = Vcell_traj + Vcell_add/(o/2);

            Ubead_traj = Ubead_traj + Ubead_add/(o/2);
            Vbead_traj = Vbead_traj + Vbead_add/(o/2);
        end

        Ucell_temp = Ucell_traj - U_initial; Vcell_temp = Vcell_traj -V_initial;
        Ubead_temp = Ubead_traj - U_bead_initial; Vbead_temp = Vbead_traj -V_bead_initial;
        Ucell(:,:,t) = griddata(U_initial,V_initial,Ucell_temp,X_pos,Y_pos,'cubic');
        Vcell(:,:,t) = griddata(U_initial,V_initial,Vcell_temp,X_pos,Y_pos,'cubic');
        Ubead(:,:,t) = griddata(U_bead_initial,V_bead_initial,Ubead_temp,X_pos,Y_pos,'cubic');
        Vbead(:,:,t) = griddata(U_bead_initial,V_bead_initial,Vbead_temp,X_pos,Y_pos,'cubic');

        if ti<10
            num_ti = ['000',num2str(ti)];
        elseif ti<100
            num_ti = ['00',num2str(ti)];
        elseif ti<1000
            num_ti = ['0',num2str(ti)];
        else
            num_ti = num2str(ti);
        end

        %% File loading for plottings
        % if i <= 8
        %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %         ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % 
        % 
        %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %         ,char(Sample_name(i)), ' Max intensity/Ch2 bead/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % 
        %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        %         ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        % else
        %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %         ,char(Sample_name(i)), ' Max intensity/Ch2 nucleus/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % 
        % 
        %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %         ,char(Sample_name(i)), ' Max intensity/Ch3 bead/Max_',char(Sample_name(i)),'_C3_t',num_ti,'.tif'];
        % 
        % 
        %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        %         ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        % end
        % 
        % % % For PEG data
        % % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        % %     ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % % 
        % % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/2 '...
        % %     ,char(Sample_name(i)), ' Max intensity/Ch2 bead/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % % 
        % % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        % %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        % % 
        % % if i <= 5
        % %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        % %         ' nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        % %         ' bead/Max_',char(Sample_name(i)),'_C4_t',num_ti,'.tif'];
        % %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        % %         ' Mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % % else
        % %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        % %         ' ch2/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        % %         ' ch4/Max_',char(Sample_name(i)),'_C4_t',num_ti,'.tif'];
        % %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        % %         ' mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % % end

        Mask_time = imread(Mask_time_name);
        boundaries = bwboundaries(Mask_time);
        % change the domain to 1-0 scale
        Mask_time = double(Mask_time);
        Mask_time = Mask_time/max(Mask_time(:));
        stats = regionprops(Mask_time, 'Centroid');
        % Extract the centroid coordinates
        centroid = stats.Centroid*0.65; % This will be a 1x2 vector [x, y]

        ratio = 1;
        X_pos_prin = downsample(downsample(X_pos,ratio)',ratio)';
        Y_pos_prin = downsample(downsample(Y_pos,ratio)',ratio)';
        Umap_prin = downsample(downsample(Umap_do(:,:,ti),ratio)',ratio)';
        Vmap_prin = downsample(downsample(Vmap_do(:,:,ti),ratio)',ratio)';
        Ucell_prin = downsample(downsample(Ucell_optic(:,:,ti),ratio)',ratio)';
        Vcell_prin = downsample(downsample(Vcell_optic(:,:,ti),ratio)',ratio)';
        REF_x = X_pos_prin - centroid(1);
        REF_y = Y_pos_prin - centroid(2);
        Distance_cent = sqrt(REF_x.^2+REF_y.^2);
        ROI = Distance_cent<250;
        ROI_do = double(ROI);
        ROI_do(ROI_do==0) = NaN;

        angle = atan2(REF_y,REF_x);
        R_cell_prin = Ucell_prin.*cos(angle) + Vcell_prin.*sin(angle);
        A_cell_prin = Ucell_prin.*cos(angle+pi/2) + Vcell_prin.*sin(angle+pi/2);
        R_map_prin = Umap_prin.*cos(angle) + Vmap_prin.*sin(angle);
        A_map_prin = Umap_prin.*cos(angle+pi/2) + Vmap_prin.*sin(angle+pi/2);


        image = double(imread(Nuc_name));
        im = image/max(image(:));

        idx = randperm(size(X_pos(:),1),ceil(size(X_pos(:),1)/spacing)); % random vector position
        Xi = X_pos(idx);
        Yi = Y_pos(idx);
        Umap = Ucell_prin; Vmap = Vcell_prin;
        U = Ucell_prin(idx); V = Vcell_prin(idx);
        W = zeros(size(Xi)); Z = zeros(size(Xi));
        R_val = A_cell_prin(idx);

        U_t = Umap_prin(idx); V_t = Vmap_prin(idx);
        R_val_t = R_map_prin(idx);


        % Compute magnitude of displacement

        mag = sqrt(U.^2 + V.^2 + W.^2);

        % dimag_max = max(abs(A_cell_prin), [], 'all');
        dimag_max = 10;

        %dimag_max_t = max(abs(R_map_prin),[],'all');
        dimag_max_t = 1;

        mag = mag(:);
        R_val = R_val(:);
        R_val_t = R_val_t(:);
        f = figure('Visible','on');

        IM = dimag_max_t.*2*im - dimag_max_t-0.4;
        [m,n] = size(IM);
        x = 1:n;
        y = 1:m;
        [x_im,y_im] = meshgrid(x,y);
        x_im = x_im*0.65; y_im = y_im.*0.65;
        [~, hC] = contourf(x_im,y_im,IM, 128);
        set(hC, 'LineStyle', 'none');
        set(gcf,'Visible', 'on');
        colormap(tt);  
        freezeColors
        freezeColors(colorbar(gca, 'south', 'Visible','off'))
        axis image
        axis ij
        axis off

        % hc = coneplot(Xi, Yi, Z, U, V, W, 0.02,'nointerp');
        % clim([-dimag_max, dimag_max]);
        % fvc = repmat(R_val.', [42 1]);
        % set(hc, 'FaceColor', 'flat', 'FaceVertexCData', fvc(:));
        % hc.EdgeColor = 'none';
        % hc.AmbientStrength = 0.6;
        % hc.DiffuseStrength = 0.75;
        % hc.SpecularStrength = 0.4;
        % hold on
        % 
        % colormap(Cmap_A);
        % freezeColors
        % h_color = jicolorbar('vshort');

        % set(h_color, 'ylim', [-dimag_max,dimag_max]);
        % set(h_color, 'FontSize', 12);
        % axis image
        % axis ij
        % lighting phong;


        hc2 = coneplot(Xi, Yi, Z, U_t, V_t, W, 0.025,'nointerp');
        clim([-dimag_max_t, 1]);
        fvc2 = repmat(R_val_t.', [42 1]);
        set(hc2, 'FaceColor', 'flat', 'FaceVertexCData', fvc2(:));
        hc2.EdgeColor = 'none';
        hc2.AmbientStrength = 0.6;
        hc2.DiffuseStrength = 0.75;
        hc2.SpecularStrength = 0.4;
        hold on

        colormap(gradient_colormap);
        freezeColors

        h_color = jicolorbar('vshort');

        set(h_color, 'ylim', [-dimag_max_t, dimag_max_t]);
        set(h_color, 'FontSize', 12);
        axis image
        axis ij
        lighting phong;
        set(gcf,'Visible', 'off');

        for k = 1:length(boundaries)
            boundary = boundaries{k};
            plot(boundary(:,2)*0.65, boundary(:,1)*0.65, 'y-', 'LineWidth', 1.5);
        end


        
        % Figure limits for the spheroid
        xmin = centroid(1) - 250;
        xmax = centroid(1) + 250;
        ymin = centroid(2) - 250;
        ymax = centroid(2) + 250;
        xlim([xmin, xmax])
        ylim([ymin, ymax])
        ylabel(h_color, 'Radial Displacement');
        hold off

        set(gcf,'Visible', 'off');

        title({['Frame: ' num2str(ti)]},'interpreter','none', 'fontsize', 18)
        exportgraphics(f, append(fol2, 'QDIC_spacing4_displacement_Frame', num2str(ti),'.png'),Resolution=150);

        close;
        
    end
end
