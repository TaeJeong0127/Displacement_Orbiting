clear; close all; clc;

%% This code implements the methods described in:
% Jiwon Kim, Hyuntae Jeong, Carles Falc´o, Alex M. Hruska, W. Duncan 
% Martinson, Alejandro Marzoratti, Mauricio Araiza, Haiqian Yang,
% Christian Franck, Jos´e A. Carrillo, Ming Guo, and Ian Y. Wong,
% "Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids"
% The original code is written by Hyuntae Jeong, Brown Univ., Ian Wong Lab,
% 2024 and Mauricio Araiza, Christian Frank's Lab
% The quiver plots for the displacement fields are generated with this code.

%% sample names
% Sample_name = {'G11'; 'G13'; 'G14'; 'G15'; 'G16'; 'M1'; 'M2'; 'M3'; ...
%     'M4'; 'M5'; 'M6'; 'M7'; 'M8'; 'M9'; 'M10'; 'M11'; 'M12'};
% Sample_name = {'C4A'; 'C4B'; 'C5A'; 'C5B'; 'C5C'; 'C7C'; 'D4E'; 'D6B'; ...
%    'D6C'; 'D7B'; 'D8A'; 'E4A'; 'E7C'; 'E4A2'; 'E6A2'; 'E7C2'};
% Sample_name = {'N2'; 'N6'; 'N9'; 'N11'; 'N15'; 'I1'; 'I2'; 'I3'; ...
%    'I4'; 'I5'; 'I6';'I7';'I8';'I9';'I10';'I12';'I13'};
% For PEG before sample
%Sample_name = {'E001';'E002';'E041';'E042';'E043';'E044';'E045';'E046';'E401';'E402';'E403';'E404';'E405';'E406';'E407';'E408';'E441';'E442';'E443'};
%% For PEG after sample
% Sample_name = {'F001';'F002';'F041';'F042';'F043';'F044';'F045';'F046';'F401';'F402';'F403';'F404';'F405';'F406';'F407';'F408';'F441';'F442';'F443'};


Sample_num = length(Sample_name);
spacing = 5; % spacing for vector maps

for i = 11:17
    % savenameheader = 'Leader_test';
    % foldername = ['DART_plot_franklab_deformation_1015_realscale',char(Sample_name(i))];
    % mkdir(foldername)
    % fol2 = [foldername,'/'];
    file_name = ['Results_compliation_',char(Sample_name(i)),'_new.mat'];
    load(file_name);

    time_length = size(Ucell_optic,3);
    o = 4;
    m = 4; % for frame/hour
    for t = 1+m:time_length-(m+1)
        ti = t;

        Ucell_traj = X_pos(~isnan(Ucell_optic(:,:,ti))); Vcell_traj = Y_pos(~isnan(Vcell_optic(:,:,ti)));
        Ubead_traj = X_pos(~isnan(Udis_do(:,:,ti))); Vbead_traj = Y_pos(~isnan(Vdis_do(:,:,ti)));
        U_initial = Ucell_traj; V_initial = Vcell_traj;
        U_bead_initial = Ubead_traj; V_bead_initial = Vbead_traj;

        for t1 = 1:o
            Ucell_add = griddata(X_pos, Y_pos, Ucell_optic(:,:,ti+(t1-m)), Ucell_traj, Vcell_traj);
            Vcell_add = griddata(X_pos, Y_pos, Vcell_optic(:,:,ti+(t1-m)), Ucell_traj, Vcell_traj);
            Ubead_add = griddata(X_pos, Y_pos, Udis_do(:,:,ti+(t1-m)), Ubead_traj, Vbead_traj);
            Vbead_add = griddata(X_pos, Y_pos, Vdis_do(:,:,ti+(t1-m)), Ubead_traj, Vbead_traj);

            Ucell_add(isnan(Ucell_add)) = 0; Vcell_add(isnan(Vcell_add)) = 0;

            Ucell_traj = Ucell_traj + Ucell_add/m;
            Vcell_traj = Vcell_traj + Vcell_add/m;

            Ubead_traj = Ubead_traj + Ubead_add/m;
            Vbead_traj = Vbead_traj + Vbead_add/m;
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

        %% For loading the mask, image, bead files
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


        %% For PEG data
        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/2 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch2 bead/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];

        %% for GM data
        % if i <= 5
        %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        %         ' nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        %         ' bead/Max_',char(Sample_name(i)),'_C4_t',num_ti,'.tif'];
        %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        %         ' Mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % else
        %     Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        %         ' ch2/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        %     bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        %         ' ch4/Max_',char(Sample_name(i)),'_C4_t',num_ti,'.tif'];
        %     Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        %         ' mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % end
        
        %% For PEG before data (1223)

        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        %     ,char(Sample_name(i)), ' nucleus/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        %     ,char(Sample_name(i)), ' bead/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        % ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        % 
        % %% For PEG after data (1225)
        % 
        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241225 peg media change/'...
        %     ,char(Sample_name(i)), ' nucleus/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241225 peg media change/'...
        %     ,char(Sample_name(i)), ' bead/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241225 peg media change/tiff2/'...
        % ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];


        Mask_time = imread(Mask_time_name);
        boundaries = bwboundaries(Mask_time);
        % change the domain to 1-0 scale
        Mask_time = double(Mask_time);
        Mask_time = Mask_time/max(Mask_time(:));
        stats = regionprops(Mask_time, 'Centroid');
        % Extract the centroid coordinates
        centroid = stats.Centroid*0.65; % This will be a 1x2 vector [x, y]
        % Compute magnitude of displacement
        xc = centroid(1);
        yc = centroid(2);
        [thetas,~] = cart2pol(X_pos-xc, Y_pos-yc);
        edges = -pi:pi/36:pi;
        IDX = discretize(thetas, edges);
        DARTvalues = zeros(12, 4);
        ratio = 1;
        X_pos_prin = downsample(downsample(X_pos,ratio)',ratio)';
        Y_pos_prin = downsample(downsample(Y_pos,ratio)',ratio)';
        Umap_prin = downsample(downsample(Umap_do(:,:,ti),ratio)',ratio)';
        Vmap_prin = downsample(downsample(Vmap_do(:,:,ti),ratio)',ratio)';
        Ubead_prin = downsample(downsample(Ubead(:,:,t),ratio)',ratio)';
        Vbead_prin = downsample(downsample(Vbead(:,:,t),ratio)',ratio)';
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

        R_bead_prin = Ubead_prin.*cos(angle) + Vbead_prin.*sin(angle);
        A_bead_prin = Ubead_prin.*cos(angle+pi/2) + Vbead_prin.*sin(angle+pi/2);
        
        Mag_U = sqrt(Umap_prin.^2 + Vmap_prin.^2);
        percentile_95(i,t) = prctile(Mag_U(:), 95, 'all', 'Method','approximate');
        percent25(i,t) = prctile(Mag_U(:), 25, 'all', 'Method','approximate');
        percent75(i,t) = prctile(Mag_U(:), 75, 'all', 'Method','approximate');


        Mag_U_bead = sqrt(Ubead_prin.^2 + Vbead_prin.^2);
        percentile_95_bead(i,t) = prctile(Mag_U_bead(:), 95, 'all', 'Method','approximate');
        percent25_bead(i,t) = prctile(Mag_U_bead(:), 25, 'all', 'Method','approximate');
        percent75_bead(i,t) = prctile(Mag_U_bead(:), 75, 'all', 'Method','approximate');


        for ii = 1:72
            tempBin = R_map_prin(IDX == ii);
            R_bead_target = R_bead_prin(IDX == ii);
            A_bead_target = A_bead_prin(IDX == ii);
            R_bead_mean(ii) = mean(R_bead_target,'all','omitnan');
            A_bead_mean(ii) = mean(abs(A_bead_target),'all','omitnan');

            % Classify vectors as protrusive or contractile
            ur_prot = tempBin(tempBin > 0);
            ur_cont = tempBin(tempBin < 0);

            % Classify vectors as counterclockwise or clockwise
            tempBin2 = A_map_prin(IDX == ii);
            ut_clock = tempBin2(tempBin2 < 0);
            ut_count = tempBin2(tempBin2 > 0);

            % Take the mean of each variable
            ur_prot_mean = mean(ur_prot, 'all', 'omitnan');
            ur_cont_mean = mean(ur_cont, 'all', 'omitnan');
            ur_cont_mean = abs(ur_cont_mean);

            ut_clock_mean = mean(ut_clock, 'all', 'omitnan');
            ut_clock_mean = abs(ut_clock_mean);
            ut_count_mean = mean(ut_count, 'all', 'omitnan');

            DARTvalues(ii, 1) = ur_prot_mean;
            DARTvalues(ii, 2) = ur_cont_mean;

            DARTvalues(ii, 3) = ut_count_mean;
            DARTvalues(ii, 4) = ut_clock_mean;
        end

        DARTvalues(isnan(DARTvalues)) = 0;
        R_bead_total(:,t) = R_bead_mean;
        data_pro = DARTvalues(:,2)';
        Angle_data = edges;

        fs = 1;

        smooth_data = smoothdata(data_pro,"gaussian",20);
        Y = fft(smooth_data);
        L = length(smooth_data);


        gradient_signal = diff(smooth_data);
        radient_signal = [gradient_signal, gradient_signal(end)];
        edges = edges(2:end);

        extended_signal = [smooth_data, smooth_data, smooth_data];
        R_extended = [R_bead_mean, R_bead_mean, R_bead_mean];
        A_extended = [A_bead_mean, A_bead_mean, A_bead_mean];

        [ext_pks, ext_locs] = findpeaks(extended_signal);
        valid_peak2 = ext_locs(ext_locs > length(smooth_data) & ext_locs <= 2*length(smooth_data)) - length(smooth_data);
        valid_peak = ext_locs(ext_locs > length(smooth_data) & ext_locs <= 2*length(smooth_data));
        peak_values = smooth_data(valid_peak2);

        extended_signal = [smooth_data, smooth_data, smooth_data];
        [ext_pks, ext_locs] = findpeaks(-extended_signal);
        valid_troughs2 = ext_locs(ext_locs > length(smooth_data) & ext_locs <= 2*length(smooth_data)) - length(smooth_data);
        valid_troughs = ext_locs(ext_locs > length(smooth_data) & ext_locs <= 2*length(smooth_data));
        trough_values = smooth_data(valid_troughs2);

         offsets = [-1, 0 ,1];
         Valid_peak = unique(valid_peak + offsets');
         Valid_through = unique(valid_troughs + offsets');

         Radial_temp_peak(i,t) = mean(R_extended(Valid_peak));
         Angular_temp_peak(i,t) = mean(A_extended(Valid_peak));
         Radial_temp_trough(i,t) = mean(R_extended(Valid_through));
         Angular_temp_trough(i,t) = mean(A_extended(Valid_through));


        % figure;
        % plot(edges, data_pro, 'k--', 'LineWidth',2 ,'DisplayName', 'Original Signal'); hold on; % 원 신호
        % plot(edges, smooth_data, 'r-', 'LineWidth', 2, 'DisplayName', 'Smoothed Signal'); % 스무딩 신호
        % plot(edges(valid_peak2), peak_values, 'ko', 'MarkerSize', 15,'MarkerFaceColor','k','MarkerEdgeColor','k', 'DisplayName', 'Peaks'); % 마루
        % plot(edges(valid_troughs2), trough_values, 'bo', 'MarkerSize',15,'MarkerFaceColor','b','MarkerEdgeColor','b','DisplayName', 'Troughs'); % 골
        % xlabel('Time');
        % ylabel('Amplitude');
        % title('Gradient-Based Peak and Trough Detection');




   end
end