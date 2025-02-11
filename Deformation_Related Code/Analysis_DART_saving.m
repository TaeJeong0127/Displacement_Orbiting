clear; close all; clc;
%% This code implements the methods described in:
% Jiwon Kim, Hyuntae Jeong, Carles Falc´o, Alex M. Hruska, W. Duncan 
% Martinson, Alejandro Marzoratti, Mauricio Araiza, Haiqian Yang,
% Christian Franck, Jos´e A. Carrillo, Ming Guo, and Ian Y. Wong,
% "Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids"
% The original code is written by Hyuntae Jeong, Brown Univ., Ian Wong Lab,
% 2024 and Mauricio Araiza, Christian Frank's Lab
% The quiver plots for the displacement fields are generated with this code.

%% plot_tractions from the Christian Frank Lab DART plots.
Sample_name = {'G11'; 'G13'; 'G14'; 'G15'; 'G16'; 'M1'; 'M2'; 'M3'; ...
   'M4'; 'M5'; 'M6'; 'M7'; 'M8'; 'M9'; 'M10'; 'M11'; 'M12'};
% Sample_name = {'C4A'; 'C4B'; 'C5A'; 'C5B'; 'C5C'; 'C7C'; 'D4E'; 'D6B'; ...
%      'D6C'; 'D7B'; 'D8A'; 'E4A2'; 'E6A2'; 'E7C2'};
%Sample_name = {'E4A'; 'E6A'; 'E7C'};

% Sample_name = {'N2'; 'N6'; 'N9'; 'N11'; 'N15'; 'I1'; 'I2'; 'I3';
%     'I4'; 'I5'; 'I6';'I7';'I8';'I9';'I10';'I11';'I12';'I13'};
Sample_num = length(Sample_name);
spacing = 5; % spacing for vector maps

for i = 1:Sample_num
    % savenameheader = 'Leader_test';
    % foldername = ['DART_plot_franklab_temporal_new',char(Sample_name(i))];
    % mkdir(foldername)
    % fol2 = [foldername,'/'];
    file_name = ['Results_compliation_',char(Sample_name(i)),'_new.mat'];
    load(file_name);

    time_length = size(Ucell_optic,3);
    o = 4;
    % for the M series = /2, G series =/4
    m = 2;



savearray = [];

    for t = 1+m:time_length-(m+1)

        ti = t;

        U_cell_track_line = []; V_cell_track_line = [];
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

    % G,M sample
        if i <= 8
            Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
                ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];


            bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
                ,char(Sample_name(i)), ' Max intensity/Ch2 bead/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];

            Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
                ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        else
            Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
                ,char(Sample_name(i)), ' Max intensity/Ch2 nucleus/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];


            bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
                ,char(Sample_name(i)), ' Max intensity/Ch3 bead/Max_',char(Sample_name(i)),'_C3_t',num_ti,'.tif'];


            Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
                ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];
        end

        % 
        % % For PEG data
        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/2 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch2 bead/Max_',char(Sample_name(i)),'_C2_t',num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];

        % N,I sample
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
        edges = -pi:pi/18:pi;
        IDX = discretize(thetas, edges);
        DARTvalues = zeros(12, 4);
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
        
        Mag_U = sqrt(Umap_prin.^2 + Vmap_prin.^2);
        percentile_95 = prctile(Mag_U(:), 95, 'all', 'Method','approximate');

        for ii = 1:36
            tempBin = R_map_prin(IDX == ii);


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
        
        header = ti*ones(36,1);
        savearraytemp = [header, DARTvalues];
        savearray = [savearray;savearraytemp];

        Linear_DART(t,:) = DARTvalues(:,2);
    end


    savetable = array2table(savearray, 'VariableNames',{'slide','protrusion','contraction','counterclockwise','clockwise'});
    saveexcelname = [char(Sample_name(i)),'DART_deformation_1216'];
    writetable(savetable, saveexcelname)

    
    Mean_value = mean(Linear_DART,2);
    Angle = linspace(5,355,36);
    Angle_2 = Angle*2;
    Angle_2(Angle_2>=360)= Angle_2(Angle_2>=360)-360;

    for t = 1+m:time_length-(m+1)
      Threshold_val = prctile(Linear_DART(t,:),80); 
      target_angle = Angle_2(Linear_DART(t,:)>Mean_value(t));
      radian = deg2rad(target_angle);
      num_target = length(target_angle);
      sum_sin = sum(sind(target_angle)); sum_cos = sum(cosd(target_angle));
      Y_a = sum_sin/num_target; X_a = sum_cos/num_target;
      r_a = sqrt(X_a.^2+Y_a.^2);
      sin_ref = Y_a/r_a; cos_ref = X_a/r_a;
      Theta(t) = atand(sin_ref/cos_ref);
      [pval,z] = circ_rtest(radian);
      P_value(t) = pval;

    end



end

figure()
plot(Theta);
ylim([-90 90]);
figure()
plot(P_value)
ylim([0 1]);