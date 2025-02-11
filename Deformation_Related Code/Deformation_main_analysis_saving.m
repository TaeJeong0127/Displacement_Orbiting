clear; close all; clc;
%% This code implements the methods described in:
% Jiwon Kim, Hyuntae Jeong, Carles Falc´o, Alex M. Hruska, W. Duncan 
% Martinson, Alejandro Marzoratti, Mauricio Araiza, Haiqian Yang,
% Christian Franck, Jos´e A. Carrillo, Ming Guo, and Ian Y. Wong,
% "Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids"
% The original code is written by Hyuntae Jeong, Brown Univ., Ian Wong Lab, 2024
% It calculated various type of displacement that originated from 
% the optical flow by implementing the trajectory calculation. 
 
hf = figure;
pos = get(hf,'Position');
close(hf);

%% Taking the name of sample for analysis
% Sample_name = {'G11'; 'G13'; 'G14'; 'G15'; 'G16'; 'M1'; 'M2'; 'M3'; ...
%     'M4'; 'M5'; 'M6'; 'M7'; 'M8'; 'M9'; 'M10'; 'M11'; 'M12'};
% Sample_name = {'C4A'; 'C4B'; 'C5A'; 'C5B'; 'C5C'; 'C7C'; 'D4E'; 'D6B'; ...
%     'D6C'; 'D7B'; 'D8A'; 'E4A2'; 'E6A2'; 'E7C2'; 'E4A';'E6A';'E7A';'E7B';'E7C'};
% Sample_name = {'N2'; 'N6'; 'N9'; 'N11'; 'N15'; 'I1'; 'I2'; 'I3'; ...
%     'I4'; 'I5'; 'I6'; 'I7'; 'I8'; 'I9'; 'I10'; 'I11'; 'I12'; 'I13'; 'I14'; 'I15'; 'I16'; 'I17'; 'I18'};
%% For PEG before sample
% Sample_name = {'E001';'E002';'E041';'E042';'E043';'E044';'E045';'E046';'E401';'E402';'E403';'E404';'E405';'E406';'E407';'E408';'E441';'E442';'E443'};
% For PEG after sample
% Sample_name = {'F001';'F002';'F041';'F042';'F043';'F044';'F045';'F046';'F401';'F402';'F403';'F404';'F405';'F406';'F407';'F408';'F441';'F442';'F443'};
Sample_num = length(Sample_name);

%% color map allocating
cmap_R = flipud(slanCM('RdYlBu',256));
cmap_A = slanCM('PRGn',256);
cmap_A = flipud(cmap_A);
cmap_RedCH = slanCM('Reds',256);
cmap_Deform_M = slanCM('heat',256);
cmap_meanstress = flipud(slanCM('RdYlGn',256));

%% color map data
data1 = load('blue2red.mat');
Cmap_R = data1.colmap;
data2 = load('green2purple.mat');
Cmap_A = data2.colmap;

%% main loop for analyzing the deformation fields

for i = 1:Sample_num
    
    %% Directory for the initial mask file to choose the ROI for deformation
    % Dir_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
    %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask0000.tif'];
    % Mask_initial = imread(Dir_name);
    %% for N,I samples
    % if i <= 5
    %     Dir_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
    %         ' Mask/',char(Sample_name(i)),' Mask0000.tif'];
    % else
    %     Dir_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
    %         ' mask/',char(Sample_name(i)),' Mask0000.tif'];
    % end
    % Mask_initial = imread(Dir_name);

    %     %% for PEG before sample
    % Dir_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
    %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask0000.tif'];
    % Mask_initial = imread(Dir_name);

    % %% for PEG after sample
    % Dir_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241225 peg media change/tiff2/'...
    %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask0000.tif'];
    % Mask_initial = imread(Dir_name);


    % change the domain to 1-0 scale
    Mask_initial = double(Mask_initial);
    Mask_initial = Mask_initial/max(Mask_initial(:));
    Mask_logic = Mask_initial > 0.1;
    Mask_initial = double(Mask_logic);

    % get a centroid of mask
    stats = regionprops(Mask_initial, 'Centroid');
    % Extract the centroid coordinates
    centroid = stats.Centroid*0.65; % This will be a 1x2 vector [x, y]
    % Dilation mask (For getting an enclosing boundary)
    SE = strel("disk",350);
    Mask_dil = imdilate(Mask_initial,SE,8);
    Mask_dropoff = Mask_dil.*~Mask_initial;
    
    %% For Cell displacement results
    Dataname_cell = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Optical flow data_nucleus_Farneback/'...
        ,char(Sample_name(i)),'_nucleus_Farneback_tiff_10um.xlsx'];
    Data_cell = readtable(Dataname_cell);
    time_cell = Data_cell.t;
    X_cell_pos = Data_cell.x; Y_cell_pos = Data_cell.y;
    U_cell = Data_cell.u; V_cell = Data_cell.v;
    time_point_cell = unique(time_cell);
    l_time_cell = length(time_point_cell);
    Time_iter_cell = time_point_cell(2) - time_point_cell(1);
    Ucell_scale = U_cell*Time_iter_cell; Vcell_scale = V_cell*Time_iter_cell;
    

    %% For beads displacement results
    folderPath = '/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Optical flow data_bead_Farneback/';
    filePattern = [char(Sample_name(i)), '_bead_Farneback_tiff_*']; % 패턴 설정
    fileList = dir(fullfile(folderPath, [filePattern, '.xlsx'])); % 패턴을 이용해 파일 검색
    % debugging to prevent when there is no file
    if ~isempty(fileList)
        % choose the first file if name is matched
        Dataname = fullfile(folderPath, fileList(1).name);
        disp(['Selected file: ', Dataname]);
        Data = readtable(Dataname);
    else
        error('No file matching the pattern was found.');
    end
    %% Feature extractions
    time = Data.t;
    Xposition = Data.x; Yposition = Data.y;
    Uvalue = Data.u; Vvalue = Data.v;
    time_point = unique(time);
    l_time = length(time_point);
    Time_iteraction = time_point(2) - time_point(1); 
    Uval_scale = Uvalue*Time_iteraction; Vval_scale = Vvalue*Time_iteraction;
    U_mag_scale = sqrt(Uval_scale.^2 + Vval_scale.^2);
    % step for removing negligible data field to improve the trajectory
    % tracking
    Xposition(U_mag_scale<0.05) = [];  Yposition(U_mag_scale<0.05) = [];
    Uval_scale(U_mag_scale<0.05) = [];  Vval_scale(U_mag_scale<0.05) = [];
    time(U_mag_scale<0.05) = [];

    % Image size
    [m_m,m_n] = size(Mask_dropoff);
    % making a mesh grid for the decreasing computation
    % 0.65 : scale for the um/pixel
    X_max = ceil(m_n*0.65); Y_max = ceil(m_m*0.65);
    X_array = 1:5:X_max; Y_array = 1:5:Y_max;
    [X_pos,Y_pos] = meshgrid(X_array,Y_array);
    % to get how much downregulate the matrix size
    x_pos = X_pos; y_pos = Y_pos; % get a position as a mesh grid
    [m, n] = size(x_pos);
    im_domain = Mask_dropoff;
    [M, N] = size(im_domain);
    Ratio = floor(M/m);
    % down sizing the domain to matching the matirix size with PIV results
    domain_edit = downsample(im_domain,Ratio);
    domain_edit = downsample(domain_edit',Ratio)';
    % to match the domain and PIV results correctly, the center point shoulb be
    % same between domain_edit and PIV_x
    [P,Q] = size(domain_edit);
    offset_y = floor((P-m)/2)+1; offset_x = floor((Q-n)/2)+1;
    if offset_y<0 ||offset_x < 0
        domain = zeros(m,n);
        Ab_off_y = abs(offset_y); Ab_off_x = abs(offset_x);
        domain(Ab_off_y:Ab_off_y+P-1,Ab_off_x:Ab_off_x+Q-1) = domain_edit;
    else
        domain = domain_edit(offset_y:offset_y+m-1,offset_x+1:offset_x+n);
    end

    time_length = l_time; % total number of time frame
    num_traj = sum(domain(:));
    IDX = logical(domain);

    traj_x = nan*zeros(num_traj,time_length);
    traj_y = nan*zeros(num_traj,time_length);
    % First set of trajectories are given by points inside domain (ie, where IDX==1)
    traj_x(:,1) = x_pos(IDX);
    traj_y(:,1) = y_pos(IDX);

    for k = 2:time_length
        ID = find(time==time_point(k));
        Xval = Xposition(ID); Yval = Yposition(ID);
        Uval = Uval_scale(ID); Vval = Vval_scale(ID);
        % Interpolate k-th cell displacements to gridpoints from (k-1)th
        % timepoint
        displ_x(:,k) = griddata(Xval,Yval,Uval,traj_x(:,k-1),traj_y(:,k-1),'cubic');
        displ_y(:,k) = griddata(Xval,Yval,Vval,traj_x(:,k-1),traj_y(:,k-1),'cubic');
        
        displ_x(isnan(displ_x)) = 0; displ_y(isnan(displ_y)) = 0;

        % Add to trajectory arrays
        traj_x(:,k) = traj_x(:,k-1) + displ_x(:,k);
        traj_y(:,k) = traj_y(:,k-1) + displ_y(:,k);

        % For the strain fields
        U_displ(:,k) = traj_x(:,k) - traj_x(:,1);
        V_displ(:,k) = traj_y(:,k) - traj_y(:,1);

        Umap(:,:,k) = griddata(traj_x(:,k),traj_y(:,k),U_displ(:,k),x_pos,y_pos,'cubic');
        Vmap(:,:,k) = griddata(traj_x(:,k),traj_y(:,k),V_displ(:,k),x_pos,y_pos,'cubic');
    end



    for ti = 1:l_time
        if ti<10
            num_ti = ['000',num2str(ti)];
        elseif ti<100
            num_ti = ['00',num2str(ti)];
        elseif ti<1000
            num_ti = ['0',num2str(ti)];
        else
            num_ti = num2str(ti);
        end
        
        %% Loading the cell mask, image, and beads files for analysis and figure plotting
        %% for mosaic samples
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
        %% For PEG data(1st)
        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch1 nucleus/Max_',char(Sample_name(i)),'_C1_t',num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/3 Max Intensity/3 '...
        %     ,char(Sample_name(i)), ' Max intensity/Ch3 bead/Max_',char(Sample_name(i)),'_C3_t',num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/4 Mask/4 '...
        %     ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];

        %% For N~I dat

        % if i <= 5
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240602 mosaic/',char(Sample_name(i)),...
        %     ' Mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % else
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/240724 MMP inhibitor/',char(Sample_name(i)),...
        %     ' mask/',char(Sample_name(i)),' Mask',num_ti,'.tif'];
        % end
        
        %% For PEG before data (1223)
        % 
        % Nuc_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        %     ,char(Sample_name(i)), ' nucleus/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % bead_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        %     ,char(Sample_name(i)), ' bead/',char(Sample_name(i)),num_ti,'.tif'];
        % 
        % Mask_time_name = ['/Volumes/LRSResearch/ENG_BBCancer_Shared/group/Jiwon/Experimental data/241223 peg/'...
        % ,char(Sample_name(i)), ' Mask/' ,char(Sample_name(i)), ' Mask',num_ti,'.tif'];

        %% For PEG after data (1225)
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
        % Image size
        [m_m,m_n] = size(Mask_time);

        X_max = ceil(m_n*0.65); Y_max = ceil(m_m*0.65);
        X_array = 1:5:X_max; Y_array = 1:5:Y_max;
        [X_pos,Y_pos] = meshgrid(X_array,Y_array);

        % to get how much downregulate the matrix size
        x_pos = X_pos; y_pos = Y_pos; % get a position as a mesh grid
        [m, n] = size(x_pos);
        Ratio = floor(m_m/m);
        % down sizing the domain t

        stats = regionprops(Mask_time, 'Centroid');
        % Extract the centroid coordinates
        centroid = stats.Centroid*0.65; % This will be a 1x2 vector [x, y]

        % down sizing the domain to matching the matirix size with PIV results
        domain_edit = downsample(Mask_time,Ratio);
        domain_edit = downsample(domain_edit',Ratio)';
        % to match the domain and PIV results correctly, the center point shoulb be
        % same between domain_edit and PIV_x
        [P,Q] = size(domain_edit);

        offset_y = floor((P-m)/2)+1; offset_x = floor((Q-n)/2)+1;
        if offset_y<0 ||offset_x < 0
            domain = zeros(m,n);
            Ab_off_y = abs(offset_y); Ab_off_x = abs(offset_x);
            domain(Ab_off_y:Ab_off_y+P-1,Ab_off_x:Ab_off_x+Q-1) = domain_edit;
        else
            domain = domain_edit(offset_y:offset_y+m-1,offset_x+1:offset_x+n);
        end

        ID = find(time == time_point(ti));
        Xval = Xposition(ID); Yval = Yposition(ID);
        Uval = Uval_scale(ID); Vval = Vval_scale(ID);

        Umap_optic = griddata(Xval,Yval,Uval,X_pos,Y_pos,'cubic');
        Vmap_optic = griddata(Xval,Yval,Vval,X_pos,Y_pos,'cubic');

        ID_cell = find(time_cell==time_point_cell(ti));
        Xval_cell = X_cell_pos(ID_cell); Yval_cell = Y_cell_pos(ID_cell);
        Ucell = Ucell_scale(ID_cell); Vcell = Vcell_scale(ID_cell);

        Ucell_optic(:,:,ti) = griddata(Xval_cell,Yval_cell,Ucell,X_pos,Y_pos,'cubic');
        Vcell_optic(:,:,ti) = griddata(Xval_cell,Yval_cell,Vcell,X_pos,Y_pos,'cubic');


        REF_x = x_pos - centroid(1);
        REF_y = y_pos - centroid(2);

        angle = atan2(REF_y,REF_x);

        R_U(:,:,ti) = Umap_optic.*cos(angle) + Vmap_optic.*sin(angle);
        A_U(:,:,ti) = Umap_optic.*cos(angle+pi/2) + Vmap_optic.*sin(angle+pi/2);
        R_U_cell(:,:,ti) = Ucell_optic(:,:,ti).*cos(angle) + Vcell_optic(:,:,ti).*sin(angle);
        A_U_cell(:,:,ti) = Ucell_optic(:,:,ti).*cos(angle+pi/2) + Vcell_optic(:,:,ti).*sin(angle+pi/2);

        R_Umap(:,:,ti) = Umap(:,:,ti).*cos(angle) + Vmap(:,:,ti).*sin(angle);
        A_Umap(:,:,ti) = Umap(:,:,ti).*cos(angle+pi/2) + Vmap(:,:,ti).*sin(angle+pi/2);


        Rad = R_U(:,:,ti).*~domain;
        Ang = A_U(:,:,ti).*~domain;
        Umap_do(:,:,ti) = Umap(:,:,ti).*~domain; Vmap_do(:,:,ti) = Vmap(:,:,ti).*~domain;
        Udis_do(:,:,ti) = Umap_optic(:,:).*~domain; Vdis_do(:,:,ti) = Vmap_optic(:,:).*~domain;
        domain_out = double(~domain);
        domain_out(domain_out==0) = NaN;
        domain(domain==0) = NaN;

        % im_k = imread(Nuc_name);
        % [M,N] = size(im_k);
        % 
        % im_b = imread(bead_name);
        % pix_size = 0.65;
        Mag_displ = sqrt(Umap_optic.^2+Vmap_optic.^2);
        Mag_deform = sqrt(Umap(:,:,ti).^2 + Vmap(:,:,ti).^2);
        ANG_deform = atan2(Umap(:,:,ti),Vmap(:,:,ti)).*domain_out;
        %   ANG_deform = atan2(Umap_optic,Vmap_optic).*domain_out;


    end
    
    result_name = ['Results_compliation_',char(Sample_name(i)),'_new0120.mat'];
    save(result_name,"X_pos","Y_pos","Ucell_optic","Vcell_optic","Umap_do","Vmap_do","Udis_do","Vdis_do","traj_x","traj_y",'-mat')

    clear var displ_x displ_y ANG_deform Mag_deform Mag_displ A_U A_U_cell A_Umap Ang R_Umap R_U_cell ...
        R_U Rad Ang Ucell_optic Vcell_optic U_displ V_displ Umap Vmap Umap_do Vmap_do
end
