%% Jiwon Kim, Ph.D. Brown University; Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids

close all;
clear all;
numImages = 161; % Update this with the number of images in your sequence.
% imageFolder = ['']; % Define the folder where your TIFF files are stored.
fileName_prefix = ['Nucleus']; % Update this with your filen name accordingly. 
savename = ['OpticalFlow_Farneback'];

opticFlow = opticalFlowFarneback(FilterSize=15);

h = figure;
movegui(h);
hViewPanel = uipanel(h, 'Position', [0 0 1 1], 'Title', 'Plot of Optical Flow Vectors');
hPlot = axes(hViewPanel);

Vx_cell = cell(numImages, 1);
Vy_cell = cell(numImages, 1);
time = zeros(numImages, 1);

fs = opticFlow.FilterSize;
fs = double(fs);

for i = 1:numImages
    fileName = [fileName_prefix,sprintf('%04d',i-1),'.tif'];     % 
    % fullFileName = fullfile(imageFolder, fileName);
    fullFileName = fileName;
    if ~isfile(fullFileName)
        warning('File does not exist: %s', fullFileName);
        continue;
    end
    frameRGB = imread(fullFileName);
    frameGray = im2gray(frameRGB);
    flow = estimateFlow(opticFlow, frameGray);

    imshow(frameRGB)
    hold on
    plot(flow, 'DecimationFactor', [5 5], 'ScaleFactor', 1, 'Parent', hPlot);

    pause(0.01)
    hold off
    Vx_cell{i, 1} = flow.Vx;
    Vy_cell{i, 1} = flow.Vy;
    time(i) = i;
end

save(savename, 'time', 'Vx_cell', 'Vy_cell', 'fs', '-v7.3')
