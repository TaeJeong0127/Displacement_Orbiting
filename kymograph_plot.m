%% Jiwon Kim, Ph.D. Brown University; Collective Transitions from Orbiting to Matrix Invasion in 3D Multicellular Spheroids


clear all; close all;

timestep = 161;
dt = 0.25;
theta_origin = 93.5;

savename = 'Kymograph_plot';
kymodata = 'Kymograph_data_slope';

load(kymodata)
kymo1 = slope(:,1:timestep);
kymo1 = circshift(kymo1,-round(theta_origin/2)+90);

kymo1_output = locallapfilt(kymo1, 0.9, 0.1);

[h, w] = size(kymo1_output);

grey = gray;
uicmap2 = (grey);

q = imagesc(kymo1_output);
colormap(uicmap2)
q.AlphaData = 0.9;

xratio = 1.5;
xtick = 12;
yspace = linspace(1, h, 5);
yticks(yspace);
yticklabels({'180','90','0','-90','-180'});     % deg
% yticklabels({'{\pi}','1/2{\pi}','0','-1/2{\pi}','-{\pi}'});   % rad
ylabel('Angle,\theta ({\circ})')

xspace = 0:(xtick)/dt:w;
xticks(xspace);
xticknum = length(xspace);
xtickstring = string(0:xtick:w*dt); 
xticklabels(xtickstring)


xlabel('Time (hour)')
xtickangle(0)
axis tight

fontsize(23,'point');
fontname('Helvetica Neue')
set(gcf,'PaperPositionMode','auto');
set(gcf,'color', 'none');
set(gca,'color', 'none');
set(gcf,'InvertHardcopy','off');
pbaspect([3 1 1])

savekymo = [savename, '.tif'];
exportgraphics(gcf,savekymo,'Resolution',300)
