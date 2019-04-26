function fig_6_a(rootfolder)
%% load/refine data

NO_sens = [1 3 5 10];
time_end = [20 20 30 50];
clear WF
dt = 0.01;
LineColor = {[1 0 0],[0 0 1]};

for ii = 1:length(NO_sens)
    time{ii} = [0:dt:time_end(ii)];
    hold_data = importdata(['1.2NO_1NOd_' num2str(NO_sens(ii)) 'NOsens_50bGC_0Hz_10VD_6sNOkernel_Oscillations.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    WF{1,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')'; %concentration
    WF{1,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')'; %dilation
    
     time{ii} = [0:dt:time_end(ii)];
    hold_data = importdata(['1.2NO_1NOd_' num2str(NO_sens(ii)) 'NOsens_50bGC_0Hz_10VD_6sNOkernel_OscillationsStatic.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    WF{2,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')'; %concentration
    WF{2,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')'; %dilation
end

%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 600; ht = 400; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*16/20)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Figure 6 a: Coupling NO degradation to vasodilation yields oscillations if the vessel is sufficiently senstive to NO','NumberTitle','off',...
    'position',pos,'color','w');

kk = 1;
for ii = 1:length(NO_sens)
    hold_conc = movmean(WF{kk,1,ii}(3:end),25);
    subplot(2,1,1), hold on, plot(time{ii}(3:end),hold_conc,'Color',[LineColor{kk} (ii)/4],'LineWidth',2)
    title('Concentration [nM]'),xlabel('time [s]'), ylabel('NO in smooth muscle'),%axis([6 30 5 20])
    subplot(2,1,2), hold on, plot(time{ii},WF{kk,2,ii},'Color',[LineColor{kk} (ii)/4],'LineWidth',2),%axis([6 30 -5 40])
    title('Vessel Response'),xlabel('time [s]'), ylabel('\DeltaVessel Diameter [%]')
end
subplot(211),legend('x1','x3','x5','x10')
subplot(212),legend('x1','x3','x5','x10')


kk = 2;
ax1 = axes('Position',[.6 .7 .2 .2]);
box off
hold(ax1,'on')
ax2 = axes('Position',[.6 .3 .2 .2]);
box off
hold(ax2,'on')
for ii = 1:length(NO_sens)
    hold_conc = movmean(WF{kk,1,ii}(3:end),25);
    hold on, plot(ax1,time{ii}(3:end),hold_conc,'Color',[LineColor{kk} (ii)/4],'LineWidth',2);
    title(ax1,'static RBC core'),xlabel(ax1,'time [s]'), ylabel(ax1,'NO in smooth muscle');xlim([0 20])
    hold on, plot(ax2,time{ii},WF{kk,2,ii},'Color',[LineColor{kk} (ii)/4],'LineWidth',2);axis([0 20 0 20])
    title('static RBC core'),xlabel('time [s]'), ylabel('\DeltaVessel Diameter [%]');
end

    legend(ax1,'x1','x3','x5','x10');
    legend(ax2,'x1','x3','x5','x10');


end

function [output] = getCOMSOLsimulation(input)

output = input; %carry over variables

stim_duration = {'0.7','10'};
t_end = [20 30];
ii_end = 2;

for ii = 1:ii_end
time{ii} = [0:0.01:t_end(ii)];
hold_data = importdata(['1.4NO_' stim_duration{ii} 'NOd_3NOsens_50bGC_0Hz_10VD_6sNOkernel_static_PNAS.csv']);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
WF{1,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')';
WF{1,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')';

hold_data = importdata(['1.4NO_' stim_duration{ii} 'NOd_3NOsens_50bGC_0Hz_10VD_6sNOkernel_PNAS.csv']);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
WF{2,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')';
WF{2,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')';
end

offset  = 6.5; %stimulus starts 6.5 seconds into simulation to allow time for 6 sec vessel response kernel

for ii = 1:ii_end
time{ii} = [0:0.01:t_end(ii)];
hold_data = importdata(['1.4NO_' stim_duration{ii} 'NOd_3NOsens_50bGC_0Hz_10VD_6sNOkernel_static_PNAS.csv']);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
WF{1,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')';
WF{1,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')';

hold_data = importdata(['1.4NO_' stim_duration{ii} 'NOd_3NOsens_50bGC_0Hz_10VD_6sNOkernel_PNAS.csv']);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
WF{2,1,ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')';
WF{2,2,ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')';
end

output.comsol.static.onePuff.time = time{1}-offset;
output.comsol.static.tenSecPuff.time = time{2}-offset;
output.comsol.dynamic.onePuff.time = time{1}-offset;
output.comsol.dynamic.tenSecPuff.time = time{2}-offset;

output.comsol.static.onePuff.conc = WF{1,1,1};
output.comsol.static.onePuff.dilation = WF{1,2,1};
output.comsol.static.tenSecPuff.conc = WF{1,1,2};
output.comsol.static.tenSecPuff.dilation = WF{1,2,2};

output.comsol.dynamic.onePuff.conc = WF{2,1,1};
output.comsol.dynamic.onePuff.dilation = WF{2,2,1};
output.comsol.dynamic.tenSecPuff.conc = WF{2,1,2};
output.comsol.dynamic.tenSecPuff.dilation = WF{2,2,2};

end