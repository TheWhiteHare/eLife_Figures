function fig_4_e_f(data)
%% load data


conc = data.fig_4.conc;
dilation = data.fig_4.dilation;
time = data.fig_4.time;

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

figure('name','Figure 4 e,f: Coupling NO concentration to vessel size and NO degradation reproduces a post stimulus undershoot from a purely dilatory stimulus','NumberTitle','off',...
    'position',pos,'color','w');

LineColor = {[1 0 1],[0 0 0]};

for ii = 1:3
    
    subplot(2,1,1), hold on, plot(time,movmean(conc{ii},20),'Color',[LineColor{2} ii/3],'LineWidth',2)
    title('Concentration [nM]'),xlabel('time [s]'), ylabel('NO in smooth muscle'), xlim([5 15])
    subplot(2,1,2), hold on, plot(time,dilation{ii},'Color',[LineColor{2} ii/3],'LineWidth',2)
    title('Vessel Response'),xlabel('time [s]'), ylabel('\DeltaVessel Diameter [%]'), axis([5 15 -5 25])
end

    ax1 = axes('Position',[.7 .25 .2 .2]);
    box off
    hold(ax1,'on')
for ii = 1:3
    plot(ax1,time,dilation{ii}/max(dilation{ii}),'Color',[LineColor{2} ii/3],'LineWidth',3)
    title([{'Relative Undershoot'}]),xlabel('time [s]'), xlim([5 15]), axis square
end
end


