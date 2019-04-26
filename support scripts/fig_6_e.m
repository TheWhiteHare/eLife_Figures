function fig_6_e(data)
%% load/refine data
%[oscillate] = get_OscillationsWithVaryingBaselineGC();

oscillate = data.fig_6.e;

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

figure('name','Figure 6 d: NO consumption coupled to dilation is a linear system','NumberTitle','off',...
    'position',pos,'color','w'); hold on

hold on
plot(oscillate.baselineGC,oscillate.amplitude,'k','LineWidth',2)
scatter(oscillate.baselineGC,oscillate.amplitude,'k','fill')
title('oscillation amplitude is maximized when baseline GC activity is low'),xlabel('baseline GC (%)'), ylabel('relative amplitude of oscillations (a.u.)'),%axis([6 30 5 20])
axis square

end

function [output] = get_OscillationsWithVaryingBaselineGC()

dt = 0.01;
bGC = [1.8591 4.1810 8.0808 16.868 33.121 55.206 75.530 88.566 95.110];
time = [15:dt:30]; %get 15 seconds into trial to avoid edge effects

for ii = 1:length(bGC)
hold_data = importdata(['2NO_0NOd_3NOsens_' num2str(bGC(ii)) 'bGC_0.16Hz_10VD_6sNOkernel_GCBaselines.csv']);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
WF{1,ii} = interp1(hold_data(:,1),hold_data(:,2),time')'; %concentration
WF{2,ii} = interp1(hold_data(:,1),hold_data(:,3),time')'; %dilation
end

MAX = cellfun(@max, WF); %maximum size of vessel
MIN = cellfun(@min, WF); %minumum size of vessel

Amp = (MAX - MIN)/2; %amplitude of vessel oscillations

output.amplitude = Amp(2,:)/max(Amp(2,:)); %normalize to 1
output.baselineGC = bGC; %remember baseline GC activity of smooth muscle
end