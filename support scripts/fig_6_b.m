function fig_6_b(data)
%% load/refine data


% Hz = [0.05 0.1 0.15 0.2 0.25 0.5 0.75 1]; %frequencies tested
% [Amp] = calculateAmplitude(Hz);

Amp = data.fig_6.b;

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

figure('name','Figure 6 b: Coupling NO degradation to vasodilation yields oscillations if the vessel is sufficiently senstive to NO','NumberTitle','off',...
    'position',pos,'color','w');



hold on, grid on
title('Frequency Response Dyanamic RBC Core'),xlabel('frequency [Hz]'), ylabel('Amplitude'),xlim([0 1])
scatter(Amp.Hz,Amp.dynamic,'ro','fill')
h1 = plot(Amp.Hz_pchip,Amp.dynamic_pchip,'r-');

scatter(Amp.Hz,Amp.static,'bo','fill')
h2 = plot(Amp.Hz_pchip,Amp.static_pchip,'b-');

legend([h1 h2],'Dynamic RBC Core','Static RBC Core')

end

function [output] = calculateAmplitude()

time_end = [60 30 30 30 30 30 30 15]; %duration of simulation
ii_end = length(Hz);
LineColor = {[1 0 0],[0 0 1]}; %red = dynamic, blue = static
start_time = 15; %seconds; begin analysis 15 seconds into simultion
dt = 0.005;
condition = {'','static_'};

for con = 1:length(condition) %dynamic
for ii = 1:ii_end
    time = [0:dt:time_end(ii)];
    hold_data = importdata(['2NO_1NOd_3NOsens_50bGC_' num2str(Hz(ii)) 'Hz_10VD_6sNOkernel_' condition{con} 'positive.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    WF{1,ii} = interp1(hold_data(:,1),hold_data(:,2),time')'; %concentration
    WF{2,ii} = interp1(hold_data(:,1),hold_data(:,3),time')'; %dilation
    
    Amp(ii) = (max(WF{2,ii}(start_time/dt+1:end)) - min(WF{2,ii}(start_time/dt+1:end)))/2; %calculate amplitude of vessel oscillations
end
    if con == 1 %dynamic RBC core case
        output.dynamic = Amp;
        Hz_smooth = [0:0.01:1];
        output.dynamic_pchip = pchip(Hz,Amp,Hz_smooth);
    else %static RBC core case
        output.static = Amp;
        Hz_smooth = [0:0.01:1];
        output.static_pchip = pchip(Hz,Amp,Hz_smooth);
    end
    
    output.Hz_pchip = Hz_smooth;
end

end