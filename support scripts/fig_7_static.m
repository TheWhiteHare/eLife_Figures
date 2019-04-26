function fig_7_static(data)
%________________________________________________________________________________________________________________________
% Written by W. Davis Haselden
% M.D./Ph.D. Candidate, Department of Neuroscience
% The Pennsylvania State University & Herhsey College of Medicine
%________________________________________________________________________________________________________________________
%
%   Purpose: plot how the changing the sensitivity of a arteriole to NO
%   alters vasodynamics (slope, m = 1/3/5)
%________________________________________________________________________________________________________________________
%
%   Inputs: data.fig_7 - contains NO production rate, NO in smooth muscle,
%   and vasodynamics of 1500 second long simulation of a penetrating arteriole
%   responding to a white gaussian noise production rate of NO low pass
%   filtered below 2 Hz.
%
%   Outputs: figure reporting example of m = 1/3/5 for GC activation:vessel
%   diameter relationship, variance change, HRF, Frequency spectrum of
%   vasculat response, example (20s long snapshot) of vessel response to
%   the same white noise stimulus with m = 1/3/5.
%________________________________________________________________________________________________________________________
%%
color_index = {[0 0 1],[0 0 1],[166,189,219]./256,[54,144,192]./256,[1,100,80]./256}; %blue to dark green with increasing alpha
data_index = {'static_m_1','static','m_1','m_3','m_5'};
dt = 1/6;

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

figure('name','Figure 7: Sensitivity to NO alters vasodynamics','NumberTitle','off',...
    'position',pos,'color','w');



%% static HRF
subplot(2,2,3), hold on
time = [0:dt:15];
for ii = 1:2
    h1{ii} = plot(time,data.fig_7.(data_index{ii}).HRF,'Color',[color_index{ii} ii/3],'LineWidth',3);
end
axis([0 10 -20 40])
%legend([h1{1} h1{2} h1{3}],{'m = 1','m = 3','m = 5'})
xlabel('time (s)')
ylabel('a.u.')
title('hemodynamics response function')

%% static frequency response
subplot(2,2,4), hold on
Hz = data.fig_7.(data_index{1}).spectrumHz;

for ii = 1:2
    h2{ii} = plot(Hz,data.fig_7.(data_index{ii}).spectrumPower,'Color',[color_index{ii} ii/3],'LineWidth',3);
    
    high = data.fig_7.(data_index{ii}).Serr{1}(1,:) %high CI - jack-knife resampling p = 0.05
    low = data.fig_7.(data_index{ii}).Serr{1}(2,:) %low CI - jack-knife resampling p = 0.05
    
    f = fill([Hz flip(Hz)],[high flip(low)],'r','Linestyle','none');
    set(f,'facea',1/5);
    set(f,'facec',color_index{ii})
end
%legend([h2{1} h2{2} h2{3}],{'m = 1','m = 3','m = 5'})
set(gca,'XScale','log','YScale','linear')
xlim([0.025 0.5])
%ylim([10^-8 10])
xlabel('frequency (Hz)')
ylabel('Power')
title('vasodynamics')
yt=[0.01:0.01:0.09 0.1:0.1:0.9 1];
set(gca,'xtick',yt)
set(gca,'xminorgrid','on')

%% dynamic HRF
subplot(2,2,1), hold on
time = [0:dt:15];
for ii = 3:5
    h1{ii} = plot(time,data.fig_7.(data_index{ii}).HRF,'Color',[color_index{ii}],'LineWidth',3);
end
axis([0 10 -20 40])
%legend([h1{1} h1{2} h1{3}],{'m = 1','m = 3','m = 5'})
xlabel('time (s)')
ylabel('a.u.')
title('hemodynamics response function')

%% dynamic frequency response
subplot(2,2,2), hold on
Hz = data.fig_7.(data_index{1}).spectrumHz;

for ii = 3:5
    h2{ii} = plot(Hz,data.fig_7.(data_index{ii}).spectrumPower,'Color',[color_index{ii}],'LineWidth',3);
    
    high = data.fig_7.(data_index{ii}).Serr{1}(1,:) %high CI - jack-knife resampling p = 0.05
    low = data.fig_7.(data_index{ii}).Serr{1}(2,:) %low CI - jack-knife resampling p = 0.05
    
    f = fill([Hz flip(Hz)],[high flip(low)],'r','Linestyle','none');
    set(f,'facea',1/5);
    set(f,'facec',color_index{ii})
end
%legend([h2{1} h2{2} h2{3}],{'m = 1','m = 3','m = 5'})
set(gca,'XScale','log','YScale','linear')
xlim([0.025 0.5])
%ylim([10^-8 10])
xlabel('frequency (Hz)')
ylabel('Power')
title('vasodynamics')
yt=[0.01:0.01:0.09 0.1:0.1:0.9 1];
set(gca,'xtick',yt)
set(gca,'xminorgrid','on')

end