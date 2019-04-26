function fig_7_s(data)
%________________________________________________________________________________________________________________________
% Written by W. Davis Haselden
% M.D./Ph.D. Candidate, Department of Neuroscience
% The Pennsylvania State University & Herhsey College of Medicine
%________________________________________________________________________________________________________________________
%
%   Purpose: plot how the changing the plasma free hemoglobin in the cell
%   free layer alters vasodynamics, baseline vessel diameter, the
%   hemodynamic response function, and perivascular NO concentrations 
%   pfHgb = 1, 20, 40 uM
%________________________________________________________________________________________________________________________
%
%   Inputs: data.fig_7 - contains NO production rate, NO in smooth muscle,
%   and vasodynamics of 1500 second long simulation of a penetrating arteriole
%   responding to a white gaussian noise production rate of NO low pass
%   filtered below 2 Hz.
%
%   Outputs: figure reporting example vasodynamics of pfHgb = 1, 20, 40 uM.
%   baseline vessel diameter from changing pfHgb. perivascular NO
%   concentrations, HRF, and power spectrum of vasodynamics from changing
%   plasma free hemoglobin.
%________________________________________________________________________________________________________________________
%%
color_index = {[0,0,256]./256,[256,0,0]./256}; %orange to dark red with increasing alpha
data_index = {'static','m_3'};
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

figure('name','Figure 7: vasomotion','NumberTitle','off',...
    'position',pos,'color','w');
%%
subplot(2,3,1)
%%
subplot(2,3,2)
%% variance change in vasodynamics from pfHgb = 1, 20, 40
subplot(2,3,3), hold on
for ii = 1:length(data_index)
notBoxPlot([data.fig_7.(data_index{ii}).variance],[ii]);
end
set(gca,'XTick',[1 2 3])
set(gca,'XTickLabels',{'1','20','40'})
ylabel('variance')
title('baseline vasodynamics')
ylim([0 2])
xlim([0.5 3.5])

%% example vasodynamics
subplot(2,3,4), hold on
time_all = [15:dt:1500];
for ii = 1:length(data_index)
    plot(time_all,data.fig_7.(data_index{ii}).dilation,'Color',[color_index{ii}],'LineWidth',2)
end
xlabel('time (s)')
ylabel('\Deltavessel diameter')
title('example trial vasodynamics')
%legend({'NOx1','NOx3','NOx5'})
xlim([220 250])
ylim([-5 5])

%% hemodynamic Response Function when pfHgb = 1, 20, 40
subplot(2,3,5), hold on
time = [0:dt:15];
for ii = 1:length(data_index)
    h1{ii} = plot(time,data.fig_7.(data_index{ii}).HRF,'Color',[color_index{ii}],'LineWidth',3);
end
axis([0 10 -20 40])
legend([h1{1} h1{2}],{'static','dynamic'})
xlabel('time (s)')
ylabel('a.u.')
title('hemodynamics response function')

%% Frequency response of the system when pfHgb = 1, 20, 40
subplot(2,3,6), hold on
Hz = data.fig_7.(data_index{1}).spectrumHz;

for ii = 1:length(data_index)
    h2{ii} = plot(Hz,data.fig_7.(data_index{ii}).spectrumPower,'Color',[color_index{ii}],'LineWidth',3);
    
    high = data.fig_7.(data_index{ii}).Serr{1}(1,:); %high CI - jack-knife resampling p = 0.05
    low = data.fig_7.(data_index{ii}).Serr{1}(2,:); %low CI - jack-knife resampling p = 0.05
    
    f = fill([Hz flip(Hz)],[high flip(low)],'r','Linestyle','none');
    set(f,'facea',1/5);
    set(f,'facec',color_index{ii});
end
legend([h2{1} h2{2}],{'static','dynamic'})
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




























