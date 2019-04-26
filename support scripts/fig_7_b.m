function fig_7_b(data)
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
%color_index = {[252,146,114]./256,[239,59,44]./256,[153,0,13]./256}; %white to dark red with increasing alpha
color_index = {[254,178,76]./256,[252,78,4]./256,[177,0,38]./256}; %orange to dark red with increasing alpha
data_index = {'m_3','Hgb_20','Hgb_40'};
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

figure('name','Figure 7: Impact of plasma free hemoglobin on NO dynamics','NumberTitle','off',...
    'position',pos,'color','w');

%% plot baseline vasodilation
subplot(2,3,1), hold on
plot(data.fig_7.baseline.Hgb(:,1),data.fig_7.baseline.Hgb(:,2),'-sk','LineWidth',3,'MarkerFaceColor',[0 0 0])
plot(data.fig_7.baseline.Hgb_feedback(:,1),data.fig_7.baseline.Hgb_feedback(:,2)./3,'-sr','LineWidth',3,'MarkerFaceColor',[0 0 0])
ylabel('\DeltaGC activation (%)')
xlabel('plasma free Hemoglobin (\muM)')
title('increased pfHgb decreases vessel diameter')

a = -15; b = 5;
ylim([a b])

yyaxis right
plot(data.fig_7.baseline.Hgb(:,1),data.fig_7.baseline.Hgb(:,2).*3,'k','LineWidth',1)
plot(data.fig_7.baseline.Hgb_feedback(:,1),data.fig_7.baseline.Hgb_feedback(:,2),'r','LineWidth',1)
ylabel('\Deltadilation (m = 3)')
ylim([a*3 b*3])

%% plot perivascular NO concentrations
subplot(2,3,2), hold on
for ii = 1:length(data_index)
    plot(data.fig_7.(data_index{ii}).perivascular(:,1),data.fig_7.(data_index{ii}).perivascular(:,2),'Color',[color_index{ii}],'LineWidth',2)
end
kk = 1;
for cfl = [2.8 2.8 2.8]
f = fill([10 10-cfl 10-cfl 10],[0 0 14 14],'r','Linestyle','none');
set(f,'facea',1/10)
set(f,'facec',color_index{kk})
kk = kk + 1;
end

f2 = fill([13 11 11 13],[0 0 14 14],'r','Linestyle','none'); %smooth muscle
set(f2,'facea',1/5)
set(f2,'facec',[241,163,64]./256)

axis([0 50 0 14])
ylabel('NO (nM)')
xlabel('radius (\mum)')
title('Perivascular NO')
legend('1 \muM','20 \muM','40 \muM')

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
legend([h1{1} h1{2} h1{3}],{'1','20','40'})
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
legend([h2{1} h2{2} h2{3}],{'1','20','40'})
set(gca,'XScale','log','YScale','linear')
xlim([0.025 0.5])
%ylim([10^-8 10])
xlabel('frequency (Hz)')
ylabel('Power')
title('vasodynamics')
yt=[0.01:0.01:0.09 0.1:0.1:0.9 1];
set(gca,'xtick',yt)
set(gca,'xminorgrid','on')

%% figure
figure, hold on
MAP = 120;
plot(data.fig_7.baseline.Hgb(:,1),(MAP-((100+data.fig_7.baseline.Hgb(:,2).*3)/100).^4.*MAP)./MAP*100,'-sk','LineWidth',3,'MarkerFaceColor',[0 0 0])
plot(data.fig_7.baseline.Hgb_feedback(:,1),(MAP-((100+data.fig_7.baseline.Hgb_feedback(:,2))/100).^4.*MAP)./MAP*100,'-sr','LineWidth',3,'MarkerFaceColor',[0 0 0])
ylabel('\DeltaMAP (%)')
xlabel('\DeltaHemoglobin (\muM)')
title('\DeltaMAP from increased NO absorption')
% xlim([-10 40])
% ylim([-20 10])
grid on


end
