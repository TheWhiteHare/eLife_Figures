function fig_1(rootfolder)
%% load data


%% refine data


%% Generate figure for illustrator: only plot lines and fill this later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 200; ht = 400; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Predicting perivascular NO concentrations and smooth muscle signaling',...
    'position',pos,'color','w');
subplot(211)
% pt1 = shadedErrorBar(CBV.LTA.time, nanmean(CBV.LTA.RP),nanstd(CBV.LTA.RP),...
%    'g',1);
% hold on;
% pt2 = shadedErrorBar(CBV.LTA.time, nanmean(CBV.LTA.RF),nanstd(CBV.LTA.RF),...
%   'c',1);
legend([pt1.mainLine, pt2.mainLine],{'FL/HL (n = 11)', 'FC (n = 11)'});
ylabel('\DeltaR/R_0 (%)');
xlabel('Time from locomotion onset (s)');
axis([-3 5 -7 2]);
subplot(212)
%pt1 = plot()
hold on;
%pt2 = plot()
legend([pt1.mainLine, pt2.mainLine],{'FL/HL (n = 5)', 'FC (n = 5)'});
ylabel('CBF (%)');
xlabel('Time from locomotion onset (s)');
axis([-3 5 -7 20]);
end
