function fig_6_c(data)
%% load/refine data
%[WGN] = get_mtspectrumFromWGN();

WGN = data.fig_6.c;

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

figure('name','Figure 6 c: increased NO degradation during vasodilation produces resonant frequencies similar to vasomotion','NumberTitle','off',...
    'position',pos,'color','w');

LineColor = {[1 0 0],[0 0 1]};

hold on
for ii = 1:length(WGN.dynamic.f)
loglog(WGN.dynamic.fr,WGN.dynamic.f{ii},'Color',[LineColor{1} 1/10],'Linewidth',2);
loglog(WGN.static.fr,WGN.static.f{ii},'Color',[LineColor{2} 1/10],'Linewidth',2);   
end

h1 = loglog(WGN.dynamic.fr,WGN.dynamic.avg_f,'Color',[LineColor{1}],'Linewidth',2);
h2 = loglog(WGN.static.fr,WGN.static.avg_f,'Color',[LineColor{2}],'Linewidth',2);
set(gca,'XScale','linear','YScale','linear'), xlim([0 1])
legend([h1 h2],{'WGN dynamic','WGN static'})
title('vessel response from gamma band')
xlabel('frequency [Hz]')
ylabel('Power')

end

function [output] = get_mtspectrumFromWGN()
dt = 0.01;
LineColor = {[1 0 0],[0 0 1]};
RBC_core = {'','_static'}
figure, hold on
for l = 1:2 %iterate over dynamic and static cases
    for g = 1:10 %iterate over number of 5 minute trials
        time = [15:dt:300];
        hold_data = importdata(['2NO_0NOd_3NOsens_50bGC_10' num2str(g) 'Hz_10VD_6sNOkernel_GammaBand' RBC_core{l} '.csv']);
        [a index] = unique(hold_data(:,1));
        hold_data = hold_data(index,:); %don't take replicate values
        WF{l,1,g} = interp1(hold_data(:,1),hold_data(:,2),time')'; %concentration
        WF{l,2,g} = interp1(hold_data(:,1),hold_data(:,3),time')'; %dilation
        
        dilation_hold = detrend(WF{l,2,g}); %mean subtract
        
        params.Fs = 1/dt;
        params.tapers = [5 9];
        params.fpass = [0.025 3];
        
        if l == 1
            output.dynamic.dilation{g} = dilation_hold;
            [output.dynamic.f{g} output.dynamic.fr] = mtspectrumc(dilation_hold,params);
            output.dynamic.f{g} = output.dynamic.f{g};
        else
            output.static.dilation{g} = dilation_hold;
            [output.static.f{g} output.static.fr] = mtspectrumc(dilation_hold,params);
            output.static.f{g} = output.static.f{g};
        end
        
    end
end

output.dynamic.avg_f = mean(cell2mat(output.dynamic.f),2);
output.static.avg_f = mean(cell2mat(output.static.f),2)
end
