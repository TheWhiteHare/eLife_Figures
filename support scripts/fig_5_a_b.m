function fig_5_a_b(data)
%% load/refine data

% [model] = getImpulseResponse();
% [model] = getCOMSOLsimulation(model);

model = data.fig_5.model;

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

figure('name','Figure 5 a,b: Increased NO consumption during vasodilation replicates the indershoot/overshoot during NVC events','NumberTitle','off',...
    'position',pos,'color','w');

subplot(2,1,1), hold on, 
plot(model.measured.onePuff.time,model.measured.kernel_norm,'r','LineWidth',2)
plot(model.fit.time,model.fit.kernel_norm,'m','LineWidth',2)
title(['R^{2} = ' num2str(model.fit.rsq) ' fit between 1-4 sec'])
legend('Impulse Response','fitted Kernel')

LineColor = {[1 0 0],[0 0 1],[0 0 0]};
subplot(2,1,2), hold on
plot(model.measured.onePuff.time,model.measured.onePuff.dilation,'Color',[LineColor{3} 1/2],'LineWidth',3)
plot(model.comsol.static.onePuff.time,model.comsol.static.onePuff.dilation/100,'Color',[LineColor{2} 1/2],'LineWidth',3)
plot(model.comsol.dynamic.onePuff.time,model.comsol.dynamic.onePuff.dilation/100,'Color',[LineColor{1} 1/2],'LineWidth',3)

plot(model.measured.tenSecPuff.time,model.measured.tenSecPuff.dilation,'Color',[LineColor{3}],'LineWidth',3)
plot(model.comsol.static.tenSecPuff.time,model.comsol.static.tenSecPuff.dilation/100,'Color',[LineColor{2}],'LineWidth',3)
plot(model.comsol.dynamic.tenSecPuff.time,model.comsol.dynamic.tenSecPuff.dilation/100,'Color',[LineColor{1}],'LineWidth',3)


title('Vessel Response'),xlabel('time [s]'), ylabel('\DeltaVessel Diameter [%]')
legend('0.7s true','0.7s static','0.7s dynamic','10s true','10s static','10s dynamic')

end

function [output] = getImpulseResponse()
%get data from PNAS paper

hold_figure = open('Art_vein_diamters.fig')
h = gcf; %current figure handle
axesObjs = get(h, 'Children');  %axes handles
dataObjs = get(axesObjs, 'Children');
xdata = get(dataObjs{2,1}, 'XData');  %data from low-level grahics objects
ydata = get(dataObjs{2,1}, 'YData');
zdata = get(dataObjs{2,1}, 'ZData');
close(hold_figure)

% plot PNAS data
ii = 15; %15th row is average HRF from a single puff
measured_HRF = ydata{ii,1}/max(ydata{ii,1});
        %plot(xdata{ii,1},measured_HRF_ROI,'Color',C,'LineWidth',3)

% fit HRF to data between 1-4 seconds after stimulus to avoid pre-stim
% undershoot (artifact) and post-stim undershoot (to see if the model can predict it)
start = find(xdata{ii,1}==1);
stop = find(xdata{ii,1}==4);

clear rsq
measured_HRF_ROI = measured_HRF(start:stop);
A = 6.5;
res = 0.01;
disp('optimizing vessel response kernel fit...')
for a1 = res:res:10;
    for b1 = res:res:10;
        
        x = linspace(1,4,length(measured_HRF_ROI));
        
        kernel_g = A.*((x).^(a1-1).*b1.^(a1).*exp(-b1.*(x))./gamma(a1));
        kernel = kernel_g/sum(kernel_g);
        fit = kernel/max(kernel);
        %figure(10), plot(x,fit)
        
        yresid = measured_HRF_ROI - fit;
        ssresid = sum(yresid.^2);
        sstotal = (length(measured_HRF_ROI-1)) * var(measured_HRF_ROI);
        rsq{round(a1/res),round(b1/res)} = 1 - ssresid/sstotal;
    end
end

[a b] = find(abs(cell2mat(rsq)-1)<= 0.008,1);

x = linspace(0,6,100);
a1 = res*a;
b1 = res*b;
kernel_g = A.*((x).^(a1-1).*b1.^(a1).*exp(-b1.*(x))./gamma(a1));
kernel = kernel_g/sum(kernel_g);
fit = kernel/max(kernel);

output.fit.kernel = kernel;
output.fit.time = x;
output.fit.kernel_norm = fit;
output.fit.rsq = round(rsq{a,b},4);

output.measured.onePuff.dilation = ydata{15,1};
output.measured.tenSecPuff.dilation = ydata{12,1};
output.measured.onePuff.time = xdata{15,1};
output.measured.tenSecPuff.time = xdata{12,1};
output.measured.kernel_norm = measured_HRF;
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