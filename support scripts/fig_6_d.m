function fig_6_d(data)
%% load/refine data
% [WGN] = get_mtspectrumFromWGN();
% [WGN] = get_WGN_neural(WGN);
% [WGN] = analyze_WGN_neural(WGN);

WGN = data.fig_6.d;

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

LineColor = {[1 0 0],[0 0 1]};

subplot(3,2,1), hold on
h1 = plot(WGN.dynamic.HRF_time,WGN.dynamic.HRF_avg,'r');
h2 = plot(WGN.static.HRF_time,WGN.static.HRF_avg,'b');
title('HRF kernel deconvolved from white noise NO production rate. input = 5 trials, 300s duration')
legend([h1 h2],'dynamic','static')
xlabel('time (s)')

subplot(3,2,2), hold on
dynamic_r2 = cell2mat(WGN.dynamic.r2);
static_r2 = cell2mat(WGN.static.r2);
bar(1,mean(dynamic_r2),0.5,'r')
scatter(1.*ones(1,5),dynamic_r2,'filled','MarkerFaceColor',[.502 0 0.125])
scatter(1,dynamic_r2(1),'filled','r')
scatter(1,dynamic_r2(1),'k')


bar(2,mean(static_r2),0.5,'b')
scatter(2.*ones(1,5),static_r2,'filled','MarkerFaceColor',[0 0 0.502])
scatter(2,static_r2(1),'filled','b')
scatter(2,static_r2(1),'k')

axis([0.5 2.5 0.9 1])
title('predictive value of HRF to 5 independent 300s trials')
xticklabels({'dynamic','static'})
xticks([1 2])
ylabel('R^{2}')

subplot(3,2,[3 4]), hold on
d1 = plot(WGN.time,WGN.dynamic.dilation{1},'Color',[0 0 0],'LineWidth',2);
d2 = plot(WGN.time,WGN.dynamic.prediction{1},'Color',[1 0 0],'LineWidth',2);
legend([d1 d2],'true','prediction from kernel')
title('dynamic RBC core example trial')
xlabel('time (s)'), ylabel('\Delta vessel diameter (%)')
axis([15 300 -2 3])

subplot(3,2,[5 6]), hold on
s1 = plot(WGN.time,WGN.static.dilation{1},'Color',[0 0 0],'LineWidth',2);
s2 = plot(WGN.time,WGN.static.prediction{1},'Color',[0 0 1],'LineWidth',2);
legend([s1 s2],'true','prediction from kernel')
title('static RBC core example trial')
xlabel('time (s)'), ylabel('\Delta vessel diameter (%)')
axis([15 300 -2 3])

end

function [output] = get_mtspectrumFromWGN()
dt = 0.01;
LineColor = {[1 0 0],[0 0 1]};
RBC_core = {'','_static'};

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

output.time = time;
output.dynamic.avg_f = mean(cell2mat(output.dynamic.f),2);
output.static.avg_f = mean(cell2mat(output.static.f),2)
end

function [output] = get_WGN_neural(input);
%import gamma band power (simulated with white gaussian noise)
output = input;
fs = 30;
t = input.time;
for ii = 1:10 %iterate across all WGN
    if ii == 10
        Gam{ii} = importdata(['GammaBandPower_110.csv'])';
        time = 1/fs:1/fs:length(Gam{ii})/fs;
        interp_Gam{ii} = interp1(time,Gam{ii}',t,'spline','extrap');
    else
        Gam{ii} = importdata(['GammaBandPower_10' num2str(ii) '.csv'])';
        time = 1/fs:1/fs:length(Gam{ii})/fs;
        interp_Gam{ii} = interp1(time,Gam{ii}',t,'spline','extrap');
    end
end

output.neural = interp_Gam;
end

function [output] = analyze_WGN_neural(input);

output = input;
%figure, hold on
Fs = round(1/(input.time(2)-input.time(1)));
kernel_length = 10*Fs; %get a 10s kernel estimate
for ii = 6:10
    disp(['calculating kernel for 5 min trial #' num2str(ii-5) ' out of 5...'])
    From = detrend(input.neural{ii});
    To = detrend(input.dynamic.dilation{ii});
    [HRF_dynamic{ii}] = deconvolveHRF(From,To,kernel_length);
    %HRF{ii} = ifft(fft(To)./fft(From)); %inferior but faster TF calculation
    %plot(input.time(1:length(HRF_dynamic{1}))-15,HRF_dynamic{ii},'Color',[1 0 0 1/5]),axis([0 6 -0.1 0.35])
    
    
    From = detrend(input.neural{ii});
    To = detrend(input.static.dilation{ii});
    [HRF_static{ii}] = deconvolveHRF(From,To,kernel_length);
    %TF_static{ii} = ifft(fft(To)./fft(From)); %inferior but faster TF calculation
    %plot(input.time(1:length(HRF_static{1}))-15,HRF_static{ii},'Color',[0 0 1 1/5]),axis([0 6 -0.1 0.35])
    
end
dynamic_avg = mean(cell2mat(cellfun(@(x) x',HRF_dynamic,'UniformOutput',0)'));
static_avg = mean(cell2mat(cellfun(@(x) x',HRF_static,'UniformOutput',0)'));


%apply TF to remaining 5 trials

%figure, hold on
HRF_dur = 10;
for ii = 1:5
    From = detrend(input.neural{ii});
    actual = detrend(input.dynamic.dilation{ii});
    prediction{ii} = conv(From,dynamic_avg(1:HRF_dur*Fs));
    prediction{ii} = prediction{ii}(1:length(actual)); %remove prediction beyond actual
    t = input.time-15;
    %plot(t,actual(1:length(t)),'k')
    %plot(t,prediction(1:length(t)),'m')
    
    r2_dynamic{ii} = CalculateRsquared(prediction{ii},actual);
end

% figure, hold on
HRF_dur = 10;
for ii = 1:5
    From = detrend(input.neural{ii});
    actual = detrend(input.static.dilation{ii});
    prediction{ii} = conv(From,static_avg(1:HRF_dur*Fs));
    prediction{ii} = prediction{ii}(1:length(actual)); %remove prediction beyond actual
%     t = input.time-15;
%     plot(t,actual(1:length(t)),'k')
%     plot(t,prediction(1:length(t)),'m')
    
    r2_static{ii} = CalculateRsquared(prediction{ii},actual);
end


output.static.HRF_time = [0:1/Fs:kernel_length/Fs];
output.static.r2 = r2_static;
output.static.HRF = HRF_static;
output.static.HRF_avg = static_avg;
output.static.prediction = prediction;

output.dynamic.HRF_time = [0:1/Fs:kernel_length/Fs];
output.dynamic.r2 = r2_dynamic;
output.dynamic.HRF = HRF_dynamic;
output.dynamic.HRF_avg = dynamic_avg;
output.dynamic.prediction = prediction;

output.HRF_dur = HRF_dur;
end

function [HRF] = deconvolveHRF(input,output,kernel_length)

K = kernel_length;
N = length(input);
Toeplitz = toeplitz([input zeros(1, K-1)], [input(1) zeros(1, K-1)]);
Toeplitz = [ones(size(Toeplitz,1), 1) Toeplitz];
output = [output zeros(1, size(Toeplitz, 1)-length(output))];
HRF = Toeplitz\[output]';

end

function r2 = CalculateRsquared(pred,act)
%   function r2 = CalculateRsquared(pred,act)
%
%   Author: Aaron Winder
%   Affiliation: Engineering Science and Mechanics, Penn State University
%   https://github.com/awinde
%
%   DESCRIPTION: Calculates the coefficient of determination between two
%   signals. Formula was obtained from 
%       http://en.wikipedia.org/wiki/Coefficient_of_determination
%   and is 1 - [(sum of squares of residuals)/(total sum of squares)]
%
%_______________________________________________________________
%   PARAMETERS:
%                   pred - [array] column vector of the model which 
%                   predicts the actual signal
%
%                   act - [array] the actual signal as a column
%                   vector
%_______________________________________________________________
%   RETURN:
%                   r2 - [double] the coefficient of determination
%_______________________________________________________________

% Error Variance
SSE = sum((act-pred).^2);

% Total Variance
SST = sum((act-(ones(size(act,1),1)*mean(act))).^2);

% Check that the sum of the residuals is small compared to SSE + SSR
r2 = ones(size(SSE))-SSE./SST;
end