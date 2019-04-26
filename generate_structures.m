function [data] = generate_structures(rootfolder)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global f;
%% Figure 2
f = waitbar(1/100,'Analyzing data for figure 2...');
data.fig_2 = fig2();
%% Figure 3
waitbar(2/100,f,'Analyzing data for figure 3...');
data.fig_3 = fig3();
%% Figure 4
waitbar(3/100,f,'Analyzing data for figure 4...');
data.fig_4 = fig4();
%% Figure 5
waitbar(5/100,f,'Analyzing data for figure 5...');
data.fig_5 = fig5();
%% Figure 6
waitbar(25/100,f,'Analyzing data for figure 6...');
data.fig_6 = fig6();

%% Figure 7
waitbar(80/100,f,{'Analyzing data for figure 7...','Multi-taper analysis'});
data.fig_7 = fig7();

%% completion notification
waitbar(6/6,f,'Analysis Complete');
pause(2)
close(f)
end

%% figure 2 functions
function [output] = fig2()

output.proximal.SM.GC_raw = importdata('Proximal_ParametricSweep_GCActivation.csv');
output.proximal.tissue.NO_raw = importdata('Proximal_ParametricSweep_PerivascularNO.txt');
output.proximal.tissue.dr = 1; %interpolate mesh to 1 um

output.regional.SM.GC_raw = importdata('Regional_ParametricSweep_GCActivation.csv');
output.regional.tissue.NO_raw = importdata('Regional_ParametricSweep_PerivascularNO.txt');
output.regional.tissue.dr = 1; %interpolate mesh to 1 um

output.uniform.SM.GC_raw = importdata('Uniform_ParametricSweep_GCActivation.csv');
output.uniform.tissue.NO_raw = importdata('Uniform_ParametricSweep_PerivascularNO.txt');
output.uniform.tissue.dr = 1; %interpolate mesh to 1 um

%% refine data

[output.proximal] = createGCandNOmatrix(output.proximal);
[output.regional] = createGCandNOmatrix(output.regional);
[output.uniform] = createGCandNOmatrix(output.uniform);

end

function [output] = createGCandNOmatrix(input)

%make GC activation in SM matrix
vessel_size = round(unique(input.SM.GC_raw.data(:,5)).*10^6,2); %um
NO_prod = round(unique(input.SM.GC_raw.data(:,6)),2); %10^NO_prod (M/s)
GC_raw = input.SM.GC_raw.data(:,7); % percent GC activation in SM

m = length(NO_prod);
n = length(vessel_size);

index = [1:m:length(GC_raw) length(GC_raw)+1];
for ii = 1:length(index)-1
    GC(:,ii) = flip(GC_raw(index(ii):index(ii+1)-1));
end

output.SM.GC_raw = input.SM.GC_raw;
output.SM.GC = GC;
output.NO_prod = NO_prod;
output.vessel_size = vessel_size;
output.dimensions = [m n];

%make perivascular NO concentration matrix

[index] = [find(input.tissue.NO_raw(:,1)==0)' length(input.tissue.NO_raw(:,1))];
radius = [0:input.tissue.dr:100];

for jj = 1:n %vessel sizes
for ii = 1:m %NO productions
    radius_raw{ii,jj} = input.tissue.NO_raw(index(ii+(jj-1)*m):index(ii+(jj-1)*m+1)-1,1);
    conc_raw{ii,jj} = input.tissue.NO_raw(index(ii+(jj-1)*m):index(ii+(jj-1)*m+1)-1,2);
    
    conc{ii,jj} = interp1(radius_raw{ii,jj},conc_raw{ii,jj},radius); %interpolate concentration with dr
end
end

output.tissue.dr = input.tissue.dr;
output.tissue.NO_raw = input.tissue.NO_raw;
output.tissue.NO.radius_raw = radius_raw;
output.tissue.NO.conc_raw = conc_raw;
output.tissue.NO.conc = conc;
output.tissue.NO.radius = radius;
end

%% figure 3 functions
function [output] = fig3()
proximal.SM.GC_raw = importdata('Proximal_ParametricSweep_GCActivation.csv');
proximal.tissue.NO_raw = importdata('Proximal_ParametricSweep_PerivascularNO.txt');
proximal.tissue.dr = 0.1; %interpolate mesh to 1 um

regional.SM.GC_raw = importdata('Regional_ParametricSweep_GCActivation.csv');
regional.tissue.NO_raw = importdata('Regional_ParametricSweep_PerivascularNO.txt');
regional.tissue.dr = 0.1; %interpolate mesh to 1 um

uniform.SM.GC_raw = importdata('Uniform_ParametricSweep_GCActivation.csv');
uniform.tissue.NO_raw = importdata('Uniform_ParametricSweep_PerivascularNO.txt');
uniform.tissue.dr = 0.1; %interpolate mesh to 1 um

[proximal] = createGCandNOmatrix(proximal);
[regional] = createGCandNOmatrix(regional);
[uniform] = createGCandNOmatrix(uniform);

[O2_mmHg] = getO2(proximal.vessel_size, proximal.tissue.dr);

proximal = calculate_CcO_inhibition(proximal,O2_mmHg);
regional = calculate_CcO_inhibition(regional,O2_mmHg);
uniform = calculate_CcO_inhibition(uniform,O2_mmHg);

CcO_limit = 12; % below 12% CcO activity is considered non-physiological 
[proximal] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(proximal,CcO_limit);
[regional] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(regional,CcO_limit);
[uniform] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(uniform,CcO_limit);

[proximal] = normalize_to_percent_dilation(proximal);
[regional] = normalize_to_percent_dilation(regional);
[uniform] = normalize_to_percent_dilation(uniform);

output.proximal = proximal;
output.regional = regional;
output.uniform = uniform;
output.O2_mmHg = O2_mmHg;
end

function [O2] = getO2(vessel_size, dr)
%perivascular oxygen profile fitted to Devor 2011, with an arterial O2 of
%35 mmHg from Lyons 2016
O2_in_mmHg = 35;
p1 = 0.03052; %[1/um]
p2 = 0.0002471; %[1/um]
r = [0:dr:100];
O2_hold = (O2_in_mmHg/(53.44+14.13))*(53.44*exp(-p1*(r))+14.13*exp(p2*(r)));
for ii = 1:length(vessel_size)
    O2{ii} = [ones(1,length(0:dr:vessel_size(ii))).*O2_in_mmHg O2_hold];
    O2{ii} = O2{ii}(1:length(r));
end
end
function [output] = calculate_CcO_inhibition(input,O2_mmHg)

NO = input.tissue.NO.conc; %grab NO concentration in the tissue

for m = 1:input.dimensions(1)
    for n = 1:input.dimensions(2)
        
        O2_uM{n} = O2_mmHg{n}./(760.*0.21).*215.919; %convert O2 from mmHg to uM
        Ko2 = 0.21; %uM
        Kno = 0.225; %nM
        CcO{m,n} = O2_uM{n}./(O2_uM{n}+Ko2.*(1+(NO{m,n})./Kno)); %calculate CcO inhibition (competitive inhibition)
        
    end
end

output = input; %carry variables from input to output
output.tissue.CcO = CcO; %add CcO inhibition to output
end
function [output] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(input,CcO_limit)

    for m = 1:input.dimensions(1)
        for n = 1:input.dimensions(2)
            toxic_region = length(find(input.tissue.CcO{m,n}<=CcO_limit/100));
            total_region = length(input.tissue.CcO{m,n});
            toxic_fraction(m,n) = toxic_region/total_region;
        end
    end

    output = input;
    output.tissue.CcO_fraction = flip(toxic_fraction);
end
function [output] = normalize_to_percent_dilation(input)

GC = input.SM.GC;
CcO_fraction = input.tissue.CcO_fraction;

xi = [5:5:95];

    for n = 1:input.dimensions(2)
        x = GC(:,n);
        y = CcO_fraction(:,n);
        yi = interp1(x,y,xi,'linear');
        CcO_norm2dilation(:,n) = yi;
    end
output = input;
output.tissue.CcO_norm2dilation.CcO_fraction = CcO_norm2dilation;
output.tissue.CcO_norm2dilation.dilation = xi;
end

%% figure 4 functions
function [output] = fig4()

NO_sens = {1, 2, 3};
stim_duration = {'0p5','1','2','4'};
ii_end = length(NO_sens);
time = [0:0.005:15];

for ii = 1:ii_end
    hold_data = importdata(['2NO_1NOd_' num2str(NO_sens{ii}) 'NOsens_50bGC_0Hz_10VD_6sNOkernel.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    conc{ii} = interp1(hold_data(:,1),hold_data(:,2),time')'; %concentration
    dilation{ii} = interp1(hold_data(:,1),hold_data(:,3),time')'; %dilation
end

output.conc = conc;
output.dilation = dilation;
output.time = time;

end

%% figure 5 functions
function [output] = fig5()

[model] = getImpulseResponse();
[model] = getCOMSOLsimulation(model);

output.model = model;
end

function [output] = getImpulseResponse()
%get data from PNAS paper
global f;

hold_figure = open('Art_vein_diamters.fig');
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
%disp('optimizing vessel response kernel fit...')
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
    waitbar(5/100+20/100*a1/10,f,{'Analyzing data for figure 5...' 'optimizing kernel'});
end


%%

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

%% figure 6 functions
function [output] = fig6()
global f
[output.a] = fig6a();
[output.b] = calculateAmplitude();
[output.c] = get_mtspectrumFromWGN();

waitbar(32.5/100,f,'Analyzing data for figure 6...');

[WGN] = get_WGN_neural(output.c);
[WGN] = analyze_WGN_neural(WGN);
[output.d] = WGN;

[output.e] = get_OscillationsWithVaryingBaselineGC();
end

function [output] = fig6a()

NO_sens = [1 3 5 10];
time_end = [20 20 30 50];
dt = 0.01;

for ii = 1:length(NO_sens)
    time{ii} = [0.05:dt:time_end(ii)];
    hold_data = importdata(['1.2NO_1NOd_' num2str(NO_sens(ii)) 'NOsens_50bGC_0Hz_10VD_6sNOkernel_Oscillations.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    dynamic.conc{ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')'; %concentration
    dynamic.dilation{ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')'; %dilation
    
    hold_data = importdata(['1.2NO_1NOd_' num2str(NO_sens(ii)) 'NOsens_50bGC_0Hz_10VD_6sNOkernel_OscillationsStatic.csv']);
    [a index] = unique(hold_data(:,1));
    hold_data = hold_data(index,:); %don't take replicate values
    static.conc{ii} = interp1(hold_data(:,1),hold_data(:,2),time{ii}')'; %concentration
    static.dilation{ii} = interp1(hold_data(:,1),hold_data(:,3),time{ii}')'; %dilation
end

output.dynamic = dynamic; 
output.static = static;
output.time = time;
end
function [output] = calculateAmplitude()
Hz = [0.05 0.1 0.15 0.2 0.25 0.5 0.75 1];
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
    
    output.Hz = Hz;
    output.Hz_pchip = Hz_smooth;
end

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
output.static.avg_f = mean(cell2mat(output.static.f),2);
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
global f
output = input;
%figure, hold on
Fs = round(1/(input.time(2)-input.time(1)));
kernel_length = 10*Fs; %get a 10s kernel estimate
for ii = 6:10
    waitbar(25/100+50/100*(ii-5)/6,f,{'Analyzing data for figure 6...',['calculating kernel for 5 min trial ' num2str(ii-5) '/5']});
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

%% figure 7 functions
function [output] = fig7()
% unless specfically stated otherwise: slope, m = 3; plasma free Hemoglobin = 1 uM; Hematocrit = 45%; dynamically moving RBC core

f_name_label = {
    'm_1',
    'm_3',
    'm_5',
    'Hgb_20',
    'Hgb_40',
    'Hct20',
    'Hct60',
    'static',
    'static_m_1'};

f_name = {
    '5NO_0NOd_1NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_45_Hct_RawGamma_dynamic_v24.csv', %slope, m = 1
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_45_Hct_RawGamma_dynamic_v24.csv', %slope, m = 3
    '5NO_0NOd_5NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_45_Hct_RawGamma_dynamic_v24.csv', %slope, m = 5
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith20Hgb_45_Hct_RawGamma_dynamic_v24.csv', %plasma free Hemoglobin = 20 uM
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith40Hgb_45_Hct_RawGamma_dynamic_v24.csv', %plasma free Hemoglobin = 40 uM
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_20_Hct_RawGamma_dynamic_v24.csv', %Hematocrit = 20 uM
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_60_Hct_RawGamma_dynamic_v24.csv', %Hematocrit = 60 uM
    '5NO_0NOd_3NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_45_Hct_RawGamma_static_v24.csv', %static red blood cell (RBC) core, m = 3
    '5NO_0NOd_1NOsens_50bGC_1500Hz_20VD_6sNOkernel_GammaWith1Hgb_45_Hct_RawGamma_dynamic_v25.csv'}; %static red blood cell (RBC) core, m = 1
tapers = [101 201]; %time-bandwidth product determines number of tapers. 
% for a 1500s trial and a resolution of 0.1 Hz a maximum of 150 tapers are available for use.
% using 101 instead for a smaller degree of frequency smearing.

for ii = 1:length(f_name)
    [output.(f_name_label{ii})] = ExtractVasodynamics_v24(f_name{ii}, tapers);
end

% get perivascular NO concentrations
f_name_2 = {
    'perivascular_Hgb_1_Hct_45.csv',
    'perivascular_Hgb_20_Hct_45.csv',
    'perivascular_Hgb_40_Hct_45.csv',
    'perivascular_Hgb_1_Hct_20.csv',
    'perivascular_Hgb_1_Hct_60.csv'};
kk = 1;
for ii = [2 4 5 6 7]
    [output.(f_name_label{ii}).perivascular] = importdata(f_name_2{kk})
    kk = kk + 1;
end

% get (%) change in baseline vessel diameter from Hgb and Hct changes
hold_hgb = importdata('baseline_Hgb.csv');
GC_activation = 100./((10^0.95./hold_hgb.data(:,12)).^0.8+1);
GC_activation_baseline = 100./((10^0.95./[hold_hgb.data(2,12)]).^0.8+1);;%normalize to 1 uM pfHgb
dilation = GC_activation - GC_activation_baseline;
%[data.fig_7.baseline.Hgb] = [hold_hgb.data(:,8) dilation hold_hgb.data(:,12)];
[output.baseline.Hgb] = [hold_hgb.data(:,8) dilation hold_hgb.data(:,12)];

hold_hct = importdata('baseline_Hct_v3.csv');
GC_activation = 100./((10^0.95./[hold_hct.data(:,11)]).^0.8+1);
GC_activation_baseline = 100./((10^0.95./[hold_hct.data(10,11)]).^0.8+1); %normalize to 45% Hct
dilation = GC_activation - GC_activation_baseline;
%[data.fig_7.baseline.Hct] = [hold_hct.data(:,9) dilation hold_hct.data(:,11)];
[output.baseline.Hct] = [hold_hct.data(:,9) dilation hold_hct.data(:,11)];

% get (%) change in baseline vessel diameter with feedback from a changing vessel
hold_hgb = importdata('baseline_Hgb_feedback.csv');
[time a] = unique(hold_hgb(:,1));
dilation = hold_hgb(a,3)-30;
Hgb = [0 1 5:5:40];
%[data.fig_7.baseline.Hgb_feedback] = [Hgb' interp1(time,dilation,[27:15:170])'];
[output.baseline.Hgb_feedback] = [Hgb' interp1(time,dilation,[27:15:170])'];

hold_hct = importdata('baseline_Hct_feedback.csv');
[time a] = unique(hold_hct(:,1));
dilation = hold_hct(a,3)-5;
Hct = [0.45 0.6 0.75 0.9 1 0.3 0.15 0];
Hct = [0 0.15 0.3 0.45 0.6 0.75 0.9 1];
[160 140 120 20 40 60 80 100]-5
%[data.fig_7.baseline.Hct_feedback] = [Hct' interp1(time,dilation,[160 140 120 20 40 60 80 100]-5)'];
[output.baseline.Hct_feedback] = [Hct' interp1(time,dilation,[160 140 120 20 40 60 80 100]-5)'];

end

function [output] = ExtractVasodynamics_v24(f_name, tapers)
% this subfunction takes a csv file containing time, NO in smooth muscle,
% and vasodilation for a 1500s long trial and calculates:
%   variance of dilation
%   NO production spectrum Power (from WGN_1500.csv)
%   vasodynamics spectrum Power
%   hemodynamic response function (calculated from Toeplitz)
%
% also reports:
%   NO concentration
%   vasodilation
%   NO production rate (white gaussian noise low pass filtered below 2Hz)

dt = 1/6; %time step resolution
time = [15:dt:1500]; %trials are 1500 seconds long - ignore first 15 seconds

params.Fs = 1/dt;
params.tapers = tapers;
params.fpass = [0.025 3];
params.err = [2 0.05]; %use jack-knife resampling confidence intervals p = 0.05 
% (could also try multitaper frequency domain bootstrapping (MFDB) for less noisy CI - but jackknife used by so many people in NVC)

file_number = [1500];

hold_data = importdata(f_name);
[a index] = unique(hold_data(:,1));
hold_data = hold_data(index,:); %don't take replicate values
ii = 1; %only 1 trail - don't need to average across multiple shorter simulations

% get concentrations and dilations predicted by the model_______________________________________________________________
output.concentration(ii,:) = detrend(interp1(hold_data(:,1),hold_data(:,2),time')'); %concentration
output.dilation(ii,:) = detrend(interp1(hold_data(:,1),hold_data(:,3),time')'); %dilation

% get the variance of the vessel dynamics_______________________________________________________________________________
output.variance(ii,:) = var(output.dilation(ii,:));

% get the power spectrum of the vessel dynamics_________________________________________________________________________
[output.spectrumPower(ii,:), output.spectrumHz(ii,:), output.Serr{ii}] = mtspectrumc(output.dilation(ii,:),params);

% deconvolve the hemodynamics response function from the white
% gaussian noise input and the vasodilation output from COMSOL__________________________________________________________

fs = 30; %get white gaussian noise NO production given to COMSOL
Gam = importdata(['WGN_' num2str(file_number(ii)) '.csv'])';

time_Gam = 1/fs:1/fs:length(Gam)/fs;
output.NO_production(ii,:) = interp1(time_Gam,Gam',time,'spline','extrap');
[output.NOspectrumPower(ii,:), output.NOspectrumHz(ii,:)] = mtspectrumc(output.NO_production(ii,:),params);


kernel_length = 15/dt; %get a 15s kernel estimate
From = detrend(output.NO_production(ii,:)); %remove DC component of NO production
To = detrend(output.dilation(ii,:)); %remove DC component of vasodilation
[output.HRF(ii,:)] = deconvolveHRF(From,To,kernel_length); %calculate HRF with Toeplitz matrix
end
%%





















