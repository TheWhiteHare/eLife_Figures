function fig_3_a_b_c(data)
%% load data

% proximal.SM.GC_raw = importdata('Proximal_ParametricSweep_GCActivation.csv');
% proximal.tissue.NO_raw = importdata('Proximal_ParametricSweep_PerivascularNO.txt');
% proximal.tissue.dr = 0.1; %interpolate mesh to 1 um
% 
% regional.SM.GC_raw = importdata('Regional_ParametricSweep_GCActivation.csv');
% regional.tissue.NO_raw = importdata('Regional_ParametricSweep_PerivascularNO.txt');
% regional.tissue.dr = 0.1; %interpolate mesh to 1 um
% 
% uniform.SM.GC_raw = importdata('Uniform_ParametricSweep_GCActivation.csv');
% uniform.tissue.NO_raw = importdata('Uniform_ParametricSweep_PerivascularNO.txt');
% uniform.tissue.dr = 0.1; %interpolate mesh to 1 um
% 
% %% refine data
% 
% [proximal] = createGCandNOmatrix(proximal);
% [regional] = createGCandNOmatrix(regional);
% [uniform] = createGCandNOmatrix(uniform);
% 
% [O2_mmHg] = getO2(proximal.vessel_size, proximal.tissue.dr);
% 
% proximal = calculate_CcO_inhibition(proximal,O2_mmHg);
% regional = calculate_CcO_inhibition(regional,O2_mmHg);
% uniform = calculate_CcO_inhibition(uniform,O2_mmHg);

proximal = data.fig_3.proximal;
regional = data.fig_3.regional;
uniform = data.fig_3.uniform;
O2_mmHg = data.fig_3.O2_mmHg;

example_production = 35;
yheight = max(proximal.tissue.NO.conc{example_production,1});
yheight = round(yheight*1.2, 2, 'significant');

%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 600; ht = 200; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*16/20)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Figure 3 a,b,c Perivascular CcO inhibition','NumberTitle','off',...
    'position',pos,'color','w');

model_index = [1 9];

subplot(131), hold on
for ii = 1:length(model_index)
yyaxis right, axis([0 100 0 yheight]), title('NO in the smooth muscle [nM]')
h1{model_index(ii)} = plot(proximal.tissue.NO.radius,proximal.tissue.NO.conc{example_production,model_index(ii)},'k-','Color',[0 0 0 ii/length(model_index)]);
yyaxis left, axis([0 100 0 100])
plot(proximal.tissue.NO.radius,proximal.tissue.CcO{example_production,ii}.*100,'r-','Color',[1 0 0 ii/length(model_index)])
plot(proximal.tissue.NO.radius,O2_mmHg{model_index(ii)},'b-','Color',[0 0 1 ii/length(model_index)])
end
ylabel('CcO activity (%)'),xlabel('radius (\mum)'),title('proximal')
legend([h1{1} h1{9}],{[num2str(round(proximal.SM.GC(72-example_production,1))) '% GC'],[num2str(round(proximal.SM.GC(72-example_production,9))) '% GC']})
axis square, %axis([0 100 0 100])
hold on

subplot(132), hold on
for ii = 1:length(model_index)
yyaxis right, axis([0 100 0 yheight]), title('NO in the smooth muscle [nM]')
h2{model_index(ii)} = plot(regional.tissue.NO.radius,regional.tissue.NO.conc{example_production,model_index(ii)},'k-','Color',[0 0 0 ii/length(model_index)]);
yyaxis left, axis([0 100 0 100])
plot(regional.tissue.NO.radius,regional.tissue.CcO{example_production,model_index(ii)}.*100,'r-','Color',[1 0 0 ii/length(model_index)])
plot(regional.tissue.NO.radius,O2_mmHg{model_index(ii)},'b-','Color',[0 0 1 ii/length(model_index)])
end
ylabel('CcO activity (%)'),xlabel('radius (\mum)'),title('regional')
legend([h2{1} h2{9}],{[num2str(round(regional.SM.GC(72-example_production,1))) '% GC'],[num2str(round(regional.SM.GC(72-example_production,9))) '% GC']})
axis square, %axis([0 100 0 100])
hold on

subplot(133), hold on
for ii = 1:length(model_index)
yyaxis right, axis([0 100 0 yheight]), title('NO in the smooth muscle [nM]')
h3{model_index(ii)} = plot(uniform.tissue.NO.radius,uniform.tissue.NO.conc{example_production,model_index(ii)},'k-','Color',[0 0 0 ii/length(model_index)]);
yyaxis left, axis([0 100 0 100])
plot(uniform.tissue.NO.radius,uniform.tissue.CcO{example_production,model_index(ii)}.*100,'r-','Color',[1 0 0 ii/length(model_index)])
plot(uniform.tissue.NO.radius,O2_mmHg{model_index(ii)},'b-','Color',[0 0 1 ii/length(model_index)])
end
ylabel('CcO activity (%)'),xlabel('radius (\mum)'),title('uniform')
legend([h3{1} h3{9}],{[num2str(round(uniform.SM.GC(72-example_production,1))) '% GC'],[num2str(round(uniform.SM.GC(72-example_production,9))) '% GC']})
axis square, %axis([0 100 0 100])
hold on


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




