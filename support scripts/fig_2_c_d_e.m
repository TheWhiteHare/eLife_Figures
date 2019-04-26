function fig_2_c_d_e(data)
%% load data

% proximal.SM.GC_raw = importdata('Proximal_ParametricSweep_GCActivation.csv');
% proximal.tissue.NO_raw = importdata('Proximal_ParametricSweep_PerivascularNO.txt');
% proximal.tissue.dr = 1; %interpolate mesh to 1 um
% 
% regional.SM.GC_raw = importdata('Regional_ParametricSweep_GCActivation.csv');
% regional.tissue.NO_raw = importdata('Regional_ParametricSweep_PerivascularNO.txt');
% regional.tissue.dr = 1; %interpolate mesh to 1 um
% 
% uniform.SM.GC_raw = importdata('Uniform_ParametricSweep_GCActivation.csv');
% uniform.tissue.NO_raw = importdata('Uniform_ParametricSweep_PerivascularNO.txt');
% uniform.tissue.dr = 1; %interpolate mesh to 1 um
% 
% %% refine data
% 
% [proximal] = createGCandNOmatrix(proximal);
% [regional] = createGCandNOmatrix(regional);
% [uniform] = createGCandNOmatrix(uniform);

proximal = data.fig_2.proximal;
regional = data.fig_2.regional;
uniform = data.fig_2.uniform;

%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 600; ht = 200; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [(0)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*9/20)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Figure 2 c,d,e:Predicting perivascular NO concentrations and smooth muscle signaling','NumberTitle','off',...
    'position',pos,'color','w');
subplot(131)
pt1 = imagesc(proximal.vessel_size.*2,proximal.NO_prod,proximal.SM.GC);
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('proximal')
axis square, caxis([0 100])
hold on
contour(uniform.vessel_size.*2,uniform.NO_prod,proximal.SM.GC,[10 10 50 50 90 90],'w','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})

subplot(132)
pt2 = imagesc(regional.vessel_size.*2,regional.NO_prod,regional.SM.GC);
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('regional')
axis square, caxis([0 100])
hold on
contour(uniform.vessel_size.*2,uniform.NO_prod,regional.SM.GC,[10 10 50 50 90 90],'w','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})

subplot(133)
pt3 = imagesc(uniform.vessel_size.*2,uniform.NO_prod,uniform.SM.GC);
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('uniform')
axis square, caxis([0 100]), %colorbar
colormap(summer)
hold on
contour(uniform.vessel_size.*2,uniform.NO_prod,uniform.SM.GC,[10 10 50 50 90 90],'w','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})
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














