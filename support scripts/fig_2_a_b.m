function fig_2_a_b(data)
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
w = 600; ht = 300; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [(0)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*9/10)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

index50 = [43 46 47];

figure('name','Figure 2 a,b: Perivascular NO concentrations','NumberTitle','off',...
    'position',pos,'color','w');
subplot(121), hold on
pt1 = plot(proximal.tissue.NO.radius,proximal.tissue.NO.conc{index50(1),3},'k','LineWidth',3);
pt2 = plot(regional.tissue.NO.radius,regional.tissue.NO.conc{index50(2),3},'r','LineWidth',3);
pt3 = plot(uniform.tissue.NO.radius,uniform.tissue.NO.conc{index50(3),3},'m','LineWidth',3);
ylabel('NO concentration (nM)'),xlabel('radius (\mum)'),title('perivascular NO')
legend([pt3 pt2 pt1],{'uniform','regional','proximal'},'Location','southeast')
axis square

subplot(122)
c = 10.^[-2:0.1:4];
GC_activation = 100./((10^0.95./c).^0.8+1);
semilogx(c,GC_activation,'k','LineWidth',3)
ylabel('NO concentration (nM)'),xlabel('GC activation (%)'),title('NO signaling of GC')
axis([10^-2 10^4 0 100]), axis square

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














