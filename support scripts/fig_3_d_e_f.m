function fig_3_d_e_f(data)
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
% 
% 
% CcO_limit = 12; % 12%
% [proximal] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(proximal,CcO_limit);
% [regional] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(regional,CcO_limit);
% [uniform] = calculate_fraction_of_tissue_with_toxic_CcO_inhibition(uniform,CcO_limit);
% 
% [proximal] = normalize_to_percent_dilation(proximal);
% [regional] = normalize_to_percent_dilation(regional);
% [uniform] = normalize_to_percent_dilation(uniform);

proximal = data.fig_3.proximal;
regional = data.fig_3.regional;
uniform = data.fig_3.uniform;

%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 900; ht = 200; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)/2)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Figure 3 d,e,f: Perivascular CcO inhibition','NumberTitle','off',...
    'position',pos,'color','w');
subplot(131), hold on
pt1 = imagesc(proximal.vessel_size.*2,proximal.NO_prod,flip(proximal.tissue.CcO_fraction));
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('proximal: toxic fraction')
axis square, caxis([0 1])
hold on
contour(proximal.vessel_size.*2,proximal.NO_prod,flip(proximal.SM.GC),[10 10 50 50 90 90],'k','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})

subplot(132), hold on
pt1 = imagesc(regional.vessel_size.*2,proximal.NO_prod,flip(regional.tissue.CcO_fraction));
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('regional: toxic fraction')
axis square, caxis([0 1])
hold on
contour(regional.vessel_size.*2,regional.NO_prod,flip(regional.SM.GC),[10 10 50 50 90 90],'k','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})


subplot(133), hold on
pt1 = imagesc(uniform.vessel_size.*2,proximal.NO_prod,flip(uniform.tissue.CcO_fraction));
ylabel('NO production (M/s)'),xlabel('vessel diameter (\mum)'),title('uniform: toxic fraction')
axis square, caxis([0 1])
hold on
contour(uniform.vessel_size.*2,proximal.NO_prod,flip(uniform.SM.GC),[10 10 50 50 90 90],'k','ShowText','on')
yticklabels({'10^{1}','10^{0}','10^{-1}','10^{-2}','10^{-3}','10^{-4}','10^{-5}','10^{-6}'})
colormap(flip(hot))
colorbar
%%
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 900; ht = 200; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)/20)+monitorPos(mon,2)/2,...
    w ht];

figure('name','Figure 3 g,h,i: Perivascular CcO inhibition','NumberTitle','off',...
    'position',pos,'color','w');
subplot(131), hold on
pt4 = imagesc(proximal.vessel_size.*2,proximal.tissue.CcO_norm2dilation.dilation,proximal.tissue.CcO_norm2dilation.CcO_fraction);
ylabel('GC activation (%)'),xlabel('vessel diameter (\mum)'),title('proximal: toxic fraction')
axis square, caxis([0 1]),axis([10 50 5 95])
yticks([10:10:90])
hold on

subplot(132)
pt5 = imagesc(regional.vessel_size.*2,regional.tissue.CcO_norm2dilation.dilation,flip(regional.tissue.CcO_norm2dilation.CcO_fraction));
ylabel('GC activation (%)'),xlabel('vessel diameter (\mum)'),title('regional: toxic fraction')
axis square, caxis([0 1]),axis([10 50 5 95])
yticklabels(flip({'10','20','30','40','50','60','70','80','90'}))
yticks([10:10:90])
hold on


subplot(133), hold on
pt6 = imagesc(uniform.vessel_size.*2,uniform.tissue.CcO_norm2dilation.dilation,uniform.tissue.CcO_norm2dilation.CcO_fraction);
ylabel('GC activation (%)'),xlabel('vessel diameter (\mum)'),title('uniform: toxic fraction')
axis square, caxis([0 1]),axis([10 50 5 95])
yticklabels({'10','20','30','40','50','60','70','80','90'})
yticks([10:10:90])
hold on
colormap(flip(flip(fliplr(hot))))
colorbar
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
        yi = interp1(x,y,xi,'linear','extrap');
        CcO_norm2dilation(:,n) = yi;
    end
output = input;
output.tissue.CcO_norm2dilation.CcO_fraction = CcO_norm2dilation;
output.tissue.CcO_norm2dilation.dilation = xi;
end







