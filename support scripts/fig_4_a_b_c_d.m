function fig_4_a_b_c_d()


%% Generate figure for illustrator: only plot data and refine later in the illustrator
monitorPos = get(0,'MonitorPositions'); % [x y w ht] X # of monitors
w = 1200; ht = 200; % Define window size: width(w) and height(ht)
mon = 1; % Only one monitor

% Determine positioning for main figure
pos = [((monitorPos(mon,3)-w)/2)+monitorPos(mon,1)/2,...
    ((monitorPos(mon,4)-ht)*16/20)+monitorPos(mon,2)/2,...
    w ht];

% turn off warning to prevent MATLAB throwing warning
warning('off','MATLAB:legend:PlotEmpty');
warning('off','MATLAB:legend:IgnoringExtraEntries');

figure('name','Figure 4 a,b,c,d:Coupling NO concentration to vessel size and NO degradation reproduces a post stimulus undershoot from a purely dilatory stimulus','NumberTitle','off',...
    'position',pos,'color','w');

subplot(1,4,1)
t = [0:0.01:15];
y = rectangularPulse(6,7,t)+1;
plot(t,y,'k','LineWidth',2)
xlim([5 11])
xlabel('time [s]')
ylabel('Fold increased production [NO]')
title('nNOS neuron NO production')
axis([5 15 0.5 2.5])
axis square

delay = 6;

    A = 6.5; %doesn't matter
    a1 = 8.91;
    b1 = 3.7;

    x = linspace(0,delay,length([-delay:0.01:0]));
    
    kernel_g = A.*((x).^(a1-1).*b1.^(a1).*exp(-b1.*(x))./gamma(a1));
    kernel = kernel_g/sum(kernel_g);

subplot(1,4,2)
plot([-6:0.01:0],[flip(kernel)],'k','LineWidth',3)
xlabel('time [s]')
ylabel('Weight')
title('NO in SM time kernel')
axis square
xlim([-6 0])

subplot(1,4,3)
past_c = 10.^[-3:0.01:5];
GC_activation = 100./((10^0.95./past_c).^0.8+1);
semilogx(past_c,GC_activation,'k','LineWidth',3);
xlabel('NO in SM [nM]')
ylabel('GC Activation [%]')
title('GC activation curve')
xlim([10^-3 10^4])
axis square

subplot(1,4,4), hold on
title('GC activation to vasodilation')
axis([-30 30 -30 30])
GC = [-30:0.1:30];
VDx1 = GC;
VDx2 = GC.*2;
VDx3 = GC.*3;
plot(GC,VDx1,'Color',[0 0 0 1/3],'LineWidth',3)
plot(GC,VDx2,'Color',[0 0 0 2/3],'LineWidth',3)
plot(GC,VDx3,'Color',[0 0 0 1],'LineWidth',3)
xlabel('GC activation [%]')
ylabel('\Deltavessel diameter [%]')
axis square

end


