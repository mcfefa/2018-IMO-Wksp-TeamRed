close all
clear all


%%% Initial values of state variables
% Escape through highly proliferating immunogenic cells, relapse
TV = 1.4 * 10e10;
HI0 = 0.5 * TV;
HN0 = 0.0 * TV;
LI0 = 0.25 * TV;
LN0 = 0.25 * TV;

% Escape through resistance to immune system
% HN0 = 0.01 * TV;
% LI0 = 0.245 * TV;
% LN0 = 0.245 * TV;

% Cured
% HN0 = 0.0 * TV;
% LI0 = 0.3 * TV;
% LN0 = 0.2 * TV;

A0 = 0;
S0 = 0;

x0 = [HI0, HN0, LI0, LN0, A0, S0];


% define schedule of treatment

BCGtimes = xlsread('IL10data.xlsx','Responders','A1:A27');
BCGtimes = 7.*(BCGtimes - 1); 
BCGtimes2 = xlsread('IL10data.xlsx','Responders','G1:G27');
BCGtimes2 = 7.*(BCGtimes2 - 1); 

schedule_BCG(:,1) = BCGtimes';
schedule_BCG(:,2) = ones(size(BCGtimes'));
BCGdose = 4e6; % each treatment (usually 2.2e6 to 6.4e6 each week)
schedule_BCG(:,2) = BCGdose.*schedule_BCG(:,2);

schedule_BCG2(:,1) = BCGtimes2';
schedule_BCG2(:,2) = ones(size(BCGtimes2'));
%BCGdose = 4e6; % each treatment (usually 2.2e6 to 6.4e6 each week)
schedule_BCG2(:,2) = BCGdose.*schedule_BCG2(:,2);


% define growth parameters
rHI = 0.4; % PI producers
rHN = 0.4; % X2 free riders
rLI = rHI * 0.1; % Xr resistant
rLN = rHI * 0.1;
K = 10e13; % carrying capacity of tumor

delta0 = 2;
rhoA = 1;
rhoS = 1;
lambdaA = 0.5;
lambdaS = 0.5;
beta = 1;
eta = 1e3;
gamma = 0*0.02;

% BCG half life
BCGhl = 0.3;

par = [rHI, rHN, rLI, rLN, K, delta0, rhoA, rhoS, lambdaA, lambdaS, beta, eta, gamma];

tmax = 300;
tspan = 0:0.1:tmax;
tscat = [1 21 42]'; 

BCG = zeros(size(tspan));
for i = 1:length(tspan)-1
    tnow = tspan(i);
    for j = 1:length(schedule_BCG(:,1))
        if tnow == schedule_BCG(j,1)
            %fprintf('\n \t At %g day (when BCG at %1.1e), injected %1.1e BCG. \n', tnow, BCG(i), schedule_BCG(j,2))
            BCG(i+1) = BCG(i) + schedule_BCG(j,2);
            BCG0 = BCG(i+1);
            t0 = tnow;
            %fprintf('\t Now there is %1.1e BCG in system after treatment. \n', BCG0)
        end
    end
    BCG(i+1) = BCG0*exp(-(tnow-t0)*BCGhl);
end
    
clear BCG0 t0 tnow
BCG2 = zeros(size(tspan));
for i = 1:length(tspan)-1
    tnow = tspan(i);
    for j = 1:length(schedule_BCG2(:,1))
        if tnow == schedule_BCG2(j,1)
            %fprintf('\n \t At %g day (when BCG at %1.1e), injected %1.1e BCG. \n', tnow, BCG(i), schedule_BCG(j,2))
            BCG2(i+1) = BCG2(i) + schedule_BCG2(j,2);
            BCG0 = BCG2(i+1);
            t0 = tnow;
            %fprintf('\t Now there is %1.1e BCG in system after treatment. \n', BCG0)
        end
    end
    BCG2(i+1) = BCG0*exp(-(tnow-t0)*BCGhl);
end
    



TTimes = [];
for i = 1:45
    newTT = [(i-1)*28, (i-1)*28+7, (i-1)*28+14];
    TTimes = [TTimes, newTT];
end

TT = zeros(size(tspan));
treatment00 = [BCG; TT];
treatment20 = [BCG2; TT];
TTstart = 200;
TTdose = 60;
TThl = 1.59;
clear tnow t0;
TT0 = 0;
t0 = 1;
for i = 1:length(tspan)-1
    tnow = tspan(i);
    for j = 1:length(TTimes)
        if tnow == TTimes(j) + TTstart
            %fprintf('\n \t At %g day (when BCG at %1.1e), injected %1.1e BCG. \n', tnow, BCG(i), schedule_BCG(j,2))
            TT(i+1) = TT(i) + TTdose;
            TT0 = TT(i+1);
            t0 = tnow;
            %fprintf('\t Now there is %1.1e BCG in system after treatment. \n', BCG0)
        end
    end
    TT(i+1) = TT0*exp(-(tnow-t0)*TThl);
end



%for i = 1:length(TT)
 %   TT(i) = TT(i) + TTstart;
%end

treatment = [BCG; TT];

figure;
plot(tspan, BCG); hold on;
plot(tspan, 10000.*TT,'r')




IL10_no = zeros(size(BCGtimes));
IL10_no(1:15) = xlsread('IL10data.xlsx','Non-Responders','B1:B15');
IL10_yay = xlsread('IL10data.xlsx','Responders','B1:B27');

sol0 = ode23( @bladder_ODE_ana, tspan, x0, [], par, treatment00 );
sol = ode23( @bladder_ODE_ana, tspan, x0, [], par, treatment );
Ys = deval(sol, tspan)';
Ys0 = deval(sol0, tspan)';
HIsol = Ys(:,1);
HNsol = Ys(:,2);
LIsol = Ys(:,3);
LNsol = Ys(:,4);
Asol = Ys(:,5);
Ssol = Ys(:,6);
totalcells = HIsol + HNsol + LIsol + LNsol;
totalcells0 = Ys0(:,1) + Ys0(:,2) + Ys0(:,3) + Ys0(:,4);
%%
figure;
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])
subplot(3,1,1);
plot(tspan, totalcells0, 'k','LineWidth', 1.5); hold on;
plot(tspan, HIsol, 'b','LineWidth', 2.5); hold on;
plot(tspan, HNsol, 'b--','LineWidth', 2.5); hold on;
plot(tspan, LIsol, 'r','LineWidth', 2.5); hold on;
plot(tspan, LNsol, 'r--','LineWidth', 2.5); hold on;
title('Cell numbers')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'linewidth',1.5);
%l = legend('high-proliferating - immunogenic', 'high-proliferating - low-immunogenic', ...
  %  'low-proliferating - immunogenic', 'low-proliferating - low-immunogenic','fontsize',6,'location','best')
%set(l,'Fontsize',20);

subplot(3,1,2);
plot(tspan,Asol,'c','LineWidth',2.5); hold on;
plot(tspan,Ssol,'m','LineWidth',2.5); hold on;
title('Cytokines')
%legend('inflammatory', 'anti-inflamatory') 
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'linewidth',1.5);
% 
% subplot(3,1,2);
% plot(tspan,Ssol,'k'); hold on;
% plot(tscat, IL10,'ko');
% plot(BCGtimes,IL10_yay,'bo');
% plot(BCGtimes,IL10_no,'rx');
% title('Immunosupressor (IL-10)');
% set(gca,'Xlim',[-10 tmax+5]);
% 
subplot(3,1,3);
plot(tspan, BCG,'k','LineWidth',2.5); hold on;
plot(tspan, BCG2,'g','LineWidth',2.5); hold on;
%plot(tspan, 10000.*TT,'r','LineWidth',2); hold on;
title('BCG and Everolimus (MTOR inhibitor)');
xlabel('Time (days)')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'linewidth',1.5);
% set(gca,'Ylim',[-10 tmax+5]);








figure;
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])
subplot(2,2,1);
plot(tspan, totalcells0, 'k','LineWidth', 1.5); hold on;
title('BCG Only')
ylabel('Cells')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'Ylim',[0 2e12],'linewidth',1.5);
%l = legend('high-proliferating - immunogenic', 'high-proliferating - low-immunogenic', ...
  %  'low-proliferating - immunogenic', 'low-proliferating - low-immunogenic','fontsize',6,'location','best')
%set(l,'Fontsize',20);

subplot(2,2,2);
plot(tspan, totalcells, 'k','LineWidth', 1.5); hold on;
title('BCG and Everolimus (MTOR inhibitor)')
ylabel('Cells')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'Ylim',[0 2e12],'linewidth',1.5);

subplot(2,2,3);
plot(tspan, BCG,'k','LineWidth',2.5); hold on;
%plot(tspan, 10000.*TT,'r','LineWidth',2); hold on;
ylabel('BCG')
xlabel('Time (days)')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'linewidth',1.5);



subplot(2,2,4);
plot(tspan, BCG2,'k','LineWidth',2.5); hold on;
plot(tspan, 10000.*TT,'r','LineWidth',2); hold on;
%plot(tspan, 10000.*TT,'r','LineWidth',2); hold on;
ylabel('BCG and Everolimus')
xlabel('Time (days)')
set(gca,'Xlim',[-10 tmax+5],'fontsize',26,'linewidth',1.5);
% set(gca,'Ylim',[-10 tmax+5]);


