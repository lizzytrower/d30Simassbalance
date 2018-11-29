%%
%Silicon isotope mass balance Monte Carlo model
%Elizabeth Trower, November 2018
%this code was designed with Matlab 2017b

clear

load('Precambriand30Sidata.mat');
load('clayd30Sidata');
Clay = clayd30Sidata;
load('BIFabund.mat');
load('Tatzeletal2017data.mat');
spicularcht = Tatzeletal2017spicularchert;

n = 10000;
tstep = 0.001; %time step in Ga

BSE = -0.29;
BSE_1sd = 0.04;
pdfBSE = makedist('Normal','mu',BSE,'sigma',BSE_1sd);

z00 =rand(1,n);
icdfBSE = icdf(pdfBSE,z00);

pdfClay = fitdist(Clay,'Kernel');

z0 = rand(1,n);
icdfClay = icdf(pdfClay,z0);

xfall = min(spicularcht(:,3)):tstep:max(PrChtall_sort(:,3));
Chall = cat(1,spicularcht,silcarball_sort,PrChtall_sort);

[Chfit,Chgof] = fit(Chall(:,3),Chall(:,1),...
    'SmoothingSpline','SmoothingParam',0.9);
yCh = Chfit(xfall);

Chsims = Chgof.rmse.*randn(n,length(xfall)) + yCh'.*ones(n,length(xfall));

trange = min(xfall):tstep:max(xfall);

fClaymat = zeros(n,length(trange));
fChmat = zeros(n,length(trange));

percs = 10:10:90;
fClaypercs = zeros(length(percs),length(trange));
fChpercs = zeros(length(percs),length(trange));

for counter0 = 1:length(trange)

        fClay = (icdfBSE' - Chsims(:,counter0))./(icdfClay' - ...
            Chsims(:,counter0));
        fCh = 1 - fClay;
        
        fClay(fClay > 1) = NaN;
        fClay(fClay < 0) = NaN;
        fClay(fCh > 1) = NaN;
        fClay(fCh < 0) = NaN;
        
        fCh(fCh < 0) = NaN;
        fCh(fCh > 1) = NaN;
        fCh(isnan(fClay) == 1) = NaN;
        
        fClaypercs(:,counter0) = prctile(fClay,percs);
        fChpercs(:,counter0) = prctile(fCh,percs);
        fClaymat(:,counter0) = fClay;
        fChmat(:,counter0) = fCh;
    
end

%%
%count how many NaN's are in each vector
fChnans = sum(isnan(fChmat));
fClaynans = sum(isnan(fClaymat));

figure
plot(trange,fChnans)
hold on
plot(trange,fClaynans)

%%
fig1 = figure;
tmat = trange.*ones(n,length(trange));
h1 = histogram2(tmat,Chsims,'DisplayStyle','tile');
h1.XBinLimits = [.52 3.7];
h1.YBinLimits = [-4 6];
h1.NumBins = [64 100];
Chsims_histvals = h1.Values;
normfactor = sum(Chsims_histvals,2);
Chsims_histnorm = Chsims_histvals./normfactor;

x_chsims = 0.525:0.05:3.675;
x_chsims(1) = 0.52;
x_chsims(end) = 3.7;
y_chsims = -3.95:.1:5.95;
y_chsims(1) = -4;
y_chsims(end) = 6;

p1 = pcolor(x_chsims,y_chsims,Chsims_histnorm');
hold on
p1.EdgeColor = 'none';
xwnsilcarb = awgn(silcarball_sort(:,3),40);
scatter(xwnsilcarb,silcarball_sort(:,1),'.k')
xwncht = awgn(PrChtall_sort(:,3),40);
scatter(xwncht,PrChtall_sort(:,1),'.k')
xwnspicularcht = awgn(spicularcht(:,3),40);
scatter(xwnspicularcht,spicularcht(:,1),'.k')
ylim([-4 6])
xlim([0 4])
colorbar
caxis([0 0.05])
xlabel('time (Ga)')
ylabel('chert d^3^0Si')
fig1.Renderer = 'painters';

saveas(gcf,'chertSisims','epsc')

%%
fig2 = figure;
tmat = trange.*ones(n,length(trange));
h2 = histogram2(tmat,fChmat,'DisplayStyle','tile');
h2.XBinLimits = [.52 3.7];
h2.NumBins = [32 21];
fCh_histvals = h2.Values;
normfactor2 = sum(fCh_histvals,2);
fCh_histnorms = fCh_histvals./normfactor2;

x_fch = 0.525:0.1:3.675;
x_fch(1) = 0.52;
x_fch(end) = 3.7;

y_fch = 0.:0.05:1;

p2 = pcolor(x_fch,y_fch,fCh_histnorms');
p2.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fChpercs(5,:),'k','LineWidth',2)
plot(trange,fChpercs(1,:),'k','LineWidth',0.25)
plot(trange,fChpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_c_h_e_r_t')
colorbar
caxis([0 0.1])

fig2.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv1_fch','epsc')

%%
fig3 = figure;
h3 = histogram2(tmat,fClaymat,'DisplayStyle','tile');
h3.XBinLimits = [.52 3.7];
h3.NumBins = [32 21];
fClay_histvals = h3.Values;
normfactor3 = sum(fClay_histvals,2);
fClay_histnorms = fClay_histvals./normfactor3;

x_fclay = x_fch;
y_fclay = y_fch;

p3 = pcolor(x_fclay,y_fclay,fClay_histnorms');
p3.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fClaypercs(5,:),'k','LineWidth',2)
plot(trange,fClaypercs(1,:),'k','LineWidth',0.25)
plot(trange,fClaypercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_c_l_a_y')
colorbar
caxis([0 0.1])

fig3.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv1_fclay','epsc')