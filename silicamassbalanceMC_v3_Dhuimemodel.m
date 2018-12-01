%%
%Silicon isotope mass balance Monte Carlo model
%Elizabeth Trower, November 2018
%this code was designed with Matlab 2017b

%this version of the model includes iron formation (IF) as a silica sink
%and assumes a variable Si isotope composition for clay, here based on a
%Mesoarchean initiation of plate tectonics, following Dhuime et al. (2015).

clear

load('Precambriand30Sidata.mat');
load('clayd30Sidata');
Clay = clayd30Sidata;
load('IFabund.mat');
load('Tatzeletal2017data.mat');
spicularcht = Tatzeletal2017spicularchert;
load('Dhuime_fdata.mat');

n = 10000;
tstep = 0.001; %time step in Ga

BSE = -0.29;
BSE_1sd = 0.04;
pdfBSE = makedist('Normal','mu',BSE,'sigma',BSE_1sd);

z00 =rand(1,n);
icdfBSE = icdf(pdfBSE,z00);

feldspar = -0.2;
olivine = -0.4;

authclayfract = -1;

claytime = Age;
f_fspar = ffelsic;
f_ol = fmafic;

claycombined = f_fspar*(feldspar+authclayfract) + f_ol*(olivine...
    +authclayfract);

%%
F2 = griddedInterpolant(claytime,claycombined);

xfall = min(spicularcht(:,3)):tstep:max(PrChtall_sort(:,3));
Chall = cat(1,spicularcht,silcarball_sort,PrChtall_sort);

Claysims = 0.9*randn(n,length(xfall)) + F2(xfall);

[fIFfit,fIFgof] = fit(IFabund(:,1),IFabund(:,2),'SmoothingSpline');
yfIF = fIFfit(xfall);
yfIF(yfIF<0) = 0;

%%
fIFsims = fIFgof.rmse.*randn(n,length(xfall)) + yfIF'.*ones(n,length(xfall));
fIFsims(fIFsims<0) = 0;

[IFfit,IFgof] = fit(IFall_sort(:,3),IFall_sort(:,1),'smoothingspline',...
    'SmoothingParam',0.9);
yIF = IFfit(xfall);

[Chfit,Chgof] = fit(Chall(:,3),Chall(:,1),...
    'SmoothingSpline','SmoothingParam',0.9);
yCh = Chfit(xfall);

IFsims = IFgof.rmse.*randn(n,length(xfall)) + yIF'.*ones(n,length(xfall));
Chsims = Chgof.rmse.*randn(n,length(xfall)) + yCh'.*ones(n,length(xfall));

trange = min(xfall):tstep:max(xfall);

fClaymat = zeros(n,length(trange));
fChmat = zeros(n,length(trange));

percs = 10:10:90;
fClaypercs = zeros(length(percs),length(trange));
fChpercs = zeros(length(percs),length(trange));
fIFpercs = zeros(length(percs),length(trange));

%%
for counter0 = 1:length(trange)
    
    fClay = (icdfBSE' - fIFsims(:,counter0).*(IFsims(:,counter0) -...
        Chsims(:,counter0)) - Chsims(:,counter0))./(Claysims(:,counter0) - ...
        Chsims(:,counter0));
    fCh = 1 - fClay - fIFsims(:,counter0);

    fClay(fClay > 1) = NaN;
    fClay(fClay < 0) = NaN;
    fClay(fCh > 1) = NaN;
    fClay(fCh < 0) = NaN;

    fCh(fCh < 0) = NaN;
    fCh(fCh > 1) = NaN;
    fCh(isnan(fClay) == 1) = NaN;

    fClaypercs(:,counter0) = prctile(fClay,percs);
    fChpercs(:,counter0) = prctile(fCh,percs);
    fIFpercs(:,counter0) = prctile(fIFsims(:,counter0),percs);
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
tmat = trange.*ones(n,length(trange));

fig1 = figure;
h1 = histogram2(tmat,fChmat,'DisplayStyle','tile');
h1.XBinLimits = [.52 3.7];
h1.NumBins = [32 21];
fCh_histvals = h1.Values;
normfactor1 = sum(fCh_histvals,2);
fCh_histnorms = fCh_histvals./normfactor1;

x_fch = 0.525:0.1:3.675;
x_fch(1) = 0.52;
x_fch(end) = 3.7;

y_fch = 0.:0.05:1;

p1 = pcolor(x_fch,y_fch,fCh_histnorms');
p1.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fChpercs(5,:),'k','LineWidth',2)
plot(trange,fChpercs(1,:),'k','LineWidth',.25)
plot(trange,fChpercs(9,:),'k','LineWidth',.25)
xlabel('age (Ga)')
ylabel('f_c_h_e_r_t')
colorbar
caxis([0 .1])

fig1.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv3_dhuime_fch','epsc')

%%
fig2 = figure;
h2 = histogram2(tmat,fClaymat,'DisplayStyle','tile');
h2.XBinLimits = [.52 3.7];
h2.NumBins = [32 21];
fClay_histvals = h2.Values;
normfactor2 = sum(fClay_histvals,2);
fClay_histnorms = fClay_histvals./normfactor2;

x_fclay = x_fch;
y_fclay = y_fch;

p2 = pcolor(x_fclay,y_fclay,fClay_histnorms');
p2.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fClaypercs(5,:),'k','LineWidth',2)
plot(trange,fClaypercs(1,:),'k','LineWidth',0.25)
plot(trange,fClaypercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_c_l_a_y')
colorbar
caxis([0 .1])

fig2.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv3_dhuime_fclay','epsc')

%%
tmat = trange.*ones(n,length(trange));

fig3 = figure;
h3 = histogram2(tmat,Claysims,'DisplayStyle','tile');
h3.XBinLimits = [.52 3.7];
h3.YBinLimits = [-6 2];
h3.NumBins = [32 100];

dClay_histvals = h3.Values;
normfactor3 = sum(dClay_histvals,2);
dClay_histnorms = dClay_histvals./normfactor3;

x_dClay = x_fch;
y_dClay = linspace(-6,2,100);

p3 = pcolor(x_dClay,y_dClay,dClay_histnorms');
p3.EdgeColor = 'none';
xlim([0 4])
ylim([-5 1])
xlabel('age (Ga)')
ylabel('d^3^0Si_c_l_a_y')
colorbar

fig3.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv3_dhuime_dClay','epsc')
