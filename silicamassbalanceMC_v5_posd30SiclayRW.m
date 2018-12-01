%%
%Silicon isotope mass balance Monte Carlo model
%Elizabeth Trower, November 2018
%this code was designed with Matlab 2017b

%this version of the model includes iron formation (IF) as a silica sink
%and applies f_clay based on the estimate of the reverse weathering
%authigenic clay sink from Isson and Planavsky (2018). Unlike previous
%model versions, this model assumes that d30Si of clays formed via reverse
%weathering is > d30Si_BSE and includes an additional, separate sink for
%authigenic clays formed via regular silicate weathering.

clear

load('Precambriand30Sidata.mat');
load('clayd30Sidata');
Clay = clayd30Sidata;
load('IFabund.mat');
load('Tatzeletal2017data.mat');
spicularcht = Tatzeletal2017spicularchert;

n = 10000;
tstep = 0.001; %time step in Ga

BSE = -0.29;
BSE_1sd = 0.04;
pdfBSE = makedist('Normal','mu',BSE,'sigma',BSE_1sd);

z00 =rand(1,n);
icdfBSE = icdf(pdfBSE,z00);

pdfClay_sw = fitdist(Clay,'Kernel');
z0 = rand(1,n);
icdfClay_sw = icdf(pdfClay_sw,z0);

xfall = min(spicularcht(:,3)):tstep:max(PrChtall_sort(:,3));
Chall = cat(1,spicularcht,silcarball_sort,PrChtall_sort);

fclay_isson = -.0126.*xfall.^2 - 0.0234*xfall + 0.3552;
fclay_issonsims = abs(0.1*randn(n,length(xfall)) + fclay_isson);

pdfClay_rw = makedist('Normal','mu',1,'sigma',1);
zrw = rand(1,n);
icdfClay_rw = icdf(pdfClay_rw,zrw);

[fIFfit,fIFgof] = fit(IFabund(:,1),IFabund(:,2),'SmoothingSpline');
yfIF = fIFfit(xfall);
yfIF(yfIF<0) = 0;

%%
fIFsims = abs(fIFgof.rmse.*randn(n,length(xfall)) + ...
    yfIF'.*ones(n,length(xfall)));

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
dClaymat = zeros(n,length(trange));
fChmat = zeros(n,length(trange));

percs = 10:10:90;
fClaypercs = zeros(length(percs),length(trange));
fChpercs = zeros(length(percs),length(trange));
fIFpercs = zeros(length(percs),length(trange));

%%
for counter0 = 1:length(trange)
    
    fClay_sw = (Chsims(:,counter0).*(1 - fclay_issonsims(:,counter0) - ...
        fIFsims(:,counter0)) + icdfClay_rw'.*fclay_issonsims(:,counter0) +...
        fIFsims(:,counter0).*IFsims(:,counter0) - icdfBSE')./...
        (Chsims(:,counter0) - icdfClay_rw');

    fCh = 1 - fclay_issonsims(:,counter0) - fClay_sw - fIFsims(:,counter0);

    fCh(fCh < 0) = NaN;
    fCh(fCh > 1) = NaN;
    fClay_sw(fClay_sw < 0) = NaN;
    fClay_sw(fClay_sw > 1) = NaN;
    fClay_sw(isnan(fCh) == 1) = NaN;
    fCh(isnan(fClay_sw) == 1) = NaN;

    fChpercs(:,counter0) = prctile(fCh,percs);
    fIFpercs(:,counter0) = prctile(fIFsims(:,counter0),percs);
    fClaypercs(:,counter0) = prctile(fClay_sw,percs);
    fClaymat(:,counter0) = fClay_sw;
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
h1.XBinLimits = [.55 3.7];
h1.NumBins = [32 21];
fCh_histvals = h1.Values;
normfactor1 = sum(fCh_histvals,2);
fCh_histnorms = fCh_histvals./normfactor1;

x_fch = 0.575:0.1:3.675;
x_fch(1) = 0.57;
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
caxis([0 .2])

fig1.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv5_posd30SiclayRW_fch','epsc')

%%
fig2 = figure;
h2 = histogram2(tmat,fClaymat,'DisplayStyle','tile');
h2.XBinLimits = [.55 3.7];
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
ylabel('f_c_l_a_y_,_s_w')
colorbar
caxis([0 .2])

fig2.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv5_posd30SiclayRW_fclay','epsc')
