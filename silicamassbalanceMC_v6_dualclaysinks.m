%%
%Silicon isotope mass balance Monte Carlo model
%Elizabeth Trower, February 2019
%this code was designed with Matlab 2018b

%this version of the model includes iron formation (IF) as a silica sink
%and applies f_RW based on the estimate of the reverse weathering
%authigenic clay sink from Isson and Planavsky (2018) and a fractionation
%factor between chert and reverse weathering-type clays of -2 permil
%based on Ehlert et al. (2016).

clear

load('Precambriand30Sidata.mat');
load('IFabund.mat');
load('Tatzeletal2017data.mat');
spicularcht = Tatzeletal2017spicularchert;
load('clayd30Sidata');
kaol = clayd30Sidata;

fractionation_chert_RW = -2;

n = 10000;
tstep = 0.001; %time step in Ga

BSE = -0.29;
BSE_1sd = 0.04;
pdfBSE = makedist('Normal','mu',BSE,'sigma',BSE_1sd);

z00 =rand(1,n);
icdfBSE = icdf(pdfBSE,z00);

pdfkaol = fitdist(kaol,'Kernel');

z0 = rand(1,n);
icdfkaol = icdf(pdfkaol,z0);

xfall = min(spicularcht(:,3)):tstep:max(PrChtall_sort(:,3));
Chall = cat(1,spicularcht,silcarball_sort,PrChtall_sort);

fRW_isson = -0.0126.*xfall.^2 - 0.0234*xfall + 0.3552;
fRW_issonsims = abs(0.1*randn(n,length(xfall)) + fRW_isson);

[fIFfit,fIFgof] = fit(IFabund(:,1),IFabund(:,2),'SmoothingSpline');
yfIF = fIFfit(xfall);
yfIF(yfIF<0) = 0;

%%
fIFsims = abs(fIFgof.rmse.*randn(n,length(xfall)) + ...
    yfIF'.*ones(n,length(xfall)));
fChsims = 1 - fRW_issonsims - fIFsims;

[IFfit,IFgof] = fit(IFall_sort(:,3),IFall_sort(:,1),'smoothingspline',...
    'SmoothingParam',0.9);
yIF = IFfit(xfall);

[Chfit,Chgof] = fit(Chall(:,3),Chall(:,1),...
    'SmoothingSpline','SmoothingParam',0.9);
yCh = Chfit(xfall);

IFsims = IFgof.rmse.*randn(n,length(xfall)) + yIF'.*ones(n,length(xfall));
Chsims = Chgof.rmse.*randn(n,length(xfall)) + yCh'.*ones(n,length(xfall));
dRWsims = Chsims + fractionation_chert_RW;

trange = min(xfall):tstep:max(xfall);

fRWmat = zeros(n,length(trange));
fIFmat = zeros(n,length(trange));
fChmat = zeros(n,length(trange));
fkaolmat = zeros(n,length(trange));
fclaytotmat = zeros(n,length(trange));
dclaytotmat = zeros(n,length(trange));

percs = 10:10:90;
dRWpercs = zeros(length(percs),length(trange));
fChpercs = zeros(length(percs),length(trange));
fIFpercs = zeros(length(percs),length(trange));
fRWpercs = zeros(length(percs),length(trange));
fkaolpercs = zeros(length(percs),length(trange));
fclaytotpercs = zeros(length(percs),length(trange));
dclaytotpercs = zeros(length(percs),length(trange));

%%
for counter0 = 1:length(trange)
    
    dmarine = fIFsims(:,counter0).*IFsims(:,counter0) + ...
        dRWsims(:,counter0).*fRW_issonsims(:,counter0) + ...
        Chsims(:,counter0).*fChsims(:,counter0);
    
    fkaol = (icdfBSE' - dmarine)./(icdfkaol' - dmarine);
    
    fmarine = 1 - fkaol;
    fIF = fmarine.*fIFsims(:,counter0);
    fRW = fmarine.*fRW_issonsims(:,counter0);
    fCh = fmarine.*fChsims(:,counter0);
    fIF(fkaol<0) = NaN;
    fRW(fkaol<0) = NaN;
    fCh(fkaol<0) = NaN;
    fIF(fkaol>1) = NaN;
    fRW(fkaol>1) = NaN;
    fCh(fkaol>1) = NaN;
    fkaol(fkaol<0) = NaN;
    fkaol(fkaol>1) = NaN;
    fkaolmat(:,counter0) = fkaol;
    fIFmat(:,counter0) = fIF;
    fRWmat(:,counter0) = fRW;
    fChmat(:,counter0) = fCh;
    fclaytotmat(:,counter0) = fkaol + fRW;
    dclaytotmat(:,counter0) = fkaol.*icdfkaol' + fRW.*fRW_issonsims(:,counter0);
    
    fkaolpercs(:,counter0) = prctile(fkaol,percs);
    fChpercs(:,counter0) = prctile(fCh,percs);
    fRWpercs(:,counter0) = prctile(fRW,percs);
    fIFpercs(:,counter0) = prctile(fIF,percs);
    dRWpercs(:,counter0) = prctile(dRWsims(:,counter0),percs);
    fclaytotpercs(:,counter0) = prctile(fkaol+fRW,percs);

    
end

tots = fkaolmat + fChmat + fIFmat + fRWmat;
mismatch = icdfBSE' - fkaolmat.*icdfkaol' - fChmat.*Chsims - ...
    fIFmat.*IFsims - fRWmat.*dRWsims;

%%
%count how many NaN's are in each vector
fChnans = sum(isnan(fChmat));

figure
plot(trange,fChnans)

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
caxis([0 .12])

fig1.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_fch','epsc')

%%
fig2 = figure;
h2 = histogram2(tmat,fkaolmat,'DisplayStyle','tile');
h2.XBinLimits = [.52 3.7];
h2.NumBins = [32 21];
fkaol_histvals = h2.Values;
normfactor2 = sum(fkaol_histvals,2);
fkaol_histnorms = fkaol_histvals./normfactor2;

x_fclay = x_fch;
y_fclay = y_fch;

p2 = pcolor(x_fclay,y_fclay,fkaol_histnorms');
p2.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fkaolpercs(5,:),'k','LineWidth',2)
plot(trange,fkaolpercs(1,:),'k','LineWidth',0.25)
plot(trange,fkaolpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_k_a_o_l')
colorbar
caxis([0 .12])

fig2.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_fkaol','epsc')
%%
fig3 = figure;
h3 = histogram2(tmat,fRWmat,'DisplayStyle','tile');
h3.XBinLimits = [.52 3.7];
h3.NumBins = [32 21];
fRW_histvals = h3.Values;
normfactor3 = sum(fRW_histvals,2);
fRW_histnorms = fRW_histvals./normfactor3;

x_fclay = x_fch;
y_fclay = y_fch;

p3 = pcolor(x_fclay,y_fclay,fRW_histnorms');
p3.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fRWpercs(5,:),'k','LineWidth',2)
plot(trange,fRWpercs(1,:),'k','LineWidth',0.25)
plot(trange,fRWpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_R_W')
colorbar
caxis([0 .4])

fig3.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_fRW','epsc')

%%
fig4 = figure;
h4 = histogram2(tmat,fIFmat,'DisplayStyle','tile');
h4.XBinLimits = [.52 3.7];
h4.NumBins = [32 21];
fIF_histvals = h4.Values;
normfactor4 = sum(fIF_histvals,2);
fIF_histnorms = fIF_histvals./normfactor4;

x_fIF = x_fch;
y_fIF = y_fch;

p4 = pcolor(x_fIF,y_fIF,fIF_histnorms');
p4.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fIFpercs(5,:),'k','LineWidth',2)
plot(trange,fIFpercs(1,:),'k','LineWidth',0.25)
plot(trange,fIFpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_I_F')
colorbar
caxis([0 .7])

fig4.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_fIF','epsc')

%%
fig5 = figure;
h5 = histogram2(tmat,dRWsims,'DisplayStyle','tile');
h5.XBinLimits = [.52 3.7];
h5.YBinLimits = [-4 3];
h5.NumBins = [32 29];
dRWhistvals = h5.Values;
normfactor5 = sum(dRWhistvals,2);
dRW_histnorms = dRWhistvals./normfactor5;

y_dclay = -4:.25:3;

p5 = pcolor(x_fclay,y_dclay,dRW_histnorms');
p5.EdgeColor = 'none';
xlim([0 4])
ylim([-4 3])
hold on
plot(trange,dRWpercs(5,:),'k','LineWidth',2)
plot(trange,dRWpercs(1,:),'k','LineWidth',0.25)
plot(trange,dRWpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('d^3^0Si_R_W')
colorbar
caxis([0 .12])

fig5.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_dRW','epsc')

%%
fig6 = figure;
h6 = histogram2(tmat,fclaytotmat,'DisplayStyle','tile');
h6.XBinLimits = [.52 3.7];
h6.NumBins = [32 21];
fclaytot_histvals = h6.Values;
normfactor6 = sum(fclaytot_histvals,2);
fclaytot_histnorms = fclaytot_histvals./normfactor6;

x_fclay = x_fch;
y_fclay = y_fch;

p6 = pcolor(x_fclay,y_fclay,fclaytot_histnorms');
p6.EdgeColor = 'none';
xlim([0 4])
ylim([0 1])
hold on
plot(trange,fclaytotpercs(5,:),'k','LineWidth',2)
plot(trange,fclaytotpercs(1,:),'k','LineWidth',0.25)
plot(trange,fclaytotpercs(9,:),'k','LineWidth',0.25)
xlabel('age (Ga)')
ylabel('f_c_l_a_y_t_o_t')
colorbar
caxis([0 .12])

fig6.Renderer = 'painters';
saveas(gcf,'Troweretalmodelv6_revision_fclaytot','epsc')