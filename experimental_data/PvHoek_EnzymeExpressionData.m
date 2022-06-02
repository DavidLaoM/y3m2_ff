% Interpolation of the correlation in enzyme concentration. Using the
% PvHoek dataset.
    % The if-statement aims to speed up the process when the assays do not
    % involve different growth rates and are fixed at d = 0.1. Interpolation
    % can become time consuming in out dynamic simulations.

if d == 0.1 
    p.HXK_ExprsCor  = 1;
    p.PGI_ExprsCor  = 1;
    p.PFK_ExprsCor  = 1;
    p.FBA_ExprsCor  = 1;
    p.TPI_ExprsCor  = 1;
    p.GAPDH_ExprsCor= 1;
    p.PGK_ExprsCor  = 1;
    p.PGM_ExprsCor  = 1;
    p.ENO_ExprsCor  = 1;
    p.PYK_ExprsCor  = 1;
    p.PDC_ExprsCor  = 1;
    p.ADH_ExprsCor  = 1;
else
    PvHoek.D_HXK=[0.048535	0.065089	0.079577	0.094065	0.112703	0.135481	0.154118	0.176879	0.199657	0.222436	0.245214	0.270084	0.303244	0.332263	0.34886	0.375839];
    PvHoek.HXK=[2.612394 2.455598 2.330085 2.204571 2.078798 1.921615 1.795843 1.576355 1.419172 1.261989 1.104806 1.009798 0.883121 0.787855 0.78682 0.816292];

    PvHoek.D_PGI=[0.052684	0.069281	0.090028	0.10455	0.123222	0.146044	0.168865	0.187537	0.214507	0.235254	0.260149	0.28712	0.31409	0.336912	0.359733	0.374264];
    PvHoek.PGI=[2.643322	2.642287	2.640995	2.64009	2.638927	2.637505	2.636083	2.634919	2.633239	2.631946	2.630395	2.628715	2.627034	2.625612	2.62419	2.654438];

    PvHoek.D_PFK=[0.055866	0.082835	0.116028	0.149221	0.178265	0.209348	0.238392	0.26951	0.29233	0.319334	0.338039	0.348498	0.358991	0.377817];
    PvHoek.PFK=[0.313703	0.313367	0.312953	0.31254	0.312178	0.30556	0.305198	0.30481	0.304525	0.31042	0.316417	0.331864	0.353541	0.381345];

    PvHoek.D_FBA=[0.049436	0.078411	0.111552	0.146716	0.169467	0.192287	0.217147	0.239933	0.268942	0.293888	0.308531	0.323225	0.333804	0.350572	0.359043	0.371714];
    PvHoek.FBA=[1.114507	1.063375	1.024458	0.948171	0.897348	0.896214	0.870132	0.844154	0.817866	0.853896	0.94013	1.063632	1.212185	1.335584	1.459395	1.620268];

    PvHoek.D_TPI=[0.049587	0.070248	0.086777	0.115702	0.154959	0.183884	0.216942	0.247934	0.280992	0.307851	0.332645	0.349174	0.376033];
    PvHoek.TPI=[53.582555	50.46729	49.221184	46.728972	44.859813	43.613707	41.121495	42.367601	47.975078	54.205607	61.05919	65.420561	72.897196];

    PvHoek.D_GAPDH=[0.051867	0.072614	0.09751	0.114108	0.134855	0.155602	0.178423	0.19917	0.221992	0.248963	0.275934	0.29668	0.319502	0.348548	0.375519];
    PvHoek.GAPDH=[5.981308	5.358255	4.672897	4.23676	3.862928	3.489097	3.239875	3.05296	2.990654	3.115265	3.613707	3.925234	4.423676	5.109034	6.292835];

    PvHoek.D_PGK=[0.052523	0.072215	0.095844	0.119474	0.135221	0.154929	0.178566	0.202218	0.229816	0.255468	0.287065	0.306811	0.32855	0.34636	0.364177	0.382034];
    PvHoek.PGK=[7.979653	7.444776	6.791089	6.137403	5.662165	5.245629	4.651114	4.17494	3.639128	3.340233	3.158978	3.038297	3.035725	3.211133	3.445711	3.976145];

    PvHoek.D_PGM=[0.050367	0.068099	0.085838	0.101601	0.119333	0.139056	0.160741	0.190345	0.22591	0.253578	0.283268	0.307077	0.321059	0.339017	0.353015	0.377043];
    PvHoek.PGM=[6.63844	6.220898	5.862704	5.504744	5.087203	4.788121	4.370111	4.069856	4.006287	4.003003	4.355567	5.064918	6.190873	7.4944	8.739048	11.110121];

    PvHoek.D_ENO=[0.051383	0.063241	0.08498	0.108696	0.134387	0.158103	0.179842	0.211462	0.243083	0.262846	0.286561	0.302372	0.322134	0.337945	0.351779	0.373518];
    PvHoek.ENO=[0.66568	0.606509	0.532544	0.473373	0.399408	0.325444	0.295858	0.281065	0.295858	0.295858	0.340237	0.414201	0.56213	0.724852	0.843195	1.050296];

    PvHoek.D_PYK=[0.05336	0.075099	0.096838	0.12253	0.146245	0.167984	0.193676	0.221344	0.250988	0.27668	0.296443	0.310277	0.324111	0.335968	0.351779	0.359684	0.37747];
    PvHoek.PYK=[2.751479	2.633136	2.514793	2.426036	2.366864	2.278107	2.189349	2.189349	2.278107	2.514793	2.95858	3.431953	3.934911	4.497041	5.147929	5.621302	6.508876];

    PvHoek.D_PDC=[0.047059	0.058824	0.068627	0.086275	0.1	0.115686	0.139216	0.170588	0.201961	0.229412	0.256863	0.276471	0.288235	0.301961	0.309804	0.319608	0.32549	0.333333	0.341176	0.35098	0.358824	0.376471];
    PvHoek.PDC=[0.792899	0.751479	0.704142	0.650888	0.60355	0.568047	0.544379	0.52071	0.514793	0.514793	0.532544	0.579882	0.650888	0.751479	0.840237	0.952663	1.04142	1.147929	1.254438	1.378698	1.443787	1.609467];

    PvHoek.D_ADH=[0.051181	0.076772	0.098425	0.124016	0.149606	0.173228	0.204724	0.222441	0.257874	0.281496	0.311024	0.324803	0.344488	0.375984];
    PvHoek.ADH=[9.881657	9.053254	8.343195	7.514793	6.745562	6.094675	5.088757	4.615385	3.727811	3.254438	2.781065	2.662722	2.95858	4.260355];

    % d = 0.1; %h^-1
    p.HXK_ExprsCor=interp1(PvHoek.D_HXK,PvHoek.HXK,d,'pchip','extrap')./interp1(PvHoek.D_HXK,PvHoek.HXK,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,1),plot(PvHoek.D_HXK,PvHoek.HXK./mean(PvHoek.HXK),'k*')
    % hold on
    % plot(d,p.HXK_ExprsCor,'r*')
    % title('HXK')

    p.PGI_ExprsCor=interp1(PvHoek.D_PGI,PvHoek.PGI,d,'pchip','extrap')./interp1(PvHoek.D_PGI,PvHoek.PGI,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,2),plot(PvHoek.D_PGI,PvHoek.PGI./mean(PvHoek.PGI),'k*')
    % hold on
    % plot(d,p.PGI_ExprsCor,'r*')
    % title('PGI')

    p.PFK_ExprsCor=interp1(PvHoek.D_PFK,PvHoek.PFK,d,'pchip','extrap')./interp1(PvHoek.D_PFK,PvHoek.PFK,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,3),plot(PvHoek.D_PFK,PvHoek.PFK./mean(PvHoek.PFK),'k*')
    % hold on
    % plot(d,p.PFK_ExprsCor,'r*')
    % title('PFK')

    p.FBA_ExprsCor=interp1(PvHoek.D_FBA,PvHoek.FBA,d,'pchip','extrap')./interp1(PvHoek.D_FBA,PvHoek.FBA,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,4),plot(PvHoek.D_FBA,PvHoek.FBA./mean(PvHoek.FBA),'k*')
    % hold on
    % plot(d,p.FBA_ExprsCor,'r*')
    % title('FBA')

    p.TPI_ExprsCor=interp1(PvHoek.D_TPI,PvHoek.TPI,d,'pchip','extrap')./interp1(PvHoek.D_TPI,PvHoek.TPI,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,5),plot(PvHoek.D_TPI,PvHoek.TPI./mean(PvHoek.TPI),'k*')
    % hold on
    % plot(d,p.TPI_ExprsCor,'r*')
    % title('TPI')

    p.GAPDH_ExprsCor=interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,d,'pchip','extrap')./interp1(PvHoek.D_GAPDH,PvHoek.GAPDH,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,6),plot(PvHoek.D_GAPDH,PvHoek.GAPDH./mean(PvHoek.GAPDH),'k*')
    % hold on
    % plot(d,p.GAPDH_ExprsCor,'r*')
    % title('GAPDH')

    p.PGK_ExprsCor=interp1(PvHoek.D_PGK,PvHoek.PGK,d,'pchip','extrap')./interp1(PvHoek.D_PGK,PvHoek.PGK,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,7),plot(PvHoek.D_PGK,PvHoek.PGK./mean(PvHoek.PGK),'k*')
    % hold on
    % plot(d,p.PGK_ExprsCor,'r*')
    % title('PGK')

    p.PGM_ExprsCor=interp1(PvHoek.D_PGM,PvHoek.PGM,d,'pchip','extrap')./interp1(PvHoek.D_PGM,PvHoek.PGM,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,8),plot(PvHoek.D_PGM,PvHoek.PGM./mean(PvHoek.PGM),'k*')
    % hold on
    % plot(d,p.PGM_ExprsCor,'r*')
    % title('PGM')

    p.ENO_ExprsCor=interp1(PvHoek.D_ENO,PvHoek.ENO,d,'pchip','extrap')./interp1(PvHoek.D_ENO,PvHoek.ENO,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,9),plot(PvHoek.D_ENO,PvHoek.ENO./mean(PvHoek.ENO),'k*')
    % hold on
    % plot(d,p.ENO_ExprsCor,'r*')
    % title('ENO')

    p.PYK_ExprsCor=interp1(PvHoek.D_PYK,PvHoek.PYK,d,'pchip','extrap')./interp1(PvHoek.D_PYK,PvHoek.PYK,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,10),plot(PvHoek.D_PYK,PvHoek.PYK./mean(PvHoek.PYK),'k*')
    % hold on
    % plot(d,p.PYK_ExprsCor,'r*')
    % title('PYK')

    p.PDC_ExprsCor=interp1(PvHoek.D_PDC,PvHoek.PDC,d,'pchip','extrap')./interp1(PvHoek.D_PDC,PvHoek.PDC,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,11),plot(PvHoek.D_PDC,PvHoek.PDC./mean(PvHoek.PDC),'k*')
    % hold on
    % plot(d,p.PDC_ExprsCor,'r*')
    % title('PDC')

    p.ADH_ExprsCor=interp1(PvHoek.D_ADH,PvHoek.ADH,d,'pchip','extrap')./interp1(PvHoek.D_ADH,PvHoek.ADH,0.1,'pchip','extrap');
    % figure(3),
    % subplot(4,3,12),plot(PvHoek.D_ADH,PvHoek.ADH./mean(PvHoek.ADH),'k*')
    % hold on
    % plot(d,p.ADH_ExprsCor,'r*')
    % title('ADH')  
end





% B = struct2cell(PvHoek);
% figure(4)
% for i = 1:12
% j = 2*i;
% o = j-1;
% subplot(3,4,i)
% plot(B{o}, B{j})
% end
% 
% % Loading van Heerden trehalose data file
% load('vHeerden_trehalose_data_micromolgdw.mat')
% figure(5)
% B1 = zeros(size(data.metabolites(:,1)));
% B2 = data.metabolites(:,2:end);
% data.metabolites2 = [B1 B2];
% plot(data.time_metabolites, data.metabolites2)    % taking out glucose. Too high concentration
% legend(data.legenda_metabolites)
% figure(6)
% plot(data.time_metabolites(:,1), data.metabolites(:,1))    % taking out glucose. Too high concentration
% legend(data.legenda_metabolites{1})
% figure(7)
% plot(data.time_nucleotides, data.nucleotides)
% legend(data.legenda_nucleotides)
% figure(8)
% plot(data.time_fluxes, data.fluxes)
% legend(data.legenda_fluxes)
% figure(9)
% plot(data.time_totP, data.totP)

% % @mu = 0.1 h^{-1}
%     p.HXK_ExprsCor  = 1;
%     p.PGI_ExprsCor  = 1;
%     p.PFK_ExprsCor  = 1;
%     p.FBA_ExprsCor  = 1;
%     p.TPI_ExprsCor  = 1;
%     p.GAPDH_ExprsCor= 1;
%     p.PGK_ExprsCor  = 1;
%     p.PGM_ExprsCor  = 1;
%     p.ENO_ExprsCor  = 1;
%     p.PYK_ExprsCor  = 1;
%     p.PDC_ExprsCor  = 1;
%     p.ADH_ExprsCor  = 1;
