% % PAREST1_REPLICATE_INITIAL.m
% In this code the tre cycle and glt, hxk are optimized as it happen before
% to obtain the parameter set x32.

% structure is:
    % early sampling
    % parameter estimation with paraemter and cost function combinations
    

%% (0) Initial setup.

% startup
clear, %close all
set_paths;
dbstop if error

% assays to run
setup.runSScanelas = 0;
setup.runGSvanHeerden = 0;
setup.plotResultsMode = 0; %30; %20; %10;
setup.conditionsSS = 0;
setup.conditionsGP = 0;

% model structure used
setup.biomass = 0; % biomass
setup.clamp10.Pi = 0;%1; 
setup.clamp10.TRE = 0;%1; 
setup.clamp10.NADX = 0;%1; 
setup.clamp10.AXP = 0;%1; 
setup.clamp10.IXP = 1; 
setup.clamp10.UDP_GLC = 1; % UDP_GLC_exp in AGT1 rate equation
setup.clamp = setup.clamp10;

% select parameters to study:
parsNTH1 = [122, 123];
parsATH1 = [149, 150];
parsAGT1 = [151:154, 158];  % include par158 (Ki for UDPG)
% parsvacuoleT = [155, 156];
parsvacuoleT = [155, 156, 157];
parsGLT = [35:36];
parsGLK = [28, 29, 30, 31, 32, 33, 34]; % 33 is k_T6P
parsPGI = [57:60];
parsPFK = [43:56,87];
parsPGM1 = [83:86];
parsUGP = [144:148];
parsTPS1 = [124:128];
parsTPS2 = [119:121];
parsPFK_sweep = [43:45,47,51,52,87];
parsPFK_selected = [51, 52];
parCombs = [33, parsNTH1, parsATH1, parsAGT1, parsvacuoleT, parsPGM1, parsUGP, parsTPS1, parsTPS2];

% get datasets (others needed)
[legenda] = legendaFull;
loadData_Y3M1;
% get GPFF datasets 
load('TUDdata.mat'); % loads the datasets
    dataset.FF01.time_mets = [0;5;10;15;20;25;30;60;90;120;150;180;220;250;300;350;400];
    dataset.FF01.timeECgluc= [0;5;11;15;20;   30;60;90;    150;180;220;250;300;350;400];
    dataset.FF03.time_mets = [0;5;10;20;30;40;60;90;120;150;200;250;300;400;550;700;800;900;1000;1200;1400;1600;1700;1803];
    dataset.FF04.time_mets = [0;5;10;15;20;30;60;90;120;150;180;220;250;300;350;398];
reorganiseTUDdata; % puts it in order for easily coding in the next section


%% (1) simulations (GP suarez)
% simulation FF01
setup.GPdataset.GP400WT = 1;
setup.GPdataset.GP1800WT = 0;
setup.GPdataset.GP400M = 0;
setup.GPdataset.Fructose  = 0;
setup.GPdataset.Sucrose = 0;

su = 0;
fn = fieldnames(setup.GPdataset);
for i = 1:numel(fn)
    su = su+ setup.GPdataset.(fn{i});
end
if su > 1
    disp('too many datasets used')
end

z = 2; %l/gcdw
if setup.GPdataset.GP400WT        == 1
    InitCond_GP_TUD;
elseif setup.GPdataset.GP1800WT 	== 1
    InitCond_GP_TUD_1800;
elseif setup.GPdataset.GP400M     == 1
    InitCond_GP_TUD_mutant;
elseif setup.GPdataset.Fructose == 1
    FFFruc = readFructose(z);
    dataset.FFFruc = FFFruc; 
    reorganiseFrucData;
    IC = ICFrucData(z, FFFruc);
elseif setup.GPdataset.Sucrose == 1
    FFSuc = readSucrose(z);
    dataset.FFSuc = FFSuc;
    reorganiseSucData;
    InitConditionsSucrose;
end
NumberCycles = 5;
plotflag = 0;


%% updated parameter set from pE2
load('pset_pE2res.mat');
x1 = x_pE2_start;
x2 = x_pE2_end;
%% glycerol_SynthDeg active
setup.glycSynthDeg = 1; % this adds a reaction to the mass balance (ODEs), reaction calculation in rate equations and the starting parameter value
x1(159) = -4; % -10;
x2(159) = -4; %-10;
%% extension glycogen metabolism
% load parameters
load('pset_pE6_x1res.mat','x_pE6_x1_end'); x1 = x_pE6_x1_end;
load('pset_pE6_x2res.mat','x_pE6_x2_end'); x2 = x_pE6_x2_end;
%% adding a more detailed glycogen metabolism
setup.glycSynthDeg_separate = 1;
%% recall last case
x = x1;
x(162) = 10;
setup.ATHinhibitionT6P = 1;
load('pset_pE7_xres.mat','x_pE7_start','x_pE7_end');
clear x x1 x2
x = x_pE7_end;
%% separating ath1_ec and ath1_vac
x(163:164) = x(149:150); % km_tre, k_cat
x(165) = x(162); % ki_t6p
setup.ATH_separate_EC_VAC = 1;
%% setting the option
setup.glycogenReactionsSink = 1;
setup.dataset = dataset;
%
load('pset_pE10_xres.mat','x_pE10_start','x_pE10_end'); x = x_pE10_end;
%% latest setup
% % % % setup.updated_bmf_Cx_ATH1ec = 1;
setup.TREec_brothOut_OFF = 1;
setup.updated_bmf_Cx_ATH1ec = 1;


%% parameter estimation setup
% latest added
setup.clamp_GLCec = 1;
NumberCycles = 5;

% %% testing the change in the mass balance
[legenda] = legendaFull; %legenda for the names needed
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
% plotflag = 2; % variables by iteration
choosedataset

% % previous setup
% setup.changing_Keq_glt = 1;
% tempValues = -[0 0.1 0.5 1 2 3];
% setup.Keq_glt_inc = tempValues(i);
% setup.clamp_GLCec = 1;

setup.changing_Keq_glt = 2;
x(166) = 0; % % <== change x166 if needed here.
x(167) = 0; % % <== change x167 if needed here.


%%
% blank and constant setup
blankWeight = zeros(1,88); % 85+3 for glycerol
    warray = blankWeight;
lambdalist = 0;
setup.parEst.lambda = lambdalist(1); lam = setup.parEst.lambda;
% % % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx});
% % % options = optimoptions('lsqnonlin','Display','iter','OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% % %     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% % options = optimoptions('lsqnonlin','Display','iter',...
% %     'OutputFcn',{@saveIterationsMain},'PlotFcn',{@optimplotx},...
% %     'FiniteDifferenceStepSize',0.2);
% options = optimoptions('lsqnonlin','Display','iter',...
%     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4,...
%     'DiffMinChange',0.1);
options = optimoptions('lsqnonlin','Display','off',...
    'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4,...
    'DiffMinChange',0.1);
%     'OutputFcn',{@saveIterationsMain},...
%     'FunctionTolerance', 1e-4, 'OptimalityTolerance', 1e-4);
% setup.parEst.costfun = 10015; 

% names_cell
ntests = 1;
% names_cell
names_cell = cell(ntests,1);
names_cell{1} = 'FF_pE12g_nMS%d.mat';
% parComb_cell
parsGLT = [35 36 38 166 167];
parsGLK = 28:34; % hxk
pC4 = [parsGLT, parsGLK];
% 
parComb_cell = cell(ntests,1);
parComb_cell{1} = pC4;

% number of multi starts
nMS = 1000;
xres_cell = cell(nMS,1);
simData_cell = cell(nMS,1);
simRes_cell = cell(nMS,1);

% array of parameters
xArray = zeros(nMS, length(x));
for i = 2:length(x)
    rng(i), xArray(i,pC4) = -3 + 6 * rand(1,length(pC4));
    rng(i+nMS), xArray(i,166:167) = -6 + 12 * rand(1,2);
end

% preparing coise for GLCin as well
A = 0.0001; B = 0.5;
LA = log10(A); LB = log10(B);
GLCin_data = 10.^(LA + (LB-LA) * rand(1,nMS));
%     figure, hist(log10(temp))


%% PARAMETER ESTIMATION: MS
selPars = parComb_cell{ntests};
tempName = names_cell{ntests};
%
for o = 1:nMS % SPLIT FOR SOME
    % saveName
    saveName = sprintf(tempName,o);

    % core options:
    plength = length(selPars);
    x_temp = xArray(o,selPars);
    setup.GLCin_data = GLCin_data(o);

    % boundaries
    lb = -3*ones(1,plength); %lb(end-12:end) = -5 * ones;
    ub = 3*ones(1,plength); %ub(end-12:end) = 5 * ones;
        % specific for glt
        lb(4:5) = -6*ones;
        ub(4:5) = 6*ones;

    % parameter estimation
    sprintf('MSnum = %d.', o)
    tic
    [xres,resnorm,residual,exitflag] = lsqnonlin(@costfunGLKsmallScale,x_temp,lb,ub,options,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
    t = toc; 
%     disp(xres);

    % simulation data
    [error, simData] = costfunGLKsmallScale_simData(x_temp,canelas_SS,setup,x,data,dataset,NumberCycles,IC,selPars,warray);
    
    % recall results
    xres_cell{o} = xres;
    simData_cell{o} = simData;
    
end
save('pEsts_12g.mat','xres_cell','simData_cell'); % SPLIT FOR SOME


%% see results
load('pEsts_12g.mat','xres_cell','simData_cell'); % SPLIT FOR SOME

% reorganise parameters
xAll = zeros(nMS,5);
for i = 1:nMS
    xAll(i,:) = xres_cell{i}(1:5);
end

% plot parameter estimates
% clf(100)
figure(100)
for j = 1:5
    subplot(2,3,j)
    histogram(xAll(:,j),50)%,'BinCounts',20)
    title(sprintf('pNum = %d',parsGLT(j)))
    if j < 4, xlim([-3 3]), else, xlim([-6 6]), end
end

% reorganise simData
v_GLT_All = zeros(nMS,5);
GLCec_All = zeros(nMS,5);
GLCi_All = zeros(nMS,1);
for i = 1:nMS
    xAll(i,:) = xres_cell{i}(1:5);
    v_GLT_All(i,:) = simData_cell{i}.v_GLT; % = zeros(nMS,5);
    GLCec_All(i,:) = simData_cell{i}.GLCec; % = zeros(nMS,5);
    GLCi_All(i,:) = simData_cell{i}.GLCi; % = zeros(nMS,1);
end

%% plotting in simulation view
% clf(101)
figure(101)
plot(dataset.FF01.fluxes_times,v_GLT_All(1:1000,:),'-','color',[.5 .5 .5])
hold on
plot(dataset.FF01.fluxes_times, dataset.FF01.fluxes{1},'ko','MarkerFaceColor','r')
ylim([0 1])
% %%
% temp1 = 1:5;
% % randArr = 
% rng(1), temp2 = rand(4,5);
% figure, plot(temp1, temp2, 'o-')
% rng(1), temp3 = rand(5,4);
% figure, plot(temp1, temp3, 'o-')

% creating color array (based on GLCi)
b = [0 0 1]; % Blue
r = [1 0 0]; % Red
% x = 0:0.01:1; % Values to interpolate on (This would be your actual vector of values)
x = GLCi_All;
xOk = [min(x);max(x)];
% Interpolate over all of the values in x
y = interp1(xOk,[b;r],x); % Here Blue corresponds to xOk == 0 and Red to xOk == 1
% figure
% scatter(x,x,[],y) % Scatter x vs. x but color based on the interpolated color.

% plotting with color accordign to GLCi
% clf(102)
figure(102)
for i = 750:1000%nMS
    plot(dataset.FF01.fluxes_times,v_GLT_All(i,:),'-','color',y(i,:))
    hold on
end
plot(dataset.FF01.fluxes_times, dataset.FF01.fluxes{1},'ko','MarkerFaceColor','r')
ylim([0 1])


%% plotting the scatted plot from before
clf(11)
figure(11), 
scatter(v_GLT_All(:,2), v_GLT_All(:,5))
xlim([0 50])
ylim([0 50])
% xlim([0 15])
% ylim([0 3.5])

%%
clf(10)
figure(10)
% range in area plot
x = 0:50;
y1 = x * 0.8139/3.787;
y2 = x * 3.132/3.315;
s2 = patch([x fliplr(x)], [y1 fliplr(y2)], [.9 .9 .9]);
s2.EdgeColor = [1 1 1];

% scatter data
hold on
s = scatter(v_GLT_All(:,2), v_GLT_All(:,5));
s.LineWidth = 0.6;
s.SizeData = 15;
s.MarkerEdgeColor = 'k';
s.MarkerFaceColor = 'r';
xlabel('v_{GLT,15}')
ylabel('v_{GLT,end}')
xlim([0 15])%50])
ylim([0 3.5])
    % focus
%     xlim([0 20])%50])
%     ylim([0 20])

% specific data point
hold on
plot(dataset.FF01.fluxes{1}(2), dataset.FF01.fluxes{1}(5),'ko','MarkerFaceColor','b')
legend('MPSA range', 'MPSA data points','experimental')


%% memoryDump

% % %%
% for i = 1:15
%     fid = fopen('parEst12f_estGLTGLK_selStart.m','rt');
%     X = fread(fid);
%     fclose(fid);
%     X = char(X.');
% %     str2rep = sprintf('vals2run = %d; % <--',i);
%     str2rep = sprintf('selCase = %d; % % <== changed already',i);
%     Y = strrep(X,'selCase = 1; % % <== HERE WE CHANGE FOR THE FAKE PARALLEL CALLS',str2rep);
%     tempName = sprintf('pE12f_%d.m',i);
%     fid2 = fopen(tempName,'wt') ;
%     fwrite(fid2,Y) ;
%     fclose(fid2);
% end

