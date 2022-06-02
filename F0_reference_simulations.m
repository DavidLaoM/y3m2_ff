% % F0_REFERENCE_SIMULATIONS.m
% This code runs to reproduce the feast famine simulations.
% David Lao-Martil, 2022-06-01

% initializing
set_paths;

% experimental data
[legenda] = legendaFull;
metNames = legenda.metabolites;
reactNames = legenda.fluxes;
legendaMetabolites_addEnrichment;
loadData_Y3M1; % SS and GP data. Only used when FF data is not available.
load('TUDdata.mat'); % FF data
    dataset.FF01.time_mets = [0;5;10;15;20;25;30;60;90;120;150;180;220;250;300;350;400];
    dataset.FF01.timeECgluc= [0;5;11;15;20;   30;60;90;    150;180;220;250;300;350;400];
    dataset.FF03.time_mets = [0;5;10;20;30;40;60;90;120;150;200;250;300;400;550;700;800;900;1000;1200;1400;1600;1700;1803];
    dataset.FF04.time_mets = [0;5;10;15;20;30;60;90;120;150;180;220;250;300;350;398];
reorganiseTUDdata;
load('datasetEnrich.mat'); % enrichment data
reorganiseEnrichData;

% select simulation to perform
setup.GPdataset.GP400WT = 1; % wild type strain, 400 seconds FF cycle
setup.GPdataset.GP1800WT = 0; % wild type strain, 1800 seconds FF cycle
setup.GPdataset.GP400M = 0; % TPS1 gene deletion mutant, 400 seconds FF cycle
if setup.GPdataset.GP400WT        == 1
    InitCond_GP_TUD;
elseif setup.GPdataset.GP1800WT 	== 1
    InitCond_GP_TUD_1800;
elseif setup.GPdataset.GP400M     == 1
    InitCond_GP_TUD_mutant;
end
choosedataset
setup.dataset = dataset; % required for some locations in the pipeline

% setup model structure
setup.biomass = 0; % 1 == include biomass. Else fixed at 0.1 h-1.
setup.clamp.Pi = 0; % 0 == cytsolic inorganic phosphate as dynamic variable. 1 == fixed at 10 mM.
setup.clamp.TRE = 0; % 0 == trehalose cycle as dynamic variable. 1 == fixed at experimental data values. 2 == fixed only at steady state data values
setup.clamp.NADX = 0; % 0 == mitochondrial NADH recycle. 1 == no mitochondiral recycle + NAD(H) concentrations fixed. 
setup.clamp.AXP = 0; % 0 == AXP cofactor metabolism modelled. 1 == not modelled + A(X)P concentrations fixed. 
setup.clamp.IXP = 1; % 0 == inosine salvage pathway modelled. 1 == not modelled + I(X)P concentrations fixed. 

% simulation and visualization setup
plotflag = 0; % 1 = visualize simulations for each cycle, 2 = visualize only concentrations and fluxes at the end of each cycle.
NumberCycles = 5;
setup.clamp_GLCec = 0; % 0 = extracellular glucose as dynamic variable, 1 = extracellular glucose concentration fixed to experimental data.
setup.csmin = 0.094; % 0 = no baseline concentration, 0.094 = most accurate baseline concentration
setup.ratio_decrease_sinkPYR = 0.1; % decrease in the sink of pyruvate, set to fit experimental data. This accounts for the lower fluxes (to TCA) for FF, compared to SS and GP.

% parameter set
load('pset_Y3M2.mat','x16d_initial','x16d_final')
x = x16d_final;

% simulation
[T_FF01_1,Y_FF01_1,V_FF01_1] = simulate_FF_enrichment(x,canelas_SS,data,dataset,setup,NumberCycles,plotflag,IC);
selResCell{1}.T_FF01 = T_FF01_1;
selResCell{1}.Y_FF01 = Y_FF01_1;
selResCell{1}.V_FF01 = V_FF01_1;

% % visualization (standard style)
% plotMode = 2; % multiple simulations
% referencePlotSimulations_enrichment
% plotMode = 0;

% visualization (manuscript style)
referencePlotSimulations_forManuscript




