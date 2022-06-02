function [error]=costfunSystemY3M1_FF2(x_temp,canelas_SS,setup,xarray,data,dataset,n,IC0)

% preliminary
% xref = xarray;
xarray(setup.caseStudy.parameters) = x_temp;
lambda = setup.parEst.lambda; 
w = setup.w;

% simulation
plotflag = 0; 
[TgpFF,YgpFF,VgpFF] = simulate_FF(xarray,canelas_SS,data,dataset,setup,n,plotflag,IC0);
simResults.TgpFF = TgpFF;
simResults.YgpFF = YgpFF;
simResults.VgpFF = VgpFF;


% error calculation 
% [EYss, EVss, EYgp, EVgp, EYgpff, EVgpff] = calcErrorY3M1_FF(simResults,canelas_SS,data,dataset,setup);
[EYgpff, EVgpff] = ErrorCalc_FF(simResults,canelas_SS,data,dataset,setup);

% Metabolites

errorG6P = EYgpff.expG6P - EYgpff.simG6P;
errorF6P = EYgpff.expF6P - EYgpff.simF6P;
errorFBP = EYgpff.expFBP - EYgpff.simFBP;
errorGAP = EYgpff.expGAP - EYgpff.simGAP;
errorDHAP = EYgpff.expDHAP - EYgpff.simDHAP;
error2PG = EYgpff.exp2PG - EYgpff.sim2PG;
errorP3G = EYgpff.expPG3 - EYgpff.simPG3;
errorPEP = EYgpff.expPEP - EYgpff.simPEP;

errorG1P = EYgpff.expG1P - EYgpff.simG1P;
errorUDP_GLC = EYgpff.expUDP_GLC - EYgpff.simUDP_GLC;
errorT6P = EYgpff.expT6P - EYgpff.simT6P;
errorTREic = EYgpff.expTREic - EYgpff.simTREic;
errorG3P = EYgpff.expG3P - EYgpff.simG3P;

errorTREic = EYgpff.expTREic - EYgpff.simTREic;
errorGLCic = EYgpff.expGLCic - EYgpff.simGLCic;
% errorGLYC = EYgpff.expGLYC - EYgpff.sim
errorTREec = EYgpff.expTREec - EYgpff.simTREec;
errorGLCec = EYgpff.expGLCec - EYgpff.simGLCec;

if ((setup.GPdataset.Sucrose == 1))
    errorPYR = 0;
    errorATP = 0;
    errorADP = 0;
    errorAMP = 0;
    errorPI = 0;
    
else
    errorPYR = EYgpff.expPYR - EYgpff.simPYR;
    errorATP = EYgpff.expATP - EYgpff.simATP;
    errorADP = EYgpff.expADP - EYgpff.simADP;
    errorAMP = EYgpff.expAMP - EYgpff.simAMP;
    errorPI = 10 - EYgpff.simAMP(end);
end

% Fluxes

errorPGM1 = EVgpff.expPGM1 - EVgpff.simPGM1;
errorTPS1 = EVgpff.expTPS1 - EVgpff.simTPS1;
errorTPS2 = EVgpff.expTPS2 - EVgpff.simTPS2;
errorNTH1 = EVgpff.expNTH1 - EVgpff.simNTH1;
errorUGP = EVgpff.expUGP - EVgpff.simUGP;
errorAGT1 = EVgpff.expAGT1 - EVgpff.simAGT1;
errorATH1ec = EVgpff.expATH1ec - EVgpff.simATH1ec;
errorATH1vac = EVgpff.expATH1vac - EVgpff.simATH1vac;
% errorvacuoleT = EVgpff.expvacuoleT - EVgpff.simvacuoleT;

errorGLT = EVgpff.expGLT - EVgpff.simGLT;
errorGLK = EVgpff.expGLK - EVgpff.simGLK;
errorPGI = EVgpff.expPGI - EVgpff.simPGI;
errorPFK = EVgpff.expPFK - EVgpff.simPFK;

errorRegularization = lambda * x_temp; %lambda*x_temp';


error = [...
    w(5).*errorG6P; % = EYgpff.expG6P - EYgpff.simG6P;
    w(4).*errorF6P; % = EYgpff.expF6P - EYgpff.simF6P;
    w(3).*errorFBP; % = EYgpff.expFBP - EYgpff.simFBP;
%     w(14).*errorGAP; % = EYgpff.expGAP - EYgpff.simGAP;
%     w(17).*errorDHAP; % = EYgpff.expDHAP - EYgpff.simDHAP;
%     w(10).*error2PG; % = EYgpff.exp2PG - EYgpff.sim2PG;
%     w(11).*errorP3G; % = EYgpff.expP3G - EYgpff.simP3G;
%     w(12).*errorPEP; % = EYgpff.expPEP - EYgpff.simPEP;
%     w(13).*errorPYR; % = EYgpff.expPYR - EYgpff.simPYR;
    w(21).*errorG1P; % = EYgpff.expG1P - EYgpff.simG1P;
    w(24).*errorUDP_GLC; % = EYgpff.expUDP_GLC - EYgpff.simUDP_GLC;
    w(26).*errorT6P; % = EYgpff.expT6P - EYgpff.simT6P;
    w(25).*errorTREic; % = EYgpff.expTRE - EYgpff.simTRE;
%     w(18).*errorG3P; % = EYgpff.expG3P - EYgpff.simG3P;
%     w(9).*errorATP; % = EYgpff.expATP - EYgpff.simATP;
%     w(15).*errorADP; % = EYgpff.expADP - EYgpff.simADP;
%     w(16).*errorAMP; % = EYgpff.expAMP - EYgpff.simAMP;
%     w(27).*errorPI; % = 10 - EYgpff.simAMP(end);
%     w(6).*errorGLCic; % = EYgpff.exp_GLCic - EYgpff.simGLCic;
    w(37).*errorTREec; % = EYgpff.exp_TREec - EYgpff.simTREec;
    w(36).*errorGLCec; % = EYgpff.exp_GLCec - EYgpff.simGLCec;
    
    % Fluxes:
    
    w(55).*errorPGM1; % = EVgpff.expPGM1 - EVgpff.simPGM1;
    w(59).*errorTPS1; % = EVgpff.expTPS1 - EVgpff.simTPS1;
    w(57).*errorTPS2; % = EVgpff.expTPS2 - EVgpff.simTPS2;
    w(58).*errorNTH1; % = EVgpff.expNTH1 - EVgpff.simNTH1;
    w(56).*errorUGP; % = EVgpff.expUGP - EVgpff.simUGP;
    w(85).*errorAGT1; % = EVgpff.expAGT1 - EVgpff.simAGT1;
    w(83).*errorATH1ec; % = EVgpff.expATH1ec - EVgpff.simATH1ec;
    w(84).*errorATH1vac; % = EVgpff.expATH1vac - EVgpff.simATH1vac;
    w(39).*errorGLT; % = EVgpff.expGLT - EVgpff.simGLT;
    w(40).*errorGLK; % = EVgpff.expGLK - EVgpff.simGLK;
    w(41).*errorPGI; % = EVgpff.expPGI - EVgpff.simPGI;
    w(42).*errorPFK; % = EVgpff.expPFK - EVgpff.simPFK;
%     %     errorRegularization'
    ];

if isfield(setup, 'visualizeCostFun')
    if setup.visualizeCostFun == 1
        % visualize error
        time_Y = dataset.FF01.time_mets;
        time_Y_glcEC = dataset.FF01.timeECgluc;
        time_V = dataset.FF01.fluxes_times;
        figure
        subplot(1,3,1), plot(time_Y, EYgpff.simGLCic', '.-'), hold on, plot(time_Y, EYgpff.expGLCic', 'ko'), title('GLCic');
        subplot(1,3,2), plot(time_Y, EYgpff.simTREec', '.-'), hold on, plot(time_Y, EYgpff.expTREec', 'ko'), title('TREec');
        subplot(1,3,3), plot(time_Y, EYgpff.simGLCec', '.-'), hold on, plot(time_Y, EYgpff.expGLCec', 'ko'), title('GLCec');
        
    end
    
    
end



end

