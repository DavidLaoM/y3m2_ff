function [error]=costfunSystemY3M1_FF_pRun(x_temp,canelas_SS,setup,xarray,data,dataset,n,IC0,selPars,warray)
% recall from parallel
setup.caseStudy.parameters = selPars;
setup.w = warray;

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

% % Metabolites
% 
% errorG6P = EYgpff.expG6P - EYgpff.simG6P;
% errorF6P = EYgpff.expF6P - EYgpff.simF6P;
% errorFBP = EYgpff.expFBP - EYgpff.simFBP;
% errorGAP = EYgpff.expGAP - EYgpff.simGAP;
% errorDHAP = EYgpff.expDHAP - EYgpff.simDHAP;
% error2PG = EYgpff.exp2PG - EYgpff.sim2PG;
% errorP3G = EYgpff.expPG3 - EYgpff.simPG3;
% errorPEP = EYgpff.expPEP - EYgpff.simPEP;
% 
% errorG1P = EYgpff.expG1P - EYgpff.simG1P;
% errorUDP_GLC = EYgpff.expUDP_GLC - EYgpff.simUDP_GLC;
% errorT6P = EYgpff.expT6P - EYgpff.simT6P;
% errorTREic = EYgpff.expTREic - EYgpff.simTREic;
% errorG3P = EYgpff.expG3P - EYgpff.simG3P;
% 
% errorGLCic = EYgpff.expGLCic - EYgpff.simGLCic;
% % errorGLYC = EYgpff.expGLYC - EYgpff.sim
% errorTREec = EYgpff.expTREec - EYgpff.simTREec;
% errorGLCec = EYgpff.expGLCec - EYgpff.simGLCec;
% 
% errorTREvac = EYgpff.expTREvac - EYgpff.simTREvac;

% Metabolites (Normalized to individual experimental error)

errorG6P = ( EYgpff.expG6P - EYgpff.simG6P ) ./ EYgpff.expG6P;
errorF6P = ( EYgpff.expF6P - EYgpff.simF6P ) ./ EYgpff.expF6P;
errorFBP = ( EYgpff.expFBP - EYgpff.simFBP ) ./ EYgpff.expFBP;
errorGAP = ( EYgpff.expGAP - EYgpff.simGAP ) ./ EYgpff.expGAP;
errorDHAP = ( EYgpff.expDHAP - EYgpff.simDHAP ) ./ EYgpff.expDHAP;
error2PG = ( EYgpff.exp2PG - EYgpff.sim2PG ) ./ EYgpff.exp2PG;
errorP3G = ( EYgpff.expPG3 - EYgpff.simPG3 ) ./ EYgpff.expPG3;
errorPEP = ( EYgpff.expPEP - EYgpff.simPEP ) ./ EYgpff.expPEP;

errorG1P = ( EYgpff.expG1P - EYgpff.simG1P ) ./ EYgpff.expG1P;
errorUDP_GLC = ( EYgpff.expUDP_GLC - EYgpff.simUDP_GLC ) ./ EYgpff.expUDP_GLC;
errorT6P = ( EYgpff.expT6P - EYgpff.simT6P ) ./ EYgpff.expT6P;
errorTREic = ( EYgpff.expTREic - EYgpff.simTREic ) ./ EYgpff.expTREic;
errorG3P = ( EYgpff.expG3P - EYgpff.simG3P ) ./ EYgpff.expG3P;

errorGLCic = ( EYgpff.expGLCic - EYgpff.simGLCic ) ./ EYgpff.expGLCic;
% errorGLYC = ( EYgpff.expGLYC - EYgpff.sim
errorTREec = ( EYgpff.expTREec - EYgpff.simTREec ) ./ EYgpff.expTREec;
errorGLCec = ( EYgpff.expGLCec - EYgpff.simGLCec ) ./ EYgpff.expGLCec;

errorTREvac = ( EYgpff.expTREvac - EYgpff.simTREvac ) ./ EYgpff.expTREvac;

%

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

% % Fluxes
% 
% errorPGM1 = EVgpff.expPGM1 - EVgpff.simPGM1;
% errorTPS1 = EVgpff.expTPS1 - EVgpff.simTPS1;
% errorTPS2 = EVgpff.expTPS2 - EVgpff.simTPS2;
% errorNTH1 = EVgpff.expNTH1 - EVgpff.simNTH1;
% errorUGP = EVgpff.expUGP - EVgpff.simUGP;
% errorAGT1 = EVgpff.expAGT1 - EVgpff.simAGT1;
% errorATH1ec = EVgpff.expATH1ec - EVgpff.simATH1ec;
% errorATH1vac = EVgpff.expATH1vac - EVgpff.simATH1vac;
% % errorvacuoleT = EVgpff.expvacuoleT - EVgpff.simvacuoleT;
% 
% errorGLT = EVgpff.expGLT - EVgpff.simGLT;
% errorGLK = EVgpff.expGLK - EVgpff.simGLK;
% errorPGI = EVgpff.expPGI - EVgpff.simPGI;
% errorPFK = EVgpff.expPFK - EVgpff.simPFK;

if isfield(setup,'fakeATHcostFuns')
    if setup.fakeATHcostFuns == 1
        EVgpff.expATH1ec = 0.75 * 1E-2 * [0.94 0.75 0.19 0.82 0.94]';
        EVgpff.expATH1vac = 1E-3 * [0.94 0.75 0.19 0.82 0.94]';
    else
    end
else
end

% Fluxes (Normalized to maximum experimental error)
errorPGM1 = ( EVgpff.expPGM1 - EVgpff.simPGM1 ) ./ EVgpff.expPGM1 ;
errorTPS1 = ( EVgpff.expTPS1 - EVgpff.simTPS1 ) ./ EVgpff.expTPS1([2 2 3 4 4]); %EVgpff.expTPS1 ;
errorTPS2 = ( EVgpff.expTPS2 - EVgpff.simTPS2 ) ./ EVgpff.expTPS2 ;
errorNTH1 = ( EVgpff.expNTH1 - EVgpff.simNTH1 ) ./ EVgpff.expNTH1([1 2 3 5 5]) ; %EVgpff.expNTH1 ;
errorUGP = ( EVgpff.expUGP - EVgpff.simUGP ) ./ EVgpff.expUGP ;
errorAGT1 = ( EVgpff.expAGT1 - EVgpff.simAGT1 ) ./ EVgpff.expAGT1 ;
errorATH1ec = ( EVgpff.expATH1ec - EVgpff.simATH1ec ) ./ EVgpff.expATH1ec ;
errorATH1vac = ( EVgpff.expATH1vac - EVgpff.simATH1vac ) ./ EVgpff.expATH1vac ;
% errorvacuoleT = ( EVgpff.expvacuoleT - EVgpff.simvacuoleT ) ./ EVgpff.expPGM1 ;

errorGLT = ( EVgpff.expGLT - EVgpff.simGLT ) ./ EVgpff.expGLT ;
errorGLK = ( EVgpff.expGLK - EVgpff.simGLK ) ./ EVgpff.expGLK ;
errorPGI = ( EVgpff.expPGI - EVgpff.simPGI ) ./ EVgpff.expPGI ;
errorPFK = ( EVgpff.expPFK - EVgpff.simPFK ) ./ EVgpff.expPFK ;

% errorr Glycogen
errorGlycogen = ( EYgpff.expGlycogen - EYgpff.simGlycogen) ./ EYgpff.expGlycogen;
errorGlycSynth = ( EVgpff.expGlySynth - EVgpff.simGlySynth ) ./ EVgpff.expGlySynth;
errorGlycDeg = ( EVgpff.expGlyDeg - EVgpff.simGlyDeg ) ./ EVgpff.expGlySynth;
%

errorRegularization = lambda * x_temp; %lambda*x_temp';

% % % % % if isfield(setup,'updated_bmf_Cx_ATH1ec')
% % % % %     if setup.updated_bmf_Cx_ATH1ec == 1
        % errorATH1ec
        temp_errorATH1ec = errorATH1ec;
        clear errorATH1ec
        errorATH1ec = temp_errorATH1ec([1,end]);
        % errorAGT1
        temp_errorAGT1 = errorAGT1;
        clear errorAGT1
        errorAGT1 = temp_errorAGT1([1,end]);
        % errorATH1vac
        temp_errorATH1vac = errorATH1vac;
        clear errorATH1vac
        errorATH1vac = temp_errorATH1vac([1,end]);
% % % % %     else
% % % % %     end
% % % % % else
% % % % % end


error = [...
    w(5).*errorG6P; % = EYgpff.expG6P - EYgpff.simG6P;
    w(4).*errorF6P; % = EYgpff.expF6P - EYgpff.simF6P;
    w(3).*errorFBP; % = EYgpff.expFBP - EYgpff.simFBP;
%     w(14).*errorGAP; % = EYgpff.expGAP - EYgpff.simGAP;
%     w(17).*errorDHAP; % = EYgpff.expDHAP - EYgpff.simDHAP;
%     w(10).*error2PG; % = EYgpff.exp2PG - EYgpff.sim2PG;
%     w(11).*errorP3G; % = EYgpff.expP3G - EYgpff.simP3G;
    w(12).*errorPEP; % = EYgpff.expPEP - EYgpff.simPEP;
%     w(13).*errorPYR; % = EYgpff.expPYR - EYgpff.simPYR;
    w(21).*errorG1P; % = EYgpff.expG1P - EYgpff.simG1P;
    w(24).*errorUDP_GLC; % = EYgpff.expUDP_GLC - EYgpff.simUDP_GLC;
    w(26).*errorT6P; % = EYgpff.expT6P - EYgpff.simT6P;
    w(25).*errorTREic; % = EYgpff.expTRE - EYgpff.simTRE;
%     w(18).*errorG3P; % = EYgpff.expG3P - EYgpff.simG3P;
    w(9).*errorATP; % = EYgpff.expATP - EYgpff.simATP;
%     w(15).*errorADP; % = EYgpff.expADP - EYgpff.simADP;
%     w(16).*errorAMP; % = EYgpff.expAMP - EYgpff.simAMP;
%     w(27).*errorPI; % = 10 - EYgpff.simAMP(end);
    w(6).*errorGLCic; % = EYgpff.exp_GLCic - EYgpff.simGLCic;
    w(37).*errorTREec; % = EYgpff.exp_TREec - EYgpff.simTREec;
    w(36).*errorGLCec; % = EYgpff.exp_GLCec - EYgpff.simGLCec;
    w(38).*errorTREvac; % = EYgpff.expTREvac - EYgpff.simTREvac;
    
    w(86).*errorGlycogen;

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

    w(87).*errorGlycSynth;
    w(88).*errorGlycDeg;

%     %     errorRegularization'
    ];


% %% visualizing the errors in the cost function
% tempNames_concs = cell(12,1);
% %
% tempNames_concs{1} = 'errorG6P';
% tempNames_concs{2} = 'errorF6P';
% tempNames_concs{3} = 'errorFBP';
% tempNames_concs{4} = 'errorPEP';
% tempNames_concs{5} = 'errorG1P';
% tempNames_concs{6} = 'errorUDP_GLC';
% %
% tempNames_concs{7} = 'errorT6P';
% tempNames_concs{8} = 'errorTREic';
% tempNames_concs{9} = 'errorGLCic';
% tempNames_concs{10} = 'errorTREec';
% tempNames_concs{11} = 'errorGLCec';
% tempNames_concs{12} = 'errorTREvac';
% 
% tempNames_rates = cell(12,1);
% %
% tempNames_rates{1} = 'errorPGM1';
% tempNames_rates{2} = 'errorTPS1';
% tempNames_rates{3} = 'errorTPS2';
% tempNames_rates{4} = 'errorNTH1';
% tempNames_rates{5} = 'errorUGP';
% tempNames_rates{6} = 'errorAGT1';
% %
% tempNames_rates{7} = 'errorATH1ec';
% tempNames_rates{8} = 'errorATH1vac';
% tempNames_rates{9} = 'errorGLT';
% tempNames_rates{10} = 'errorGLK';
% tempNames_rates{11} = 'errorPGI';
% tempNames_rates{12} = 'errorPFK';
%     
% figure
% %
% subplot(1,2,1), title('concentrations')
% plot(errorG6P,'.-')
% hold on
% plot(errorF6P,'.-')
% plot(errorFBP,'.-')
% plot(errorPEP,'.-')
% plot(errorG1P,'.-')
% plot(errorUDP_GLC,'.-')
% plot(errorT6P,'.-')
% plot(errorTREic,'.-')
% plot(errorGLCic,'.-')
% plot(errorTREec,'.-')
% plot(errorGLCec,'.-')
% plot(errorTREvac,'.-')
% %
% legend(tempNames_concs)
% 
% %
% subplot(1,2,2), title('fluxes')
% plot(errorPGM1,'.-')
% hold on
% 
% plot(errorTPS1,'.-')
% plot(errorTPS2,'.-')
% plot(errorNTH1,'.-')
% plot(errorUGP,'.-')
% plot(errorAGT1,'.-')
% 
% plot(errorATH1ec,'.-')
% plot(errorATH1vac,'.-')
% plot(errorGLT,'.-')
% plot(errorGLK,'.-')
% plot(errorPGI,'.-')
% 
% plot(errorPFK,'.-')
% %
% legend(tempNames_rates)
% %%


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

