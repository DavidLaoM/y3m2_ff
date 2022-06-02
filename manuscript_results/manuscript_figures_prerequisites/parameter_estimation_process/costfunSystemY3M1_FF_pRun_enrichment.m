function [error]=costfunSystemY3M1_FF_pRun_enrichment(x_temp,canelas_SS,setup,xarray,data,dataset,n,IC0,selPars,warray)
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
setup.experiment = 1;
[TgpFF,YgpFF,VgpFF] = simulate_FF_enrichment(xarray,canelas_SS,data,dataset,setup,n,plotflag,IC0);
simResults.TgpFF = TgpFF;
simResults.YgpFF = YgpFF;
simResults.VgpFF = VgpFF;
    % visualize
%     [legenda] = legendaFull; %legenda for the names needed
%     metNames = legenda.metabolites;
%     reactNames = legenda.fluxes;
%     legendaMetabolites_addEnrichment;
%     choosedataset
%     ExpData = setup.ExpData;
%     plotMode = 2; 
%     selResCell{1}.T_FF01 = simResults.TgpFF;
%     selResCell{1}.Y_FF01 = simResults.YgpFF;
%     selResCell{1}.V_FF01 = simResults.VgpFF;
%     referencePlotSimulations_enrichment
%     plotMode = 0;
%     figure,
%     for i = 1:38
%         subplot(7,6,i)
%         plot(TgpFF, YgpFF(:,i), 'b.')
%         hold on
%     end
%     % diplay inputs
%     disp(xarray), 
%     disp(canelas_SS), 
%     disp(data), 
%     disp(dataset), 
%     disp(setup), 
%     disp(n), 
%     disp(plotflag), 
%     disp(IC0)
    

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

% errors added in the final regularization
if(sum(abs(w(60:66))) ~= 0)
    
    % calculate errors
    % 
    temp_time = setup.simData.T_FF01;
    temp_hxk = setup.simData.V_FF01(:,2);
    temp_glt = setup.simData.V_FF01(:,1);
    temp_pgm1 = setup.simData.V_FF01(:,17);
    temp_tps1 = setup.simData.V_FF01(:,21);
    temp_tps2 = setup.simData.V_FF01(:,19);
    temp_nth1 = setup.simData.V_FF01(:,20);
    temp_g1p = setup.simData.Y_FF01(:,21);
    temp_t6p = setup.simData.Y_FF01(:,26);
    temp_udpglc = setup.simData.Y_FF01(:,24);
    % 
    temp_timepoints = [6:5:26, 76:50:376];
    % 
    temp_HXKexp = interp1(temp_time, temp_hxk, temp_timepoints, 'pchip');
    temp_GLTexp = interp1(temp_time, temp_glt, temp_timepoints, 'pchip');
    temp_PGM1exp = interp1(temp_time, temp_pgm1, temp_timepoints, 'pchip');
    temp_TPS1exp = interp1(temp_time, temp_tps1, temp_timepoints, 'pchip');
    temp_TPS2exp = interp1(temp_time, temp_tps2, temp_timepoints, 'pchip');
    temp_NTH1exp = interp1(temp_time, temp_nth1, temp_timepoints, 'pchip');
    temp_G1Pexp = interp1(temp_time, temp_g1p, temp_timepoints, 'pchip');
    temp_T6Pexp = interp1(temp_time, temp_t6p, temp_timepoints, 'pchip');
    temp_UDPGLCexp = interp1(temp_time, temp_udpglc, temp_timepoints, 'pchip');
    % 
    temp_HXKsim = interp1(TgpFF, VgpFF(:,2), temp_timepoints, 'pchip');
    temp_GLTsim = interp1(TgpFF, VgpFF(:,1), temp_timepoints, 'pchip');
    temp_PGM1sim = interp1(TgpFF, VgpFF(:,17), temp_timepoints, 'pchip');
    temp_TPS1sim = interp1(TgpFF, VgpFF(:,21), temp_timepoints, 'pchip');
    temp_TPS2sim = interp1(TgpFF, VgpFF(:,19), temp_timepoints, 'pchip');
    temp_NTH1sim = interp1(TgpFF, VgpFF(:,20), temp_timepoints, 'pchip');
    temp_G1Psim = interp1(TgpFF, YgpFF(:,21), temp_timepoints, 'pchip');
    temp_T6Psim = interp1(TgpFF, YgpFF(:,26), temp_timepoints, 'pchip');
    temp_UDPGLCsim = interp1(TgpFF, YgpFF(:,24), temp_timepoints, 'pchip');
    % 
    error_simHXK = ( temp_HXKexp - temp_HXKsim ) ./ temp_HXKexp;
    error_simGLT = ( temp_GLTexp - temp_GLTsim ) ./ temp_GLTexp;
    error_simPGM1 = ( temp_PGM1exp - temp_PGM1sim ) ./ temp_PGM1exp;
    error_simTPS1 = ( temp_TPS1exp - temp_TPS1sim ) ./ temp_TPS1exp;
    error_simTPS2 = ( temp_TPS2exp - temp_TPS2sim ) ./ temp_TPS2exp;
    error_simNTH1 = ( temp_NTH1exp - temp_NTH1sim ) ./ temp_NTH1exp;
    error_simG1P = ( temp_G1Pexp - temp_G1Psim ) ./ temp_G1Pexp;
    error_simT6P = ( temp_T6Pexp - temp_T6Psim ) ./ temp_T6Pexp;
    error_simUDPGlc = ( temp_UDPGLCexp - temp_UDPGLCsim ) ./ temp_UDPGLCexp;
    
    % calculate errors (finalReg)
%     error_finalReg = (x_temp - setup.pset_Y3M1(setup.selPars));  
    error_finalReg = (x_temp - setup.pset_Y3M1(selPars));    
    if(isfield(setup,'regularization_Smallbone2011')&&(setup.regularization_Smallbone2011 == 1))
        % 
        clear error_finalReg
%         % 
%         params_involved = ...
%             83 %[parsPGM1, parsPGM1 = 
%             84
%             85
%             86;
%             144 parsUGP, parsUGP = 144:148
%             145
%             146
%             147
%             148
%             124 parsTPS1, parsTPS1 = 124:128
%             125
%             126 
%             127
%             128            
%             119 parsTPS2, parsTPS2 = 119:121
%             120
%             121
%             122 parsNTH1] parsNTH1 = 122:123;
%             122
        parsCurrent = [10.^xarray(83).*1/6,...% 83 %[parsPGM1, parsPGM1 = 
            10.^xarray(84).*0.023,... % 84
            10.^xarray(85).*0.05,... % 85
            10.^xarray(86).*100.*(1+0+0),... % 86; (Vmax) p.PGM1_kcat = 10.^x(86).*100.*(f.PGM1+f.PGM2+f.PGM3);
            10.^xarray(144).* 0.11,... % 144 parsUGP, parsUGP = 144:148
            10.^xarray(145).* 0.11,... % 145
            10.^xarray(146).* 0.32,... % 146
            10.^xarray(147).* 0.035,... % 147
            10.^xarray(148).* 1000.*0.00237,... % 148 (Vmax) p.UGP_kcat = 10.^x(148).* 1000.*f.UGP;
            10.^xarray(124).*3.8,... % 124 parsTPS1, parsTPS1 = 124:128
            10.^xarray(125).*0.886,... % 125
            10.^xarray(126).*1000.*0.00145,... % 126 (Vmax) p.TPS1_kcat = 10.^x(126).*1000.*f.TPS1;%
            10.^xarray(127).*1,... % 127
            10.^xarray(128).*1,... % 128            
            10.^xarray(119).*0.5,... % 119 parsTPS2, parsTPS2 = 119:121
            10.^xarray(120).*81.45.*0.00133,... % 120 (Vmax) p.TPS2_kcat = 10.^x(120).*81.45.*f.TPS2;%
            10.^xarray(121).*1,... % 121
            10.^xarray(122).*2.99,... % 122 parsNTH1] parsNTH1 = 122:123;
            10.^xarray(123).*100.*1.5.*0.00196]; % 123 (Vmax) p.NTH1_kcat = 10.^x(123).*100 .* 1.5.*f.NTH1;
        parsSmallbone2011 = [1/6,...% 83 %[parsPGM1, parsPGM1 = 
            0.023,... % 84
            0.05,... % 85
            0.3545*60,... % 86; (Vmax)
            0.11,... % 144 parsUGP, parsUGP = 144:148
            0.11,... % 145
            0.32,... % 146
            0.0035,... % 147
            36.82*60,... % 148 (Vmax)
            3.8,... % 124 parsTPS1, parsTPS1 = 124:128
            0.886,... % 125
            1.371*60,... % 126 (Vmax)
            1,... % 127
            1,... % 128            
            0.5,... % 119 parsTPS2, parsTPS2 = 119:121
            6.5*60,... % 120 (Vmax)
            1,... % 121
            2.99,... % 122 parsNTH1] parsNTH1 = 122:123;
            15.2*60]; % 123 (Vmax)
        error_finalReg = log10(parsCurrent./parsSmallbone2011); % also regularizing to zero here.
        % reducing weight Vmax vs other parameters (set to 1/5)
        error_finalReg(4) = error_finalReg(4) / 5;
        error_finalReg(9) = error_finalReg(9) / 5;
        error_finalReg(12) = error_finalReg(12) / 5;
        error_finalReg(16) = error_finalReg(16) / 5;
        error_finalReg(19) = error_finalReg(19) / 5;
    end
end


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


% Adding error for the enrichment data (legenda.metabolites([36 65]))
    %
    time_exp = setup.ExpData.metabolites{65}.time;
    exp_ratio = setup.ExpData.metabolites{65}.fraction;
    % 
    sim_GLC_ec_total = interp1(TgpFF, YgpFF(:,36), time_exp, 'pchip');
    sim_GLC_ec_13C_l = interp1(TgpFF, YgpFF(:,65), time_exp, 'pchip');
    sim_ratio = sim_GLC_ec_13C_l ./ sim_GLC_ec_total;
    % 
    error_GLCenrich = (sim_ratio - exp_ratio) ./ exp_ratio; % check result looks like it should (numbers in plot)
    error = [error; error_GLCenrich * w(89)]; % check result looks like it should (shape)

    
% errors added in the final regularization
if(sum(abs(w(60:66))) ~= 0)
%     error_simHXK = ( temp_HXKexp - temp_HXKsim ) ./ temp_HXKexp;
%     error_simGLT = ( temp_GLTexp - temp_GLTsim ) ./ temp_GLTexp;
%     error_simPGM1 = ( temp_PGM1exp - temp_PGM1sim ) ./ temp_PGM1exp;
%     error_simTPS1 = ( temp_TPS1exp - temp_TPS1sim ) ./ temp_TPS1exp;
%     error_simTPS2 = ( temp_TPS2exp - temp_TPS2sim ) ./ temp_TPS2exp;
%     error_simNTH1 = ( temp_NTH1exp - temp_NTH1sim ) ./ temp_NTH1exp;
    error = [error;...
        w(60) .* error_simHXK';...
        w(61) .* error_finalReg';...
        w(62) .* error_simGLT';...
        w(63) .* error_simPGM1';...
        w(64) .* error_simTPS1';...
        w(65) .* error_simTPS2';...
        w(66) .* error_simNTH1';...
        w(67) .* error_simG1P';...
        w(68) .* error_simT6P';...
        w(69) .* error_simUDPGlc';...
        ];
%     k2 = find(error)
    % 311:324 for error_finalReg
    % 299:end since simHXK to UDPGlc, included
    
    
    
%     %
%     figure,
%     plot(setup.ExpData.metabolites{65}.time,...
%         setup.ExpData.metabolites{65}.fraction,...
%         'r+')
%     hold on
%     plot(setup.ExpData.metabolites{65}.time,...
%         sim_ratio,...
%         'r-')
%     % [TgpFF,YgpFF,VgpFF]
%     figure,
%     for i = 1:38
%         subplot(7,6,i)
%         plot(TgpFF, YgpFF(:,i), 'b.')
%         hold on
%     end
%     disp('stop here.')
    
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

    % 
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

% (2021-10-06) Use an inreased error if there is NaN
if any(isnan(error))
    error(isnan(error)) = 100;
end


end
