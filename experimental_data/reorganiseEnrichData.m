% % reorganiseTUDdata
% gets the data loaded from the file 'TUDdata.mat' and reorders so that it
% matches the legend being used in this study.
% dataFF01.metabolites = cell(66,1);
% dataset.FF01.time_mets = [0;5;10;15;20;25;30;60;90;120;150;180;220;250;300;350;400];
% dataset.FF01.timeECgluc = [0;5;11;15;20;   30;60;90;    150;180;220;250;300;350;400];
    
    
fractionsFF01.metabolites = cell(28,1);


% 'ACE^{L}, f_{1}';
%     'BPG^{L}, f_{2}';
fractionsFF01.metabolites{3} = datasetEnrich.FF01.metabolites.ICFBP; %     'F16BP^{L}, f_{3}';
fractionsFF01.metabolites{4} = datasetEnrich.FF01.metabolites.ICF6P;%     'F6P^{L}, f_{4}';
fractionsFF01.metabolites{5} = datasetEnrich.FF01.metabolites.ICG6P;%   'G6P^{L}, f_{5}';
fractionsFF01.metabolites{6} = datasetEnrich.FF01.metabolites.ICglucose;%     'GLCi^{L}, f_{6}';
%     'P2G^{L}, f_{7}';
%     'P3G^{L}, f_{8}';
fractionsFF01.metabolites{9} = datasetEnrich.FF01.metabolites.ICPEP;%     'PEP^{L}, f_{9}';
fractionsFF01.metabolites{10} = datasetEnrich.FF01.metabolites.ICPYR;%     'PYR^{L}, f_{10}';
%     'GLYCERAL3P^{L}, f_{11}';
fractionsFF01.metabolites{12} = datasetEnrich.FF01.metabolites.ICDHAP;%     'DHAP^{L}, f_{12}';
%     'GLYC3P^{L}, f_{13}';
%     'GLYCEROL^{L}, f_{14}';
%     'ETOH^{L}, f_{15}';
fractionsFF01.metabolites{16} = datasetEnrich.FF01.metabolites.ICG1P;%     'G1P^{L}, f_{16}';
fractionsFF01.metabolites{17} = datasetEnrich.FF01.metabolites.ICUDPG;%     'UDP_{GLC}^{L}, f_{17}';
fractionsFF01.metabolites{18} = datasetEnrich.FF01.metabolites.ICtreh;%     'TRE_{IC}^{L}, f_{18}';
fractionsFF01.metabolites{19} = datasetEnrich.FF01.metabolites.ICT6P;%     'T6P^{L}, f_{19}';
%     'ETOH_{EC}^{L}, f_{20}';
%     'GLYC_{EC}^{L}, f_{21}';
%     'FRC_{IC}^{L}, f_{22}';
%     'FRC_{EC}^{L}, f_{23}';
%     'SUC_{EC}^{L}, f_{24}';
fractionsFF01.metabolites{25} = datasetEnrich.FF01.metabolites.ECglucose;%     'GLC_{EC}^{L}, f_{25}';
fractionsFF01.metabolites{26} = datasetEnrich.FF01.metabolites.ECtreh;%     'TRE_{EC}^{L}, f_{26}';
%     'TRE_{cyt}^{L}, f_{27}'};
%     'TRE_{vac}^{L}, f_{28}'};

%% % % % DATASET FF01
% dataFF01.time_mets = dataset.FF01.time_mets;
% dataFF01.timeECgluc = dataset.FF01.timeECgluc;

% 
%     METABOLITES
% %     {'ACE, x_{1}, [mM]'        }
% %     {'BPG, x_{2}, [mM]'        }
% dataFF01.metabolites{3} = dataset.FF01.metabolites.ICFBP;%     {'F16BP, x_{3}, [mM]'      }
% dataFF01.metabolites{4} = dataset.FF01.metabolites.ICF6P;%     {'F6P, x_{4}, [mM]'        }
% dataFF01.metabolites{5} = dataset.FF01.metabolites.ICG6P;%     {'G6P, x_{5}, [mM]'        }
% dataFF01.metabolites{6} = dataset.FF01.metabolites.ICglucose;%     {'GLCi, x_{6}, [mM]'       }
% %     {'NAD, x_{7}, [mM]'        }
% %     {'NADH, x_{8}, [mM]'       }
% dataFF01.metabolites{9} = dataset.FF01.metabolites.ICATP;%     {'ATP, x_{9}, [mM]'        }
% dataFF01.metabolites{10} = dataset.FF01.metabolites.IC2PG;%     {'P2G, x_{10}, [mM]'       }
% dataFF01.metabolites{11} = dataset.FF01.metabolites.IC3PG;%     {'P3G, x_{11}, [mM]'       }
% dataFF01.metabolites{12} = dataset.FF01.metabolites.ICPEP;%     {'PEP, x_{12}, [mM]'       }
% dataFF01.metabolites{13} = dataset.FF01.metabolites.ICPYR;%     {'PYR, x_{13}, [mM]'       }
% dataFF01.metabolites{14} = dataset.FF01.metabolites.ICGAP;%     {'GLYCERAL3P, x_{14}, [mM]'}
% dataFF01.metabolites{15} = dataset.FF01.metabolites.ICADP;%     {'ADP, x_{15}, [mM]'       }
% dataFF01.metabolites{16} = dataset.FF01.metabolites.ICAMP;%     {'AMP, x_{16}, [mM]'       }
% dataFF01.metabolites{17} = dataset.FF01.metabolites.ICDHAP;%     {'DHAP, x_{17}, [mM]'      }
% dataFF01.metabolites{18} = dataset.FF01.metabolites.ICG3P;%     {'GLYC3P, x_{18}, [mM]'    }
% %     {'GLYCEROL, x_{19}, [mM]'  }
% %     {'ETOH, x_{20}, [mM]'      }
% dataFF01.metabolites{21} = dataset.FF01.metabolites.ICG1P;%     {'G1P, x_{21}, [mM]'       }
% %     {'UTP, x_{22}, [mM]'       }
% %     {'UDP, x_{23}, [mM]'       }
% dataFF01.metabolites{24} = dataset.FF01.metabolites.ICUDPG;%     {'UDP_GLC, x_{24}, [mM]'   }
% 
% dataFF01.metabolites{25} = dataset.FF01.metabolites.ICtreh;%     {'TRE, x_{25}, [mM]'       }
% dataFF01.metabolites{25}.conc = 0.1.* dataset.FF01.metabolites.ICtreh.conc;
% dataFF01.metabolites{25}.stdev = 0.1.* dataset.FF01.metabolites.ICtreh.stdev;
% dataFF01.metabolites{25}.sterr = 0.1.* dataset.FF01.metabolites.ICtreh.sterr;
% 
% dataFF01.metabolites{26} = dataset.FF01.metabolites.ICT6P;%     {'T6P, x_{26}, [mM]'       }
% %     {'PI, x_{27}, [mM]'        }
% %     {'IMP, x_{28}, [mM]'       }
% %     {'INO, x_{29}, [mM]'       }
% %     {'HYP, x_{30}, [mM]'       }
% %     {'ETOH_{EC}, x_{31}, [mM]' }
% %     {'GLYC_{EC}, x_{32}, [mM]' }
% %     {'FRC_{IC}, x_{33}, [mM]' }
% %     {'FRC_{EC}, x_{34}, [mM]' }
% %     {'SUC_{EC}, x_{35}, [mM]' }
% dataFF01.metabolites{36} = dataset.FF01.metabolites.ECglucose;%     {'GLC_{EC}, x_{36}, [mM]' }
% dataFF01.metabolites{37} = dataset.FF01.metabolites.ECtreh;%     {'TRE_{EC}, x_{37}, [mM]' }
% 
% dataFF01.metabolites{38} = dataset.FF01.metabolites.ICtreh; %     {'TRE_{vac}, x_{38}, [mM]' }
% dataFF01.metabolites{38}.conc = 0.9.* dataset.FF01.metabolites.ICtreh.conc;
% dataFF01.metabolites{38}.stdev = 0.9.* dataset.FF01.metabolites.ICtreh.stdev;
% dataFF01.metabolites{38}.sterr = 0.9.* dataset.FF01.metabolites.ICtreh.sterr;

%   'ACE^{L}, x_{41}, [mM]';
%     'BPG^{L}, x_{42}, [mM]';
dataFF01.metabolites{43} = datasetEnrich.FF01.metabolites.ICFBP; % NOT THERE %     'F16BP^{L}, x_{42}, [mM]';
dataFF01.metabolites{44} = datasetEnrich.FF01.metabolites.ICF6P;%     'F6P^{L}, x_{43}, [mM]';
dataFF01.metabolites{45} = datasetEnrich.FF01.metabolites.ICG6P;%     'G6P^{L}, x_{44}, [mM]';
dataFF01.metabolites{46} = datasetEnrich.FF01.metabolites.ICglucose;%     'GLCi^{L}, x_{45}, [mM]';
%     'P2G^{L}, x_{47}, [mM]';
%     'P3G^{L}, x_{48}, [mM]'; % (MISSING) %
dataFF01.metabolites{49} = datasetEnrich.FF01.metabolites.ICPEP;%     'PEP^{L}, x_{48}, [mM]';
dataFF01.metabolites{50} = datasetEnrich.FF01.metabolites.ICPYR; % NOT THERE %     'PYR^{L}, x_{49}, [mM]';
%     'GLYCERAL3P^{L}, x_{51}, [mM]';
dataFF01.metabolites{52} = datasetEnrich.FF01.metabolites.ICDHAP;%     'DHAP^{L}, x_{51}, [mM]';
%     'GLYC3P^{L}, x_{53}, [mM]';
%     'GLYCEROL^{L}, x_{54}, [mM]';
%     'ETOH^{L}, x_{55}, [mM]';
dataFF01.metabolites{56} = datasetEnrich.FF01.metabolites.ICG1P;%     'G1P^{L}, x_{55}, [mM]';
dataFF01.metabolites{57} = datasetEnrich.FF01.metabolites.ICUDPG;%     'UDP_{GLC}^{L}, x_{56}, [mM]';
dataFF01.metabolites{58} = datasetEnrich.FF01.metabolites.ICtreh;%     'TRE_{IC}^{L}, x_{57}, [mM]';
dataFF01.metabolites{59} = datasetEnrich.FF01.metabolites.ICT6P;%     'T6P^{L}, x_{58}, [mM]';
%     'ETOH_{EC}^{L}, x_{60}, [mM]';
%     'GLYC_{EC}^{L}, x_{61}, [mM]';
%     'FRC_{IC}^{L}, x_{62}, [mM]';
%     'FRC_{EC}^{L}, x_{63}, [mM]';
%     'SUC_{EC}^{L}, x_{64}, [mM]';
dataFF01.metabolites{65} = datasetEnrich.FF01.metabolites.ECglucose;%     'GLC_{EC}^{L}, x_{64}, [mM]';
dataFF01.metabolites{66} = datasetEnrich.FF01.metabolites.ECtreh; % NOT THERE %     'TRE_{EC}^{L}, x_{65}, [mM]';
dataFF01.metabolites{67} = cell(1,1); %     'TRE_{vac}^{L}, x_{67}, [mM]'};


%% Direct visual extraction of % (MISSING) % values
% P3G
measTime = dataset.FF01.metabolites.IC3PG.time;
measConc = dataset.FF01.metabolites.IC3PG.conc;
tempExtract = [0, 0.9199690488034094;
10.810810810810807, 3.672967445973555;
16.21621621621621, 17.779251644282326;
20.180180180180187, 30.352898911181114;
30.990990990990994, 45.988780191234156;
60.18018018018016, 64.37130381915658;
89.72972972972973, 77.84557563698667;
121.8018018018018, 80.58171668601115;
150.6306306306306, 82.09335101973139;
180.18018018018017, 80.84369645719337;
221.26126126126127, 78.05145636433978;
249.72972972972974, 77.7228762504836;
302.70270270270277, 73.69452274360248;
350.990990990991, 69.66976178632619;
400, 67.4843862267175];
tempTime = tempExtract(:,1);
tempFrac = tempExtract(:,2)/100;
datasetEnrich.FF01.metabolites.ICP3G.time = tempTime;
datasetEnrich.FF01.metabolites.ICP3G.fraction = tempFrac;
datasetEnrich.FF01.metabolites.ICP3G.conc = tempFrac.*interp1(measTime,measConc,tempTime,'pchip');
dataFF01.metabolites{48} = datasetEnrich.FF01.metabolites.ICP3G; 
ExpData.metabolites{48} = dataFF01.metabolites{48};
% % PYR
% measTime = dataset.FF01.metabolites.ICPYR.time;
% measConc = dataset.FF01.metabolites.ICPYR.conc;
% tempTime = 1;
% tempFrac = 1;
% datasetEnrich.FF01.metabolites.ICPYR.time = tempTime;
% datasetEnrich.FF01.metabolites.ICPYR.fraction = tempFrac;
% datasetEnrich.FF01.metabolites.ICPYR.conc = tempFrac.*interp1(measTime,measConc,tempTime,'pchip');
% dataFF01.metabolites{50} = datasetEnrich.FF01.metabolites.ICPYR;
% % FBP
% measTime = dataset.FF01.metabolites.ICFBP.time;
% measConc = dataset.FF01.metabolites.ICFBP.conc;
% tempTime = 1;
% tempFrac = 1;
% datasetEnrich.FF01.metabolites.ICFBP.time = tempTime;
% datasetEnrich.FF01.metabolites.ICFBP.fraction = tempFrac;
% datasetEnrich.FF01.metabolites.ICFBP.conc = tempFrac.*interp1(measTime,measConc,tempTime,'pchip');
% dataFF01.metabolites{43} = datasetEnrich.FF01.metabolites.ICFBP; 
% % TREec
% measTime = dataset.FF01.metabolites.ECtreh.time;
% measConc = dataset.FF01.metabolites.ECtreh.conc;
% tempTime = 1;
% tempFrac = 1;
% datasetEnrich.FF01.metabolites.ECTRE.time = tempTime;
% datasetEnrich.FF01.metabolites.ECTRE.fraction = tempFrac;
% datasetEnrich.FF01.metabolites.ECTRE.conc = tempFrac.*interp1(measTime,measConc,tempTime,'pchip');
% dataFF01.metabolites{66} = datasetEnrich.FF01.metabolites.ECTRE; 
% TREic
% % measTime = dataset.FF01.metabolites.ICtreh.time;
% % measConc = dataset.FF01.metabolites.ICtreh.conc;
% % tempTime = 1;
% % tempFrac = 1;
% % datasetEnrich.FF01.metabolites.ICTRE.time = tempTime;
% % datasetEnrich.FF01.metabolites.ICTRE.fraction = tempFrac;
% % datasetEnrich.FF01.metabolites.ICTRE.conc = tempFrac.*interp1(measTime,measConc,tempTime,'pchip');
% % dataFF01.metabolites{68} = dataFF01.metabolites{58}; 
 

