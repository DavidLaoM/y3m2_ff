% % reorganiseTUDdata
% gets the data loaded from the file 'TUDdata.mat' and reorders so that it
% matches the legend being used in this study.
dataFF01.metabolites = cell(38,1);
dataFF03.metabolites = cell(38,1);
dataFF04.metabolites = cell(38,1);
dataFF01.fluxes = cell(48,1);
dataFF03.fluxes = cell(48,1);
dataFF04.fluxes = cell(48,1);


%% % % % DATASET FF01
dataFF01.time_mets = dataset.FF01.time_mets;
dataFF01.timeECgluc = dataset.FF01.timeECgluc;
dataFF01.fluxes_times = dataset.FF01.fluxes_times;
% 
%     METABOLITES
%     {'ACE, x_{1}, [mM]'        }
%     {'BPG, x_{2}, [mM]'        }
dataFF01.metabolites{3} = dataset.FF01.metabolites.ICFBP;%     {'F16BP, x_{3}, [mM]'      }
dataFF01.metabolites{4} = dataset.FF01.metabolites.ICF6P;%     {'F6P, x_{4}, [mM]'        }
dataFF01.metabolites{5} = dataset.FF01.metabolites.ICG6P;%     {'G6P, x_{5}, [mM]'        }
dataFF01.metabolites{6} = dataset.FF01.metabolites.ICglucose;%     {'GLCi, x_{6}, [mM]'       }
%     {'NAD, x_{7}, [mM]'        }
%     {'NADH, x_{8}, [mM]'       }
dataFF01.metabolites{9} = dataset.FF01.metabolites.ICATP;%     {'ATP, x_{9}, [mM]'        }
dataFF01.metabolites{10} = dataset.FF01.metabolites.IC2PG;%     {'P2G, x_{10}, [mM]'       }
dataFF01.metabolites{11} = dataset.FF01.metabolites.IC3PG;%     {'P3G, x_{11}, [mM]'       }
dataFF01.metabolites{12} = dataset.FF01.metabolites.ICPEP;%     {'PEP, x_{12}, [mM]'       }
dataFF01.metabolites{13} = dataset.FF01.metabolites.ICPYR;%     {'PYR, x_{13}, [mM]'       }
dataFF01.metabolites{14} = dataset.FF01.metabolites.ICGAP;%     {'GLYCERAL3P, x_{14}, [mM]'}
dataFF01.metabolites{15} = dataset.FF01.metabolites.ICADP;%     {'ADP, x_{15}, [mM]'       }
dataFF01.metabolites{16} = dataset.FF01.metabolites.ICAMP;%     {'AMP, x_{16}, [mM]'       }
dataFF01.metabolites{17} = dataset.FF01.metabolites.ICDHAP;%     {'DHAP, x_{17}, [mM]'      }
dataFF01.metabolites{18} = dataset.FF01.metabolites.ICG3P;%     {'GLYC3P, x_{18}, [mM]'    }
%     {'GLYCEROL, x_{19}, [mM]'  }
%     {'ETOH, x_{20}, [mM]'      }
dataFF01.metabolites{21} = dataset.FF01.metabolites.ICG1P;%     {'G1P, x_{21}, [mM]'       }
%     {'UTP, x_{22}, [mM]'       }
%     {'UDP, x_{23}, [mM]'       }
dataFF01.metabolites{24} = dataset.FF01.metabolites.ICUDPG;%     {'UDP_GLC, x_{24}, [mM]'   }

dataFF01.metabolites{25} = dataset.FF01.metabolites.ICtreh;%     {'TRE, x_{25}, [mM]'       }
dataFF01.metabolites{25}.conc = 0.1.* dataset.FF01.metabolites.ICtreh.conc;
dataFF01.metabolites{25}.stdev = 0.1.* dataset.FF01.metabolites.ICtreh.stdev;
dataFF01.metabolites{25}.sterr = 0.1.* dataset.FF01.metabolites.ICtreh.sterr;

dataFF01.metabolites{26} = dataset.FF01.metabolites.ICT6P;%     {'T6P, x_{26}, [mM]'       }
%     {'PI, x_{27}, [mM]'        }
%     {'IMP, x_{28}, [mM]'       }
%     {'INO, x_{29}, [mM]'       }
%     {'HYP, x_{30}, [mM]'       }
%     {'ETOH_{EC}, x_{31}, [mM]' }
%     {'GLYC_{EC}, x_{32}, [mM]' }
%     {'FRC_{IC}, x_{33}, [mM]' }
%     {'FRC_{EC}, x_{34}, [mM]' }
%     {'SUC_{EC}, x_{35}, [mM]' }
dataFF01.metabolites{36} = dataset.FF01.metabolites.ECglucose;%     {'GLC_{EC}, x_{36}, [mM]' }
dataFF01.metabolites{37} = dataset.FF01.metabolites.ECtreh;%     {'TRE_{EC}, x_{37}, [mM]' }

dataFF01.metabolites{38} = dataset.FF01.metabolites.ICtreh; %     {'TRE_{vac}, x_{38}, [mM]' }
dataFF01.metabolites{38}.conc = 0.9.* dataset.FF01.metabolites.ICtreh.conc;
dataFF01.metabolites{38}.stdev = 0.9.* dataset.FF01.metabolites.ICtreh.stdev;
dataFF01.metabolites{38}.sterr = 0.9.* dataset.FF01.metabolites.ICtreh.sterr;


%     FLUXES
for i = 1:48
    dataFF01.fluxes{i}.time = dataset.FF01.fluxes_times';
end
dataFF01.fluxes{1}.rate = dataset.FF01.fluxes{1}';%     {'v_{GLT}, r_{1}, [mM s^{-1}]'               }
dataFF01.fluxes{2}.rate = dataset.FF01.fluxes{1}';%     {'v_{GLK}, r_{2}, [mM s^{-1}]'               }
dataFF01.fluxes{3}.rate = dataset.FF01.fluxes{2}' - dataset.FF01.fluxes{3}';%     {'v_{PGI}, r_{3}, [mM s^{-1}]'               }
dataFF01.fluxes{4}.rate = dataset.FF01.fluxes{4}';%     {'v_{PFK}, r_{4}, [mM s^{-1}]'               }
dataFF01.fluxes{5}.rate = dataset.FF01.fluxes{5}';%     {'V_{ALD}, r_{5}, [mM s^{-1}]'               }
%     {'V_{G3PDH}, r_{6}, [mM s^{-1}]'             }
%     {'V_{TPI1}, r_{7}, [mM s^{-1}]'              }
dataFF01.fluxes{8}.rate = dataset.FF01.fluxes{6}' - dataset.FF01.fluxes{7}';%     {'V_{TDH1}, r_{8}, [mM s^{-1}]'              }
dataFF01.fluxes{9}.rate = dataset.FF01.fluxes{6}' - dataset.FF01.fluxes{7}';%     {'V_{PGK}, r_{9}, [mM s^{-1}]'               }
dataFF01.fluxes{10}.rate = dataset.FF01.fluxes{8}';%     {'V_{PGM}, r_{10}, [mM s^{-1}]'              }
dataFF01.fluxes{11}.rate = dataset.FF01.fluxes{8}';%     {'V_{ENO}, r_{11}, [mM s^{-1}]'              }
dataFF01.fluxes{12}.rate = dataset.FF01.fluxes{9}';%     {'V_{PYK}, r_{12}, [mM s^{-1}]'              }
dataFF01.fluxes{13}.rate = dataset.FF01.fluxes{10}';%     {'V_{PDC}, r_{13}, [mM s^{-1}]'              }
%     {'V_{ADK1}, r_{14}, [mM s^{-1}]'             }
%     {'V_{HOR2}, r_{15}, [mM s^{-1}]'             }
%     {'V_{RHR2}, r_{16}, [mM s^{-1}]'             }
dataFF01.fluxes{17}.rate = dataset.FF01.fluxes{39}' - dataset.FF01.fluxes{40}';%     {'V_{PGM1}, r_{17}, [mM s^{-1}]'             }
dataFF01.fluxes{18}.rate = dataset.FF01.fluxes{41}';%     {'V_{UGP}, r_{18}, [mM s^{-1}]'              }
dataFF01.fluxes{19}.rate = dataset.FF01.fluxes{12}' - dataset.FF01.fluxes{13}';%     {'V_{TPS2}, r_{19}, [mM s^{-1}]'             }
dataFF01.fluxes{20}.rate = dataset.FF01.fluxes{14}';%     {'V_{NTH1}, r_{20}, [mM s^{-1}]'             }
dataFF01.fluxes{21}.rate = dataset.FF01.fluxes{11}';%     {'V_{TPS1}, r_{21}, [mM s^{-1}]'             }
%     {'V_{ACE}, r_{22}, [mM s^{-1}]'              }
%     {'V_{ETOHtransport}, r_{23}, [mM s^{-1}]'    }
%     {'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]'}
%     {'V_{PDH}, r_{25}, [mM s^{-1}]'              }
%     {'V_{mitoNADH}, r_{26}, [mM s^{-1}]'         }
%     {'V_{ATPase}, r_{27}, [mM s^{-1}]'           }
%     {'V_{mito}, r_{28}, [mM s^{-1}]'             }
%     {'V_{Amd1}, r_{29}, [mM s^{-1}]'             }
%     {'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]'    }
%     {'V_{Isn1}, r_{31}, [mM s^{-1}]'             }
%     {'V_{Pnp1}, r_{32}, [mM s^{-1}]'             }
%     {'V_{Hpt1}, r_{33}, [mM s^{-1}]'             }
%     {'V_{sinkG6P}, r_{34}, [mM s^{-1}]'          }
%     {'V_{sinkF6P}, r_{35}, [mM s^{-1}]'          }
%     {'V_{sinkGAP}, r_{36}, [mM s^{-1}]'          }
%     {'V_{sinkP3G}, r_{37}, [mM s^{-1}]'          }
%     {'V_{sinkPEP}, r_{38}, [mM s^{-1}]'          }
%     {'V_{sinkPYR}, r_{39}, [mM s^{-1}]'          }
%     {'V_{sinkACE}, r_{40}, [mM s^{-1}]'          }
%     {'V_{ADH}, r_{41}, [mM s^{-1}]'              }
%     {'V_{vacuolePi}, r_{42}, [mM s^{-1}]'        }
%     {V_{FRT}, r_{43}, [mM s^{-1}]'               }
%     {V_{FRK}, r_{44}, [mM s^{-1}]'               }
%     V_{ATH1} total  %   dataFF01.fluxes{45}.rate = dataset.FF01.fluxes{34}';
% dataFF01.fluxes{45}.rate = 0.1.* dataset.FF01.fluxes{34}';%     {V_{ATH1ec}, r_{45}, [mM s^{-1}]'            }
% dataFF01.fluxes{46}.rate = 0.9.* dataset.FF01.fluxes{34}';%     {V_{ATH1vac}, r_{46}, [mM s^{-1}]'           }
dataFF01.fluxes{45}.rate = 0.9.* dataset.FF01.fluxes{34}';%     {V_{ATH1ec}, r_{45}, [mM s^{-1}]'            }
dataFF01.fluxes{46}.rate = 0.1.* dataset.FF01.fluxes{34}';%     {V_{ATH1vac}, r_{46}, [mM s^{-1}]'           }

dataFF01.fluxes{47}.rate = dataset.FF01.fluxes{32}' - dataset.FF01.fluxes{33}';%     {V_{AGT1}, r_{47}, [mM s^{-1}]'              }
%     {V_{vacuoleT}, r_{48}, [mM s^{-1}]'          }
% dataFF01.fluxes{49}.rate = dataset.FF01.fluxes{42}';%     {V_{GSY}, r_{49}, [mM s^{-1}]'          }
% dataFF01.fluxes{50}.rate = dataset.FF01.fluxes{43}' - dataset.FF01.fluxes{44}';%     {V_{GPH1}, r_{50}, [mM s^{-1}]'          }
% dataFF01.fluxes{51}.rate = dataset.FF01.fluxes{45}';%     {V_{Gly_vacdeg}, r_{51}, [mM s^{-1}]'          }

%% % % % DATASET FF03
% 
%     METABOLITES
%     {'ACE, x_{1}, [mM]'        }
%     {'BPG, x_{2}, [mM]'        }
dataFF03.metabolites{3} = dataset.FF03.metabolites.ICFBP;%     {'F16BP, x_{3}, [mM]'      }
dataFF03.metabolites{4} = dataset.FF03.metabolites.ICF6P;%     {'F6P, x_{4}, [mM]'        }
dataFF03.metabolites{5} = dataset.FF03.metabolites.ICG6P;%     {'G6P, x_{5}, [mM]'        }
dataFF03.metabolites{6} = dataset.FF03.metabolites.ICglucose;%     {'GLCi, x_{6}, [mM]'       }
%     {'NAD, x_{7}, [mM]'        }
%     {'NADH, x_{8}, [mM]'       }
dataFF03.metabolites{9} = dataset.FF03.metabolites.ICATP;%     {'ATP, x_{9}, [mM]'        }
dataFF03.metabolites{10} = dataset.FF03.metabolites.IC2PG;%     {'P2G, x_{10}, [mM]'       }
dataFF03.metabolites{11} = dataset.FF03.metabolites.IC3PG;%     {'P3G, x_{11}, [mM]'       }
dataFF03.metabolites{12} = dataset.FF03.metabolites.ICPEP;%     {'PEP, x_{12}, [mM]'       }
dataFF03.metabolites{13} = dataset.FF03.metabolites.ICPYR;%     {'PYR, x_{13}, [mM]'       }
dataFF03.metabolites{14} = dataset.FF03.metabolites.ICGAP;%     {'GLYCERAL3P, x_{14}, [mM]'}
dataFF03.metabolites{15} = dataset.FF03.metabolites.ICADP;%     {'ADP, x_{15}, [mM]'       }
dataFF03.metabolites{16} = dataset.FF03.metabolites.ICAMP;%     {'AMP, x_{16}, [mM]'       }
dataFF03.metabolites{17} = dataset.FF03.metabolites.ICDHAP;%     {'DHAP, x_{17}, [mM]'      }
dataFF03.metabolites{18} = dataset.FF03.metabolites.ICG3P;%     {'GLYC3P, x_{18}, [mM]'    }
%     {'GLYCEROL, x_{19}, [mM]'  }
%     {'ETOH, x_{20}, [mM]'      }
dataFF03.metabolites{21} = dataset.FF03.metabolites.ICG1P;%     {'G1P, x_{21}, [mM]'       }
%     {'UTP, x_{22}, [mM]'       }
%     {'UDP, x_{23}, [mM]'       }
dataFF03.metabolites{24} = dataset.FF03.metabolites.ICUDPG;%     {'UDP_GLC, x_{24}, [mM]'   }
dataFF03.metabolites{25} = dataset.FF03.metabolites.ICtreh;%     {'TRE, x_{25}, [mM]'       }
dataFF03.metabolites{26} = dataset.FF03.metabolites.ICT6P;%     {'T6P, x_{26}, [mM]'       }
%     {'PI, x_{27}, [mM]'        }
%     {'IMP, x_{28}, [mM]'       }
%     {'INO, x_{29}, [mM]'       }
%     {'HYP, x_{30}, [mM]'       }
%     {'ETOH_{EC}, x_{31}, [mM]' }
dataFF03.metabolites{36} = dataset.FF03.metabolites.ECglucose;%     {'GLYC_{EC}, x_{32}, [mM]' }

%     FLUXES
for i = 1:48
    dataFF03.fluxes{i}.time = dataset.FF03.fluxes_times';
end
dataFF03.fluxes{1}.rate = dataset.FF03.fluxes{1}';%     {'v_{GLT}, r_{1}, [mM s^{-1}]'               }
dataFF03.fluxes{2}.rate = dataset.FF03.fluxes{1}';%     {'v_{GLK}, r_{2}, [mM s^{-1}]'               }
dataFF03.fluxes{3}.rate = dataset.FF03.fluxes{2}' - dataset.FF03.fluxes{3}';%     {'v_{PGI}, r_{3}, [mM s^{-1}]'               }
dataFF03.fluxes{4}.rate = dataset.FF03.fluxes{4}';%     {'v_{PFK}, r_{4}, [mM s^{-1}]'               }
dataFF03.fluxes{5}.rate = dataset.FF03.fluxes{5}';%     {'V_{ALD}, r_{5}, [mM s^{-1}]'               }
%     {'V_{G3PDH}, r_{6}, [mM s^{-1}]'             }
%     {'V_{TPI1}, r_{7}, [mM s^{-1}]'              }
dataFF03.fluxes{8}.rate = dataset.FF03.fluxes{6}' - dataset.FF03.fluxes{7}';%     {'V_{TDH1}, r_{8}, [mM s^{-1}]'              }
dataFF03.fluxes{9}.rate = dataset.FF03.fluxes{6}' - dataset.FF03.fluxes{7}';%     {'V_{PGK}, r_{9}, [mM s^{-1}]'               }
dataFF03.fluxes{10}.rate = dataset.FF03.fluxes{8}';%     {'V_{PGM}, r_{10}, [mM s^{-1}]'              }
dataFF03.fluxes{11}.rate = dataset.FF03.fluxes{8}';%     {'V_{ENO}, r_{11}, [mM s^{-1}]'              }
dataFF03.fluxes{12}.rate = dataset.FF03.fluxes{9}';%     {'V_{PYK}, r_{12}, [mM s^{-1}]'              }
dataFF03.fluxes{13}.rate = dataset.FF03.fluxes{10}';%     {'V_{PDC}, r_{13}, [mM s^{-1}]'              }
%     {'V_{ADK1}, r_{14}, [mM s^{-1}]'             }
%     {'V_{HOR2}, r_{15}, [mM s^{-1}]'             }
%     {'V_{RHR2}, r_{16}, [mM s^{-1}]'             }
dataFF03.fluxes{17}.rate = dataset.FF03.fluxes{39}' - dataset.FF03.fluxes{40}';%     {'V_{PGM1}, r_{17}, [mM s^{-1}]'             }
%     {'V_{UGP}, r_{18}, [mM s^{-1}]'              }
dataFF03.fluxes{19}.rate = dataset.FF03.fluxes{12}' - dataset.FF03.fluxes{13}';%     {'V_{TPS2}, r_{19}, [mM s^{-1}]'             }
dataFF03.fluxes{20}.rate = dataset.FF03.fluxes{14}';%     {'V_{NTH1}, r_{20}, [mM s^{-1}]'             }
dataFF03.fluxes{21}.rate = dataset.FF03.fluxes{11}';%     {'V_{TPS1}, r_{21}, [mM s^{-1}]'             }
%     {'V_{ACE}, r_{22}, [mM s^{-1}]'              }
%     {'V_{ETOHtransport}, r_{23}, [mM s^{-1}]'    }
%     {'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]'}
%     {'V_{PDH}, r_{25}, [mM s^{-1}]'              }
%     {'V_{mitoNADH}, r_{26}, [mM s^{-1}]'         }
%     {'V_{ATPase}, r_{27}, [mM s^{-1}]'           }
%     {'V_{mito}, r_{28}, [mM s^{-1}]'             }
%     {'V_{Amd1}, r_{29}, [mM s^{-1}]'             }
%     {'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]'    }
%     {'V_{Isn1}, r_{31}, [mM s^{-1}]'             }
%     {'V_{Pnp1}, r_{32}, [mM s^{-1}]'             }
%     {'V_{Hpt1}, r_{33}, [mM s^{-1}]'             }
%     {'V_{sinkG6P}, r_{34}, [mM s^{-1}]'          }
%     {'V_{sinkF6P}, r_{35}, [mM s^{-1}]'          }
%     {'V_{sinkGAP}, r_{36}, [mM s^{-1}]'          }
%     {'V_{sinkP3G}, r_{37}, [mM s^{-1}]'          }
%     {'V_{sinkPEP}, r_{38}, [mM s^{-1}]'          }
%     {'V_{sinkPYR}, r_{39}, [mM s^{-1}]'          }
%     {'V_{sinkACE}, r_{40}, [mM s^{-1}]'          }
%     {'V_{ADH}, r_{41}, [mM s^{-1}]'              }
%     {'V_{vacuolePi}, r_{42}, [mM s^{-1}]'        }


% % % % DATASET FF04
% 
%     METABOLITES
%     {'ACE, x_{1}, [mM]'        }
%     {'BPG, x_{2}, [mM]'        }
dataFF04.metabolites{3} = dataset.FF04.metabolites.ICFBP;%     {'F16BP, x_{3}, [mM]'      }
dataFF04.metabolites{4} = dataset.FF04.metabolites.ICF6P;%     {'F6P, x_{4}, [mM]'        }
dataFF04.metabolites{5} = dataset.FF04.metabolites.ICG6P;%     {'G6P, x_{5}, [mM]'        }
dataFF04.metabolites{6} = dataset.FF04.metabolites.ICglucose;%     {'GLCi, x_{6}, [mM]'       }
%     {'NAD, x_{7}, [mM]'        }
%     {'NADH, x_{8}, [mM]'       }
dataFF04.metabolites{9} = dataset.FF04.metabolites.ICATP;%     {'ATP, x_{9}, [mM]'        }
dataFF04.metabolites{10} = dataset.FF04.metabolites.IC2PG;%     {'P2G, x_{10}, [mM]'       }
dataFF04.metabolites{11} = dataset.FF04.metabolites.IC3PG;%     {'P3G, x_{11}, [mM]'       }
dataFF04.metabolites{12} = dataset.FF04.metabolites.ICPEP;%     {'PEP, x_{12}, [mM]'       }
dataFF04.metabolites{13} = dataset.FF04.metabolites.ICPYR;%     {'PYR, x_{13}, [mM]'       }
dataFF04.metabolites{14} = dataset.FF04.metabolites.ICGAP;%     {'GLYCERAL3P, x_{14}, [mM]'}
dataFF04.metabolites{15} = dataset.FF04.metabolites.ICADP;%     {'ADP, x_{15}, [mM]'       }
dataFF04.metabolites{16} = dataset.FF04.metabolites.ICAMP;%     {'AMP, x_{16}, [mM]'       }
dataFF04.metabolites{17} = dataset.FF04.metabolites.ICDHAP;%     {'DHAP, x_{17}, [mM]'      }
dataFF04.metabolites{18} = dataset.FF04.metabolites.ICG3P;%     {'GLYC3P, x_{18}, [mM]'    }
%     {'GLYCEROL, x_{19}, [mM]'  }
%     {'ETOH, x_{20}, [mM]'      }
dataFF04.metabolites{21} = dataset.FF04.metabolites.ICG1P;%     {'G1P, x_{21}, [mM]'       }
%     {'UTP, x_{22}, [mM]'       }
%     {'UDP, x_{23}, [mM]'       }
dataFF04.metabolites{24} = dataset.FF04.metabolites.ICUDPG;%     {'UDP_GLC, x_{24}, [mM]'   }
dataFF04.metabolites{25} = dataset.FF04.metabolites.ICtreh;%     {'TRE, x_{25}, [mM]'       }
dataFF04.metabolites{26} = dataset.FF04.metabolites.ICT6P;%     {'T6P, x_{26}, [mM]'       }
%     {'PI, x_{27}, [mM]'        }
%     {'IMP, x_{28}, [mM]'       }
%     {'INO, x_{29}, [mM]'       }
%     {'HYP, x_{30}, [mM]'       }
%     {'ETOH_{EC}, x_{31}, [mM]' }
dataFF04.metabolites{36} = dataset.FF04.metabolites.ECglucose;%     {'GLYC_{EC}, x_{32}, [mM]' }

%     FLUXES
for i = 1:48
    dataFF04.fluxes{i}.time = dataset.FF04.fluxes_times';
end
% % dataFF04.fluxes{1}.rate = dataset.FF04.fluxes{1}';%     {'v_{GLT}, r_{1}, [mM s^{-1}]'               }
% % dataFF04.fluxes{2}.rate = dataset.FF04.fluxes{1}';%     {'v_{GLK}, r_{2}, [mM s^{-1}]'               }
% % dataFF04.fluxes{3}.rate = dataset.FF04.fluxes{2}' - dataset.FF04.fluxes{3}';%     {'v_{PGI}, r_{3}, [mM s^{-1}]'               }
% % dataFF04.fluxes{4}.rate = dataset.FF04.fluxes{4}';%     {'v_{PFK}, r_{4}, [mM s^{-1}]'               }
% % dataFF04.fluxes{5}.rate = dataset.FF04.fluxes{5}';%     {'V_{ALD}, r_{5}, [mM s^{-1}]'               }
%     {'V_{G3PDH}, r_{6}, [mM s^{-1}]'             }
%     {'V_{TPI1}, r_{7}, [mM s^{-1}]'              }
% % dataFF04.fluxes{8}.rate = dataset.FF04.fluxes{6}' - dataset.FF04.fluxes{7}';%     {'V_{TDH1}, r_{8}, [mM s^{-1}]'              }
% % dataFF04.fluxes{9}.rate = dataset.FF04.fluxes{6}' - dataset.FF04.fluxes{7}';%     {'V_{PGK}, r_{9}, [mM s^{-1}]'               }
% % dataFF04.fluxes{10}.rate = dataset.FF04.fluxes{8}';%     {'V_{PGM}, r_{10}, [mM s^{-1}]'              }
% % dataFF04.fluxes{11}.rate = dataset.FF04.fluxes{8}';%     {'V_{ENO}, r_{11}, [mM s^{-1}]'              }
% % dataFF04.fluxes{12}.rate = dataset.FF04.fluxes{9}';%     {'V_{PYK}, r_{12}, [mM s^{-1}]'              }
% % dataFF04.fluxes{13}.rate = dataset.FF04.fluxes{10}';%     {'V_{PDC}, r_{13}, [mM s^{-1}]'              }
%     {'V_{ADK1}, r_{14}, [mM s^{-1}]'             }
%     {'V_{HOR2}, r_{15}, [mM s^{-1}]'             }
%     {'V_{RHR2}, r_{16}, [mM s^{-1}]'             }
% % dataFF04.fluxes{17}.rate = dataset.FF04.fluxes{39}' - dataset.FF04.fluxes{40}';%     {'V_{PGM1}, r_{17}, [mM s^{-1}]'             }
%     {'V_{UGP}, r_{18}, [mM s^{-1}]'              }
% % dataFF04.fluxes{19}.rate = dataset.FF04.fluxes{12}' - dataset.FF04.fluxes{13}';%     {'V_{TPS2}, r_{19}, [mM s^{-1}]'             }
% % dataFF04.fluxes{20}.rate = dataset.FF04.fluxes{14}';%     {'V_{NTH1}, r_{20}, [mM s^{-1}]'             }
% % dataFF04.fluxes{21}.rate = dataset.FF04.fluxes{11}';%     {'V_{TPS1}, r_{21}, [mM s^{-1}]'             }
%     {'V_{ACE}, r_{22}, [mM s^{-1}]'              }
%     {'V_{ETOHtransport}, r_{23}, [mM s^{-1}]'    }
%     {'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]'}
%     {'V_{PDH}, r_{25}, [mM s^{-1}]'              }
%     {'V_{mitoNADH}, r_{26}, [mM s^{-1}]'         }
%     {'V_{ATPase}, r_{27}, [mM s^{-1}]'           }
%     {'V_{mito}, r_{28}, [mM s^{-1}]'             }
%     {'V_{Amd1}, r_{29}, [mM s^{-1}]'             }
%     {'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]'    }
%     {'V_{Isn1}, r_{31}, [mM s^{-1}]'             }
%     {'V_{Pnp1}, r_{32}, [mM s^{-1}]'             }
%     {'V_{Hpt1}, r_{33}, [mM s^{-1}]'             }
%     {'V_{sinkG6P}, r_{34}, [mM s^{-1}]'          }
%     {'V_{sinkF6P}, r_{35}, [mM s^{-1}]'          }
%     {'V_{sinkGAP}, r_{36}, [mM s^{-1}]'          }
%     {'V_{sinkP3G}, r_{37}, [mM s^{-1}]'          }
%     {'V_{sinkPEP}, r_{38}, [mM s^{-1}]'          }
%     {'V_{sinkPYR}, r_{39}, [mM s^{-1}]'          }
%     {'V_{sinkACE}, r_{40}, [mM s^{-1}]'          }
%     {'V_{ADH}, r_{41}, [mM s^{-1}]'              }
%     {'V_{vacuolePi}, r_{42}, [mM s^{-1}]'        }



%%
%     METABOLITES
%     {'ACE, x_{1}, [mM]'        }
%     {'BPG, x_{2}, [mM]'        }
%     {'F16BP, x_{3}, [mM]'      }
%     {'F6P, x_{4}, [mM]'        }
%     {'G6P, x_{5}, [mM]'        }
%     {'GLCi, x_{6}, [mM]'       }
%     {'NAD, x_{7}, [mM]'        }
%     {'NADH, x_{8}, [mM]'       }
%     {'ATP, x_{9}, [mM]'        }
%     {'P2G, x_{10}, [mM]'       }
%     {'P3G, x_{11}, [mM]'       }
%     {'PEP, x_{12}, [mM]'       }
%     {'PYR, x_{13}, [mM]'       }
%     {'GLYCERAL3P, x_{14}, [mM]'}
%     {'ADP, x_{15}, [mM]'       }
%     {'AMP, x_{16}, [mM]'       }
%     {'DHAP, x_{17}, [mM]'      }
%     {'GLYC3P, x_{18}, [mM]'    }
%     {'GLYCEROL, x_{19}, [mM]'  }
%     {'ETOH, x_{20}, [mM]'      }
%     {'G1P, x_{21}, [mM]'       }
%     {'UTP, x_{22}, [mM]'       }
%     {'UDP, x_{23}, [mM]'       }
%     {'UDP_GLC, x_{24}, [mM]'   }
%     {'TRE, x_{25}, [mM]'       }
%     {'T6P, x_{26}, [mM]'       }
%     {'PI, x_{27}, [mM]'        }
%     {'IMP, x_{28}, [mM]'       }
%     {'INO, x_{29}, [mM]'       }
%     {'HYP, x_{30}, [mM]'       }
%     {'ETOH_{EC}, x_{31}, [mM]' }
%     {'GLYC_{EC}, x_{32}, [mM]' }
% 'FRC_{IC}, x_{33}, [mM]';
%     'FRC_{EC}, x_{34}, [mM]';
%     'SUC_{EC}, x_{35}, [mM]';
%     'GLC_{EC}, x_{36}, [mM]';
%     'TRE_{EC}, x_{37}, [mM]';
%     'TRE_{vac}, x_{38}, [mM]'};

%     FLUXES
%     {'v_{GLT}, r_{1}, [mM s^{-1}]'               }
%     {'v_{GLK}, r_{2}, [mM s^{-1}]'               }
%     {'v_{PGI}, r_{3}, [mM s^{-1}]'               }
%     {'v_{PFK}, r_{4}, [mM s^{-1}]'               }
%     {'V_{ALD}, r_{5}, [mM s^{-1}]'               }
%     {'V_{G3PDH}, r_{6}, [mM s^{-1}]'             }
%     {'V_{TPI1}, r_{7}, [mM s^{-1}]'              }
%     {'V_{TDH1}, r_{8}, [mM s^{-1}]'              }
%     {'V_{PGK}, r_{9}, [mM s^{-1}]'               }
%     {'V_{PGM}, r_{10}, [mM s^{-1}]'              }
%     {'V_{ENO}, r_{11}, [mM s^{-1}]'              }
%     {'V_{PYK}, r_{12}, [mM s^{-1}]'              }
%     {'V_{PDC}, r_{13}, [mM s^{-1}]'              }
%     {'V_{ADK1}, r_{14}, [mM s^{-1}]'             }
%     {'V_{HOR2}, r_{15}, [mM s^{-1}]'             }
%     {'V_{RHR2}, r_{16}, [mM s^{-1}]'             }
%     {'V_{PGM1}, r_{17}, [mM s^{-1}]'             }
%     {'V_{UGP}, r_{18}, [mM s^{-1}]'              }
%     {'V_{TPS2}, r_{19}, [mM s^{-1}]'             }
%     {'V_{NTH1}, r_{20}, [mM s^{-1}]'             }
%     {'V_{TPS1}, r_{21}, [mM s^{-1}]'             }
%     {'V_{ACE}, r_{22}, [mM s^{-1}]'              }
%     {'V_{ETOHtransport}, r_{23}, [mM s^{-1}]'    }
%     {'V_{GLYCEROLtransport}, r_{24}, [mM s^{-1}]'}
%     {'V_{PDH}, r_{25}, [mM s^{-1}]'              }
%     {'V_{mitoNADH}, r_{26}, [mM s^{-1}]'         }
%     {'V_{ATPase}, r_{27}, [mM s^{-1}]'           }
%     {'V_{mito}, r_{28}, [mM s^{-1}]'             }
%     {'V_{Amd1}, r_{29}, [mM s^{-1}]'             }
%     {'V_{Ade13_v_Ade12}, r_{30}, [mM s^{-1}]'    }
%     {'V_{Isn1}, r_{31}, [mM s^{-1}]'             }
%     {'V_{Pnp1}, r_{32}, [mM s^{-1}]'             }
%     {'V_{Hpt1}, r_{33}, [mM s^{-1}]'             }
%     {'V_{sinkG6P}, r_{34}, [mM s^{-1}]'          }
%     {'V_{sinkF6P}, r_{35}, [mM s^{-1}]'          }
%     {'V_{sinkGAP}, r_{36}, [mM s^{-1}]'          }
%     {'V_{sinkP3G}, r_{37}, [mM s^{-1}]'          }
%     {'V_{sinkPEP}, r_{38}, [mM s^{-1}]'          }
%     {'V_{sinkPYR}, r_{39}, [mM s^{-1}]'          }
%     {'V_{sinkACE}, r_{40}, [mM s^{-1}]'          }
%     {'V_{ADH}, r_{41}, [mM s^{-1}]'              }
%     {'V_{vacuolePi}, r_{42}, [mM s^{-1}]'        }
% 'V_{FRT}, r_{43}, [mM s^{-1}]';
%     'V_{FRK}, r_{44}, [mM s^{-1}]';
%     'V_{ATH1ec}, r_{45}, [mM s^{-1}]';
%     'V_{ATH1vac}, r_{46}, [mM s^{-1}]';
%     'V_{AGT1}, r_{47}, [mM s^{-1}]';
%     'V_{vacuoleT}, r_{48}, [mM s^{-1}]'};

