if setup.GPdataset.GP400WT == 1
    ExpData = dataFF01;
elseif setup.GPdataset.GP1800WT == 1
    ExpData = dataFF03;
elseif setup.GPdataset.GP400M == 1
    ExpData = dataFF04;
elseif setup.GPdataset.Fructose == 1
    ExpData = dataFFFruc;
elseif setup.GPdataset.Sucrose == 1
    ExpData = dataFFSuc;
end