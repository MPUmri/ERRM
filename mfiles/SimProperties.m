function [simProp] = SimProperties()
    % General properties for the simulation

    % Pharmacokinetics of reference tissue
    % Units are min^{-1} for Ktrans, while EES is unitless
    % Walker-Samuel et al. (2007), PMB, 52(1), 75–89. doi:10.1088/0031-9155/52/1/006
    simProp.KtRR = 0.07;
    simProp.veRR = 0.14;
    simProp.kepRR = simProp.KtRR/simProp.veRR;
    
    simProp.sigmaC = 0:0.01:0.05; % Range of stdDev of noise in concentration
    simProp.nVox = 10000; % Number of replications for each CNR
    simProp.TRes = [1,5,10,15]; % Temporal resolutions (in seconds)
    simProp.tDuration = 10; % Duration of DCE Acquisition (in minutes)
    simProp.initTRes = 0.1; % Temporal resolution of initially simulated data

    % Pharmacokinetics of tumour tissue
    simProp.Kt = 0.25; % /min
    simProp.ve = 0.4; % unitless (fractional)
    simProp.kep = simProp.Kt/simProp.ve; % /min
    simProp.vp = [0.001 0.005 0.05:0.05:0.1]; % unitless (fractional)
    
end