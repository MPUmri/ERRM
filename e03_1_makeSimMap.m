% Simulated scenario with a range of physiological parameters.
% This sceipt simulates curves with a combination of ktrans & ve  & vp
% and saves it for future steps

% Estimated runtime: < 1 second

%% Pre-setup

clearvars
addpath('./mfiles')
tic

%% Build map of true values

% Dimensions of map
% (These are actually the other way around in the figure)
nX = 20; % vp and ve
nY = 6; % Ktrans

% Values for parameters
valKt = [0.05,0.15,0.25,0.35,0.45,0.55]; % Ktrans - 6 elements
valVe = repmat([0.1,0.2,0.3,0.4,0.5], 1, 4); % ve - 5 elements, repeat 4 times
valVp = repmat([0.005,0.01,0.05,0.1]', 1, 5)'; % vp - 4 elements, repeat 5 times
% Flatten the vp values (initially matrix because thats how they were defined)
valVp = valVp(:); 

% Produce a map of the true values
trueKt = repmat(valKt,[nX 1]);
trueVe = repmat(valVe',[1 nY]);
trueVp = repmat(valVp,[1 nY]);
trueKep = trueKt./trueVe;

%% Simulate the concentration-time curves

initTRes = 0.1; % Initial temporal resolution for simulated data

% Define time
t=initTRes:initTRes:600; % in seconds
t=t'/60; % convert to minutes

% Use literature-based AIF with bolus arrival at 60s
Cp = ParkerAif(t,t(60/initTRes));

% Initialize array
simCt = zeros(length(t),nX,nY);

% Simulate the tissue of interest for each parameter combination
for i=1:nX
    for j=1:nY
        simCt(:,i,j)=ToftsKety(Cp,[trueKt(i,j) trueKep(i,j) trueVp(i,j)],t,1);
    end
end

% Save the simulated data - it'll be used in future steps
save('./data/simMap.mat');

%%
toc
disp('Finished simulating map of perfusion parameters')