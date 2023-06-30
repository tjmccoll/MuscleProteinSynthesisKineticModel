function [t, x, x_conc, expData_conc, expData_mass, CostParameters, interventionDuration, plasmaVolume, ...
    skeletalMuscleVolume] = runModel_OIM_230221_signallingPermissive(Protocol, ExpDataFile, SubjectData, ...
    InitialConditions, kValues, EqDur, RunDur, GlucInfRate, LeuInfRate, kValues_KD)
% This runModel function script selects the appropriate input data for the
% p70S6K permissive signaling protocol and simulates the model using the 
% appropriate dose and timing of leucine ingestion.
%   [t, x, x_conc, expData_conc, expData_mass, CostParameters, 
%       interventionDuration, plasmaVolume, skeletalMuscleVolume] = 
%       runModel_OIM_230221_signallingPermissive(Protocol, ExpDataFile, 
%       SubjectData, InitialConditions, kValues, EqDur, RunDur, 
%       GlucInfRate, LeuInfRate, kValues_KD)
%   Protocol = experimental protocol to simulate. This input allows for the
%   selection of appropriate experimental data
%   ExpDataFile = comprehensive excel document containing all
%   experimentally extracted data.
%   SubjectData = subject specific data for calculating plasma and skeletal 
%   muscle volume. SubjectData = [gender, height, mass, age, bia]
%   InitialConditions = vector of initial concentrations for model species
%   KValueVector_given = inputted vector of rate constants
%   EqDur = equilibrium duration to allow for model steady state (min)
%   RunDur = duration of model simulation (min)
%   GlucoseInfRate = Glucose infusion rate required for the Sturis module
%   LeucineInfRate_given = Leucine infusion rate required to maintain
%   leucine at a steady state during the equilibrium period
%   kValues_KD = adjusted set of k-values to simulate knockdown

%% Selecting inputs for each protocol

if isequal(Protocol, 'Calibration, young adults')
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [3, 4, 5, 21, 28, 41, 42, 6, 46]; %updated numbers

end

%% Calculating values for model simulation

% Calculating Compartment Volumes
plasmaVolume = plasmaVolumeFxn(SubjectData(1), SubjectData(2), SubjectData(3));
skeletalMuscleVolume = skeletalMuscleVolumeFxn(SubjectData(1), SubjectData(2), SubjectData(4), SubjectData(5));

% Initial Conditions - conversion to moles
x0 = InitialConditions(:,1);
x0_mass = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume);

LeucineMolarMass = 1/131.17; %mol/g
AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;

% Experimental Data - conversion to mass
if isequal(Protocol ,'Signaling knockdown') || isequal(Protocol, 'Basal p70S6K Analysis') || ...
    isequal(Protocol, 'p70S6K reduced, mTORC1 kinase increased')
    expData_conc = nan;
    expData_mass = nan;
else
    ExpData = xlsread(ExpDataFile, Protocol);
    expData_conc = ExpData(:,1:4);
    expData_mass = expData_massFxn(expData_conc, CostParameters, plasmaVolume, skeletalMuscleVolume);
end 

% K-Values
KValueVector_given = kValues(:,2);
KValueVector_factor(1:length(KValueVector_given)) = 10.^(floor(log10(KValueVector_given)));

OptimizeParameters = [];
ab = [];

%% Model Simulation
% Original modelDerivative is needed for the equilibrium phase to attain 
% equivalent post-equilibrium values. Then the signallingPermissive 
% modelDerivative function is used where p70S6K is maintained at the 
% post-equilibrium value

if isequal(Protocol, 'Calibration, young adults') 
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, KValueVector_given, OptimizeParameters, ...
        ab, KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), ...
        [0 :0.1: EqDur], x0_mass);
    
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;
    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220916_signallingPermissive(t2, x2, KValueVector_given, OptimizeParameters, ...
        ab, KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), ...
        [EqDur :0.1: (EqDur+RunDur)], x0_2);
                          
    x_bolus = [x; x2];
    t_bolus = [t; t2];
    
        x_bolus(:,41) = KValueVector_given(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,42) = KValueVector_given(9)*x_bolus(:,8)./x_bolus(:,10) + KValueVector_given(11)*x_bolus(:,7); % Fm,0; insulin mediated MPB
        x_bolus(:,43) = KValueVector_given(6)*x_bolus(:,4) - KValueVector_given(7)*x_bolus(:,6); % Net Balance
        x_bolus(:,44) = KValueVector_given(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
        x_bolus(:,45) = KValueVector_given(9)*x_bolus(:,8) ./ x_bolus(:,10); % MPB, r9 (p-IR concentration)
      
        % parameter 46 sets as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_bolus(:,46) = nan;       
        
    x = x_bolus;
    t = t_bolus;
end

%% Converting moles to concentration

x_conc = x_concFxn(x, x0_mass, plasmaVolume, skeletalMuscleVolume);

end