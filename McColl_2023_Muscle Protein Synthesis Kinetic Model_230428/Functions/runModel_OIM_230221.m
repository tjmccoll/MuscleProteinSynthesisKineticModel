function [t, x, x_conc, expData_conc, expData_mass, CostParameters, interventionDuration, plasmaVolume, ...
    skeletalMuscleVolume] = runModel_OIM_230221(Protocol, ExpDataFile, SubjectData, InitialConditions, ...
    kValues, EqDur, RunDur, GlucInfRate, LeuInfRate, kValues_KD)
% This runModel function script selects the appropriate input data for the
% selected protocol and simulates the model using the appropriate dose and
% timing of leucine ingestion.
%   [t, x, x_conc, expData_conc, expData_mass, CostParameters, 
%       interventionDuration, plasmaVolume, skeletalMuscleVolume] = 
%       runModel_OIM_230221(Protocol, ExpDataFile, SubjectData, 
%       InitialConditions, kValues, EqDur, RunDur, GlucInfRate, LeuInfRate, 
%       kValues_KD)
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

% Calibration protocol -----
if isequal(Protocol, 'Calibration, young adults')
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [3, 4, 5, 21, 28, 41, 42, 6, 46]; 

% Validation protocols -----
elseif isequal(Protocol, 'Valid., Glynn (2010) 1.85g') 
    interventionDuration = 180; % min
    AminoAcidInput = 1.85; % grams
    CostParameters = [3, 4, 21, 28, 41, 42, 46]; 

elseif isequal(Protocol, 'Valid., Mitchell (2015) bolus') 
    interventionDuration = 240; % min
    AminoAcidInput = 3.59; % grams
    CostParameters = [3, 4, 21, 28, 46];

elseif isequal(Protocol, 'Valid., Mitchell (2015) pulse') 
    interventionDuration = 240; % min
    AminoAcidInput = 3.59 /4; % grams
    IngestionTimes = [0, 45, 90, 135];
    CostParameters = [3, 4, 21, 28, 46];
    
elseif isequal(Protocol, 'Valid., Dickinson (2011) cntl')
    interventionDuration = 120; % min
    AminoAcidInput = 1.8; % grams
    CostParameters = [4, 6, 28, 46];
    
elseif isequal(Protocol, 'Valid., Drummond (2010)')
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [4, 6, 41];

elseif isequal(Protocol, 'Valid., Wilkinson (2013)')
    interventionDuration = 150; % min
    AminoAcidInput = 3.42; % grams
    CostParameters = [3, 4, 6, 21, 28, 46];

% Leucine measurements (additional to validation) -----
    % 48 g whey isolate: 78% protein cotent in why (Gorissen, 2018) and 11%
    % leucine content of total protein (Gorissen, 2018). ~4.1 g leucine 
    % content
elseif isequal(Protocol, 'Atherton (2010) - bolus') % 48 g whey protein
    interventionDuration = 240; % experimental intervention was 360 min but this is beyond the model limits
    AminoAcidInput = 4.1; % grams
    CostParameters = [3, 4];

    % 25 g whey isolate: 78% protein cotent in why (Gorissen, 2018) and 11%
    % leucine content of total protein (Gorissen, 2018). ~2.1 g leucine 
    % content
elseif isequal(Protocol, 'Moore (2009)') % 25 g whey protein
    interventionDuration = 300; % min
    AminoAcidInput = 2.1; % grams
    CostParameters = [3, 4];

    % 3 whole eggs ~ 150 grams; leucine content of whole eggs: ~1091.8
    % mg/100g whole eggs (Attia et al., 2020). ~1.64 grams leu content
elseif isequal(Protocol, 'Mazzulla (2017)') % mixed macronutrient meal containing 18 g whole egg protein (~3 eggs)
    interventionDuration = 300; % min
    AminoAcidInput = 1.64; % grams
    CostParameters = [4];

% Rapamycin protocol -----
elseif isequal(Protocol, 'Valid., Dickinson (2011) rap')
    interventionDuration = 120; % min
    AminoAcidInput = 1.8; % grams
    CostParameters = [4, 6, 28, 46];

% Signaling knockdown protocol
elseif isequal(Protocol, 'Signaling knockdown') 
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [];

% Permissive signaling (p70S6K maintained at basal levels) protocol
elseif isequal(Protocol, 'Basal p70S6K Analysis') 
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [];

% Variable levels of p70S6K, variable mTORC1 kinase activity levels 
elseif isequal(Protocol, 'p70S6K reduced, mTORC1 kinase increased') 
    interventionDuration = 180; % min
    AminoAcidInput = 3.5; % grams
    CostParameters = [];     
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

if isequal(Protocol, 'Calibration, young adults') || isequal(Protocol ,'Valid., Glynn (2010) 1.85g') ||  ...
        isequal(Protocol ,'Valid., Mitchell (2015) bolus') || isequal(Protocol, 'Valid., Dickinson (2011) cntl') || ...
        isequal(Protocol,'Valid., Drummond (2010)') || isequal(Protocol, 'Valid., Wilkinson (2013)') || ...
        isequal(Protocol, 'Atherton (2010) - bolus') || isequal(Protocol, 'Moore (2009)') || ...
        isequal(Protocol, 'Mazzulla (2017)') || isequal(Protocol, 'Valid., Dickinson (2011) rap') || ...
        isequal(Protocol, "Basal p70S6K Analysis")|| isequal(Protocol, "p70S6K reduced, mTORC1 kinase increased")
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, KValueVector_given, OptimizeParameters, ...
        ab, KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), ...
        [0 :0.1: EqDur], x0_mass);
    
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;
    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220719(t2, x2, KValueVector_given, OptimizeParameters, ...
        ab, KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), ...
        [EqDur :0.1: (EqDur+RunDur)], x0_2);
                          
    x_bolus = [x; x2];
    t_bolus = [t; t2];
    
        x_bolus(:,41) = KValueVector_given(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,42) = KValueVector_given(9)*x_bolus(:,8)./x_bolus(:,10) + KValueVector_given(11)*x_bolus(:,7); % Fm,0; insulin mediated MPB
        x_bolus(:,43) = KValueVector_given(6)*x_bolus(:,4) - KValueVector_given(7)*x_bolus(:,6); % Net Balance
        x_bolus(:,44) = KValueVector_given(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
        x_bolus(:,45) = KValueVector_given(9)*x_bolus(:,8) ./ x_bolus(:,10); % MPB, r9 (p-IR concentration)
      
        % parameter 46 set as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_bolus(:,46) = nan;
        
    x = x_bolus;
    t = t_bolus;
        
elseif isequal(Protocol ,'Valid., Mitchell (2015) pulse') 
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, KValueVector_given, OptimizeParameters, ab, ...
        KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), [0 :0.01: EqDur-0.01], x0_mass);
    
    % first input
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;
    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220719(t2, x2, KValueVector_given, OptimizeParameters, ab, ...
        KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), [EqDur :0.01: (EqDur+IngestionTimes(2)-0.01)], x0_2);
    
    % second input
    x_length = size(x2,1); %vertical length (rows)
    x0_3 = x2(x_length, :);
    x0_3(1) = AminoAcidInput_mass;
    [t3,x3] = ode23s(@(t3,x3) modelDerivative_OIM_220719(t3, x3, KValueVector_given, OptimizeParameters, ab, ...
        KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), [(EqDur+IngestionTimes(2)) :0.01: (EqDur+IngestionTimes(3))-0.01], x0_3);
    
    % third input
    x_length = size(x3,1); %vertical length (rows)
    x0_4 = x3(x_length, :);
    x0_4(1) = AminoAcidInput_mass;
    [t4,x4] = ode23s(@(t4,x4) modelDerivative_OIM_220719(t4, x4, KValueVector_given, OptimizeParameters, ab, ...
        KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), [(EqDur+IngestionTimes(3)) :0.01: (EqDur+IngestionTimes(4))-0.01], x0_4);

    % fourth input
    x_length = size(x4,1); %vertical length (rows)
    x0_5 = x4(x_length, :);
    x0_5(1) = AminoAcidInput_mass;
    [t5,x5] = ode23s(@(t5,x5) modelDerivative_OIM_220719(t5, x5, KValueVector_given, OptimizeParameters, ab, ...
        KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), [(EqDur+IngestionTimes(4)) :0.01: (EqDur+RunDur)], x0_5);

    x_pulse = [x; x2; x3; x4; x5];
    t_pulse = [t; t2; t3; t4; t5];
    
        x_pulse(:,41) = KValueVector_given(6)*x_pulse(:,4); % Fm,a
        x_pulse(:,42) = KValueVector_given(9)*x_pulse(:,8)./x_pulse(:,10) + KValueVector_given(11)*x_pulse(:,7); % Fm,0; insulin mediated MPB
        x_pulse(:,43) = KValueVector_given(6)*x_pulse(:,4) - KValueVector_given(7)*x_pulse(:,6); % Net Balance
        x_pulse(:,44) = KValueVector_given(15)*x_pulse(:,6).*x_pulse(:,28); % MPS, r15
        x_pulse(:,45) = KValueVector_given(9)*x_pulse(:,8) ./ x_pulse(:,10); % MPB, r9 (p-IR concentration)
      
        % parameter 46 sets as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_pulse(:,46) = nan;
        
    x = x_pulse;
    t = t_pulse;
    
elseif isequal(Protocol ,'Signaling knockdown')
    
    % K-Values: specific for KD (run portion of ode)
    KValueVector_given_KD = kValues_KD(:,2);
    KValueVector_factor_KD(1:length(KValueVector_given_KD)) = 10.^(floor(log10(KValueVector_given_KD)));
    
    % Equilibrium period - knocked down k-values used
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, KValueVector_given_KD, OptimizeParameters, ...
        ab, KValueVector_factor, EqDur, GlucInfRate, LeuInfRate), ...
        [0 :0.1: EqDur], x0_mass);

    %Run period - adjusted k-values to simulate knockdown
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;
    
    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220719(t2, x2, KValueVector_given_KD, OptimizeParameters, ...
        ab, KValueVector_factor_KD, EqDur, GlucInfRate, LeuInfRate), ...
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