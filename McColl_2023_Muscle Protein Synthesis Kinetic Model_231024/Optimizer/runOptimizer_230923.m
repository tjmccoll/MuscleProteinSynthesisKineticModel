clear
clc
close all

FolderPath = fileparts('/Users/taylormccoll/Documents/SFU/1. Ph.D./1. Research Projects/1. Network Feedback/6. Manuscript/Manuscript Versions/230912_revision/Manuscript/MuscleProteinSynthesisKineticModel/McColl_2023_Muscle Protein Synthesis Kinetic Model_230919/');
addpath(genpath(FolderPath)) % add all subfolders within 'FolderPath'

%% File Output
FilePath = fullfile(FolderPath, 'Optimizer/Output Plots/');
NewFolder = datestr(now, 'yy-mm-dd HH-MM-SS');
    CumulativeOutput = mkdir(FilePath, NewFolder);
FilePathCumulative = fullfile(FilePath, NewFolder);

%% Inputs
% Basal infusion rates
GlucoseInfRate = 75e-3 *(1/180.16); % 75 mg/min glucose infusion; converted to units of moles/min
LeucineInfRate = 2.5e-3*(1/131.17); % 25 mg/min leucine infusion; converted to units of moles/min

% Subject specific data (used to calculate compartment volumes)
gender = 1; % 1=male, 0=female
height = 173; % cm
mass = 72.1; % kg
age = 41.9; % years
bia = 517; % Bioelectical Impedance (Ohms)
SubjectData = [gender, height, mass, age, bia];

% Intervention
AminoAcidInput = 3.5; % grams

EquilibriumDuration = 300; % burn-in duration
RunDuration = 250; % ODE run duration
interventionDuration = 180; % experimental intevention (required to calculate FSR)

% Initial concentrations
x0 = xlsread('230922_Initial Values, IOM', '3. Oscillatory insulin'); % reduced phospho-proteins

% Kinetic parameters
kineticParameters = xlsread('230922_K-values, IOM', '230923');

% Experimental data (calibration data set)
ExpData = xlsread('230221_experimental data.xlsx', 'Calibration, young adults'); 

%% Optimization Information
CostParameters = [3, 4, 5, 6, 21, 28, 41, 42, 46]; 

OptimizeParameters = [4 ... % insulin stimulation
    6 7 9 10 11 12 13 14 ... % leucine module
    15 ... %protein synthesis
    16 17 18 ... % IR_B
    19 20 21 22 23 24 ... %IRS1 / IRS1_PI3K module
    25 26 ... % PDK1
    27 28 29 30 31 32 33 34 ... %Akt
    35 36 37 ... %TSC
    38 39 43 ... %mTORC1 activation
    40 41 42 ... %p70S6K stimulation
    44 45 ... % PI3K_var
    46 47 ... % mTORC2
    69 ... % MPB (mTORC1 mediated)
    ]; 

%% Calculated values
% Compartment volumes
plasmaVolume = plasmaVolumeFxn(SubjectData(1), SubjectData(2), SubjectData(3));
skeletalMuscleVolume = skeletalMuscleVolumeFxn(SubjectData(1), SubjectData(2), SubjectData(4), SubjectData(5));

% Initial Conditions - conversion to moles
x0_mass = x0_massFxn(x0(:,1), plasmaVolume, skeletalMuscleVolume);

LeucineMolarMass = 1/131.17; %mol/g
AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;

% Experimental Data - conversion to mass
expData_conc = ExpData(:,1:4);
expData_mass = expData_massFxn(expData_conc, CostParameters, plasmaVolume, skeletalMuscleVolume);

% K-Values
KValueVector_given = kineticParameters(:,2);
KValueVector_factor(1:length(KValueVector_given)) = 10.^(floor(log10(KValueVector_given)));

% Calculating inputs for optimizer (base k-values - to a factor of 10^0)
ab0 = nan(1,length(OptimizeParameters));
for n = 1:length(OptimizeParameters)
    value = KValueVector_given(OptimizeParameters(n))/KValueVector_factor(OptimizeParameters(n));
    ab0(n) = value;
end

%% Confirming objective function operates correctly
% Same cost as model execution

% change output variables if needed
[costTotal, costIndiv, FSRintegral] = objectiveFunction_230923(KValueVector_given, KValueVector_factor, x0_mass, OptimizeParameters, ...
        ab0, expData_conc, CostParameters, AminoAcidInput_mass, EquilibriumDuration, RunDuration, interventionDuration,...
        GlucoseInfRate, LeucineInfRate, plasmaVolume, skeletalMuscleVolume)

%% Global Optimizer - Single Global Minimum
tic

options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
problem.lb = ab0/100;
problem.ub = ab0*100;

% problem.lb = ab0/10000;   %**
% problem.ub = ab0*10000;   %**

problem = createOptimProblem('fmincon', ...
    'objective', @(ab0)objectiveFunction_230923(KValueVector_given, KValueVector_factor, x0_mass, OptimizeParameters, ...
        ab0, expData_conc, CostParameters, AminoAcidInput_mass, EquilibriumDuration, RunDuration, interventionDuration, ...
        GlucoseInfRate, LeucineInfRate, plasmaVolume, skeletalMuscleVolume), ...
    'x0', ab0, ...
    'options', options, ...
    'lb', problem.lb, 'ub', problem.ub);

% problem.lb = ab0/10;
% problem.ub = ab0*10;

% gs = GlobalSearch('Display', 'iter', ...
%         'NumStageOnePoints', 2000, ...
%     'NumTrialPoints', 10000);

gs = GlobalSearch('Display', 'iter', ...
        'NumStageOnePoints', 500, ...
    'NumTrialPoints', 2500);

rng(14,'twister') %reproducible data
[ab,fval] = run(gs, problem)

toc

%% Writing updated k-values to Excel file
% Creating array with optimized k-values
optimizedKValues = zeros(1,length(KValueVector_given));

for n = 1:length(KValueVector_given)
    if any(OptimizeParameters(:) == n)
        optimizedKValues(1,n) = ab(OptimizeParameters==n) * KValueVector_factor(n);
        
    else
        optimizedKValues(1,n) = KValueVector_given(n);
    end
end
KValueVector_opt = optimizedKValues';

% optimizedKValues_output = [(1:length(KValueVector_given))'; KValueVector_given];
optimizedKValues_output = [(1:length(KValueVector_given))', optimizedKValues'];

csvFileName = sprintf('%s',datestr(now), ' optimizedParameters.csv');
csvFile = fullfile(FilePath, csvFileName);
csvwrite(csvFile, optimizedKValues_output)

%% Confirm optimizer output

[costTotal] = objectiveFunction_230923(KValueVector_opt, KValueVector_factor, x0_mass, [], ...
        [], expData_conc, CostParameters, AminoAcidInput_mass, EquilibriumDuration, RunDuration, interventionDuration,...
        GlucoseInfRate, LeucineInfRate, plasmaVolume, skeletalMuscleVolume)

