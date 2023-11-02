clear
clc
close all

FolderPath = fileparts('/Users/taylormccoll/Documents/SFU/1. Ph.D./1. Research Projects/1. Network Feedback/6. Manuscript/Manuscript Versions/230207/Manuscript/MuscleProteinSynthesisKineticModel/McColl_2023_Muscle Protein Synthesis Kinetic Model_230221/');
addpath(genpath(FolderPath)) % add all subfolders within 'FolderPath'

%% Inputs

FileName = sprintf('Parameter optimization - %s', datestr(now));
FilePath = fullfile(FolderPath, 'Output Plots/');

% Basal Infusion Rates
GlucoseInfRate = 75e-3 *(1/180.16); % 75 mg/min glucose infusion; converted to units of moles/min
LeucineInfRate = 2.5e-3*(1/131.17); % 25 mg/min leucine infusion; converted to units of moles/min

% Subject specific data (used to calculate compartment volumes)
gender = 1; % 1=male, 0=female
height = 173; % cm
mass = 72.1; % kg
age = 41.9; % years
bia = 517; % Bioelectical Impedance (Ohms)
subjectData = [gender, height, mass, age, bia];

% Intervention
AminoAcidInput = 3.5;

EquilibriumDuration = 300;
RunDuration = 250;

%Experimental Data
ExpData = xlsread('230221_experimental data.xlsx', 'Calibration, young adults'); %updated dataset

% Initial Conditions
InitialValues = xlsread('230221_Initial Values, IOM', '3. Oscillatory insulin');

% K-Values
KValueVector_given = xlsread('230221_K-values, IOM', '3. Insulin oscillatory module');

%% Optimization Information
CostParameters = [3, 4, 6, 21, 28, 41, 42]; % 41,42 updated (FSR is 46, but calculated within plots)

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
    ]; 

% %% Insulin specific opt
% CostParameters = [3]; 
% 
% OptimizeParameters = [51:67];

%% Calculated Values

% Compartment Volumes
plasmaVolume = plasmaVolumeFxn(gender, height, mass);
skeletalMuscleVolume = skeletalMuscleVolumeFxn(gender, height, age, bia);

% Initial Conditions - conversion to masses
InitialValues = InitialValues(:,1);
InitialValues_mass = x0_massFxn(InitialValues, plasmaVolume, skeletalMuscleVolume);

LeucineMolarMass = 1/131.17; %mol/g
AminoAcidInput_mass = AminoAcidInput*LeucineMolarMass;

% Experimental Data - conversion to mass
expData_mass = expData_massFxn(ExpData, CostParameters, plasmaVolume, skeletalMuscleVolume);
ExpData = expData_mass;

% K-Values
KValueVector_given = KValueVector_given(:,2);
KValueVector_factor(1:length(KValueVector_given)) = 10.^(floor(log10(KValueVector_given)));

% Calculating inputs for optimizer (base k-values - to a factor of 10^0)
ab0 = nan(1,length(OptimizeParameters));
for n = 1:length(OptimizeParameters)
    value = KValueVector_given(OptimizeParameters(n))/KValueVector_factor(OptimizeParameters(n));
    ab0(n) = value;
end

% ** why are these needed?
% A = lin_ineq_con_matrix(ab0);
% % b = lin_ineq_con_vector(ab0, 100); % range of 100x above or below initial value 
% b = lin_ineq_con_vector(ab0, 10); % range of 10x above or below initial value 

ab0;


%% Global Optimizer - Single Global Minimum

tic

options = optimoptions(@fmincon,'Algorithm','sqp','Display','off');
problem.lb = ab0/100;
problem.ub = ab0*100;

% problem.lb = ab0/10000;   %**
% problem.ub = ab0*10000;   %**


problem = createOptimProblem('fmincon', ...
    'objective', @(ab0)objectiveFunction_230331(KValueVector_given, KValueVector_factor, InitialValues_mass, OptimizeParameters, ...
        ab0, ExpData, CostParameters, AminoAcidInput_mass, EquilibriumDuration, RunDuration, GlucoseInfRate, LeucineInfRate), ...
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

% rng(14,'twister') %reproducible data
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

%% test output

x0 = InitialValues_mass;

%ODE
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_210923(t, x, KValueVector_given, OptimizeParameters, ab, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, LeucineInfRate), [0 :0.01: EquilibriumDuration], x0);
    
    x_length = size(x,1); %vertical length (rows)
    x0_2 = x(x_length, :);
    x0_2(1) = AminoAcidInput_mass;
    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_210923(t2, x2, KValueVector_given, OptimizeParameters, ab, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, LeucineInfRate), [EquilibriumDuration :0.01: (EquilibriumDuration+RunDuration)], x0_2);
                          
    x_bolus = [x; x2];
    t_bolus = [t; t2];
    
        x_bolus(:,40) = KValueVector_given(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,41) = KValueVector_given(9)*x_bolus(:,8) + KValueVector_given(11)*x_bolus(:,7); % Fm,0
    
    x = x_bolus;
    t = t_bolus;
    
% Plot
clf
FileName = sprintf('1, %s - moles %s', 'Calibration', datestr(now));
% FilePath = FilePathCumulative;
figNum=1;
interventionDuration = 240;
Plot_210714(figNum, t, x, expData_mass, CostParameters, RunDuration, ...
    EquilibriumDuration, FilePath, FileName, 'Calibration', KValueVector_given, interventionDuration)




