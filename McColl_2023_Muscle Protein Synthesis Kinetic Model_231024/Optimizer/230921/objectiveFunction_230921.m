function CostTotal = objectiveFunction_230921(KValueVector_given, KValueVector_factor, InitialValues_mass, ...
    OptimizeParameters, OptimizeParameters_initial, ExpData, CostParameters, AminoAcidInput, ...
    EquilibriumDuration, runDuration, interventionDuration, GlucoseInfRate, LeucineInfRate, plasmaVolume, ...
    skeletalMuscleVolume)
    %   KValueVector_given = the full string of *initial* k-values
    %   KValueVector_factor = required to calculate K-values 
    %     -> optimization is completed on a number without the order of
    %     magnitude (i.e., m*10^n; m = OptimizeParameters_initial, 10^n = KValueVector_factor)
    %   InitialValues = initial conditions 
    %   OptimizeParameters = k-values to be optimized
    %   OptimizeParameters_initial = initial guess of k-values to optimize
    %       (iteratively changed in optimizer)
    %   ExpData = Experimental data set (calibration set)
    %   CostParameters = Experimental data set(s) to fit the model against
    %   AminoAcidInput = amino acid input at t=0 (calibration)
    %   EquilibriumDuration = duration of equilibrium period
    %   RunDuration = duration of run period
    %   interventionDuration = duration of intervention from experimental
    %   data (needed for FSR calculation)
    %   LeucineInfRate = Leucine infusion rate. Required to balance leucine
    %       oxidation during equilibrium period. Infusion is set to 0 when
    %       intervention begins

%% inital values
    x0 = InitialValues_mass;
    
%% ODE
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_230922(t, x, KValueVector_given, OptimizeParameters, ...
        OptimizeParameters_initial, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, ...
        LeucineInfRate), ...
        [0 : 0.1 : EquilibriumDuration], x0);
    
            x_length = size(x,1); %vertical length (rows)
            x0_2 = x(x_length, :);
            x0_2(1) = AminoAcidInput; %Amino Acid Ingestion

    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_230922(t2, x2, KValueVector_given, OptimizeParameters, ...
        OptimizeParameters_initial, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, LeucineInfRate), ...
        [EquilibriumDuration : 0.1 : (EquilibriumDuration+runDuration)], x0_2);
                          
    x_bolus = [x; x2];
    t_bolus = [t; t2];
    
%% Updating kValue vector to the accurate k-values (i.e., combine base value and factor)     
    KValue_opt = zeros(1,length(KValueVector_given));

    for n = 1:length(KValueVector_given)
        if any(OptimizeParameters(:) == n)
            KValue_opt(1,n) = OptimizeParameters_initial(OptimizeParameters==n) * KValueVector_factor(n);
        else
            KValue_opt(1,n) = KValueVector_given(n);
        end
    end

%% 3-pool parameters   
        x_bolus(:,41) = KValue_opt(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,42) = KValue_opt(9)*x_bolus(:,8)./x_bolus(:,20) + ...
            KValue_opt(69)*x_bolus(:,8)./x_bolus(:,26) + ...
            KValue_opt(11)*x_bolus(:,7) ; % Fm,0; insulin-mediated MPB & insulin-independent MPB
        x_bolus(:,43) = KValue_opt(6)*x_bolus(:,4) - ...
            KValue_opt(7)*x_bolus(:,6); % Net Balance
        x_bolus(:,44) = KValue_opt(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
        x_bolus(:,45) = KValue_opt(9)*x_bolus(:,8) ./ x_bolus(:,20) + ...
            KValue_opt(69)*x_bolus(:,8) ./ x_bolus(:,26); % MPB: r9 (p-Akt(T) concentration) + r69 (mTORC1 mediated)
      
        % parameter 46 set as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_bolus(:,46) = nan;
        
    x = x_bolus; % in units of moles
    t = t_bolus;
        
%% Converting moles to concentration
x_conc = x_concFxn(x, InitialValues_mass, plasmaVolume, skeletalMuscleVolume);

%% 3-pool parameters; FSR calculation

F_ma = KValue_opt(6).*x(:,4); % r6 = F_m,a
r9 = KValue_opt(9).*x(:,8)./x(:,20); % F_m,0 (1); MPB (p-Akt(T) mediated)
r11 = KValue_opt(11).*x(:,7); % F_m,0 (2)
r69 = KValue_opt(69).*x(:,8) ./ x(:,26);
F_m0 = r9+r11+r69; 

% FSR/MPS; integral for FSR comparison
t_endEq = find(t==EquilibriumDuration);
mps = KValue_opt(15)*x(:,6).*x(:,28); % FSR, r15
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% Updating simulation data - Fm,a and Fm,0 added with correct units (mol/min)
    %allows for correct cost calculation (experimental data in mol/min)
x_conc_update = x_conc;
x_conc_update(:,41) = F_ma;
x_conc_update(:,42) = F_m0;
x_conc_update(:,44) = mps;

%% Cost calculation

Akt_pS_pos = [21,22]; % both Akt species with phospho serine
p70S6KT389_pos = 28;
    
[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230920(t, x_conc_update, Akt_pS_pos, p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);
   
end

