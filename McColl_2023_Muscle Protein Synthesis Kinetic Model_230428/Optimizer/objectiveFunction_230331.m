function CostTotal = objectiveFunction_230331(KValueVector_given, KValueVector_factor, InitialValues, ...
    OptimizeParameters, OptimizeParameters_initial, ExpData, CostParameters, AminoAcidInput, ...
    EquilibriumDuration, interventionDuration, GlucoseInfRate, LeucineInfRate)
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
    %   LeucineInfRate = Leucine infusion rate. Required to balance leucine
    %       oxidation during equilibrium period. Infusion is set to 0 when
    %       intervention begins

%% inital values
    x0 = InitialValues;
    
%% ODE
    [t,x] = ode23s(@(t,x) modelDerivative_OIM_220719(t, x, KValueVector_given, OptimizeParameters, ...
        OptimizeParameters_initial, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, LeucineInfRate), ...
        [0 EquilibriumDuration], x0);
    
            x_length = size(x,1); %vertical length (rows)
            x0_2 = x(x_length, :);
            x0_2(1) = AminoAcidInput; %Amino Acid Ingestion

    [t2,x2] = ode23s(@(t2,x2) modelDerivative_OIM_220719(t2, x2, KValueVector_given, OptimizeParameters, ...
        OptimizeParameters_initial, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, LeucineInfRate), ...
        [EquilibriumDuration (EquilibriumDuration+interventionDuration)], x0_2);
                          
    x_bolus = [x; x2];
    t_bolus = [t; t2];
    
%% 3-pool parameters            
    KValue_opt = zeros(1,length(KValueVector_given));

    for n = 1:length(KValueVector_given)
        if any(OptimizeParameters(:) == n)
            KValue_opt(1,n) = OptimizeParameters_initial(OptimizeParameters==n) * KValueVector_factor(n);

        else
            KValue_opt(1,n) = KValueVector_given(n);

        end
    end

        x_bolus(:,41) = KValue_opt(6)*x_bolus(:,4); % Fm,a
        x_bolus(:,42) = KValue_opt(9)*x_bolus(:,8) + KValue_opt(11)*x_bolus(:,7); % Fm,0
        x_bolus(:,43) = KValue_opt(6)*x_bolus(:,4) - KValue_opt(7)*x_bolus(:,6); % Net Balance
        x_bolus(:,44) = KValue_opt(15)*x_bolus(:,6).*x_bolus(:,28); % MPS, r15
        x_bolus(:,45) = KValue_opt(9)*x_bolus(:,8) ./ x_bolus(:,10); % MPB, r9 (p-IR concentration)
      
        % parameter 46 set as nan's to allow for cost calculation to work
        % (the x_bolus array needs to contain values for each cost parameter to be calculated)
        x_bolus(:,46) = nan;
        
    x = x_bolus;
    t = t_bolus;
        
%% Leucine synthesized in muscle

% FSR/MPS; integral for FSR comparison
t_endEq = find(t==EquilibriumDuration);
mps = x_bolus(:,44);
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% Cost calculation

AktS473_pos = 21;
p70S6KT389_pos = 28;
    
[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230221(t, x, AktS473_pos, p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);
    
end

