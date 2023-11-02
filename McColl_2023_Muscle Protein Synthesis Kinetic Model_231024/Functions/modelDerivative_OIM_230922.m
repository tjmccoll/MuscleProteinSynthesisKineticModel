function [dx] = modelDerivative_OIM_230922(t, x, KValueVector_given, OptimizeParameters, ...
    OptimizeParameters_initial, KValueVector_factor, EquilibriumDuration, GlucoseInfRate, ...
    LeucineInfRate_given)
% System of ordinary differential equations that the ODE function solves.
%   [dx] = modelDerivative_OIM_220719(t, x, KValueVector_given,
%       OptimizeParamters, OptimizeParameters_initial, KValueVector_factor, 
%       EquilibriumDuration, GlucoseInfRate, LeucineInfRate_given)
%   t = scalar t as required by the ode solver
%   x = column vector x of all species as required by the ode solver
%   KValueVector_given = inputted K-values
%   OptimizeParameters = K-values parameters that were optimized. String of numbers
%   corresponding to k-values.
%   OptimizeParameters_initial = initial guess of k-values to optimize
%   (iteratively changed in optimizer)
%   KValueVector_factor = required to calculate K-values. Model 
%   optimization is completed on a number without the order of magnitude 
%   (i.e., m*10^n: m = KValueVector_given, n = KValueVector_factor)
%   EquilibriumDuration = duration of equilibrium period
%   GlucoseInfRate = Glucose infusion rate required for the Sturis module
%   LeucineInfRate_given = Leucine infusion rate required to maintain
%   leucine at a steady state during the equilibrium period
    
%% Concentrations
Stomach = x(1); Gut = x(2);
Leu_plasma = x(4); KIC_plasma = x(5);
Leu_in = x(6); KIC_in = x(7);
Protein = x(8);
IR_beta = x(9); IR_beta_pY1146 = x(10); IR_beta_refractory = x(11);
IRS1 = x(12); IRS1_pY = x(13); IRS1_pS636 = x(14);
PI3K = x(15); IRS1_PI3K = x(16);
PDK1 = x(17); PDK1_p = x(18);
AKT = x(19); AKT_pT308 = x(20); Akt_pS473 = x(21); AKT_pT308_pS473 = x(22);
TSC_clx = x(23); TSC_p_clx = x(24);
mTORC1_inactive = x(25); mTORC1_active = x(26);
p70S6K = x(27); p70S6K_pT389 = x(28);
PI3K_variant = x(29); PI3K_variant_p = x(30);
mTORC2_0 = x(31); mTORC2_pS2481_0 = x(32);

% Sturis (1991) - oscillatory insulin secretion module
InsulinPlasma = x(3);
InsulinRemote = x(35);
Glucose = x(36);
x1 = x(37);
x2 = x(38);
x3 = x(39);

%% K-Values
% for loop to update k-values that were changed during the model
% optimization

KValue = zeros(1,length(KValueVector_given));

for n = 1:length(KValueVector_given)
    if any(OptimizeParameters(:) == n)
        KValue(1,n) = OptimizeParameters_initial(OptimizeParameters==n) * KValueVector_factor(n);
        
    else
        KValue(1,n) = KValueVector_given(n);
        
    end
end

%% Reactions

r1 = KValue(1)*Stomach;
r2 = KValue(2)*Gut;
r3 = KValue(3)*Gut;
r4 = KValue(4)*Leu_plasma;

r6 = KValue(6)*Leu_plasma;
r7 = KValue(7) * (Leu_in - Leu_plasma); % active transport simulated

r9 = KValue(9)*Protein / AKT_pT308; % insulin-mediated MPB (p-Akt(T); Akt/FoxO interaction)

r10 = KValue(10)*Leu_in;
r11 = KValue(11)*KIC_in;
r12 = KValue(12)*KIC_plasma;
r13 = KValue(13)*KIC_in;
r14 = KValue(14)*KIC_in;
r15 = KValue(15)*Leu_in*p70S6K_pT389;
r16 = KValue(16)*IR_beta*InsulinPlasma;   
r17 = KValue(17)*IR_beta_pY1146;
r18 = KValue(18)*IR_beta_refractory;
r19 = KValue(19)*IRS1*IR_beta_pY1146;
r20 = KValue(20)*IRS1_pY;

r21 = KValue(21)*IRS1*p70S6K_pT389;
r22 = KValue(22)*IRS1_pS636;
r23 = KValue(23)*IRS1_pY*PI3K;
r24 = KValue(24)*IRS1_PI3K;
r25 = KValue(25)*PDK1*IRS1_PI3K;
r26 = KValue(26)*PDK1_p;
r27 = KValue(27)*AKT*PDK1_p;
r28 = KValue(28)*AKT_pT308;
r29 = KValue(29)*AKT*mTORC2_pS2481_0;
r30 = KValue(30)*Akt_pS473;

r31 = KValue(31)*AKT_pT308*mTORC2_pS2481_0;
r32 = KValue(32)*AKT_pT308_pS473;
r33 = KValue(33)*Akt_pS473*PDK1_p;
r34 = KValue(34)*AKT_pT308_pS473;
r35 = KValue(35)*TSC_clx*AKT_pT308_pS473;
r36 = KValue(36)*TSC_clx*AKT_pT308;
r37 = KValue(37)*TSC_p_clx;
r38 = KValue(38)*mTORC1_active*TSC_clx;
r39 = KValue(39)*mTORC1_inactive*Leu_in;
r40 = KValue(40)*p70S6K*mTORC1_active;

r41 = KValue(41)*p70S6K*PDK1_p;
r42 = KValue(42)*p70S6K_pT389;
r43 = KValue(43)*mTORC1_active*p70S6K_pT389;
r44 = KValue(44)*PI3K_variant*IR_beta_pY1146;
r45 = KValue(45)*PI3K_variant_p;
r46 = KValue(46)*mTORC2_0*PI3K_variant_p;
r47 = KValue(47)*mTORC2_pS2481_0;

r68 = KValue(68)*Gut;
r69 = KValue(69)*Protein / mTORC1_active; % insulin-indepedent MPB (mTORC1/ULK interaction)

%% Sturis (1991) - insulin oscillations

V_p = KValue(48); V_i = KValue(49); V_g = KValue(50);
E = KValue(51); 
t_p = KValue(52); t_i = KValue(53); t_d = KValue(54);
Rm = KValue(55); a1 = KValue(56); C1 = KValue(57); 
Ub = KValue(58); C2 = KValue(59); C3 = KValue(60);
U_0 = KValue(61); U_m = KValue(62); Beta = KValue(63); 
C4 = KValue(64); Rg = KValue(65); Alpha = KValue(66);
C5 = KValue(67);

% Calculated Values
f1 = Rm/(1+exp( (C1-(Glucose/V_g)) /a1));
f2 = Ub*(1-exp(-Glucose/(C2*V_g)));
f3 = Glucose/(C3*V_g);
f4 = U_0 + ((U_m-U_0) / (1+exp(-Beta*log(InsulinRemote/C4 * (1/V_i+1/(E*t_i))))));
f5 = Rg / (1+exp(Alpha*(x3/V_p - C5)));


%% ODE's
% Leucine infusion during equilibrium period to achieve steady state.
% Following which the infusion is removed and leucine dynamics are
% determined by leucine feeding

if t<EquilibriumDuration
    LeucineInfRate = LeucineInfRate_given;
else
    LeucineInfRate = 0;
end


dx = zeros(40,1);

%Digestive Tract
dx(1) = -r1;
dx(2) = r1 - r2 - r3 - r68; %splanchnic extraction included

%Blood Plasma
dx(4) = LeucineInfRate + r3 + r7 - r6;
dx(5) = r13 - r12; 

%Skeletal Muscle
dx(6) = r6 + r9 + r11 + r69 - r7 - r10 - r15; 
dx(7) = r10 + r12 - r11 - r13 - r14;
dx(8) = r15 - r9 - r69; 

dx(9) = r18 - r16;
dx(10) = r16 - r17;
dx(11) = r17 - r18;
dx(12) = r20 + r22 - r19 - r21;
dx(13) = r19 + r24 - r20 - r23;
dx(14) = r21 - r22;
dx(15) = r24 - r23; 
dx(16) = r23 - r24;
dx(17) = r26 - r25;
dx(18) = r25 - r26;

dx(19) = r28 + r30 - r27 - r29;
dx(20) = r27 + r32 - r28 - r31;
dx(21) = r29 + r34 - r30 - r33;
dx(22) = r31 + r33 - r32 - r34;

dx(23) = r37 - r35 - r36;
dx(24) = r35 + r36 - r37;
dx(25) = r38 + r43 - r39;
dx(26) = r39 - r38 - r43;
dx(27) = r42 - r40 - r41;
dx(28) = r40 + r41 - r42;

dx(29) = r45 - r44;
dx(30) = r44 - r45;
dx(31) = r47 - r46;
dx(32) = r46 - r47;

% Sturis (1991)
dx(3) = f1 - E*(InsulinPlasma/V_p - InsulinRemote/V_i) - InsulinPlasma/t_p + r4;
dx(35) = E*(InsulinPlasma/V_p - InsulinRemote/V_i) - InsulinRemote/t_i;
dx(36) = GlucoseInfRate - f2 - f3*f4 + f5;
dx(37) = 3*(InsulinPlasma-x1)/t_d;
dx(38) = 3*(x1-x2)/t_d;
dx(39) = 3*(x2-x3)/t_d;

% Excretion (not considered ODE's)
dx(33) = r2;    % excretion from gut
dx(34) = r14;   % oxidation of Intracellular KIC
dx(40) = r68;   % splanchnic extraction
