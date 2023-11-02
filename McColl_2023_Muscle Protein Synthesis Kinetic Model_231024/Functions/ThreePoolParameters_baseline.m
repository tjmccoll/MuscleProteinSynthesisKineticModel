function [Parameters, Title, ThreePool_baseline] = ThreePoolParameters_baseline(x, t, KValue, ...
    EquilibriumDuration, ThreePool_expData)
% Comparison of experimentally measure and simulated 3-pool parameter
% values in basal conditions (i.e., prior to leucine administration in the 
% model).
%   [Parameters, Title, ThreePool_baseline] = 
%       ThreePoolParameters_baseline(x, t, KValue, EquilibriumDuration, 
%       ThreePool_expData)
%   x = column vector of all species across the simulated duration in units
%   of mol/L as outputted from the ODE solver
%   t = scalar time outputted from the ODE solver
%   KValue = k-value vector use to simulate model
%   EquilibriumDuration = duration of the equilibrium period to allow model
%   to acheive steady state prior to stimulation (minutes)
%   ThreePool_expData = experimentally measured 3-pool parameter data

%%
t_endEq = find(t==EquilibriumDuration,1);

%% Calculating 3-pool parameters - Model specific
    F_ma = KValue(6)*x(t_endEq, 4);
    F_m0 = KValue(9)*x(t_endEq, 8) + KValue(11)*x(t_endEq, 7);
    F_vm = KValue(7)*x(t_endEq, 6);
    F_0m = KValue(8)*x(t_endEq,4) + KValue(10)*x(t_endEq, 4);
    KIC_ox = KValue(14)*x(t_endEq,7);
     
    ThreePool_model = [F_ma, F_m0, F_vm, F_0m, KIC_ox];
    
%% 
    Title={'Parameter' 'Model Baseline'};

    Parameters = {'F_{m,a}' 'F_{m,0}' 'F_{v,m}' 'F_{0,m}' 'KIC_{ox}**'}';
    
    ThreePool_baseline = [ThreePool_expData, ThreePool_model'];
