function [Title, RMS_CostIndividual, RMS_CostTotal, RMS_norm_CostIndividual, RMS_norm_CostTotal, RSS_CostIndividual, RSS_CostTotal] = ...
    CostCalculation_230920(t, x, Akt_pS, p70S6K_pT, EqDur, CostParameters, ExpData, FSRsim)
%   Calculating the cost value for each model component with experimental
%   data from the model simulation. I.e., calculating the difference from
%   the simulated value and the experimentally measured value.
%       [Title, ~, ~, RMS_norm_costInd, RMS_norm_costTotal, ~, ~] = 
%           CostCalculation(t, x, Akt_S_T, p70S6K_T, EqDur,
%           CostParameters, ExpData)
%       t = time vector
%       x = simualted data vector
%       Akt_S = positions within data set of Akt with phosphorylated
%       serine. The vector should include the position of the single
%       phospho site first - the first entry is referenced in the
%       ExpDataMatrix for loop
%       p70S6K_T = position within data set
%       EqDur = Equilibrium Duration
%       CostParameters = Parameters with corresponding experimental data 
%       used to calculate RMSE
%       ExpData = Experimental Data
%       FSRsim = simulated grams of leucine synthesized (i.e., integrated MPS curve)

% Find the end of the equilibrium period
t_endEq = find(t==EqDur);

% Find the values of Akt and p70S6K at the end of the eq period
    % needed to determine initial value -> factor change
Akt_S_iv = x(t_endEq(1), Akt_pS); %summing IV both Akt species with phospho serine
Akt_S_iv_sum = sum(Akt_S_iv);
p70S6K_T_iv = x(t_endEq(1), p70S6K_pT);

% Setting up the Matrix to clean exp data
ExpDataMatrix = nan(length(CostParameters), 13, 4);

ExpDataMatrix(:,1,1) = CostParameters';

for n = 1:length(CostParameters)
    if ExpDataMatrix(n,1,1) == Akt_pS(1) % 
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataMatrix(n, 1:length(rows),2) = ExpData(rows, 2) + EqDur; %Adjusted to Equilibrium period
        ExpDataMatrix(n, 1:length(rows),3) = ExpData(rows, 3) * Akt_S_iv_sum;
        ExpDataMatrix(n, 1:length(rows),4) = sqrt(Akt_S_iv_sum^2 * ExpData(rows, 4).^2); % Transforming SE data (factoring in IV)
            
    elseif ExpDataMatrix(n,1,1) == p70S6K_pT
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataMatrix(n, 1:length(rows),2) = ExpData(rows, 2) + EqDur; %Adjusted to Equilibrium period
        ExpDataMatrix(n, 1:length(rows),3) = ExpData(rows, 3) * p70S6K_T_iv;
        ExpDataMatrix(n, 1:length(rows),4) = sqrt(p70S6K_T_iv^2 * ExpData(rows, 4).^2); % Transforming SE data (factoring in IV)
        
    else
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataMatrix(n, 1:length(rows),2) = ExpData(rows, 2) + EqDur; %Adjusted to Equilibrium period
        ExpDataMatrix(n, 1:length(rows),3) = ExpData(rows, 3);
        ExpDataMatrix(n, 1:length(rows),4) = ExpData(rows, 4);
    end
end

% empty strings for the for loop
RMS_CostOutput = nan(1, length(CostParameters));    % Root Mean Square
RSS_CostOutput = nan(1, length(CostParameters));    % Residual Sum of Squares
RMS_norm_CostOutput = nan(1, length(CostParameters));   % RMS normalized
RSS_norm_CostOutput = nan(1, length(CostParameters));   % RSS normalized

% cleaning ExpDataMatrix for input into the SLS cost function (i.e., removing NaN entries)
% SLS cost calculated (output = total cost for that parameter)
    for n = 1:length(CostParameters)
        % Cleaning the ExpDataMatrix variable (removing NaNs)
        data_t = ExpDataMatrix(n,:,2);
        data_t_clean = data_t(~isnan(data_t));
        data_value = ExpDataMatrix(n,:,3);
        data_value_clean = data_value(~isnan(data_value));
        data_err = ExpDataMatrix(n,:,4);
        data_err_clean = data_err(~isnan(data_err));

        % calculating goodness of fit for each cost parameter    
        if CostParameters(n) == Akt_pS(1)
            % summing the serine phosphorylated Akt species for cost calculation
            CostParameter_Akt = Akt_pS;
            [RMS_CostOutput(n), RSS_CostOutput(n), RMS_norm_CostOutput(n), RSS_norm_CostOutput(n)] = ...
                RootMeanSquareCalc_230221(CostParameters(n), x(:,CostParameter_Akt), data_value_clean, data_err_clean, ...
                t, data_t_clean, FSRsim);

        else
            [RMS_CostOutput(n), RSS_CostOutput(n), RMS_norm_CostOutput(n), RSS_norm_CostOutput(n)] = ...
                RootMeanSquareCalc_230221(CostParameters(n), x(:,CostParameters(n)), data_value_clean, data_err_clean, ...
                t, data_t_clean, FSRsim);
        end
    end

% Root Mean Square
Title={'Parameter' 'RMS'};
RMS_CostIndividual = [CostParameters', RMS_CostOutput'];
RMS_CostTotal = sum(RMS_CostOutput);

% Normalized Standard Least Square (normalized to number of data points)
RMS_norm_CostIndividual = [CostParameters', RMS_norm_CostOutput'];
RMS_norm_CostTotal = sum(RMS_norm_CostOutput);

% Residual Sum of Squares
RSS_CostIndividual =[CostParameters', RSS_CostOutput'] ;
RSS_CostTotal = sum(RSS_CostOutput);

% % Normalized Residual Sum of Squares (normalized to number of data points)
% RMS_norm_CostIndividual = [CostParameters', RSS_norm_CostOutput'];
% RMS_norm_CostTotal = sum(RMS_norm_CostIndividual);

end