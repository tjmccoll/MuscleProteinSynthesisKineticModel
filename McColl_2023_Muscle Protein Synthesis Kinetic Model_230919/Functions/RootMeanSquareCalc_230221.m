function [RoomMeanSquare, ResidualSumSquares, RoomMeanSquare_norm, ResidualSumSquares_norm] = ...
    RootMeanSquareCalc_230221(costParameter, y, yhat, yError, t, t_point, FSRsim)
%   Calculating the Root Mean Square (RMS), Residual Sum of Squares
%   (RSS), normalized RMS, and normalized RSS for individual components of
%   the model simulation with corresponding experimental data.
%       [RootMeanSquare, ResidualSumSquares, RootMeanSquare_norm,
%           ResidualSumSquares_norm] = StandardLeastSquare_220308(costParameter, 
%           y, yhat, yError, t, t_point, FSRsim)
%       costParameter = parameter currently being evaluated
%       y = model prediction (specified to column)
%       yhat = experimental data
%       yError = error associated with the experimental data
%       t = vector of times outputed by ODE
%       t_point = vector of times to find
%       FSRsim = simulated grams of leucine synthesized (i.e., integrated MPS
%           curve)
    
    %%
    % Empty vectors of length equal to the length of experimental data
    y_test = zeros(1,length(t_point)); % simulated data at experimental data time
    SLS_ind = zeros(1,length(t_point));
    RSS_ind = zeros(1,length(t_point));
    
    % finding model simulated value at time of experimental data
    for i = 1:length(t_point) % i = number of data points for parameter

        if t_point(i) == 0      %if time 0 is a data point, it occurs at t(1)
            t_find = 1;
        else
            % locate position in t vector that is closest to experimental
            % data point
            t_find = find(t_point(i)<t & t<t_point(i)*10, 1);
        end
       
            % manually entering FSRsim value that was calculated outside function
            if costParameter == 46
                [y_test(i)] = FSRsim;      
            
            % p-Akt(S): combining the concentrations of both serine species 
            elseif costParameter == 21
                [y_test(i)] = sum(y(t_find,:));
    
            else
                % adding parameter values at same time point as the experimental data to test 
                [y_test(i)] = y(t_find);   
            end
        
        % Calculating SLS and RMS for individual data points
        [SLS_ind(i)] = (y_test(i) - yhat(i)).^2 / yError(i).^2;
        [RSS_ind(i)] = (y_test(i) - yhat(i)).^2;
    
    end
        
    % RMS: Square root of the sum of individual SLS
    RoomMeanSquare = sqrt(sum(SLS_ind));
    % RSS: Sum of the individual RSS
    ResidualSumSquares = sum(RSS_ind);
        
    % Normalizing SLS and RSS to number of data points of the specific
    % species
    RoomMeanSquare_norm = RoomMeanSquare/length(t_point);
    ResidualSumSquares_norm = ResidualSumSquares/length(t_point);
    
end
