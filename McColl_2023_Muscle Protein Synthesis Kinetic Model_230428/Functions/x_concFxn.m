function [x_conc] = x_concFxn(x, x0, plasmaVolume, skeletalMuscleVolume)
% Converts values in units of moles from the model simulation to units of
% mol/L for plotting
%   [x_conc] = x_concFxn(x, x0, plasmaVolume, skeletalMuscleVolume)
%   x = simulated values in units of mol
%   x0 = initial values - used to format the array structure
%   plasmaVolume = liters of plasma volume as calculated from the
%   'plasmaVolumeFxn' function
%   skeletalMuscleVolume = liters of skeletal muscle volume as
%   calculated from the 'skeletalMuscleVolumeFxn' function

x_conc = nan(length(x), length(x0));

for n = 1:size(x,2)
    if any(n == 1:2)        % Digestive Tract
        x_conc(:,n) = x(:,n);
        
    elseif any(n == 3:5)    % Blood Plasma
        x_conc(:,n) = x(:,n)/plasmaVolume;
        
%     elseif any(n == 6:length(x0))   % Skeletal Muscle 
    elseif any(n == 6:size(x,2))   % Skeletal Muscle 
        x_conc(:,n) = x(:,n)/skeletalMuscleVolume;
    end

end

    


