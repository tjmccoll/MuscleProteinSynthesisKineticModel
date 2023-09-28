function [x0_mass] = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume)
% Converts initial values of model species from units of mol/L to moles
%   [x0_mass] = x0_massFxn(x0, plasmaVolume, skeletalMuscleVolume)
%   x0 = string of model species initial values in units of mol/L
%   plasmaVolume = liters of plasma volume as calculated from the
%   'plasmaVolumeFxn' function
%   skeletalMuscleVolume = liters of skeletal muscle volume as
%   calculated from the 'skeletalMuscleVolumeFxn' function

x0_mass = nan(length(x0),1);

for n = 1:length(x0)        % Digestive Tract
    if any(n == 1:2)
        x0_mass(n) = x0(n);
        
    elseif any(n == 3:5)    % Blood Plasma
        x0_mass(n) = x0(n)*plasmaVolume;
        
    elseif any(n == 6:length(x0))   % Skeletal Muscle 
        x0_mass(n) = x0(n)*skeletalMuscleVolume;
    end
end
