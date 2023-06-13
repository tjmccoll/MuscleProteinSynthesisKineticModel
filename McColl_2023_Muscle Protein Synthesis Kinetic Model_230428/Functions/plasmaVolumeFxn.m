function [plasmaVolume] = plasmaVolumeFxn(gender, height_cm, mass)
% Estimating the plasma volume using the Nadler et al. (1962) formula.
%   [plasmaVolume] = plasmaVolumeFxn(gender, height_cm, mass)
%   gender: male = 1, female = 0
%   height_cm: height of the individual in centimeters
%   mass: mass of the individual in kg

height_m = height_cm/100; %convert height from cm to m

if gender == 1      %male
    plasmaVolume = (0.3669 * height_m^3) + (0.03219 * mass) + 0.6041;
    
elseif gender == 0  %female
    plasmaVolume = (0.3561 * height_m^3) + (0.03308 * mass) + 0.6041;
    
end

    


