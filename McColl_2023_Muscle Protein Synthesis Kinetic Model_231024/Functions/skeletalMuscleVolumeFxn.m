function [skeletalMuscleVolume] = skeletalMuscleVolumeFxn(gender, height, age, bia)
% Calculation of total skeletal muscle volume using the equation to
% estimate total skeletal muscle mass developed by Janssen et al. (2000)
% based on bioelectrical impedance analysis.
%   [skeletalMuscleVolume] = skeletalMuscleVolumeFxn(gender, height, age, bia)
%   Gender: male = 1, female = 0
%   Height: in units of cm
%   Age: in units of years
%   bia: Bioelectical Impedance (in units of Ohms)
%   Output: skeletalMuscleVolume in units of L

skeletalMuscleMass = ((height^2/bia * 0.401) + (gender * 3.825) + (age * -0.071)) + 5.102; % units of kg

skeletalMuscleDensity = 1.112; %g/mL (Ward et al., 2005); equivalent units of kg/L

% conversion of mass to volumen using the skeletal muscle density
skeletalMuscleVolume = skeletalMuscleMass/skeletalMuscleDensity; % kg/(kg/L) -> units of L

end