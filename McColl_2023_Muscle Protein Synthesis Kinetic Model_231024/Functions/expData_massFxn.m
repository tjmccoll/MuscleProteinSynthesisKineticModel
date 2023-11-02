function [expData_mass] = expData_massFxn(expData, ExpParameters, plasmaVolume, skeletalMuscleVolume)
% Converting experimental data from units of mol/L to units of mol.
% Phospho-data was extracted in units of fold-change, which were converted
% to physiological units in the 'plot_...' functions. FSR data was
% extracted and converted to grams in the experimental data excel file. The
% FSR data is converted to moles in the 'plot_...' functions.
%   [expData_mass] = expData_massFxn(expData, ExpParameters, plasmaVolume,
%       skeletalMuscleVolume)
%   expData = extracted experimental data. See experimental data
%   excel file for details.
%   ExpParameters = parameter numbers that have corresponding experimental
%   data 
%   plasmaVolume = liters of plasma volume as calculated from the
%   'plasmaVolumeFxn' function
%   skeletalMuscleVolume = liters of skeletal muscle volume as calculated 
%   from the 'skeletalMuscleVolumeFxn' function

expData_mass = nan(length(expData), 4);

for n = 1:length(ExpParameters)       
    if any(ExpParameters(n) == 3:5)         %Blood Plasma - units of mol/L, converted to mol
        rows = find(expData(:,1) == ExpParameters(n));
        expData_mass(rows, 1) = ExpParameters(n);
        expData_mass(rows, 2) = expData(rows, 2);
        expData_mass(rows, 3) = expData(rows, 3) * plasmaVolume;
        expData_mass(rows, 4) = expData(rows, 4) * plasmaVolume;
        
    elseif any(ExpParameters(n) == 6:8)    %Skeletal Muscle Volume - units of mol/L, converted to mol
        rows = find(expData(:,1) == ExpParameters(n));
        expData_mass(rows, 1) = ExpParameters(n);
        expData_mass(rows, 2) = expData(rows, 2);
        expData_mass(rows, 3) = expData(rows, 3) * skeletalMuscleVolume;
        expData_mass(rows, 4) = expData(rows, 4) * skeletalMuscleVolume;
        
    elseif any(ExpParameters(n) == [21,28])    %Phospho-data - units of fold change
        rows = find(expData(:,1) == ExpParameters(n));
        expData_mass(rows, 1) = ExpParameters(n);
        expData_mass(rows, 2) = expData(rows, 2);
        expData_mass(rows, 3) = expData(rows, 3);
        expData_mass(rows, 4) = expData(rows, 4);
        
   elseif any(ExpParameters(n) == [40,41,42])    %Fm,0 / Fm,a / NB - units of mol/min
        rows = find(expData(:,1) == ExpParameters(n));
        expData_mass(rows, 1) = ExpParameters(n);
        expData_mass(rows, 2) = expData(rows, 2);
        expData_mass(rows, 3) = expData(rows, 3);
        expData_mass(rows, 4) = expData(rows, 4);      
        
    elseif any(ExpParameters(n) == 46)          %FSR - units of grams
        rows = find(expData(:,1) == ExpParameters(n));
        expData_mass(rows, 1) = ExpParameters(n);
        expData_mass(rows, 2) = expData(rows, 2);
        expData_mass(rows, 3) = expData(rows, 3);
        expData_mass(rows, 4) = expData(rows, 4); 
    end

end
