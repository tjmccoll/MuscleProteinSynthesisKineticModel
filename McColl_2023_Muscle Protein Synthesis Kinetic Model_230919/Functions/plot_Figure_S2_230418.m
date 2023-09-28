function [CostIndividual, CostTotal] = plot_Figure_S2_230418(figNum, t, x_conc, x_mass, ExpData, ...
    CostParameters, RunDuration, EquilibriumDuration, Protocol, kValues, interventionDuration, ...
    LeucineDose, i)
% Plotting function that creates Supplementary Figure 3. This function
% simulates the plasma leucine time-course for each dataset with measure
% plasma leucine data and integrates them into a single 3x3 figure.
%   [CostIndividual, CostTotal] = plot_Figure_S3_230221(figNum, t, x_conc, 
%       x_mass, ExpData, CostParameters, RunDuration, EquilibriumDuration, 
%       Protocol, kValues, interventionDuration, LeucineDose, i)
%   figNum = figure number for the MATLAB plot
%   t = scalar time outputted from the ODE solver
%   x_conc = column vector of all species across the simulated duration in 
%   units of mol/L as outputted from the ODE solver
%   x_mass = column vector of all species across the simulated duration in 
%   units of moles as outputted from the ODE solver
%   ExpData = extracte experimental data that is plotted with the
%   corresponding simualted time-course
%   CostParameters = Parameters with corresponding experimental data 
%   used to calculate RMSE
%   RunDuration = simulated run duration (minutes)
%   EquilibriumDuration = duration of the equilibrium period to allow model
%   to acheive steady state prior to stimulation (minutes)
%   Protocol = study name and year to add as title to each subplot
%   kValues = k-value vector use to simulate model
%   interventionDuration = the length of the intervention duration from the
%   extracted experimental data protocol, which is used to select the
%   plotted x-axis (time, minutes)
%   LeucineDose = quantity of leucine administered in the experimental
%   intervention, which is added as text to the plot. 
%   i = iteration in for loop; used to specify location of subplot

%% Simulation duration to plot
simDuration = 200; 
    
%% Setting up subTightPlot function (allows for better plotting)

gap = [0.08/1 0.0765/1]; % gap between subplots [vert, horiz] - **with second ylabel**
marg_h = [0.1 0.08];
marg_w = [0.1 0.01*3];

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    
%% Data cleaning   
% Position of species in data frame
AktS473_pos = 21;
Akt_pos = 19;
AktT308_pos = 20;
AktS473T308_pos = 22;
p70S6K_pos = 27;
p70S6KT389_pos = 28;
LeucinePlasma_pos = 4;
LeucineInt_pos = 6;
Insulin_pos = 3;
Protein_pos = 8;
kicPlasma_pos = 5;
kicInt_pos = 7;
F_ma_pos = 41;
F_m0_pos = 42;
MPS_pos = 44; %rate of MPS
FSR_pos = 45; %total synthesized leucine (calculated from FSR)

t_endEq = find(t==EquilibriumDuration);
% p-Akt(S) incorporates both Akt species with phospho serine residues
AktS473_iv = sum(x_conc(t_endEq(1), [AktS473_pos, AktS473T308_pos]));
p70S6KT389_iv = x_conc(t_endEq(1), p70S6KT389_pos);

%% Sorting Experimental Data
ExpDataPlot = nan(length(CostParameters), 13, 4);
ExpDataPlot(:,1,1) = CostParameters';

for n = 1:length(CostParameters)
    if ExpDataPlot(n,1,1) == AktS473_pos
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3) *AktS473_iv;
        ExpDataPlot(n, 1:length(rows),4) = 1.96 * sqrt(AktS473_iv^2 * ExpData(rows, 4).^2); % Transforming SE data
            
    elseif ExpDataPlot(n,1,1) == p70S6KT389_pos
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3) *p70S6KT389_iv;
        ExpDataPlot(n, 1:length(rows),4) = 1.96 * sqrt(p70S6KT389_iv^2 * ExpData(rows, 4).^2); % Transforming SE data
    
    else
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3);
        ExpDataPlot(n, 1:length(rows),4) = ExpData(rows, 4);
    end
end

%% FSR calculation
% FSR/MPS; integral for FSR comparison
mps = x_mass(:,44);
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% x_conc vector updated to include mole/min values for 3-pool and muscle balance parameters
% required for the 'x' parameter inputted in the costCalculation function
x_conc_moleMin = x_conc;
x_conc_moleMin(:,41) = x_mass(:,41); % F_ma in units of mole/min
x_conc_moleMin(:,42) = x_mass(:,42); % F_m0 in units of mole/min
x_conc_moleMin(:,43) = nan; % NB: not correct calculation
x_conc_moleMin(:,44) = x_mass(:,44); % MPS in units of mole/min
x_conc_moleMin(:,45) = x_mass(:,45); % MPB in units of mole/min

%% Root mean squares, normalized
[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230920(t, x_conc_moleMin, [AktS473_pos, AktS473T308_pos], p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);

%% Simulation plot

x_lim = [EquilibriumDuration-5, EquilibriumDuration+simDuration];
minuteSpacing=60;
xTicksSeq = 0:minuteSpacing:RunDuration;
xTicks = EquilibriumDuration+xTicksSeq;

opacity = 0.25;  
Angle = 0;
jitter = 1.5;

length_x = 3;
height_y = 3;

fontSize_title = 15; %19.5; %20.5;
fontSize_axes = 13; %16;
fontSize_units = 10; %13.5;
fontSize_inPlotText = 12; %14;

figure(figNum(1));

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8.5 11]); %x_width=8.5" y_width=11"

startLength_x = [1, 2, 3, 4, 5, 6, 7, 8, 9];
endLength_x = [1, 2, 3, 4, 5, 6, 7, 8, 9];

%% Leucine
subplotTight(height_y, length_x, startLength_x(i):endLength_x(i))
    plot(t,(x_conc(:,LeucinePlasma_pos)), '-', 'LineWidth',1.5, 'Color',[0 0 1])
    grid on 
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.4e-3])
    yticks([0 : 0.3e-3 : 1.4e-3])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title(Protocol, 'fontSize', fontSize_title)
    xtickangle(Angle)

    hold on 
    plotRow = find(ExpDataPlot(:,1,1) == LeucinePlasma_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1.25)
    
    hold on
    if any(CostIndividual(:,1)==LeucinePlasma_pos)
        textRow = find(CostIndividual(:,1) == LeucinePlasma_pos);
        textOutput = sprintf('Quantitative\n%.2f', CostIndividual(textRow,2));
        text(0.99*(EquilibriumDuration+simDuration), 0.98*1.4e-3, textOutput, ...
            'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
            'FontSize', fontSize_inPlotText)
    end

    LeuPlasma_simPeriod = x_conc(t_endEq(1):t_endRun(1),LeucinePlasma_pos);
    LeuPlasma_max = max(LeuPlasma_simPeriod);

    textLeucineDose = sprintf('%.2f g leu', LeucineDose);
    text((EquilibriumDuration), LeuPlasma_max+0.2e-3, textLeucineDose, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left',...
        'FontSize', fontSize_inPlotText)

end


