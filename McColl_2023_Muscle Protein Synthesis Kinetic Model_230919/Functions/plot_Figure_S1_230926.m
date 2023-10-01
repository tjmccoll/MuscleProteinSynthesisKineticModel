function [CostIndividual, CostTotal, mps_integral_grams] = plot_Figure_S1_230926(figNum, t, x_conc, x_mass, ExpData, ...
    CostParameters, RunDuration, EquilibriumDuration, FilePath, FileName, Protocol, kValues, interventionDuration)
% Plotting function that creates Supplementary Figure S2. The function
% creates the comprehensive figure panel for each individual validation
% dataset, which includes panels for leucine, plasma insulin, Akt, p70S6K,
% 3-pool model parameters, and protein balance.
%   [CostIndividual, CostTotal, mps_integral_grams] = 
%       plot_Figure_S2_230221(figNum, t, x_conc, x_mass, ExpData, 
%       CostParameters, RunDuration, EquilibriumDuration, FilePath, 
%       FileName, Protocol, kValues, interventionDuration)
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
%   FilePath = file directory where figure will be saved
%   FileName = name of figure to be saved in the FilePath 
%   kValues = k-value vector use to simulate model
%   interventionDuration = the length of the intervention duration from the
%   extracted experimental data protocol, which is used to select the
%   plotted x-axis (time, minutes)

%% Simulation duration to plot
simDuration = interventionDuration*1.1; 

%% Setting up subTightPlot function (allows for better plotting)

% within plot x and y labels removed 
gap = [0.0775 0.05]; % gap between subplots [vert, horiz] - **with second ylabel**
marg_h = [0.05 0.032];
marg_w = [0.0455 0.001];

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    % marg_h(x,y) = adjusts upper and lower external boundaries of plot. x = lower border,
    % y = upper border
    % marg_w(x,y) = adjusts lateral boundaries of plot. x = left border,
    % y = right border 
    
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
MPB_pos = 45;
FSR_pos = 46; %total synthesized leucine (calculated from FSR)

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
        ExpDataPlot(n, 1:length(rows),4) = sqrt(AktS473_iv^2 * ExpData(rows, 4).^2); % Transforming SE data 
            
    elseif ExpDataPlot(n,1,1) == p70S6KT389_pos
        rows = find(ExpData(:,1) == CostParameters(n));
        ExpDataPlot(n, 1:length(rows),2) = ExpData(rows, 2) + EquilibriumDuration; %Adjusted to Equilibrium period
        ExpDataPlot(n, 1:length(rows),3) = ExpData(rows, 3) *p70S6KT389_iv;
        ExpDataPlot(n, 1:length(rows),4) = sqrt(p70S6KT389_iv^2 * ExpData(rows, 4).^2); % Transforming SE data
    
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

%% RMS, normalized

[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230920(t, x_conc_moleMin, [AktS473_pos, AktS473T308_pos], p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);

%% Simulation plot

x_lim = [EquilibriumDuration-5, EquilibriumDuration+simDuration];
minuteSpacing=60;
xTicksSeq = 0:minuteSpacing:RunDuration;
xTicks = EquilibriumDuration+xTicksSeq;

opacity = 0.25;  
Angle = 0;%45;
jitter = 1.5;

length_x = 8;
height_y = 3;

fontSize_title = 15; 
fontSize_axes = 10; 
fontSize_units = 9; 
fontSize_inPlotText = 9;

figure(figNum(1));

set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 8.5 7.333]); %x_width=10cm y_width=15cm


%Leucine
subplotTight(height_y, length_x, 1:4)
lines = LeucinePlasma_pos;
    plot (t,x_conc(:,lines), '-', 'LineWidth', 1.25, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.4e-3])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('')
    yticks([0 : 0.3e-3 : 1.4e-3])
    title('Leucine', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on
    lines = LeucineInt_pos;
    plot(t, x_conc(:,lines), '-', 'LineWidth',1.25, 'Color',[1 0 0])

hold on
if any(ExpDataPlot(:,1,1)==LeucinePlasma_pos)
    plotRow = find(ExpDataPlot(:,1,1) == LeucinePlasma_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1)
end   

hold on
if any(CostIndividual(:,1)==LeucineInt_pos)
    plotRow = find(ExpDataPlot(:,1,1) == LeucineInt_pos);
    errorbar(ExpDataPlot(plotRow,:,2)-jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[1 0.5 0.5],'LineWidth',1)
end

hold on 
%---
if any(CostIndividual(:,1)==LeucinePlasma_pos) && any(CostIndividual(:,1)==LeucineInt_pos) % text if Fm,a and Fm,0 exist
    textRow_Leu_pl = find(CostIndividual(:,1) == LeucinePlasma_pos);
    textRow_Leu_int = find(CostIndividual(:,1) == LeucineInt_pos);
    textOutput = sprintf('Plasma leucine, quantitative: %.2f\nIntracellular leucine, tracer data: %.2f', ...
        CostIndividual(textRow_Leu_pl,2),CostIndividual(textRow_Leu_int,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.82*1.4e-3, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
%---
elseif any(CostIndividual(:,1)==LeucinePlasma_pos)
    textRow = find(CostIndividual(:,1) == LeucinePlasma_pos);
    textOutput = sprintf('Quantitative\n%.2f', CostIndividual(textRow,2)); % change '%.2f' to '%.2e' for scientific notation
    text(0.99*(EquilibriumDuration+simDuration), 0.82*1.4e-3, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end
hold off
legend({'Plasma', 'Intracellular'},'Location','northeast', 'NumColumns',2);

%% Insulin
subplotTight(height_y, length_x, 5:8)
    plot(t,(x_conc(:,Insulin_pos)), '-', 'LineWidth',1.25, 'Color',[0 0 1])
    grid on 
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.2e-10])
    ylabel('')
    xlabel('')
    title('Plasma Insulin', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
if any(ExpDataPlot(:,1,1)==Insulin_pos)
    plotRow = find(ExpDataPlot(:,1,1) == Insulin_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1)
end
hold on
if any(CostIndividual(:,1)==Insulin_pos)
    textRow = find(CostIndividual(:,1) == Insulin_pos);
    textOutput = sprintf('Quantitative\n%.2f', CostIndividual(textRow,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.98*1.2e-10, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end
hold off

%% 3-pool parameters (Fm,a & Fm,0)
subplotTight(height_y, length_x, 17:20)
    plot (t, x_mass(:,F_ma_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1]) % F_m,a
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([-1E-5 18E-5])
    ylabel('Moles/min', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('3-Pool Model', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on
    plot (t, x_mass(:,F_m0_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0]) % F_m,0
hold on 
if any(ExpDataPlot(:,1,1)==F_ma_pos) % plotting F_m,a exp data
    plotRow = find(ExpDataPlot(:,1,1) == F_ma_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth', 1)
end   
hold on 
if any(ExpDataPlot(:,1,1)==F_m0_pos) % plotting F_m,0 exp data
    plotRow = find(ExpDataPlot(:,1,1) == F_m0_pos);
    errorbar(ExpDataPlot(plotRow,:,2)-jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[1 0.5 0.5],'LineWidth', 1)
end
%---
if any(CostIndividual(:,1)==F_ma_pos) && any(CostIndividual(:,1)==F_m0_pos) % text if Fm,a and Fm,0 exist
    textRow_Fma = find(CostIndividual(:,1) == F_ma_pos);
    textRow_Fm0 = find(CostIndividual(:,1) == F_m0_pos);
    textOutput = sprintf('Tracer\nF_{m,a}: %.2f\nF_{m,0}: %.2f', ...
        CostIndividual(textRow_Fma,2),CostIndividual(textRow_Fm0,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.81*18E-5, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
    
elseif any(CostIndividual(:,1)==F_ma_pos) % text for F_m,a cost
    textRow = find(CostIndividual(:,1) == F_ma_pos);
    textOutput = sprintf('Tracer data\nF_{m,a}: %.2f', CostIndividual(textRow,2));
    text(0.985*(EquilibriumDuration+simDuration), 0.75*18E-5, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
    
elseif any(CostIndividual(:,1)==F_m0_pos) % text for F_m,0 cost
    textRow = find(CostIndividual(:,1) == F_m0_pos);
    textOutput = sprintf('F_{m,0}: %.2f', CostIndividual(textRow,2));
    text(EquilibriumDuration+RunDuration-RunDuration/2, 0.55*18E-5, textOutput)
end
legend({'F_{m,a}', 'F_{m,0}'},'NumColumns',2);
hold off

%% Akt
subplotTight(height_y, length_x, 9:12)
Akt_total = x_conc(:,Akt_pos) + x_conc(:,AktT308_pos) + x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
    plot(t,Akt_total, '-', 'LineWidth', 1.25, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
        ax.YAxis.Exponent = -9; % y-axis unit in units of 10^-9
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 Akt_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('')
    yticks([0 : 2e-9 : Akt_total(1,1)*1.1])
    title('Akt', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
% plotting both Akt species with phosphorylated serine
% [p-Akt(S) & p-Akt(S,T)]
Akt_pS = x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
    plot (t,Akt_pS, '-', 'LineWidth', 1.5, 'Color',[0 0 1]);
hold on 
if any(ExpDataPlot(:,1,1)==AktS473_pos)
    plotRow = find(ExpDataPlot(:,1,1) == AktS473_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth', 1)
end
hold on
if any(CostIndividual(:,1)==AktS473_pos)
    textRow = find(CostIndividual(:,1) == AktS473_pos);
    textOutput = sprintf('Immunoblots\n%.2f ', CostIndividual(textRow,2));
    text(0.9925*(EquilibriumDuration+simDuration), 0.73*Akt_total(1,1)*1.1, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end
legend({'Total Akt', 'p-Akt^{S473}'});
hold off

%% p70S6K
subplotTight(height_y, length_x, 13:16)
p70S6K_total = x_conc(:,p70S6K_pos) + x_conc(:,p70S6KT389_pos);
lines = p70S6KT389_pos;
    plot(t, p70S6K_total, '-', 'LineWidth', 1.25, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
        ax.YAxis.Exponent = -9; % y-axis unit in units of 10^-9
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        % ifelse loop to adjust for p70S6K ylim in Glynn and Dickinson
        if Protocol == "Valid., Glynn (2010) 1.85g"
            ylim_max = 17e-9; %* Glynn et al.: 1.85 g validation dataset ylim value
            ylim([0 ylim_max])
            yticks([0 : 3e-9 : ylim_max])
        elseif Protocol == "Valid., Dickinson (2011) cntl"
            ylim_max = 20e-9; %* Dickinson et al.: validation ylim value
            ylim([0 ylim_max])
            yticks([0 : 3e-9 : ylim_max])
        else
            ylim_max = p70S6K_total(1,1)*1.1;
            ylim([0 ylim_max]) % ylim for all other validation datasets
            yticks([0 : 2e-9 : ylim_max])
        end
    ylabel('')
    xlabel('')
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot (t,x_conc(:,lines), '-', 'LineWidth',1.25, 'Color',[0 0 1]);

hold on
if any(ExpDataPlot(:,1,1)==p70S6KT389_pos)
    plotRow = find(ExpDataPlot(:,1,1) == p70S6KT389_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1)
end
hold on
if any(CostIndividual(:,1)==p70S6KT389_pos)
    textRow = find(CostIndividual(:,1) == p70S6KT389_pos);
    textOutput = sprintf('Immunoblots\n%.2f ', CostIndividual(textRow,2));
        % ifelse loop to adjust for p70S6K ylim in Glynn and Dickinson
        if Protocol == "Valid., Glynn (2010) 1.85g"
            text(0.9925*(EquilibriumDuration+simDuration), 0.73*ylim_max, textOutput, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
                'FontSize', fontSize_inPlotText)
        elseif Protocol == "Valid., Dickinson (2011) cntl"
            text(0.9925*(EquilibriumDuration+simDuration), 0.73*ylim_max, textOutput, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
                'FontSize', fontSize_inPlotText)
        else
            text(0.9925*(EquilibriumDuration+simDuration), 0.73*p70S6K_total(1,1)*1.1, textOutput, ...
                'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
                'FontSize', fontSize_inPlotText)
        end
end

legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northeast');
hold off

%% FSR
subplotTight(height_y, length_x, 21:24)
    plot (t, mps, '-', 'LineWidth',1.25, 'Color',[0 0 1]) 
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
    ylim([-0.2e-5 ymax])
    ylabel('')
%     ylabel('Moles/min', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Protein Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)

% plot MPB and NB
hold on
    plot(t, x_mass(:, MPB_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])

    NetBalance = x_mass(:,44) - x_mass(:,45);
    plot(t, NetBalance,  '-', 'LineWidth',1.25, 'Color',[0 1 0])

% MPS AUC shading
hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1])
    alpha(0.05) % opacity
    
hold on
if any(CostIndividual(:,1)==MPS_pos)
    textRow = find(CostIndividual(:,1) == MPS_pos);
    textOutput = sprintf('Drummond (2010): MPS rate\n%.2f ', CostIndividual(textRow,2));
    text(0.9925*(EquilibriumDuration+simDuration), 0.98*ymax, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
    
elseif any(ExpDataPlot(:,1,1)==FSR_pos) %plotting experimental data for total leucine synthesized to protein (FSR)
    plotRow = find(ExpDataPlot(:,1,1) == FSR_pos);
    expFSR = ExpDataPlot(plotRow,1,3);
    expFSR_se = ExpDataPlot(plotRow,1,4);
    
    textOutput = sprintf('FSR \n%.2f \x00B1 %.2f g leucine', expFSR, expFSR_se);
    text(0.9925*(EquilibriumDuration+simDuration), 0.82*ymax, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end

% FSR cost
if any(CostIndividual(:,1)==FSR_pos)
    textRow = find(CostIndividual(:,1) == FSR_pos);
    textOutput = sprintf('%.2f', CostIndividual(textRow,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.63*ymax, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end

% simulated FSR value (grams)
MPS_simPeriod = mps(t_endEq(1):t_endRun(1));
mps_max = max(MPS_simPeriod);
mps_max_t = find(mps == mps_max);    

xScale = 0.15;
yScale = 0.37e-5;
rectangle('Position', [t(mps_max_t)-xScale*simDuration, 0.5*mps_max-yScale, 2*xScale*simDuration, 2*yScale], ...
    'FaceColor',[1 1 1 0.8], 'EdgeColor',[1 1 1])

textOutput = sprintf('Simulated FSR:\n%.2f g leucine', mps_integral_grams);
text(t(mps_max_t), 0.5*mps_max, textOutput, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center',...
        'FontSize', 11) %fontSize_inPlotText)

legend({'MPS', 'MPB', 'NB'}, 'Location','northeast', 'NumColumns',3)
hold off

%% Save plot

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figure1 = figure(figNum(1));
    figfile1 = fullfile(FilePath, FileName);
        saveas(figure1, figfile1, 'epsc');

else
end


