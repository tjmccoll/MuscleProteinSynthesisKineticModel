function [CostIndividual, CostTotal, mps_integral_grams] = plot_Figure_2bc_230424(figNum, t, ...
    x_conc, x_mass, ExpData, CostParameters, RunDuration, EquilibriumDuration, FilePath, ...
    FileName, kValues, interventionDuration)
% Plotting function that creates the Figure 2 b and c panels. Each panel
% plots the simulated plasma leucine, plasma insulin, p70S6K, and Protein
% Balance plots with associated experimental data and costs.
%   [CostIndividual, CostTotal, mps_integral_grams] = 
%       plot_Figure_2bc_230221(figNum, t, x_conc, x_mass, ExpData, 
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

gap = [0.16 0.045]; % gap between subplots [vert, horiz]
marg_h = [0.1 0.08];
marg_w = [0.1 0.01];

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
MPB_pos = 45;
FSR_pos = 46; %total synthesized leucine (calculated from FSR)

t_endEq = find(t==EquilibriumDuration);
AktS473_iv = x_conc(t_endEq(1), AktS473_pos);
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

%% 3-pool parameters; FSR calculation

F_ma = kValues(6)*x_mass(:,4); % r6 = F_m,a
r9 = kValues(9)*x_mass(:,8)./x_mass(:,10); % F_m,0 (1); MPB (IR_b mediated)
r11 = kValues(11)*x_mass(:,7); % F_m,0 (2)
F_m0 = r9+r11; 

% FSR/MPS; integral for FSR comparison
mps = kValues(15)*x_mass(:,6).*x_mass(:,28); % FSR, r15
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% Updating simulation data - Fm,a and Fm,0 added with correct units (mol/min)
    %allows for correct cost calculation (experimental data in mol/min)
x_conc_update = x_conc;
x_conc_update(:,41) = F_ma;
x_conc_update(:,42) = F_m0;
x_conc_update(:,44) = mps;

%% RMS, normalized

[Title, ~, ~, CostIndividual, CostTotal] = CostCalculation_230221(t, x_conc_update, AktS473_pos, p70S6KT389_pos, ...
    EquilibriumDuration, CostParameters, ExpData, mps_integral_grams);

%% Simulation plot

x_lim = [EquilibriumDuration-5, EquilibriumDuration+simDuration];
minuteSpacing=60;
xTicksSeq = 0:minuteSpacing:RunDuration;
xTicks = EquilibriumDuration+xTicksSeq;

opacity = 0.25;  
Angle = 0;%45;
jitter = 1.5;

length_x = 12;
height_y = 2;

expand_x = 1.4;
expand_y = 1.1;

fontSize_title = 19; %20.5;
fontSize_axes = 16;
fontSize_units = 13.5;
fontSize_inPlotText = 14;

figure(figNum(1));

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 30*expand_x 15*expand_y]); %x_width=10cm y_width=15cm

%% Leucine
subplotTight(height_y, length_x, 1:3)
lines = LeucinePlasma_pos;
    plot (t,x_conc(:,lines), '-', 'LineWidth', 1.5, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.4e-3])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks([0 : 0.3e-3 : 1.4e-3])
    title('Plasma Leucine', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on
if any(ExpDataPlot(:,1,1)==LeucinePlasma_pos)
    plotRow = find(ExpDataPlot(:,1,1) == LeucinePlasma_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1.25)
end   
hold on 
if any(CostIndividual(:,1)==LeucineInt_pos)
    plot(t,x_conc(:,LeucineInt_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
    
    plotRow = find(ExpDataPlot(:,1,1) == LeucineInt_pos);
    errorbar(ExpDataPlot(plotRow,:,2)-jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[1 0.5 0.5],'LineWidth',1.25)
end
hold on 
%---
if any(CostIndividual(:,1)==LeucinePlasma_pos) && any(CostIndividual(:,1)==LeucineInt_pos) % text if Fm,a and Fm,0 exist
    textRow_Leu_pl = find(CostIndividual(:,1) == LeucinePlasma_pos);
    textRow_Leu_int = find(CostIndividual(:,1) == LeucineInt_pos);
    textOutput = sprintf('Plasma, quantitative: %.2f\nIntracellular , tracer data: %.2f', ...
        CostIndividual(textRow_Leu_pl,2),CostIndividual(textRow_Leu_int,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.98*1.4e-3, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
%---
elseif any(CostIndividual(:,1)==LeucinePlasma_pos)
    textRow = find(CostIndividual(:,1) == LeucinePlasma_pos);
%     textOutput = sprintf('V(p) = %.2f', CostIndividual(textRow,2)); % change '%.2f' to '%.2e' for scientific notation
    textOutput = sprintf('Quantitative\n%.2f', CostIndividual(textRow,2)); % change '%.2f' to '%.2e' for scientific notation
    text(0.99*(EquilibriumDuration+simDuration), 0.98*1.4e-3, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end
hold off
% legend({'Plasma', 'Intracellular'});

%% Insulin
subplotTight(height_y, length_x, 4:6)
    plot(t,(x_conc(:,Insulin_pos)), '-', 'LineWidth',1.5, 'Color',[0 0 1])
    grid on 
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1.2e-10])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Plasma Insulin', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
if any(ExpDataPlot(:,1,1)==Insulin_pos)
    plotRow = find(ExpDataPlot(:,1,1) == Insulin_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1.25)
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

%% p70S6K
subplotTight(height_y, length_x, 7:9)
p70S6K_total = x_conc(:,p70S6K_pos) + x_conc(:,p70S6KT389_pos);
lines = p70S6KT389_pos;
    plot(t, p70S6K_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 p70S6K_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
plot (t,x_conc(:,lines), '-', 'LineWidth',1.5, 'Color',[0 0 1]);
hold on
if any(ExpDataPlot(:,1,1)==p70S6KT389_pos)
    plotRow = find(ExpDataPlot(:,1,1) == p70S6KT389_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1.25)
end
hold on
if any(CostIndividual(:,1)==p70S6KT389_pos)
    textRow = find(CostIndividual(:,1) == p70S6KT389_pos);
    textOutput = sprintf('Immunoblots\n%.2f ', CostIndividual(textRow,2));
    text(0.9925*(EquilibriumDuration+simDuration), 0.73*p70S6K_total(1,1)*1.1, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end
hold off
legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northeast');

%% FSR
subplotTight(height_y, length_x, 10:12)
    plot (t, mps, '-', 'LineWidth',1.5, 'Color',[0 0 1]) 
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
    ylim([-.1e-5 ymax])
    ylabel('Moles/min', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Protein Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)
    
hold on
if any(ExpDataPlot(:,1,1)==MPS_pos) %plotting MPS rates* if data exists (Drummond)
    plotRow = find(ExpDataPlot(:,1,1) == MPS_pos);
    errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
        ExpDataPlot(plotRow,:,3), ...
        ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth',1.25)
end

% plot MPB
hold on
    plot(t, x_mass(:, MPB_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])                         
% plot NB
hold on
    NetBalance = x_mass(:,44) - x_mass(:,45);
    plot(t, NetBalance,  '-', 'LineWidth',1.5, 'Color',[0 1 0])

% MPS AUC shading
hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1])
    alpha(0.05) % opacity

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
    text(0.9925*(EquilibriumDuration+simDuration), 0.85*ymax, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end

% FSR cost
if any(CostIndividual(:,1)==FSR_pos)
    textRow = find(CostIndividual(:,1) == FSR_pos);
    textOutput = sprintf('%.2f', CostIndividual(textRow,2));
    text(0.99*(EquilibriumDuration+simDuration), 0.625*ymax, textOutput, ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'right',...
        'FontSize', fontSize_inPlotText)
end

% simulated FSR value (grams) (value + background rectangle)
MPS_simPeriod = mps(t_endEq(1):t_endRun(1));
mps_max = max(MPS_simPeriod);
mps_max_t = find(mps == mps_max);    

xScale = 0.175;
yScale = 0.37e-5;
rectangle('Position', [t(mps_max_t)-xScale*simDuration, 0.4*mps_max-yScale, 2*xScale*simDuration, 2*yScale], ...
    'FaceColor',[1 1 1 0.8], 'EdgeColor',[1 1 1])

textOutput = sprintf('Simulated FSR:\n%.2f g leucine', mps_integral_grams);
text(t(mps_max_t), 0.4*mps_max, textOutput, ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', 'center',...
        'FontSize', 11) %fontSize_inPlotText)

legend({'MPS', 'MPB', 'NB'}, 'Location','northeast', 'NumColumns',3)

%% Save plot

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figure1 = figure(figNum(1));
    figfile1 = fullfile(FilePath, FileName);
    saveas(figure1, figfile1, 'epsc' );
else
end
