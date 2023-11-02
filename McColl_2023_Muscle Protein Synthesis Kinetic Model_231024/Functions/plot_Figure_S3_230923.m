function [CostIndividual, CostTotal, mps_integral_grams] = plot_Figure_S3_230923(figNum, t, x_conc, x_mass, ExpData, ...
    CostParameters, RunDuration, EquilibriumDuration, FilePath, FileName, kValues, interventionDuration)
% Plotting function that creates Supplementary Figure S1. The plot includes
% all model species that do not include associated experimental data.
%   [CostIndividual, CostTotal, mps_integral_grams] = 
%       plot_Figure_S1_230221(figNum, t, x_conc, x_mass, ExpData, 
%       CostParameters, RunDuration, EquilibriumDuration, FilePath, 
%       FileName, kValues, interventionDuration)
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
simDuration = 200; %200
    
%% Setting up subTightPlot function (allows for better plotting)

gap = [0.16*.3 0.045*0.75]; % gap between subplots [vert, horiz]
marg_h = [0.1*.4 0.08*.3];
marg_w = [0.1*.45 0.01*.9];

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    % marg_h(x,y) = adjusts upper and lower external boundaries of plot. x = lower border,
    % y = upper border
    % marg_w(x,y) = adjusts lateral boundaries of plot. x = left border,
    % y = right border 
    
%% Data cleaning   
% Position of species in data frame
Stomach_pos = 1;
Gut_pos = 2;
Insulin_pos = 3;
LeucinePlasma_pos = 4;
kicPlasma_pos = 5;
LeucineInt_pos = 6;
kicInt_pos = 7;
Protein_pos = 8;
IR_B_pos = 9;
IR_B_p_pos = 10;
IR_B_ref_pos = 11;
IRS1_pos = 12;
IRS1_pY_pos = 13;
IRS1_pS_pos = 14;
PI3K_pos = 15;
IRS1_PI3K_pos = 16;
PDK1_pos = 17;
PDK1_p_pos = 18;
Akt_pos = 19;
AktT308_pos = 20;
AktS473_pos = 21;
AktS473T308_pos = 22;
TSC_pos = 23;
TSC_p_pos = 24;
mTORC1_inactive_pos = 25;
mTORC1_active_pos = 26;
p70S6K_pos = 27;
p70S6KT389_pos = 28;
PI3K_var_pos = 29;
PI3K_var_p_pos = 30;
mTORC2_pos = 31;
mTORC2_p_pos = 32;

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

length_x = 9;
height_y = 4;

factor = 0.6;
fontSize_title = 20.5*factor;
fontSize_axes = 16*factor;
fontSize_units = 13.5*factor;
fontSize_inPlotText = 14*factor;

figure(figNum(1));

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 16.5*1.3 19.5*1.3]); %x_width=10cm y_width=15cm

%% Digestive Tract
subplotTight(height_y, length_x, 1:3)
lines = Stomach_pos;
    plot(t,x_conc(:,lines), '-', 'LineWidth', 1.5, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
        ax.YAxis.Exponent = -2;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 3e-2])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold') %, 'fontSize', fontSize_axes)
    yticks([0 : 1e-2 : 3e-2])
    title('Digestive Tract', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:,Gut_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'Stomach', 'Gut'}, 'location', 'northeast');

%% KIC
subplotTight(height_y, length_x, 4:6)
lines = kicPlasma_pos;
    plot(t,x_conc(:,lines), '-', 'LineWidth', 1.5, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 1e-4])
    title('KIC', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:,kicInt_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold on
% if any(ExpDataPlot(:,1,1)==kicPlasma_pos)
%     plotRow = find(ExpDataPlot(:,1,1) == kicPlasma_pos);
%     errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
%         ExpDataPlot(plotRow,:,3), ...
%         ExpDataPlot(plotRow,:,4), 'o', 'Color',[0.5 0.5 1],'LineWidth', 1.0)
% end
hold off
legend({'Plasma', 'Intracellular'}, 'location', 'northeast');

%% Insulin Receptor
subplotTight(height_y, length_x, 7:9)
IR_total = x_conc(:,IR_B_pos) + x_conc(:,IR_B_ref_pos) + x_conc(:,IR_B_p_pos);
    plot(t,IR_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 IR_total(1,1)*1.1])
    title('Insulin Receptor', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, IR_B_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, IR_B_p_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold on 
    plot(t,x_conc(:, IR_B_ref_pos), '-', 'LineWidth',1.25, 'Color',[0 1 0])
hold off
legend({'IR-\beta total', 'IR-\beta', 'p-IR-\beta^{Y}', 'IR-\beta_{Ref}'}, 'location', 'east');

%% IRS1
subplotTight(height_y, length_x, 10:12)
IRS1_total = x_conc(:,IRS1_pos) + x_conc(:,IRS1_pS_pos) + x_conc(:,IRS1_pY_pos) + x_conc(:,IRS1_PI3K_pos);
    plot(t,IRS1_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 IRS1_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('IRS1', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, IRS1_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, IRS1_pS_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold on 
    plot(t,x_conc(:, IRS1_pY_pos), '-', 'LineWidth',1.25, 'Color',[0 1 0])
hold off
legend({'IRS1 total', 'IRS1', 'p-IRS1^S', 'p-IRS1^{Y}'}, 'location', 'northeast');

%% PI3K
subplotTight(height_y, length_x, 13:15)
PI3K_total = x_conc(:,PI3K_pos) + x_conc(:,IRS1_PI3K_pos);
    plot(t,PI3K_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 PI3K_total(1,1)*1.1])
    title('PI3K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, PI3K_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, IRS1_PI3K_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'PI3K total', 'PI3K', 'PI3K p-IRS1^{Y}'}, 'location', 'northeast');

%% PDK1
subplotTight(height_y, length_x, 16:18)
PDK1_total = x_conc(:,PDK1_pos) + x_conc(:,PDK1_p_pos);
    plot(t,PDK1_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 PDK1_total(1,1)*1.1])
    title('PDK1', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, PDK1_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, PDK1_p_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'PDK1 total', 'PDK1', 'p-PDK1'}, 'location', 'northeast');

%% Akt
subplotTight(height_y, length_x, 19:21)
Akt_total = x_conc(:,Akt_pos) + x_conc(:,AktT308_pos) + x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
    plot(t,Akt_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 Akt_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks([0 : 0.4e-8 : Akt_total(1,1)*1.1])
    title('Akt', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:,Akt_pos), '-', 'LineWidth', 1.25, 'Color',[57 106 177]./255); %[0.41 0.51 0.22]);
hold on
    plot(t,x_conc(:,AktT308_pos), '-', 'LineWidth',1.25, 'Color', [204 37 41]./255); %'#c04e01');
hold on
    plot(t,x_conc(:,AktS473_pos), '-', 'LineWidth',1.25, 'Color', [107 76 154]./255); %'#9a0eea');
hold on
    plot(t,x_conc(:,AktS473T308_pos), '-', 'LineWidth',1.25, 'Color', [62 150 81]./255); %'#06c2ac');
hold on 
% if any(ExpDataPlot(:,1,1)==AktS473_pos)
%     plotRow = find(ExpDataPlot(:,1,1) == AktS473_pos);
%     errorbar(ExpDataPlot(plotRow,:,2)+jitter, ...
%         ExpDataPlot(plotRow,:,3), ...
%         ExpDataPlot(plotRow,:,4), 'o', 'Color',[107 76 154]./255,'LineWidth', 1.0)
% end
hold off
legend({'Total Akt', 'Akt', 'p-Akt^{T308}', 'p-Akt^{S473}', 'p-Akt^{S,T}'});

%% TSC
subplotTight(height_y, length_x, 22:24)
TSC_total = x_conc(:,TSC_pos) + x_conc(:,TSC_p_pos);
    plot(t,TSC_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 TSC_total(1,1)*1.1])
    title('TSC', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, TSC_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, TSC_p_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'TSC total', 'TSC', 'p-TSC'}, 'location', 'northeast');

%% mTORC1
subplotTight(height_y, length_x, 25:27)
mTORC1_total = x_conc(:,mTORC1_inactive_pos) + x_conc(:,mTORC1_active_pos);
    plot(t,mTORC1_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 mTORC1_total(1,1)*1.1])
    title('mTORC1', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, mTORC1_inactive_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, mTORC1_active_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'mTORC1 total', 'mTORC1_{inactive}', 'mTORC1_{active}'}, 'location', 'east');

%% PI3K variant
subplotTight(height_y, length_x, 28:30)
PI3K_var_total = x_conc(:,PI3K_var_pos) + x_conc(:,PI3K_var_p_pos);
    plot(t,PI3K_var_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 PI3K_var_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('PI3K Variant', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, PI3K_var_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, PI3K_var_p_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'PI3K_{var} total', 'PI3K_{var}', 'p-PI3K_{var}'}, 'location', 'east');

%% mTORC2
subplotTight(height_y, length_x, 31:33)
mTORC2_total = x_conc(:,mTORC2_pos) + x_conc(:,mTORC2_p_pos);
    plot(t,mTORC2_total, '-', 'LineWidth', 1.5, 'Color',[0 0 0])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 mTORC2_total(1,1)*1.1])
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('mTORC2', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot(t,x_conc(:, mTORC2_pos), '-', 'LineWidth',1.25, 'Color',[0 0 1])
hold on 
    plot(t,x_conc(:, mTORC2_p_pos), '-', 'LineWidth',1.25, 'Color',[1 0 0])
hold off
legend({'mTORC2 total', 'mTORC2', 'p-mTORC2^{S2481}'}, 'location', 'east');

%% leucine bound to protein
subplotTight(height_y, length_x, 34:36)
    plot(t,x_conc(:,Protein_pos), '-', 'LineWidth', 1.5, 'Color',[0 0 1])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units; 
        ax.YAxis.Exponent = 0;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
%     ylim([0.133965 0.13408])
    ylim([0.13396 0.13408])
%     yticks([0.133965: 5.7500e-05 : 0.13408])
    yticks([0.1338: 0.00004 : 0.13408])
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title({'Leucine Bound to Protein'}, 'fontSize', fontSize_title)
    ytickangle(90)

%% Save plot

prompt = 'Save figure? Y/N: ';
SaveFile = input(prompt,'s');

if SaveFile == 'Y'
    figure1 = figure(figNum(1));
    figfile1 = fullfile(FilePath, FileName);
    saveas(figure1, figfile1, 'epsc' );
    
else
end
