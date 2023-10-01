function [mps_integral_grams, p_p70S6K_integral] = plot_Figure_5_230926(figNum, t, x_conc, x_mass, ...
    RunDuration, EquilibriumDuration, kValues, interventionDuration, style_variation, color_variation, ...
    Kd_Variations_legend, IV_variation_legend, lineThickness_variation)
% Plotting function that creates Figure 5a-c. The plotting function simulates
% leucine, Akt, p70S6K, and protein balance with variable levels of p70S6K
% and p-p70S6K and varaible levels of mTORC1 kinase activity. 
%   [mps_integral_grams, p_p70S6K_integral] = plot_Figure_5_230221(figNum, 
%       t, x_conc, x_mass, RunDuration, EquilibriumDuration, kValues, 
%       interventionDuration, style_variation, color_variation, 
%       Kd_Variations_legend, IV_variation_legend)
%   figNum = figure number for the MATLAB plot
%   t = scalar time outputted from the ODE solver
%   x_conc = column vector of all species across the simulated duration in 
%   units of mol/L as outputted from the ODE solver
%   x_mass = column vector of all species across the simulated duration in 
%   units of moles as outputted from the ODE solver
%   RunDuration = simulated run duration (minutes)
%   EquilibriumDuration = duration of the equilibrium period to allow model
%   to acheive steady state prior to stimulation (minutes)
%   kValues = k-value vector use to simulate model
%   interventionDuration = the length of the intervention duration from the
%   extracted experimental data protocol, which is used to select the
%   plotted x-axis (time, minutes)
%   style_variation = line style variation for each variation in p70S6K
%   levels
%   color_variation = line color variation for each variation in p70S6K
%   levels
%   Kd_variations_legend = list stating the variations in the mTORC1 kinase
%   activity. 
%   IV_variation_legend = list stating the variations in p70S6K levels

%% Simulation duration to plot
simDuration = 200; 
    
%% Setting up subTightPlot function (allows for better plotting)

gap = [0.16 0.045]; % gap between subplots [vert, horiz]
marg_h = [0.1 0.02*3];
marg_w = [0.1/2 0.01];

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    
%% Data cleaning   
% Position of species in data frame
IR_b = 9;
IR_b_Y1164 = 10;
IR_b_ref = 11;
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

%% FSR calculation
% FSR/MPS; integral for FSR comparison
mps = x_mass(:,44);
t_endRun = find(t>=EquilibriumDuration+interventionDuration & t<EquilibriumDuration+interventionDuration+20);
mps_integral = cumtrapz(t(t_endEq(1):t_endRun(1)),...
    mps(t_endEq(1):t_endRun(1)));
mps_integral=mps_integral(end);
LeucineMolarMass = 1/131.17; %mol/g
mps_integral_grams = mps_integral/LeucineMolarMass;

%% p-p70S6K AUC
p_p70S6K_integral = cumtrapz(t_endEq(1):t_endRun(1), ...
    x_conc(t_endEq(1):t_endRun(1), p70S6KT389_pos));
p_p70S6K_integral = p_p70S6K_integral(end);

%% x_conc vector updated to include mole/min values for 3-pool and muscle balance parameters
% required for the 'x' parameter inputted in the costCalculation function
x_conc_moleMin = x_conc;
x_conc_moleMin(:,41) = x_mass(:,41); % F_ma in units of mole/min
x_conc_moleMin(:,42) = x_mass(:,42); % F_m0 in units of mole/min
x_conc_moleMin(:,43) = nan; % NB: not correct calculation
x_conc_moleMin(:,44) = x_mass(:,44); % MPS in units of mole/min
x_conc_moleMin(:,45) = x_mass(:,45); % MPB in units of mole/min

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

fontSize_title = 20.5;
fontSize_axes = 16;
fontSize_units = 13.5;
fontSize_inPlotText = 14;

figure(figNum(1));

set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperPosition', [0 0 30*expand_x 15*expand_y]); %x_width=10cm y_width=15cm

%% p70S6K
subplotTight(height_y, length_x, 1:3)
p70S6K_total = x_conc(:,p70S6K_pos) + x_conc(:,p70S6KT389_pos);
lines = p70S6KT389_pos;
    plot(t, x_conc(:,lines), 'linestyle', style_variation, 'LineWidth', 1.5, 'Color', ...
        [0 0 1*color_variation], 'DisplayName', Kd_Variations_legend);
%     p1 = plot(t,(x_conc(:,LeucinePlasma_pos)), 'linestyle', style_variation,  'LineWidth', 1.5, ...
%         'Color',[0 0 1*color_variation], 'DisplayName', Kd_Variations_legend);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
        ax.YAxis.Exponent = -9; % y-axis unit in units of 10^-9
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylim_max = 8.9047e-09;
    ylim([0 ylim_max*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks(0 : 2e-9 : ylim_max*1.1)
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
%     plot (t,p70S6K_total, 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0.5 0.5 0.5]);
    p2 = plot(t,p70S6K_total, ...
        '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5], ...
        'HandleVisibility','off');
hold on

% mTORC1 kinase activity legend (only 1 legend can be plotted at a time)
leg = legend('show');
title(leg, 'mTORC1')

% legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northwest');

%% MPS and MPB
subplotTight(height_y, length_x, 4:6)
    plot (t, x_mass(:,MPS_pos), 'linestyle', style_variation,  'LineWidth', lineThickness_variation, 'Color',[0 0 1*color_variation]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
        ymin = -0.1e-5;
    ylim([ymin ymax])
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Protein Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1], 'EdgeColor',"none", ...
        'FaceAlpha',0.05, ...
        'LineStyle','none', ...
        'HandleVisibility','off')

% plot MPB
hold on
    p2 = plot(t, x_mass(:, MPB_pos), ...
        'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[1*color_variation 0 0], ...
        'HandleVisibility','off');
hold on

% FSR legend
if IV_variation_legend == "1x"
    leg = legend({'0.36', '0.46', '0.58'},'NumColumns', 1, 'Location','northeast'); % p70S6K @ 1x
elseif IV_variation_legend == "0.5x"
    leg = legend({'0.25', '0.32', '0.40'},'NumColumns', 1, 'Location','northeast'); % p70S6K @ 0.5x
elseif IV_variation_legend == "0.25x"
    leg = legend({'0.16', '0.21', '0.26', '0.38'},'NumColumns', 1, 'Location','northeast'); % p70S6K @ 0.25x
elseif IV_variation_legend == "0.125x"
    leg = legend({'0.10', '0.13', '0.16', '0.22'},'NumColumns', 1, 'Location','northeast'); % p70S6K @ 0.125x
end

title(leg, 'FSR (g leu)')

% legend({'MPS', 'MPB'}, 'Location','northwest');

%% NB
NetBalance = x_mass(:,44) - x_mass(:,45);
subplotTight(height_y, length_x, 7:9)
    plot(t, NetBalance,  ...
    'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 1*color_variation 0], ...
    'HandleVisibility','off');
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
        ymin = -0.35e-5;
    ylim([ymin ymax])
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Net Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on

end
