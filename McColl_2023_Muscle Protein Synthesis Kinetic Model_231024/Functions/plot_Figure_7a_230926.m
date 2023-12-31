function [mps_integral_grams, p_p70S6K_integral] = plot_Figure_7a_230926(figNum, t, x_conc, x_mass, RunDuration, ...
    EquilibriumDuration, kValues, interventionDuration, style_variation, color_variation, Kd_Variations_legend)
% Plotting function that creates the Figure 7a panel. The plot simulates
% leucine, IRS1, p70S6K, protein balance, and Akt (not shown in manuscript
% figure) with variying changes to the kinetic parameter controlling p70S6K
% mediated IRS1 feedback.
%   [mps_integral_grams, p_p70S6K_integral] = plot_Figure_7a_230221(figNum,
%       t, x_conc, x_mass, RunDuration, EquilibriumDuration, kValues, 
%       interventionDuration, style_variation, color_variation, 
%       Kd_Variations_legend)
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
%   style_variation = line style variation for each variation in feedback
%   color_variation = line color variation for each variation in feedback
%   Kd_variations_legend = list stating the variations for the k-value

%% Simulation duration to plot
simDuration = 200; 
    
%% Setting up subTightPlot function (allows for better plotting)

gap = [0.16*.75 0.045]; % gap between subplots [vert, horiz]
marg_h = [0.1*.8 0.02*3];
marg_w = [0.1/2 0.01];

 subplotTight = @(m,n,p) subtightplot (m, n, p, gap, marg_h, marg_w);
    % @(m,n,p): m&n = the number of subplots across and down, p = location
    % of subplot in overall plot
    
%% Data cleaning   
% Position of species in data frame
IR_b = 9;
IR_b_Y1164 = 10;
IR_b_ref = 11;
IRS1_pos = 12;
IRS1_pY_pos = 13;
IRS1_pS_pos = 14;
IRS1_PI3K_pos = 16;
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

PI3K_pos = 15;
IRS1_PI3K_pos = 16;
PDK1_pos = 17;
PDK1_p_pos = 18;
TSC_pos = 23;
TSC_p_pos = 24;
mTORC1_inactive_pos = 25;
mTORC1_active_pos = 26;
mTORC2_pos = 31;
mTORC2_p_pos = 32;

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

p_p70S6K_integral = cumtrapz(t(t_endEq(1):t_endRun(1)), ...
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
    
%% Leucine
subplotTight(height_y, length_x, 1:3)
    p1 = plot(t,(x_conc(:,LeucinePlasma_pos)), 'linestyle', style_variation,  'LineWidth', 1.5, ...
        'Color',[0 0 1*color_variation], 'DisplayName', Kd_Variations_legend);
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
    title('Leucine', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on
    p2 = plot(t,x_conc(:,LeucineInt_pos), 'linestyle', style_variation,  'LineWidth', 1.5, ...
        'Color',[1*color_variation 0 0], ...
        'HandleVisibility','off');
% ***
leg = legend('show');

% add pl vs int legend (mute above leg, run with only 1 kinase value)
% legend([p1 p2], {'Plasma', 'Intracellular'},'Location','northwest');

%% IRS1
subplotTight(height_y, length_x, 4:6)
IRS1_total = x_conc(:,IRS1_pos) + x_conc(:,IRS1_pS_pos) + x_conc(:,IRS1_pY_pos) + x_conc(:,IRS1_PI3K_pos);
    plot(t,IRS1_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;    
%         ax.YAxis.Exponent = -9; % y-axis unit in units of 10^-9
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 IRS1_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('IRS1', 'fontSize', fontSize_title)
    xtickangle(Angle)
% hold on 
%     plot(t,x_conc(:, IRS1_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
% hold on 
%     plot(t,x_conc(:, IRS1_pS_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[1*color_variation 0 0]);
% hold on 
%     plot(t,x_conc(:, IRS1_pY_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 1*color_variation 0]);
% 
% legend({'IRS1 total', 'IRS1', 'p-IRS1^S', 'p-IRS1^{Y}'}, 'location', 'northeast');

hold on 
    plot(t,x_conc(:, IRS1_pS_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
hold on 
    plot(t,x_conc(:, IRS1_pY_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[1*color_variation 0 0]);

legend({'IRS1 total', 'p-IRS1^S', 'p-IRS1^{Y}'}, 'location', 'northeast');


%% p70S6K
subplotTight(height_y, length_x, 7:9)
p70S6K_total = x_conc(:,p70S6K_pos) + x_conc(:,p70S6KT389_pos);
lines = p70S6KT389_pos;
    plot(t, p70S6K_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ylim_max = 9.3e-8;
    ylim([0 ylim_max])
    ylim([0 p70S6K_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot (t,x_conc(:,lines), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
hold on
legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northeast');

%% FSR
subplotTight(height_y, length_x, 10:12)
    p1 = plot (t, x_mass(:,MPS_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4e-5;
        ymin = -0.3e-5;
    ylim([ymin ymax])
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Protein Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1], 'EdgeColor',"none", ...
        'FaceAlpha',0.05, ...
        'HandleVisibility','off')

% plot MPB
hold on
    p2 = plot(t, x_mass(:, MPB_pos), ...
        'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[1*color_variation 0 0], ...
        'HandleVisibility','off');
% plot NB
hold on
    NetBalance = x_mass(:,44) - x_mass(:,45);
    p3 = plot(t, NetBalance,  ...
        'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 1*color_variation 0], ...
        'HandleVisibility','off');
hold on
% leg = legend({'0.36', '0.36', '0.36'},'NumColumns', 1, 'Location','northeast'); % IRS1 feedback
leg = legend({'0.36', '0.36', '0.36', '0.36', '0.37'},'NumColumns', 1, 'Location','northeast'); % IRS1 feedback
title(leg, 'FSR (g leu)');

% legend({'MPS', 'MPB', 'NB'}, 'Location','northwest') % change 'handleVisibility' to 'on;

end


