function [mps_integral_grams] = plot_Figure_3_6_230424(figNum, t, x_conc, x_mass, RunDuration, ...
    EquilibriumDuration, kValues, interventionDuration, style_variation, color_variation, ...
    Kd_Variations_legend, Kd_rxn)
% Plotting function that creates Figures 3 and 6. The plotting function
% includes subplots for leucine, plasma insulin, Akt, p70S6K, MPS, MPB, and
% net balance.
%   [mps_integral_grams] = plot_Figure_3_6_230221(figNum, t, x_conc, 
%       x_mass, RunDuration, EquilibriumDuration, kValues, 
%       interventionDuration, style_variation, color_variation, 
%       Kd_Variations_legend, Kd_rxn)
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
%   style_variation = line style variation for each variation of signalling
%   knockdown
%   color_variation = line color variation for each variation of signalling
%   knockdown
%   Kd_variations_legend = list stating the percent knockdown of signalling
%   Kd_rxn = reaction/k-value that is being knocked down

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
    p1 = plot(t,x_conc(:,LeucinePlasma_pos), 'linestyle', style_variation,  'LineWidth', 1.5, ...
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
hold on
% ***
leg = legend('show');
title(leg, 'K.D. (F.C.)')
% ***

%% Insulin
subplotTight(height_y, length_x, 4:6)
    plot(t,(x_conc(:,Insulin_pos)), 'linestyle', style_variation,  'LineWidth', 1.5, ...
        'Color',[0 0 1*color_variation]); %, 'DisplayName', Kd_Variations_legend)
    grid on 
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 1.2e-10;
    ylim([0 ymax])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Plasma Insulin', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on
% % ***
% leg = legend('show');
% title(leg, 'Percent K.D.')
% % ***
    
%% Akt
subplotTight(height_y, length_x, 7:9)
Akt_total = x_conc(:,Akt_pos) + x_conc(:,AktT308_pos) + x_conc(:,AktS473_pos) + x_conc(:,AktS473T308_pos);
lines = AktS473_pos;
    plot(t,Akt_total, '-', 'LineWidth', 1.5, 'Color',[0.5 0.5 0.5])
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    ylim([0 Akt_total(1,1)*1.1])
    ylabel('Concentration (mol/L)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    yticks([0 : 0.4e-8 : Akt_total(1,1)*1.1])
    title('Akt', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot (t,x_conc(:,lines), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
hold on
legend({'Total Akt', 'p-Akt^{S473}'});

%% p70S6K
subplotTight(height_y, length_x, 10:12)
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
    yticks([0 : 0.3e-8 : p70S6K_total(1,1)*1.1])
    title('p70S6K', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on 
    plot (t,x_conc(:,lines), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
hold on
legend({'Total p70S6K', 'p-p70S6K^{T389}'}, 'Location','northeast');

%% MPS
subplotTight(height_y, length_x, 13:16)
    plot (t, x_mass(:,MPS_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 4.2e-5;
    ylim([0 ymax])
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('MPS', 'fontSize', fontSize_title)
    xtickangle(Angle)

hold on
    area(t(t_endEq(1):t_endRun(1)), x_mass(t_endEq(1):t_endRun(1), MPS_pos), ...
        'FaceColor',[0 0 1], 'EdgeColor',"none", ...
        'FaceAlpha',0.05, ...
        'HandleVisibility','off')

hold on
if Kd_rxn == 16
    leg = legend({'0.36', '0.37', '0.38', '0.43', '0.56'},'NumColumns', 1, 'Location','best'); % Insulin KO
elseif Kd_rxn == 39
    leg = legend({'0.36', '0.32', '0.28', '0.21', '0.14'},'NumColumns', 1, 'Location','best'); % leucine KO
else
end

title(leg, 'FSR (g leu)')

%% MPB
subplotTight(height_y, length_x, 17:20)
    plot (t, x_mass(:,MPB_pos), 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
    if Kd_rxn == 16
        ylim([0 4.4e-5])
    elseif Kd_rxn == 39
        ylim([0 4.2e-5])
    else
    end
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('MPB', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on

%% NB
subplotTight(height_y, length_x, 21:24)
NB = x_mass(:,MPS_pos) - x_mass(:,MPB_pos);
    plot (t, NB, 'linestyle', style_variation,  'LineWidth', 1.5, 'Color',[0 0 1*color_variation]);
    grid on
        ax=gca;
        ax.FontSize = fontSize_units;
    xlim(x_lim)
    xticks(xTicks)
    xticklabels(xTicksSeq)
        ymax = 8e-5;
    if Kd_rxn == 16
        ylim([-3.1e-5 ymax])
    elseif Kd_rxn == 39
        ylim([-0.4e-5 ymax])
    else
    end
    ylabel('Moles/minute', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    xlabel('Time (min)', 'fontWeight', 'bold', 'fontSize', fontSize_axes)
    title('Net Balance', 'fontSize', fontSize_title)
    xtickangle(Angle)
hold on

end
