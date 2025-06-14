% This code is used to calculate the Stokes Number (Cao et al., 2003)

clc;
clear all;
close all;

% Creating hydrograph
time = linspace(0, 24, 1441); 
peak_time = 4;
peak_flow = 100;

flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);
decay_indices = time > peak_time & time <= 24*60;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2.0).^2);
flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% Defining constants
s = 0.0015;
g = 9.81;
rho = 1000;
rhos = 2650;
R = (rhos - rho) / rho;
v = 1e-6;

% Given grain size ranges (m)
D1 = 0.001:0.001:0.005;  
D2 = 0.01:0.01:0.05;     
D_all = [D1, D2];
nD = length(D_all);

% Calculation of shear velocity and water depth
u_shear = (flow * g * s / 8.3 ).^(1/3);
h = (u_shear .^2) / g / s;

% Given a/h ratio lists 
k_list = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];

% color mapping
map = flipud(mymap('plasma'));
colors_all = map(round(linspace(1, size(map, 1), nD)), :);
idx1 = linspace(1, size(map, 1) / 2, length(D1)); 
idx2 = linspace(size(map, 1) / 2 + 1, size(map, 1), length(D2)); 
colors1 = map(round(idx1), :);
colors2 = map(round(idx2), :);

% Creating figure
figure;
set(gcf, 'Position', [100, 100, 1600, 12400]);

width = 0.2;     
height = 0.20;    
gap_x = 0.06;     
gap_y = 0.07;      
x_start = 0.1;
y_start = 0.65;

for k_index = 1:length(k_list)
    k = k_list(k_index);
    a = k * h;
    z = a .* u_shear / v;

    % Calculate root-mean square velocity
    u_rootmean = u_shear .* 2.3 .* exp(-a./h) .* (1 - exp(-z ./ 10)) + 0.3 .* z .* exp(-z ./ 10);

    % Calculation of turbulent energy dissipation rate
    Epsilon = (9.8 * exp(-3 * a ./ h) .* (u_shear.^3)) ./ (((a ./ h) .^0.5) .* h);

    % Calculating Taylor microscale
    Lambda = sqrt(15 * v .* (u_rootmean.^2) ./ Epsilon);

    % Setting Subplot Position
    col = floor((k_index - 1) / 3);
    row = mod(k_index - 1, 3);
    pos_x = x_start + col * (width + gap_x);
    pos_y = y_start - row * (height + gap_y);
    axes('Position', [pos_x, pos_y, width, height]); hold on;

    % plotting
    for i = 1:nD
        D = D_all(i);
        Ws = 0.51 * (v / D) * ((D^3) * g * R / (v^2))^0.553;
        Tp = Ws / g;
        Tt = Lambda ./ u_rootmean;
        St = Tp ./ Tt;

        plot(time, St, 'LineWidth', 2.5, 'Color', colors_all(i, :));
    end

    % reference lines
    plot([0 24], [1 1], 'g--', 'LineWidth', 1.5);
    plot([0 24], [9 9], '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.5);

    % Subgraph properties
    set(gca, 'YScale', 'log');
    ax = gca;
    ax.YRuler.TickLength(1) = 0.035;  
    ax.XRuler.TickLength(1) = 0.030;  
    xlim([0 24]);
    ylim([1e-2 1e2]);
    set(gca, 'xtick', 0:4:24, 'FontSize', 16);
    xlabel('Time {\itT} (hr)', 'FontSize', 16);
    ylabel('Stokes Number {\itS_t} (-)', 'FontSize', 16);
    box on;   
end

% Setting colorbar
colormap(map);
cbar = colorbar('Position', [0.57 0.11 0.015 0.740]);
cbar.Label.String = 'Grain Size {\itD} (mm)'; 
cbar.FontSize = 16;
caxis([min([D1, D2]), max([D1, D2])]); 
cbar.Ticks = linspace(min([D1, D2]), max([D1, D2]), numel([D1, D2]));
cbar.TickLabels = arrayfun(@(x) num2str(x * 1000), [D1, D2], 'UniformOutput', false); 

%export_fig StAndTimeWITHyh.jpg -r1000 -transparent;