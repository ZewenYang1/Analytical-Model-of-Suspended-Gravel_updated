clc;
clear all;
close all;

% 1. Creating hydrograph
time = linspace(0, 24, 1441);  
peak_time = 4;     
peak_flow = 1000;  

flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);  
decay_indices = time > peak_time & time <= 24;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2).^2);  

flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% 2.Defining constants
s = 0.0015;
rho = 1000;     
rhos = 2650;    
R = (rhos - rho) / rho;
g = 9.81;       
v = 1e-6;       
karmanconstant = 0.41;

% 3.Given grain size ranges (m)
D1 = 0.001:0.001:0.005;
D2 = 0.01:0.01:0.05;
D_all = [D1, D2]; 
nD = length(D_all);

% 4.Calculation of shear velocity
u_shear = (flow * g * s / 8.3 ).^(1/3);  
tau = rho * u_shear.^2;

% 5.Creating color mappings
%colors = parula(length(D)); 
map = flipud(parula(length(D_all))); 
color_all = map(round(linspace(1, size(map, 1), nD)), :);
% Generate color index
idx1 = linspace(1, size(map, 1) / 2, length(D1)); % D1 Using colormap First Half
idx2 = linspace(size(map, 1) / 2 + 1, size(map, 1), length(D2)); % D2 Using the second half of the colormap
colors1 = map(round(idx1), :);
colors2 = map(round(idx2), :);


% 6.plotting results
% The list of a/h 
k_list = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];

% plotting
figure;
set(gcf, 'Position', [200, 150, 1200, 950]); 

% 7.Iteration a/h
for idx = 1:length(k_list)
    k_ratio = k_list(idx);  

    % Create corresponding subgraphs
    subplot(3, 2, idx); 
    hold on; box on;
    %title(['a/h = ', num2str(k_ratio)], 'FontSize', 16);

    % Iterate over all particle size values
    for i = 1:nD
        D = D_all(i); 

        % Calculation of setting velocity
        Ws = 0.51 * (v / D) * ((D^3 * g * R) / v^2)^0.553;

        % Calcualtion of water depth  h£¨m£©
        h = (u_shear.^2) / g / s;

        % Calcualtion of reference height
        a = k_ratio * h;

        % z+ parameter£¨dimensionless£© 
        z = a .* u_shear / v;

        % Calculation of u_rootmean£¨m/s£©
        u_rootmean = u_shear .* 2.3 .* exp(-a./h) .* (1 - exp(-z ./ 10)) + 0.3 .* z .* exp(-z ./ 10);

        % Calculation of turbulent energy dissipation rate ¦Å 
        Epsilon = (9.8 * exp(-3 * a ./ h) .* (u_shear.^3)) ./ (((a ./ h) .^0.5) .* h);

        % Calculation of Taylor mixing scale ¦Ë£¨m£©
        Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);

        % Calculation of the dimensionless particle size d for reference concentrations
        d = D .* ((R * g).^(1/3)) ./ (v^(2/3));

        % Calculation of reference concentration Cb(-)
        Cb = calculateCb(d);  % Calling external functions

        % Calculation of R_star
        R_star = Ws ./ u_rootmean;

        % Calculation of G_star
        G_star = g * R * Lambda .* Cb ./ (u_rootmean.^2);

        % Calculation of G_star / R_star^2
        ratio1 = G_star ./ (R_star .^2);

        % Plotting, each curve represents a particle size
        plot(time, ratio1, 'Color', color_all(i,:), 'LineWidth',2.5);
    end

    % Setting the Axis Format
    xlabel('Time {\itT} (hr)', 'FontSize', 18);
    ylabel('{\itG/R^2} (-)', 'FontSize', 18);
    xlim([0, 24]);
    ylim([0, 6]);  
    xticks(0:4:24);
    set(gca, 'FontSize', 18);
    set(gca, 'TickLength', [0.015, 0.015]);
    grid on;

    % Addition of three reference lines (corresponding to the formation of the judgment of the front)
    plot([0, 24], [2.2, 2.2], 'g--', 'LineWidth', 1.5);
    plot([0, 24], [0.5, 0.5], 'Color', [0.5 0.5 0.5],'LineStyle', '--', 'LineWidth', 1.5);                  
    plot([0, 24], [0.032, 0.032], 'y--', 'LineWidth', 1.5);    
end


% Adding colorbar
colormap(flipud(parula)); 
c = colorbar('eastoutside');
caxis([min(D_all) max(D_all)]);

c.Label.String = 'Grain size {\itD} (mm)';  
c.Label.FontSize = 18;
c.Label.FontName = 'Arial';

% Setting colorbar
tick_positions = linspace(min(D_all), max(D_all), length(D_all));
c.Ticks = tick_positions;
c.TickLabels = arrayfun(@(x) num2str(round(x * 1000)), D_all, 'UniformOutput', false); 
c.Position = [0.93 0.122 0.0187 0.803];

%export_fig GandR^2.jpg -r1000 -transparent;

function Cb = calculateCb(d)
    A = 1.3e-7;
    Rf = 0.0738 * (d.^1.021);
    Zu = 9.95 * 0.8 * (Rf.^0.882);
    Cb = (A * Zu.^5) ./ (1 + (A/0.3) * Zu.^5);
end
