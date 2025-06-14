clc;
clear all;
close all;

% Creating time series
time = linspace(0, 24, 1441); % 0 to 24 hours at 1 minute intervals

% Creating flood discharge per unit width
peak_time = 4; % Peak time is hour 4
peak_flow = 1000; % Peak discharge

% Hydrograph with with Gaussian distribution
flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2);
decay_indices = time > peak_time & time <= 24*60;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2).^2);
flow = zeros(size(time));
flow(1:numel(flow_peak)) = flow_peak;
flow(decay_indices) = flow_decay;

% Parameters
s = 0.0015; % Referring to the river gradient in the middle reaches of the Yarlung Tsangpo River (Gacha-Langxian), the river gradient is set to be 0.0015
g = 9.81; % Gravitational acceleration
u_shear = (flow * g * s / 8.3 ).^(1/3); % Shear velocity
v = 0.000001; % Kinematic viscosity of water

% Calcultions of bed shear stress
rho = 1000; % water density
rhos = 2650; % Sediment density
R = (rhos - rho) / rho; % Relative density
tau = rho * (u_shear.^2); % tau is shear stress

% Calcultions of average flow velocity
% Given (8/fc)^0.5=8.3
U = 8.3 * u_shear;

% Creating figure window
figure;
set(gcf, 'Position', [100, 100, 1200, 900]);

% Subfigure 1: Plotting hydrograph
subplot(3, 2, 1); 
plot(time, flow,'Color', [62/255 38/255 169/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Discharge per unit width {\itq} (m^2/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 1200]);
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])


% Adding auxiliary lines
line([4 4], [0 1000], 'linestyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
line([0 4], [1000 1000], 'linestyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 2);

% Subfigure 2: Plotting average flow velocity
subplot(3, 2, 2);
plot(time, U,'Color', [244/255 120/255 58/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Mean velocity {\itu} (m/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 15]); 
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])


% Subfigure 3: Plotting shear velocity
subplot(3, 2, 3); 
plot(time, u_shear,'Color', [127/255 127/255 127/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Shear velocity {\itu_*} (m/s)', 'FontSize', 14);
xlim([0 24]);
ylim([0 1.5]); 
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])
ytickformat('%.1f');


% Subfigure 4: Plotting shear stress curve
subplot(3, 2, 4); 
plot(time, tau,'Color', [12/255 172/255 232/255], 'LineWidth', 3.5); 
xlabel('Time {\itT} (hr)', 'FontSize', 14);
ylabel('Boundary  shear stress {\it\tau} (N/m^2)', 'FontSize', 14);
xlim([0 24]);
ylim([0 2000]); 
set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
set(gca, 'TickLength', [0.02 0.025])


% % Subfigure 5: Plotting the flow Reynolds number
% % calculations of water depth
% h = (u_shear .^2) / g / s;
% % Calculations of Renolds Number
% Re = U .* h ./ v;
% subplot(3, 2, 5); 
% plot(time, Re,'Color', [69/255 141/255 252/255], 'LineWidth', 3.5); 
% xlabel('Time {\itT} (hr)', 'FontSize', 14);
% ylabel('Flow Renolds number {\itRe} (-)', 'FontSize', 14);
% xlim([0 24]);
% ylim([1e2 1e10]); 
% set(gca, 'YScale', 'log'); 
% set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
% set(gca, 'TickLength', [0.025 0.025]);
% 
% 
% % Subfigure 6: Plotting the flow Froude number
% % calculations of Froude Number
% % W = 2000; % channel width
% % h = (u_shear .^2) / g / s;
% % R = h .* W ./ (2 * h + W); % hydraulic radius
% % Fr = U ./ sqrt(g .* R); 
% 
% h = (u_shear .^2) / g / s;
% Fr = U ./ sqrt(g .* h); 
% subplot(3, 2, 6); 
% plot(time, Fr,'Color', [196/255 147/255 255/255], 'LineWidth', 3.5); 
% xlabel('Time {\itT} (hr)', 'FontSize', 14);
% ylabel('Froude Number {\itFr} (-)', 'FontSize', 14);
% xlim([0 24]);
% ylim([0 1]); 
% set(gca, 'xtick', 0:4:24, 'FontSize', 14); 
% set(gca, 'TickLength', [0.02 0.025])
% 
export_fig HydraulicVariable.jpg -r1000 -transparent;




% calculations of water depth
h = (u_shear .^2) / g / s;
figure;

% Plotting water depth
plot(time, h, 'LineWidth', 4, 'Color', [69/255 141/255 252/255]); 
xlabel('Time {\itT} (hr)', 'FontSize', 28);
ylabel('Water depth {\ith} (m)', 'FontSize', 28);
xlim([0 24]);
ylim([0 120]); 
set(gca, 'xtick', 0:2:24, 'FontSize', 20); 
set(gca, 'TickLength', [0.01 0.01]);
grid on; 
hold on;

% Calculate the specified point in time
specific_times = [120, 240, 360, 480, 600, 720, 840, 960]; % minutes
specific_times_hours = specific_times / 60; % Conversion of minutes to hours

% Calculate the value of water depth and average flow velocity at the corresponding time
h_specific_waterdepth = interp1(time, h, specific_times_hours, 'linear');
U_specific_averagevelocity = interp1(time, U, specific_times_hours, 'linear');


% Printout of water depth and average flow velocity
disp('特定时间点的水深值 (m):');
disp('时间 (小时)    水深 h (m)    平均流速 U (m/s)');        
disp([specific_times_hours', h_specific_waterdepth', U_specific_averagevelocity']);
          
% Customizing the color of points
colors = [
    063/255, 007/255, 077/255; 
    068/255, 057/255, 131/255;
    047/255, 105/255, 140/255; 
    032/255, 145/255, 141/255;
    054/255, 184/255, 120/255; 
    108/255, 205/255, 90/255; 
    143/255, 215/255, 067/255; 
    242/255, 226/255, 048/255; 
];

% Plotting vertical auxiliary line and different colored intersections on a bathymetric chart
for i = 1:length(specific_times_hours)
    % Plotting vertical auxiliary line
    line([specific_times_hours(i), specific_times_hours(i)], [0, h_specific_waterdepth(i)], 'LineStyle', '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    
    % Plotting different colored intersections
    plot(specific_times_hours(i), h_specific_waterdepth(i), 'o', 'MarkerSize', 18, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', 'k'); 
end
hold off;

% Setting the figure window size
set(gcf, 'Position', [100, 100, 1500, 400]);
export_fig WaterDepth.tif -r1000 -transparent;


