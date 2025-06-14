clc;
clear all;
close all;

% Defining constants
g = 9.81; % Gravitational acceleration
K = 0.41; % Karman constant
rho = 1000; % water density
rhos = 2650; % sediment density
R = (rhos - rho) / rho; 
v = 1e-6; % Kinematic viscosity of water
Mu = 0.001; % Kinetic viscosity of water

% Given grain size ranges (m)
D = [0.001, 0.005, 0.010, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045, 0.05];

% Water depth and average flow velocity from 'Model 1'
h_values = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];
U_values = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];

% Color scheme
colors = parula(length(D)); 
colors = flipud(colors); 


figure;
map = flipud(mymap('viridis'));
colormap(gca, map);
colormap(flipud(parula)); 
set(gcf, 'Position', [100, 100, 1400, 1000]);

for i = 1:8
    h = h_values(i);
    U = U_values(i);
    
    %  Position on the vertical profile
    y = linspace(0.05, h, 100);
    y_h = y / h;
    
    % Calculations of shear velocity
    u_shear = U / 8.3;
    
    % Calculation of turbulent energy dissipation rate
    Epsilon = ((9.8 * exp(-3 * y / h)) ./ ((y / h) .^ 0.5)) .* ((u_shear ^ 3) / h);
    
    % Calculation of y1
    y1 = y * u_shear / v;
    
    % Calculations of the root-mean square velocity
    u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
    
    % Calculations of Taylor length scale
    Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);
    
    % Creating subplot
    subplot(4,2,i);
    set(gca, 'FontSize', 16);
    hold on;
    
    % Iterate over all particle sizes D
    for j = 1:length(D)
        
        % Calculation of  particle settling velocity
        Ws = 0.51 * (v / D(j)) * ((D(j)^3) * g * R / (v^2))^0.553; % Carling et al., 2020
        
        % Calculation of particle response time
        Tp = Ws / g;
        
        % Calculations of turbulence time scales
        Tk = Lambda ./ u_rootmean;
        
        % Calculations of Stokes number
        St = Tp ./ Tk;
        
        % plot curves
        plot(St, y_h, 'Color', colors(j, :), 'LineWidth', 2);
    end
    
    % Add vertical lines with St=1 and St=9 to the subplot
    plot([1, 1], [0, 1], 'g--', 'LineWidth', 2); % St = 1 
    plot([9, 9], [0, 1],'--','Color', [0.5, 0.5, 0.5], 'LineWidth', 2); % St = 9 
    
        xlabel('Stokes number {\itS_t} (-)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k'); 
    ylabel('{\ity/h} (-)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
       
    set(gca, 'YDir', 'normal', 'TickLength', [0.02 0.02]);
    set(gca, 'FontSize', 20);
    grid on;
    xlim([0 40]);
    ylim([0 1]);
    hold off;
    box on;
end

% Add a colorbar and map D values
c = colorbar('eastoutside');
caxis([min(D) max(D)]);  
c.Label.String = 'Grain size {\itD} (mm)';  
c.Label.FontSize = 20;
c.Label.FontName = 'Arial';

% setting colorbar scale
tick_positions = linspace(min(D), max(D), length(D));
c.Ticks = tick_positions; 
c.TickLabels = arrayfun(@(x) num2str(round(x * 1000)), D, 'UniformOutput', false); % unit is mm
c.Position = [0.92, 0.121, 0.018, 0.79];
c.FontSize = 20;
c.Position = [c.Position(1), c.Position(2) + 0.035, c.Position(3), c.Position(4)-0.02];
%export_fig StProfileDistribution.jpg -r1000 -transparent;


% Demonstration of the distribution of mapped u' in the bathymetric profile
figure;
set(gcf, 'Position', [100, 100, 1400, 1000]);
for i = 1:8
    h = h_values(i);
    
    % Calculations of the distribution of u_rms
    y = linspace(0.05, h, 100);
    y_h = y / h;
    u_shear = U_values(i) / 8.3;
    y1 = y * u_shear / v;
    u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
    
    % Creating subplot
    subplot(4,2,i);
    plot(u_rootmean, y_h, 'b', 'LineWidth', 4,'Color', [112/255, 48/255, 160/255]);
    set(gca, 'FontSize', 20);
    xlabel("u_{rms} (m/s)", 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ylabel('y/h', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    %title(['h = ', num2str(h), 'm'], 'FontSize', 16, 'FontName', 'Arial', 'Color', 'k');
    set(gca, 'YDir', 'normal');
    grid on;
end
