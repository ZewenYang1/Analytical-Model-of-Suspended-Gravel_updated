clc;
clear all;
close all;

% Defining constants
g = 9.81;     
rho = 1000;   
rhos = 2650;  
R = (rhos - rho) / rho; 
v = 1e-6;      
Mu = 0.001;
C_Nu = 0.09;

% Array of Karman constants
K_values = [0.41, 0.30];  % two K values

% Roughness height
z0 = 0.05 * 2.5; 

% Water depth and average flow velocity from 'Model 1'
h_values = [50.42, 99.55, 93.12, 76.25, 54.63, 34.26, 18.80, 9.03];
U_values = [7.15, 10.05, 9.72, 8.79, 7.44, 5.89, 4.37, 3.02];

% Creating figure window
figure;
set(gcf, 'Position', [100, 90, 800, 1000]);

% The color array is used to distinguish between different K
color_set = [255 145 37; 0 189 255] / 255;

% Cyclic drawing of subgraphs
for i = 1:8
    h = h_values(i);
    U = U_values(i);
    
    subplot(4, 2, i);
    hold on;

    % Storing drawing handles
    plot_handles = gobjects(length(K_values), 1);

    for j = 1:length(K_values)
        K = K_values(j);

        % Calculation of shear velocity
        u_shear = U / 8.3;

        % Calculation of TKE
        y = linspace(0.05, h, 100);
        y_h = y / h;
        Epsilon = ((9.8 * exp(-3 * y / h)) ./ ((y / h) .^ 0.5)) .* ((u_shear ^ 3) / h);
        y1 = y * u_shear / v;
        u_rootmean = 2.3 * u_shear .* exp(-y / h) .* (1 - exp(-y1 / 10)) + 0.3 * y1 .* exp(-y1 / 10);
        Lambda = sqrt(15 * v * (u_rootmean.^2) ./ Epsilon);
        TKE = Epsilon .* K * h ./ ((C_Nu^0.5)*u_shear);

        % Plotting and recording handles
        plot_handles(j) = plot(TKE, y_h, 'LineWidth', 2, 'Color', color_set(j,:));
    end

    % Setting figure 
    set(gca, 'FontSize', 20);  
    xlabel('{\itTKE} (m^2/s^2)', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ylabel('{\ity/h}', 'FontSize', 20, 'FontName', 'Arial', 'Color', 'k');
    ytickformat('%.1f');
    xlim([0, 100]);
    set(gca, 'XTick', 0:20:100); 
    set(gca, 'YDir', 'normal');
    set(gca, 'TickLength', [0.03, 0.06]); 
    box on;

    % Adding legend for first subgraph only
    if i == 1
        legend(plot_handles, {'$\beta$ = 1', '$\beta$ = 0.73'}, 'Location', 'northeast', 'FontSize', 16, 'Interpreter', 'latex');
        %  K =0.41 is corresponds to ¦Â = 1; K =0.3 is corresponds to ¦Â = 0.73
        % (K/¦Â always equal to 0.41, and ¦Â is a simple means to ¡®keep K = 0.41 in the various equations. K is Karman constant)
                
        legend boxoff;   
    end
end

% export_fig TKE_with_¦Â.jpg -r1000 -transparent;


