clc;
clear all;
close all force;

% 1.Creating hydrograph
time = linspace(0, 24, 1441);  
peak_time = 4;    
peak_flow = 1000;  

flow_peak = peak_flow * exp(-0.5 * ((time - peak_time) / 1.4).^2); 
decay_indices = time > peak_time & time <= 24;
decay_time = time(decay_indices);
flow_decay = peak_flow * exp(-0.1 * ((decay_time - peak_time) / 2.0).^2);  

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
D_values = [D1, D2];  
nD = length(D_values);

% 4.Calculation of shear velocity
u_shear = (flow * g * s / 8.3 ).^(1/3);
tau = rho * u_shear.^2;

% 5.Creating color mappings
cmap = parula(nD); 
colormap(cmap);

% 6.Setting different dimensionless height a/h values
k_list = [0.05, 0.2, 0.4, 0.6, 0.8, 0.95];

% 7.plot
figure;
set(gcf, 'Position', [200, 200, 1200, 950]); 

for idx = 1:length(k_list)
    k_ratio = k_list(idx);

    subplot(3, 2, idx); hold on; box on;
    %title(['a/h = ', num2str(k_ratio)], 'FontSize', 16);

    for i = 1:nD
        D = D_values(i); 

        % Calculation of the particles settling velocity
        Ws = 0.51 * (v / D) * ((D^3 * g * R) / v^2)^0.553;

        % water depth
        h = (u_shear.^2) / g / s;

        % reference high
        a = k_ratio * h;

        % z+ parameter
        z = a .* u_shear / v;

        % Calculation of root-mean square velocity
        u_rootmean = u_shear .* 2.3 .* exp(-a./h) .* (1 - exp(-z ./ 10)) + 0.3 .* z .* exp(-z ./ 10);

        % Calculation of R_star=
        R_star = Ws ./ u_rootmean;

        % plotting results
        plot(time, R_star, 'Color', cmap(i,:), 'LineWidth', 2.5);
    end
    xlabel('Time {\itT} (hr)', 'FontSize', 18);
    ylabel('{\itR} (-)', 'FontSize', 18);
    xlim([0, 24]);
    ylim([0, 60]);
    xticks(0:4:24);
    set(gca, 'FontSize', 18);
    set(gca, 'TickLength', [0.015, 0.015]);
    grid on;
end

hcb = colorbar('Position', [0.93 0.133 0.0185 0.791]);
caxis([1, nD]);
hcb.Ticks = linspace(1, nD, nD);
hcb.TickLabels = arrayfun(@(x) sprintf('%.d', x*1000), D_values, 'UniformOutput', false);
ylabel(hcb, 'Grain size {\itD} (mm)', 'FontSize', 18);
hold off;

%export_fig R(-).jpg -r1000 -transparent;
