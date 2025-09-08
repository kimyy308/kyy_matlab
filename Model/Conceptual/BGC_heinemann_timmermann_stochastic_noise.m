
% upwelling_N_record = [];
% upwelling_P_record = [];

run_NPZ_model


function run_NPZ_model
    % Time settings (unit: days)
    tspan = [0 500];  

    % Initial conditions [N, P, Z] in gC m^-3
    y0 = [0.002, 0.04, 0.002];  

    % Run ODE solver
    [t, y, aux] = ode23(@npz_system, tspan, y0);

    % Plot results
    figure;
    plot(t, y(:,1), 'b', 'DisplayName', 'Nutrient (N)', 'linewidth', 2);
    hold on;
    plot(t, y(:,2), 'g', 'DisplayName', 'Phytoplankton (P)', 'linewidth', 2);
    plot(t, y(:,3), 'r', 'DisplayName', 'Zooplankton (Z)', 'linewidth', 2);
    plot(t, aux.upwelling_N, 'b--', 'DisplayName', 'upwelling-N', 'linewidth', 2);
    plot(t, aux.upwelling_P, 'g--', 'DisplayName', 'upwelling-P', 'linewidth', 2);
    xlabel('Time (days)');
    ylabel('Concentration (gC m^{-3})');
    legend;
    title('N-P-Z Model Dynamics');
end

function [dydt, aux]= npz_system(t, y)
    % Unpack state variables
    N = y(1);  % Nutrient concentration (gC m^-3)
    P = y(2);  % Phytoplankton concentration (gC m^-3)
    Z = y(3);  % Zooplankton concentration (gC m^-3)

    % Ecosystem parameters from literature
    kw = 0.046;           % Light attenuation by water (m^-1)
    kc = 0.74;            % Biomass-specific attenuation (m^2 gC^-1)
    d = 1.0;              % Higher predation and natural mortality of Z (m^3 gC^-1 day^-1)
    e = 0.11;             % Half-saturation constant for N uptake (gC m^-3)
    r = 0.13;             % P respiration rate (day^-1)
    s = 0.04;             % P sinking loss rate (day^-1)
    N0 = 0.6;             % Deep nutrient concentration (gC m^-3)
    alpha = 0.25;         % Z growth efficiency
    eta = 0.33;           % Z excretion fraction
    gamma = 0.5;          % Regeneration of Z predation
    lambda = 0.6;         % Maximum Z grazing rate (day^-1)
    nu = 0.05;            % Z grazing half-saturation constant (gC m^-3)

    % Additional physical parameters (assumed values, modify as needed)
    eps = 1; % noise ampitude
%         w2 = 0.03;            % Upwelling velocity (m day^-1)
    w2 = 0.03 +rand(1).*eps;            % with stochastic noise
%     w2 = 0.03 + randn(1,eps);           % with Gaussian noise
    Hm2 = 100;             % Mixed layer depth (m)
    P0 = 0;             % Deep phytoplankton concentration (gC m^-3)

    % Biological processes
    uptake = (N / (e + N)) * (1 / (kw + kc * P)) * P;   % Nutrient uptake by P
    grazing = (lambda * P^2 / (nu^2 + P^2)) * Z;        % Grazing of P by Z
    
    upwelling_N = (w2 / Hm2) * (N0 - N);
    upwelling_P = (w2 / Hm2) * (P0 - P);
    
%     global upwelling_N_record upwelling_P_record
%     upwelling_N_record(end+1) = upwelling_N;
%     upwelling_P_record(end+1) = upwelling_P;
    aux.upwelling_N = upwelling_N;
    aux.upwelling_P = upwelling_P;

    % Differential equations
    dNdt = -uptake + r * P + eta * grazing + gamma * d * Z^2 + upwelling_N;
    dPdt = uptake - r * P - grazing - s * P + upwelling_P;
    dZdt = alpha * grazing * Z - d * Z^2;

    % Return derivatives
    dydt = [dNdt; dPdt; dZdt];
end