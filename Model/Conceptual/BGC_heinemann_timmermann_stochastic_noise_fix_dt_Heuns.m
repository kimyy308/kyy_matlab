clear all; close all;
run_NPZ_model_fixed_dt

function run_NPZ_model_fixed_dt
    % Time settings (unit: days)
    tspan = [0 1000];  
    nSteps = 50000; % total number of steps
    t = linspace(tspan(1), tspan(2), nSteps);
    dt = t(2) - t(1); % time step size

    % Noise amplitude
%     eps = 0.5;
    eps = 5;
    noise_flag = 1; % 0 = constant, 1 = stochastic, 2 = Gaussian
    
    % Define daily time grid for noise changes
    t_noise = floor(tspan(1)):1:ceil(tspan(2));
    
    % Generate noise series on daily time scale
    switch noise_flag
        case 0
            raw_noise = zeros(1, length(t_noise));
        case 1
            raw_noise = 3.4641 * (rand(1, length(t_noise)) - 0.5) * eps;
        case 2
            raw_noise = randn(1, length(t_noise)) * eps;
    end
    noise_series = raw_noise - mean(raw_noise); % zero-mean
    
    % Interpolation function, stepwise constant per day
    noise_function = @(ti) interp1(t_noise, noise_series, ti, 'previous', 'extrap');

    % Initialize state variables: [N, P, Z]
    y = zeros(nSteps, 3);
    y(1,:) = [0.13, 0.04, 0.06];

    % Preallocate arrays to record results
    upwelling_N = zeros(size(t));
    upwelling_P = zeros(size(t));
    uptake = zeros(size(t));

    % Time integration using Heun's method (RK2)
    for i = 1:nSteps-1
        [dydt1, aux1] = npz_system(t(i), y(i,:), noise_function);
        
        % Predictor step
        y_predict = y(i,:) + dt * dydt1';

        % Corrector step
        [dydt2, aux2] = npz_system(t(i) + dt, y_predict, noise_function);

        % Heun's update (average of predictor and corrector)
        y(i+1,:) = y(i,:) + (dt/2) * (dydt1' + dydt2');

        % Prevent negative concentrations
        y(i+1,:) = max(y(i+1,:), 0);

        % Record fluxes (based on current step)
        upwelling_N(i) = aux1.upwelling_N;
        upwelling_P(i) = aux1.upwelling_P;
        uptake(i) = aux1.uptake;
    end

    % Record final flux values
    [~, aux] = npz_system(t(end), y(end,:), noise_function);
    upwelling_N(end) = aux.upwelling_N;
    upwelling_P(end) = aux.upwelling_P;
    uptake(end) = aux.uptake;

    % Display mean statistics
    disp(['Mean upwelling_N = ', num2str(mean(upwelling_N))]);
    disp(['Mean upwelling_P = ', num2str(mean(upwelling_P))]);
    disp(['Mean uptake = ', num2str(mean(uptake))]);
    disp(['Mean N = ', num2str(mean(y(:,1)))]);
    disp(['Mean P = ', num2str(mean(y(:,2)))]);

    % Plot state variables
    figure;
    hold on;

%     yyaxis right
%     ylim([-0.5*1e-2, 1.5*1e-2]);
    plot(t, y(:,1), 'b', 'DisplayName', 'Nutrient (N)', 'LineWidth', 2);
    ylabel('Nutrient (N) (gC m^{-3})');
    set(gca, 'ycolor', 'b');

%     yyaxis left
    plot(t, y(:,3), 'r-', 'DisplayName', 'Zooplankton (Z)', 'LineWidth', 2);
    plot(t, y(:,2), 'g', 'DisplayName', 'Phytoplankton (P)', 'LineWidth', 2);

    % Plot mean phytoplankton line
    P_mean = mean(y(:,2));
    plot(t, P_mean * ones(size(t)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Mean P');

%     ylabel('Phytoplankton (P), Zooplankton (Z)');
    ylabel('gC/m^3');
    ylim([-0.05, 0.3]);
    set(gca, 'ycolor', 'k');

    xlabel('Time (days)');
    legend('location', 'northwest');
    set(gca,'fontsize', 15);

    % Add appropriate title
    switch noise_flag
        case 0
            title('N-P-Z Model with Constant Upwelling');
        case 1
            title('N-P-Z Model with Stochastic Upwelling');
        case 2
            title('N-P-Z Model with Gaussian Upwelling');
    end

    % Plot flux terms
    figure;
    hold on;
    plot(t, upwelling_N, 'b--', 'DisplayName', 'Upwelling N', 'LineWidth', 2);
    plot(t, upwelling_P, 'g--', 'DisplayName', 'Upwelling P', 'LineWidth', 2);
    plot(t, uptake, 'k--', 'DisplayName', 'Uptake', 'LineWidth', 2);
    ylabel('Flux terms (gC m^{-3} day^{-1})');
    ylim([-0.07, 0.07]);
    set(gca, 'ycolor', 'k');
    xlabel('Time (days)');
    legend('location', 'northwest');
    set(gca,'fontsize', 15);
end


% N-P-Z system with external noise input
function [dydt, aux] = npz_system(t, y, noise_function)
    % Unpack state variables
    N = y(1);  
    P = y(2);  
    Z = y(3);  

    % Ecosystem parameters
%     a = 0.06; % typical in equatorial Pacific alpha^B in platt et al. (2003), initial slope of the photosynthesis-irradiance curve
    a = 0.2; % Edwards 1997
    kw = 0.046;          
    kc = 0.74;           
%     d = 1.0; % original   
    d = 0.3; 
%     d = 0.25; % minimum range in Edwards 1997
    e = 0.11;             
    r = 0.13;             
    s = 0.04;             
    N0 = 0.6;             
%     alpha = 0.25; % original   
    alpha = 0.3;
%     alpha = 0.5; % maximum range in Edwards 1997
    eta = 0.33;           
    gamma = 0.5;          
%     lambda = 0.6; % original
    lambda = 1.4; % maximum range in Edwards 1997  
%     nu = 0.05;  % original
    nu = 0.02; % minimum range in Edwards 1997
    % Physical parameters
    Hm2 = 30;             
    P0 = 0;               

    % Upwelling with external noise
    w2 = 0.05 + noise_function(t);
%     w2 = 0.5 + noise_function(t);

    % Biological processes
    uptake = (N / (e + N)) * (a / (kw + kc * P)) * P;   
    grazing = (lambda * P^2 / (nu^2 + P^2)) * Z;        

    % Upwelling contributions
    upwelling_N = (w2 / Hm2) * (N0 - N);
    upwelling_P = (w2 / Hm2) * (P0 - P);

    % Differential equations
    dNdt = -uptake + r * P + eta * grazing + gamma * d * Z^2 + upwelling_N;
    dPdt = uptake - r * P - grazing - s * P + upwelling_P;
    dZdt = alpha * grazing * Z - d * Z^2;
    
    if N + dNdt < 0
        dNdt = -N;  % N should not be negative
    end
    if P + dPdt < 0
        dPdt = -P;  % N should not be negative
    end
    if Z + dZdt < 0
        dZdt = -Z;  % N should not be negative
    end

    % Return derivatives and diagnostic terms
    dydt = [dNdt; dPdt; dZdt];
    aux.upwelling_N = upwelling_N;
    aux.upwelling_P = upwelling_P;
    aux.uptake = uptake;
end
