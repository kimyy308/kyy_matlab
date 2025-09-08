clear all; close all;
run_NPZ_model

function run_NPZ_model
    % Time settings (unit: days)
    tspan = [0 1000];  
    dt_noise = 1;  % Noise 변화 간격 (1일마다 바뀜)
    t_noise = tspan(1):dt_noise:tspan(2);

    % Noise amplitude
    eps = 0.5;
    
    noise_flag = 1; % 0 = constant, 1 = stochastic, 2 = Gaussian
    
    switch noise_flag
        case 0
        % % %     Generate constant noise (zero) series (Uniform with mean 0, variance 0.25)
            raw_noise = zeros(1,length(t_noise));
            noise_series = raw_noise - mean(raw_noise); % set noise mean to 0
        case 1
        % %     % Generate stochastic noise series (Uniform with mean 0, variance 0.25)
            raw_noise = 3.4641 * (rand(1, length(t_noise)) - 0.5) * eps;
            noise_series = raw_noise - mean(raw_noise); % set noise mean to 0
        case 2
        % % %     Generate Gaussian noise series (Uniform with mean 0, variance 0.25)
            raw_noise = (randn(1, length(t_noise))) * eps;
            noise_series = raw_noise - mean(raw_noise); % set noise mean to 0
    end

    % Interpolation function for noise based on time
    noise_function = @(t) interp1(t_noise, noise_series, t, 'nearest', 'extrap');

    % Initial conditions [N, P, Z] in gC m^-3
    y0 = [0.002, 0.04, 0.002];  

    % Run ODE solver with noise passed via anonymous function
    [t, y] = ode45(@(t, y) npz_system(t, y, noise_function), tspan, y0);

    % Post-process to get upwelling terms
    upwelling_N = zeros(size(t));
    upwelling_P = zeros(size(t));
    uptake = zeros(size(t));
    noise = zeros(size(t));

    for i = 1:length(t)
        [~, aux] = npz_system(t(i), y(i,:), noise_function);
        upwelling_N(i) = aux.upwelling_N;
        upwelling_P(i) = aux.upwelling_P;
        uptake(i) = aux.uptake;
        noise(i) = aux.noise;
    end
    kurtosis(noise)
    % Display means
    disp(['Mean upwelling_N = ', num2str(mean(upwelling_N))]);
    disp(['Mean upwelling_P = ', num2str(mean(upwelling_P))]);
    disp(['Mean uptake = ', num2str(mean(uptake))]);
    disp(['Mean N = ', num2str(mean(y(:,1)))]);
    disp(['Mean P = ', num2str(mean(y(:,2)))]);
    

    % Plot results
    figure;
    % Left y-axis
    yyaxis left
    hold on;
    plot(t, y(:,2), 'g', 'DisplayName', 'Phytoplankton (P)', 'LineWidth', 2); 
    plot(t, y(:,3), 'r-', 'DisplayName', 'Zooplankton (Z)', 'LineWidth', 2);

    P_mean = mean(y(:,2));
    plot(t, P_mean * ones(size(t)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Mean P');

    ylabel('N, P, Z');
    ylim([-0.05, 0.15]);
    set(gca, 'ycolor', 'k'); % 왼쪽 축을 검정색으로

    % Right y-axis
    yyaxis right
    plot(t, y(:,1), 'b', 'DisplayName', 'Nutrient (N)', 'LineWidth', 2);
    ylabel('Nutrient (N) (gC m^{-3})');
    ylim([-0.5*1e-2, 1.5*1e-2]);
    set(gca, 'ycolor', 'b'); % 왼쪽 축을 검정색으로
    
    xlabel('Time (days)');
    legend('location', 'northwest');
    set(gca,'fontsize', 15);
    switch noise_flag
        case 0
            title('N-P-Z Model with Constant Upwelling');
        case 1
            title('N-P-Z Model with Stochastic Upwelling');
        case 2
            title('N-P-Z Model with Gaussian Upwelling');
    end

    % Plot flux results
    figure;
    hold on;
    plot(t, upwelling_N, 'b--', 'DisplayName', 'Upwelling N', 'LineWidth', 2);
    plot(t, upwelling_P, 'g--', 'DisplayName', 'Upwelling P', 'LineWidth', 2);
    plot(t, uptake, 'k--', 'DisplayName', 'Uptake', 'LineWidth', 2);
    ylabel('Flux terms');
    ylim([-0.04, 0.04]);
    set(gca, 'ycolor', 'k'); % 왼쪽 축을 검정색으로    
    xlabel('Time (days)');
    legend('location', 'northwest');
    set(gca,'fontsize', 15);
    switch noise_flag
        case 0
            title('N-P-Z Model Fluxes with Constant Upwelling');
        case 1
            title('N-P-Z Model Fluxes with Stochastic Upwelling');
        case 2
            title('N-P-Z Model Fluxes with Gaussian Upwelling');
    end
end

% Main NPZ system with noise as external input
function [dydt, aux] = npz_system(t, y, noise_function)
    % Unpack state variables
    N = y(1);  
    P = y(2);  
    Z = y(3);  

    % Ecosystem parameters
    kw = 0.046;          
    kc = 0.74;           
    d = 1.0;              
    e = 0.11;             
    r = 0.13;             
    s = 0.04;             
    N0 = 0.6;             
    alpha = 0.25;         
    eta = 0.33;           
    gamma = 0.5;          
    lambda = 0.6;         
    nu = 0.05;            
    
    Hm2 = 30;             
    P0 = 0;               

    % Upwelling with pre-generated stochastic noise
    w2 = 0.03 + noise_function(t);
    noise=noise_function(t);
    % Biological processes
    uptake = (N / (e + N)) * (1 / (kw + kc * P)) * P;   
    grazing = (lambda * P^2 / (nu^2 + P^2)) * Z;        

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

    dydt = [dNdt; dPdt; dZdt];
    aux.upwelling_N = upwelling_N;
    aux.upwelling_P = upwelling_P;
    aux.uptake = uptake;
    aux.noise = noise;
end
