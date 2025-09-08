clear all; close all;
run_NPZ_model_fixed_dt

function run_NPZ_model_fixed_dt
    % Time settings (unit: days)
    tspan = [0 1000];  
    nSteps = 50000;
    t = linspace(tspan(1), tspan(2), nSteps);
    dt = t(2) - t(1);

    % Noise amplitude
%     eps = 0.5;
    eps = 0.5;
    
    noise_flag = 1; % 0 = constant, 1 = stochastic, 2 = Gaussian

    % Noises
    switch noise_flag
        case 0
            raw_noise = zeros(1,length(t));
        case 1
            raw_noise = 3.4641 * (rand(1, length(t)) - 0.5) * eps;
        case 2
            raw_noise = (randn(1, length(t))) * eps;
    end
    noise_series = raw_noise - mean(raw_noise);
    noise_function = @(ti) interp1(t, noise_series, ti, 'nearest', 'extrap');

    % initialization [N, P, Z]
    y = zeros(nSteps, 3);
    y(1,:) = [0.002, 0.04, 0.002];

    % recording
    upwelling_N = zeros(size(t));
    upwelling_P = zeros(size(t));
    uptake = zeros(size(t));

    % integration
    for i = 1:nSteps-1
        [dydt, aux] = npz_system(t(i), y(i,:), noise_function);
        
        y(i+1,:) = y(i,:) + dt * dydt';
        
        % not negative
        y(i+1,:) = max(y(i+1,:), 0);

        % recording
        upwelling_N(i) = aux.upwelling_N;
        upwelling_P(i) = aux.upwelling_P;
        uptake(i) = aux.uptake;
    end

    % recording
    [~, aux] = npz_system(t(end), y(end,:), noise_function);
    upwelling_N(end) = aux.upwelling_N;
    upwelling_P(end) = aux.upwelling_P;
    uptake(end) = aux.uptake;

    % display
    disp(['Mean upwelling_N = ', num2str(mean(upwelling_N))]);
    disp(['Mean upwelling_P = ', num2str(mean(upwelling_P))]);
    disp(['Mean uptake = ', num2str(mean(uptake))]);
    disp(['Mean N = ', num2str(mean(y(:,1)))]);
    disp(['Mean P = ', num2str(mean(y(:,2)))]);

    % plot
    figure;
    hold on;

    yyaxis right
    plot(t, y(:,1), 'b', 'DisplayName', 'Nutrient (N)', 'LineWidth', 2);
    ylabel('Nutrient (N) (gC m^{-3})');
    ylim([-0.5*1e-2, 1.5*1e-2]);
    set(gca, 'ycolor', 'b');


    yyaxis left
    
    plot(t, y(:,3), 'r-', 'DisplayName', 'Zooplankton (Z)', 'LineWidth', 2);
    plot(t, y(:,2), 'g', 'DisplayName', 'Phytoplankton (P)', 'LineWidth', 2);

    P_mean = mean(y(:,2));
    plot(t, P_mean * ones(size(t)), '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5, 'DisplayName', 'Mean P');

    ylabel('P, Z');
    ylim([-0.05, 0.15]);
    set(gca, 'ycolor', 'k');

    

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

    figure;
    hold on;
    plot(t, upwelling_N, 'b--', 'DisplayName', 'Upwelling N', 'LineWidth', 2);
    plot(t, upwelling_P, 'g--', 'DisplayName', 'Upwelling P', 'LineWidth', 2);
    plot(t, uptake, 'k--', 'DisplayName', 'Uptake', 'LineWidth', 2);
    ylabel('Flux terms');
    ylim([-0.04, 0.04]);
    set(gca, 'ycolor', 'k');
    xlabel('Time (days)');
    legend('location', 'northwest');
    set(gca,'fontsize', 15);
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

    % Biological processes
    uptake = (N / (e + N)) * (1 / (kw + kc * P)) * P;   
    grazing = (lambda * P^2 / (nu^2 + P^2)) * Z;        

    upwelling_N = (w2 / Hm2) * (N0 - N);
    upwelling_P = (w2 / Hm2) * (P0 - P);

    % Differential equations
    dNdt = -uptake + r * P + eta * grazing + gamma * d * Z^2 + upwelling_N;
    dPdt = uptake - r * P - grazing - s * P + upwelling_P;
    dZdt = alpha * grazing * Z - d * Z^2;
    
%     if N + dNdt < 0
%         dNdt = -N;  % N should not be negative
%     end
%     if P + dPdt < 0
%         dPdt = -P;  % N should not be negative
%     end
%     if Z + dZdt < 0
%         dZdt = -Z;  % N should not be negative
%     end

    dydt = [dNdt; dPdt; dZdt];
    aux.upwelling_N = upwelling_N;
    aux.upwelling_P = upwelling_P;
    aux.uptake = uptake;
end
