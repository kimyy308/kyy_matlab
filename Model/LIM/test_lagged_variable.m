clc; clear; close all;

%% ===============================
%  1. Generate the true nonlinear system
% ===============================
n = 2000;
tau1 = 2; tau2 = 4;

Fe = filter([1 0.3],[1 -0.8],randn(n,1));    % forcing variable
Xtrue = 0.6*filter([1 0.5],[1 -0.9],Fe) + 0.3*randn(n,1);  % latent process linked to Fe
NPP = zeros(n,1);

for t = max(tau1,tau2)+1:n
    NPP(t) = 0.6*NPP(t-1) + 0.2*Fe(t) + 0.8*Xtrue(t) + 0.3*randn;
end

Fe = detrend(Fe);
Xtrue = detrend(Xtrue);
NPP = detrend(NPP);

%% ===============================
%  2. Sparse observation of X
% ===============================
mask = rand(n,1) < 0.2; % only 20% of X observed
Xobs = Xtrue;
Xobs(~mask) = NaN;

figure;
plot(1:n,Xtrue,'k','LineWidth',1.2); hold on;
plot(find(mask),Xobs(mask),'ro');
legend('True X','Observed points');
title('Sparse observation of latent variable X');
xlabel('Time step');

%% ===============================
%  3. Estimate X_hat from Fe (regression-based)
% ===============================
valid = ~isnan(Xobs);
mdl = fitlm(Fe(valid), Xobs(valid));
Xhat = predict(mdl, Fe);

fprintf('R^2 between true X and estimated Xhat = %.3f\n', corr(Xtrue,Xhat)^2);

%% ===============================
%  4. LIM setups
% ===============================
Z1 = [Fe, NPP];          % Case 1: Fe + NPP only
Z2 = [Fe, NPP, Xhat];    % Case 2: Fe + NPP + estimated Xhat
Z3 = [Fe, NPP, Xtrue];   % Case 3: Fe + NPP + true X (upper bound)

systems = {Z1, Z2, Z3};
labels = {'Fe + NPP', 'Fe + NPP + Xhat', 'Fe + NPP + true X'};

%% ===============================
%  5. LIM skill comparison
% ===============================
ntrain = 1500;
steps = 10;  % prediction lead
ntest = n - ntrain - steps;
obs = NPP(ntrain+1:ntrain+ntest);

ACC = zeros(3,1);
RMSE = zeros(3,1);
pred_all = zeros(ntest,3);

for m = 1:3
    Z = systems{m};
    Z_train = Z(1:ntrain,:);
    Z_test  = Z(ntrain:end,:);
    
    % Regularization for stability
    lambda = 1e-4;
    C0 = (Z_train(1:end-1,:)' * Z_train(1:end-1,:)) / (ntrain-1);
    C1 = (Z_train(2:end,:)'   * Z_train(1:end-1,:)) / (ntrain-1);
    L = C1 / (C0 + lambda*eye(size(C0)));
    
    % Forecast test
    pred = zeros(ntest,1);
    for t = 1:ntest
        Ztemp = Z_test(t,:)';
        for k = 1:steps
            Ztemp = L * Ztemp;
        end
        pred(t) = Ztemp(2); % predict NPP
    end
    pred_all(:,m) = pred;
    
    % Skill metrics
    ACC(m) = corr(pred, obs);
    RMSE(m) = sqrt(mean((pred - obs).^2));
    
    fprintf('\nCase %d: %s\n', m, labels{m});
    fprintf('ACC  = %.3f\n', ACC(m));
    fprintf('RMSE = %.3f\n', RMSE(m));
end

%% ===============================
%  6. Visualization
% ===============================
figure('Color','w','Position',[200 100 900 600]);
subplot(2,1,1);
plot(obs,'k','LineWidth',1.2); hold on;
plot(pred_all(:,1),'r--','LineWidth',1.2);
plot(pred_all(:,2),'b--','LineWidth',1.2);
plot(pred_all(:,3),'g--','LineWidth',1.2);
legend(['Truth', labels],'Location','best');
title('Forecasted NPP anomalies');
xlabel('Time step'); ylabel('NPP');

subplot(2,1,2);
yyaxis left
bar(ACC,'FaceColor',[0.2 0.5 0.9]);
ylabel('ACC'); ylim([0 1]);
yyaxis right
plot(RMSE,'ro-','LineWidth',1.3,'MarkerFaceColor','r');
ylabel('RMSE');
set(gca,'XTick',1:3,'XTickLabel',labels,'XTickLabelRotation',15);
title('Forecast skill comparison');
grid on;
