

clear all; close all; clc;

initial_years=1000;
%% set observation
obs(1)=0.5;
for tau=1:initial_years-1
    obs(1+tau)=obs(1+tau-1)+0.1*sin((pi/4)*tau);
end

%% set prediction
prediction_year=99; % 1~100 simulation years;
ens_size=20;
%% initialization for the first year, and prediction
predictor(1:ens_size,1,1) = obs(1)+rand(1,ens_size).*0.1;
error(1:ens_size,1,1) = predictor(1:ens_size,1,1) - obs(1);
for predic_dt=1:prediction_year
    noise=rand(1,ens_size).*0.5-0.25;
    predictor(1:ens_size,1,1+predic_dt)=predictor(1:ens_size,1,1+predic_dt-1)+0.1*sin((pi/4)*predic_dt)+noise';
    error(1:ens_size,1,1+predic_dt) = predictor(1:ens_size,1,1+predic_dt) - obs(1+predic_dt);
end

%% initialization from the second year, and prediction
for tau=1:initial_years-1
    noise=rand(1,ens_size).*0.5-0.25;
    predictor(1:ens_size,1+tau,1)=obs(1+tau)+noise';
    error(1:ens_size,1+tau,1) = noise';
    for predic_dt=1:prediction_year
        noise=rand(1,ens_size).*0.5-0.25;
        predictor(1:ens_size,1+tau,1+predic_dt)=predictor(1:ens_size,1+tau,1+predic_dt-1)+0.1*sin((pi/4)*predic_dt)+noise';
        if tau+predic_dt<=initial_years
            error(1:ens_size,1+tau,1+predic_dt) = predictor(1:ens_size,1+tau,1+predic_dt) - obs(tau+predic_dt);
        else
            error(1:ens_size,1+tau,1+predic_dt) = NaN;
        end
    end
end

%1st simulation year
plot(obs); % observation
hold on;
plot(mean(predictor(:,:,1),1)) % ensemble mean, lead year 1 (1st year simulation)
hold off;

%2nd simulation year
plot(obs(2:end));
hold on;
plot(mean(predictor(:,1:end-1,2),1))
hold off;

%10th simulation year
plot(obs(10:end));
hold on;
plot(mean(predictor(:,1:end-9,10),1))
hold off;

%% ensemble mean, 1st year initialized simulation
plot(obs(1:prediction_year+1));
hold on;
plot(squeeze(mean(predictor(:,1,:),1)))
hold off;
corrcoef(obs(1:prediction_year+1), squeeze(mean(predictor(:,1,:),1)));

%% 1st year initialized simulation, calculation of entropy
for i=1:prediction_year
    ens_mean(i)=squeeze(mean(predictor(:,1,i),1));
    ens_std(i)=std(predictor(:,1,i));
    dx=0.1;
    x=-3:dx:3;
    pdf_norm=pdf('norm',x, ens_mean(i),ens_std(i));
%     plot(x,pdf_norm)
    
    entropy(i)=-sum(pdf_norm .* log2(pdf_norm) .* dx);
end
plot(entropy)
plot(ens_mean)
plot(ens_std)

%% error of the 1st year initialized simulation, calculation of entropy using error
for i=1:prediction_year
    ens_mean(i)=squeeze(mean(error(:,1,i),1));
    ens_std(i)=std(error(:,1,i));
    dx=0.1;
    x=-3:dx:3;
    pdf_norm=pdf('norm',x, ens_mean(i),ens_std(i));
%     plot(x,pdf_norm)
    
    entropy(i)=-sum(pdf_norm .* log2(pdf_norm) .* dx);
end
plot(entropy)
plot(ens_mean)
plot(ens_std)


%% observation, calculation of entropy using obs (norm dist, std==0 -> entropy=NaN. delta func is needed)
for i=1:prediction_year
    ens_mean(i)=obs(i);
    ens_std(i)=0;
    dx=0.1;
    x=-3:dx:3;
    pdf_norm=pdf('norm',x, ens_mean(i),ens_std(i));
%     plot(x,pdf_norm)
    
    entropy(i)=-sum(pdf_norm .* log2(pdf_norm) .* dx);
end
plot(entropy)



%% joint entropy using 1st year initialized calculation, obs 
for i=1:prediction_year
    ens_mean(i)=squeeze(mean(error(:,1,i),1));
    ens_std(i)=std(error(:,1,i));
    dx=0.1;
    x=-3:dx:3;
    pdf_norm=pdf('norm',x, ens_mean(i),ens_std(i));
    pdf_norm_save(i,:)=pdf_norm;
%     plot(x,pdf_norm)
    for j=1:length(pdf_norm)
        if x(j)<= obs(i) && x(j+1) > obs(i)
            pdf_norm(j)= pdf_norm(j);
        else
            pdf_norm(j)= NaN; % 0
        end
    end
    
    entropy(i)=-sum(pdf_norm .* log2(pdf_norm) .* dx, 'omitnan');
end
plot(entropy)
plot(ens_mean)
plot(ens_std)
plot(x,pdf_norm_save(1,:));
hold on;
plot(x,pdf_norm_save(end,:));
hold off;