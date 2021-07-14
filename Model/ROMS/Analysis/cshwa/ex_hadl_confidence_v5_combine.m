clear all; close all; clc;
cd E:\OA\new_result
load ex_hadl_R_analysis_result_full
load ex_R_hadl_analysis_result_FULL_V2
save('ex_hadl_R_result_FULL_combine.mat') 

clear all; close all; clc;
cd E:\OA\new_result


clear all; close all; clc;

CI = 0; % if off this confidence 95%(CI = any number),  on : 90% (CI =1)
for i = 1:5
    
load ex_hadl_R_result_FULL_combine.mat

R_AO_88 = R_AO_88(1,i);
R_AO_98 = R_AO_98(1,i);
R_AO_16 = R_AO_16(1,i);

R_ENSO_88 = R_ENSO_88(1,i);
R_ENSO_98 = R_ENSO_98(1,i);
R_ENSO_16 = R_ENSO_16(1,i);

R_PDO_88 = R_PDO_88(1,i);
R_PDO_98 = R_PDO_98(1,i);
R_PDO_16 = R_PDO_16(1,i);

R_NINO12_88 = R_nino12_88(1,i);
R_NINO12_98 = R_nino12_98(1,i);
R_NINO12_16 = R_nino12_16(1,i);

R_NINO34_88 = R_nino34_88(1,i);
R_NINO34_98 = R_nino34_98(1,i);
R_NINO34_16 = R_nino34_16(1,i);

R_NINO3_88 = R_nino3_88(1,i);
R_NINO3_98 = R_nino3_98(1,i);
R_NINO3_16 = R_nino3_16(1,i);

R_NINO4_88 = R_nino4_88(1,i);
R_NINO4_98 = R_nino4_98(1,i);
R_NINO4_16 = R_nino4_16(1,i);

%2)rho = 0?  if( abs(t) > t_c ) means good to compute R
N_88 = 83; % data size : 1982~2016 (420 month)
N_98 = 121; % data size : 1982~2016 (420 month)
N_16 = 216; % data size : 1982~2016 (420 month)
% t = R * sqrt(N-2) / sqrt(1-R^2);    % t value, transformed from r
t_c_88 = tinv( 0.975, N_88 - 2 );           % t statistic critical value, nu=N-2, 95% level
t_c_98 = tinv( 0.975, N_98 - 2 );
t_c_16 = tinv( 0.975, N_16 - 2 );

if CI == 1    % t statistic critical value, nu=N-2, 90% level
    t_c_88 = tinv( 0.950, N_88 - 2 );           
    t_c_98 = tinv( 0.950, N_98 - 2 );
    t_c_16 = tinv( 0.950, N_16 - 2 );
end


t_ENSO_88 = R_ENSO_88 * sqrt(N_88 - 2) / sqrt(1- R_ENSO_88^2);
t_ENSO_98 = R_ENSO_98 * sqrt(N_98 - 2) / sqrt(1- R_ENSO_98^2);
t_ENSO_16 = R_ENSO_16 * sqrt(N_16 - 2) / sqrt(1- R_ENSO_16^2);

t_PDO_88 = R_PDO_88 * sqrt(N_88 - 2) / sqrt(1 - R_PDO_88^2);
t_PDO_98 = R_PDO_98 * sqrt(N_98 - 2) / sqrt(1 - R_PDO_98^2);
t_PDO_16 = R_PDO_16 * sqrt(N_16 - 2) / sqrt(1 - R_PDO_16^2);

t_AO_88 = R_AO_88 * sqrt(N_88 - 2) / sqrt(1 - R_AO_88^2);
t_AO_98 = R_AO_98 * sqrt(N_98 - 2) / sqrt(1 - R_AO_98^2);
t_AO_16 = R_AO_16 * sqrt(N_16 - 2) / sqrt(1 - R_AO_16^2);

t_NINO12_88 = R_NINO12_88 * sqrt(N_88 - 2) / sqrt(1-R_NINO12_88^2);
t_NINO12_98 = R_NINO12_98 * sqrt(N_98 - 2) / sqrt(1-R_NINO12_98^2);
t_NINO12_16 = R_NINO12_16 * sqrt(N_16 - 2) / sqrt(1-R_NINO12_16^2);

t_NINO34_88 = R_NINO34_88 * sqrt(N_88 - 2) / sqrt(1-R_NINO34_88^2);
t_NINO34_98 = R_NINO34_98 * sqrt(N_98 - 2) / sqrt(1-R_NINO34_98^2);
t_NINO34_16 = R_NINO34_16 * sqrt(N_16 - 2) / sqrt(1-R_NINO34_16^2);

t_NINO3_88 = R_NINO3_88 * sqrt(N_88 - 2) / sqrt(1-R_NINO3_88^2);
t_NINO3_98 = R_NINO3_98 * sqrt(N_98 - 2) / sqrt(1-R_NINO3_98^2);
t_NINO3_16 = R_NINO3_16 * sqrt(N_16 - 2) / sqrt(1-R_NINO3_16^2);

t_NINO4_88 = R_NINO4_88 * sqrt(N_88 - 2) / sqrt(1-R_NINO4_88^2);
t_NINO4_98 = R_NINO4_98 * sqrt(N_98 - 2) / sqrt(1-R_NINO4_98^2);
t_NINO4_16 = R_NINO4_16 * sqrt(N_16 - 2) / sqrt(1-R_NINO4_16^2);

disp(i) %disp mode    
if( abs(t_ENSO_88) > t_c_88 )
disp('can compute ENSO_88')
end
if( abs(t_ENSO_98) > t_c_98 ) 
disp('can compute ENSO_98')
end
if( abs(t_ENSO_16) > t_c_16 ) 
disp('can compute ENSO_16')
end
if( abs(t_PDO_88) > t_c_88 )
disp('can compute PDO_88')
end
if( abs(t_PDO_98) > t_c_98 ) 
disp('can compute PDO_98')
end    
if( abs(t_PDO_16) > t_c_16 )
disp('can compute PDO_16')
end    
if( abs(t_AO_88) > t_c_88 )
disp('can compute AO_88')
end    
if( abs(t_AO_98) > t_c_98 ) 
disp('can compute AO_98')
end       
if( abs(t_AO_16) > t_c_16 ) 
disp('can compute AO_16')
end
if( abs(t_NINO12_88) > t_c_88 ) 
disp('can compute NINO12_88')
end
if( abs(t_NINO12_98) > t_c_98 ) 
disp('can compute NINO12_98')
end
if( abs(t_NINO12_16) > t_c_16 ) 
disp('can compute NINO12_16')
end
if( abs(t_NINO34_88) > t_c_88 )
disp('can compute NINO34_88')
end
if( abs(t_NINO34_98) > t_c_98 ) 
disp('can compute NINO34_98')
end    
if( abs(t_NINO34_16) > t_c_16 )
disp('can compute NINO34_16')
end    
if( abs(t_NINO3_88) > t_c_88 )
disp('can compute NINO3_88')
end    
if( abs(t_NINO3_98) > t_c_98 ) 
disp('can compute NINO3_98')
end       
if( abs(t_NINO3_16) > t_c_16 ) 
disp('can compute NINO3_16')
end 
if( abs(t_NINO4_88) > t_c_88 )
disp('can compute NINO4_88')
end    
if( abs(t_NINO4_98) > t_c_98 ) 
disp('can compute NINO4_98')
end       
if( abs(t_NINO4_16) > t_c_16 ) 
disp('can compute NINO4_16')
end

%3) true rho range
z_025_88 = t_c_88;
z_025_98 = t_c_98;
z_025_16 = t_c_16;

sigma_z_88 = 1/sqrt(N_88 - 3);
sigma_z_98 = 1/sqrt(N_98 - 3);
sigma_z_16 = 1/sqrt(N_16 - 3);

z_ENSO_88 = 1/2 * log((1+R_ENSO_88)/(1-R_ENSO_88));
z_ENSO_98 = 1/2 * log((1+R_ENSO_98)/(1-R_ENSO_98));
z_ENSO_16 = 1/2 * log((1+R_ENSO_16)/(1-R_ENSO_16));

z_PDO_88 = 1/2 * log((1+R_PDO_88)/(1-R_PDO_88));
z_PDO_98 = 1/2 * log((1+R_PDO_98)/(1-R_PDO_98)); 
z_PDO_16 = 1/2 * log((1+R_PDO_16)/(1-R_PDO_16));

z_AO_88 = 1/2 * log((1+R_AO_88)/(1-R_AO_88));
z_AO_98 = 1/2 * log((1+R_AO_98)/(1-R_AO_98));
z_AO_16 = 1/2 * log((1+R_AO_16)/(1-R_AO_16));

z_NINO12_88 = 1/2 * log((1+R_NINO12_88)/(1-R_NINO12_88));
z_NINO12_98 = 1/2 * log((1+R_NINO12_98)/(1-R_NINO12_98));
z_NINO12_16 = 1/2 * log((1+R_NINO12_16)/(1-R_NINO12_16));

z_NINO34_88 = 1/2 * log((1+R_NINO34_88)/(1-R_NINO34_88));
z_NINO34_98 = 1/2 * log((1+R_NINO34_98)/(1-R_NINO34_98));
z_NINO34_16 = 1/2 * log((1+R_NINO34_16)/(1-R_NINO34_16));

z_NINO3_88= 1/2 * log((1+R_NINO3_88)/(1-R_NINO3_88));
z_NINO3_98 = 1/2 * log((1+R_NINO3_98)/(1-R_NINO3_98));
z_NINO3_16 = 1/2 * log((1+R_NINO3_16)/(1-R_NINO3_16));

z_NINO4_88 = 1/2 * log((1+R_NINO4_88)/(1-R_NINO4_88));
z_NINO4_98 = 1/2 * log((1+R_NINO4_98)/(1-R_NINO4_98));
z_NINO4_16 = 1/2 * log((1+R_NINO4_16)/(1-R_NINO4_16));


%mu_z

left_ENSO_88 = z_ENSO_88 - z_025_88 * sigma_z_88; 
right_ENSO_88 = z_ENSO_88 + z_025_88 * sigma_z_88;
left_ENSO_98 = z_ENSO_98 - z_025_98 * sigma_z_98; 
right_ENSO_98 = z_ENSO_98 + z_025_98 * sigma_z_98; 
left_ENSO_16 = z_ENSO_16 - z_025_16 * sigma_z_16;
right_ENSO_16 = z_ENSO_16 + z_025_16 * sigma_z_16; 

left_PDO_88  = z_PDO_88  - z_025_88 * sigma_z_88; 
right_PDO_88  = z_PDO_88  + z_025_88 * sigma_z_88; 
left_PDO_98 = z_PDO_98 - z_025_98 * sigma_z_98; 
right_PDO_98 = z_PDO_98 + z_025_98 * sigma_z_98; 
left_PDO_16 = z_PDO_16 - z_025_16 * sigma_z_16;
right_PDO_16 = z_PDO_16 + z_025_16 * sigma_z_16; 

left_AO_88 = z_AO_88 - z_025_88 * sigma_z_88; 
right_AO_88 = z_AO_88 + z_025_88 * sigma_z_88; 
left_AO_98 = z_AO_98 - z_025_98 * sigma_z_98; 
right_AO_98 = z_AO_98 + z_025_98 * sigma_z_98; 
left_AO_16 = z_AO_16 - z_025_16 * sigma_z_16; 
right_AO_16 = z_AO_16 + z_025_16 * sigma_z_16; 

left_NINO12_88 = z_NINO12_88 - z_025_88 * sigma_z_88; 
right_NINO12_88 = z_NINO12_88 + z_025_88 * sigma_z_88; 
left_NINO12_98 = z_NINO12_98 - z_025_98 * sigma_z_98; 
right_NINO12_98 = z_NINO12_98 + z_025_98 * sigma_z_98; 
left_NINO12_16 = z_NINO12_16 - z_025_16 * sigma_z_16; 
right_NINO12_16 = z_NINO12_16 + z_025_16 * sigma_z_16;

left_NINO34_88 = z_NINO34_88 - z_025_88 * sigma_z_88;
right_NINO34_88 = z_NINO34_88 + z_025_88 * sigma_z_88; 
left_NINO34_98 = z_NINO34_98 - z_025_98 * sigma_z_98; 
right_NINO34_98 = z_NINO34_98 + z_025_98 * sigma_z_98; 
left_NINO34_16 = z_NINO34_16 - z_025_16 * sigma_z_16; 
right_NINO34_16 = z_NINO34_16 + z_025_16 * sigma_z_16; 

left_NINO3_88 = z_NINO3_88 - z_025_88 * sigma_z_88; 
right_NINO3_88 = z_NINO3_88 + z_025_88 * sigma_z_88; 
left_NINO3_98 = z_NINO3_98 - z_025_98 * sigma_z_98;
right_NINO3_98 = z_NINO3_98 + z_025_98 * sigma_z_98; 
left_NINO3_16 = z_NINO3_16 - z_025_16 * sigma_z_16; 
right_NINO3_16 = z_NINO3_16 + z_025_16 * sigma_z_16; 

left_NINO4_88 = z_NINO4_88 - z_025_88 * sigma_z_88; 
right_NINO4_88 = z_NINO4_88 + z_025_88 * sigma_z_88; 
left_NINO4_98 = z_NINO4_98 - z_025_98 * sigma_z_98;
right_NINO4_98 = z_NINO4_98 + z_025_98 * sigma_z_98; 
left_NINO4_16 = z_NINO4_16 - z_025_16 * sigma_z_16; 
right_NINO4_16 = z_NINO4_16 + z_025_16 * sigma_z_16; 


%rho
rho_low_ENSO_88(i)= (exp(2*left_ENSO_88) - 1)/(exp(2*left_ENSO_88) + 1);
rho_high_ENSO_88(i) = (exp(2*right_ENSO_88) - 1)/(exp(2*right_ENSO_88) + 1);
rho_low_ENSO_98(i) = (exp(2*left_ENSO_98) - 1)/(exp(2*left_ENSO_98) + 1);
rho_high_ENSO_98(i) = (exp(2*right_ENSO_98) - 1)/(exp(2*right_ENSO_98) + 1);
rho_low_ENSO_16(i) = (exp(2*left_ENSO_16) - 1)/(exp(2*left_ENSO_16) + 1);
rho_high_ENSO_16(i) = (exp(2*right_ENSO_16) - 1)/(exp(2*right_ENSO_16) + 1);

rho_low_PDO_88(i) = (exp(2*left_PDO_88) - 1)/(exp(2*left_PDO_88) + 1);
rho_high_PDO_88(i) = (exp(2*right_PDO_88) - 1)/(exp(2*right_PDO_88) + 1);
rho_low_PDO_98(i) = (exp(2*left_PDO_98) - 1)/(exp(2*left_PDO_98) + 1);
rho_high_PDO_98(i) = (exp(2*right_PDO_98) - 1)/(exp(2*right_PDO_98) + 1);
rho_low_PDO_16(i) = (exp(2*left_PDO_16) - 1)/(exp(2*left_PDO_16) + 1);
rho_high_PDO_16(i) = (exp(2*right_PDO_16) - 1)/(exp(2*right_PDO_16) + 1);

rho_low_AO_88(i) = (exp(2*left_AO_88) - 1)/(exp(2*left_AO_88) + 1);
rho_high_AO_88(i) = (exp(2*right_AO_88) - 1)/(exp(2*right_AO_88) + 1);
rho_low_AO_98(i) = (exp(2*left_AO_98) - 1)/(exp(2*left_AO_98) + 1);
rho_high_AO_98(i) = (exp(2*right_AO_98) - 1)/(exp(2*right_AO_98) + 1);
rho_low_AO_16(i) = (exp(2*left_AO_16) - 1)/(exp(2*left_AO_16) + 1);
rho_high_AO_16(i) = (exp(2*right_AO_16) - 1)/(exp(2*right_AO_16) + 1);

rho_low_NINO12_88(i) = (exp(2*left_NINO12_88) - 1)/(exp(2*left_NINO12_88) + 1);
rho_high_NINO12_88(i) = (exp(2*right_NINO12_88) - 1)/(exp(2*right_NINO12_88) + 1);
rho_low_NINO12_98(i) = (exp(2*left_NINO12_98) - 1)/(exp(2*left_NINO12_98) + 1);
rho_high_NINO12_98(i) = (exp(2*right_NINO12_98) - 1)/(exp(2*right_NINO12_98) + 1);
rho_low_NINO12_16(i) = (exp(2*left_NINO12_16) - 1)/(exp(2*left_NINO12_16) + 1);
rho_high_NINO12_16(i) = (exp(2*right_NINO12_16) - 1)/(exp(2*right_NINO12_16) + 1);

rho_low_NINO34_88(i) = (exp(2*left_NINO34_88) - 1)/(exp(2*left_NINO34_88) + 1);
rho_high_NINO34_88(i) = (exp(2*right_NINO34_88) - 1)/(exp(2*right_NINO34_88) + 1);
rho_low_NINO34_98(i) = (exp(2*left_NINO34_98) - 1)/(exp(2*left_NINO34_98) + 1);
rho_high_NINO34_98(i) = (exp(2*right_NINO34_98) - 1)/(exp(2*right_NINO34_98) + 1);
rho_low_NINO34_16(i) = (exp(2*left_NINO34_16) - 1)/(exp(2*left_NINO34_16) + 1);
rho_high_NINO34_16(i) = (exp(2*right_NINO34_16) - 1)/(exp(2*right_NINO34_16) + 1);

rho_low_NINO3_88(i) = (exp(2*left_NINO3_88) - 1)/(exp(2*left_NINO3_88) + 1);
rho_high_NINO3_88(i) = (exp(2*right_NINO3_88) - 1)/(exp(2*right_NINO3_88) + 1);
rho_low_NINO3_98(i) = (exp(2*left_NINO3_98) - 1)/(exp(2*left_NINO3_98) + 1);
rho_high_NINO3_98(i) = (exp(2*right_NINO3_98) - 1)/(exp(2*right_NINO3_98) + 1);
rho_low_NINO3_16(i) = (exp(2*left_NINO3_16) - 1)/(exp(2*left_NINO3_16) + 1);
rho_high_NINO3_16(i) = (exp(2*right_NINO3_16) - 1)/(exp(2*right_NINO3_16) + 1);

rho_low_NINO4_88(i) = (exp(2*left_NINO4_88) - 1)/(exp(2*left_NINO4_88) + 1);
rho_high_NINO4_88(i) = (exp(2*right_NINO4_88) - 1)/(exp(2*right_NINO4_88) + 1);
rho_low_NINO4_98(i) = (exp(2*left_NINO4_98) - 1)/(exp(2*left_NINO4_98) + 1);
rho_high_NINO4_98(i) = (exp(2*right_NINO4_98) - 1)/(exp(2*right_NINO4_98) + 1);
rho_low_NINO4_16(i) = (exp(2*left_NINO4_16) - 1)/(exp(2*left_NINO4_16) + 1);
rho_high_NINO4_16(i) = (exp(2*right_NINO4_16) - 1)/(exp(2*right_NINO4_16) + 1);

if i == i(end)
    diary('ex_hadl_log_FULL_confidence') % save the display massage to file
end
    
end

save('ex_hadl_rho_interval_result_full','rho_*');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; clc;  close all;
load ex_hadl_R_result_FULL_combine.mat
load ex_hadl_rho_interval_result_full

%%%%%% 0 period
    % R
NINO12_0 = R_nino12_88(1:5);
NINO34_0 = R_nino34_88(1:5);
NINO3_0 = R_nino3_88(1:5);
NINO4_0 = R_nino4_88(1:5);
data_0 = [NINO12_0; NINO34_0; NINO3_0; NINO4_0];
data_0 = round(data_0,2); % cut below 10^-2 

    % rho - low
low_NINO12_0 = rho_low_NINO12_88(1:5);
low_NINO34_0 = rho_low_NINO34_88(1:5);
low_NINO3_0 = rho_low_NINO3_88(1:5);
low_NINO4_0 = rho_low_NINO4_88(1:5);
low_erro_0 = [low_NINO12_0; low_NINO34_0; low_NINO3_0; low_NINO4_0];
low_erro_0 = round(low_erro_0,2); % cut below 10^-2 

    % rho - high
high_NINO12_0 = rho_high_NINO12_88(1:5);
high_NINO34_0 = rho_high_NINO34_88(1:5);
high_NINO3_0 = rho_high_NINO3_88(1:5);
high_NINO4_0 = rho_high_NINO4_88(1:5);
high_erro_0 = [high_NINO12_0; high_NINO34_0; high_NINO3_0; high_NINO4_0];
high_erro_0 = round(high_erro_0,2); % cut behigh 10^-2 

%%%%%% 1 period
NINO12_1 = R_nino12_98(1:5);
NINO34_1 = R_nino34_98(1:5);
NINO3_1 = R_nino3_98(1:5);
NINO4_1 = R_nino4_98(1:5);
data_1 = [NINO12_1; NINO34_1; NINO3_1; NINO4_1];
data_1 = round(data_1,2); % cut below 10^-2 

% rho - low
low_NINO12_1 = rho_low_NINO12_98(1:5);
low_NINO34_1 = rho_low_NINO34_98(1:5);
low_NINO3_1 = rho_low_NINO3_98(1:5);
low_NINO4_1 = rho_low_NINO4_98(1:5);
low_erro_1 = [low_NINO12_1; low_NINO34_1; low_NINO3_1; low_NINO4_1];
low_erro_1 = round(low_erro_1,2); % cut below 10^-2 

    % rho - high
high_NINO12_1 = rho_high_NINO12_98(1:5);
high_NINO34_1 = rho_high_NINO34_98(1:5);
high_NINO3_1 = rho_high_NINO3_98(1:5);
high_NINO4_1 = rho_high_NINO4_98(1:5);
high_erro_1 = [high_NINO12_1; high_NINO34_1; high_NINO3_1;  high_NINO4_1];
high_erro_1 = round(high_erro_1,2); % cut behigh 10^-2 

% 2 period
NINO12_2 = R_nino12_16(1:5);
NINO34_2 = R_nino34_16(1:5);
NINO3_2 = R_nino3_16(1:5);
NINO4_2 = R_nino4_16(1:5);
data_2 = [NINO12_2; NINO34_2; NINO3_2; NINO4_2];
data_2 = round(data_2,2); % cut below 10^-2 

% rho - low
low_NINO12_2 = rho_low_NINO12_16(1:5);
low_NINO34_2 = rho_low_NINO34_16(1:5);
low_NINO3_2 = rho_low_NINO3_16(1:5);
low_NINO4_2 = rho_low_NINO4_16(1:5);
low_erro_2 = [low_NINO12_2; low_NINO34_2; low_NINO3_2; low_NINO4_2];
low_erro_2 = round(low_erro_2,2); % cut below 10^-2 

    % rho - high
high_NINO12_2 = rho_high_NINO12_16(1:5);
high_NINO34_2 = rho_high_NINO34_16(1:5);
high_NINO3_2 = rho_high_NINO3_16(1:5);
high_NINO4_2 = rho_high_NINO4_16(1:5);
high_erro_2 = [high_NINO12_2; high_NINO34_2; high_NINO3_2; high_NINO4_2];
high_erro_2 = round(high_erro_2,2); % cut behigh 10^-2 

myaxis_0 = linspace(0,1,2);
myaxis_1 = linspace(3,1,5);
myaxis_2 = linspace(6,1,8);

myaxis_0 = linspace(0,1,2);
myaxis_1 = linspace(3,1,5);
myaxis_2 = linspace(6,1,8);

% % % hist(myaxis_0,NINO12_P1_0)
% % bar(NINO12_P1_0,'stack');
% % NINO34_1_0 = bar(NINO34_P1_0,'stack');
% % subplot(NINO34_1_0);


figure(1)
ax1 = subplot(1,3,1); 
hb_0 = bar(ax1,data_0)
width = hb_0.BarWidth;
for i=1:length(data_0(:, 1))
    row = data_0(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
% bar(ax1,high_erro_0)
% bar(ax1,low_erro_0)
error_0 = (high_erro_0 - low_erro_0)./2 % compute error half interval 
barwitherr(error_0, data_0)
plot_line = zeros(1,6) % line plot have to match demension
plot_line(:,:) = 0.2159 % upside Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.2159 % down side Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
xlabel = {'NINO12','NINO34','NINO3','NINO4'};
set(gca,'xticklabel',xlabel);
ylim([-0.8 +0.8]);
set(gca, 'YTick', [-0.8:0.05:+0.8], 'YTickLabel',[-0.8:0.05:0.8])
ylabel('Correaltion Coefficent');
title('1982~1988 FULL region');
% title('1982~1991.07 P1 region')

ax2 = subplot(1,3,2); 
hb_1 = bar(ax2,data_1)
width = hb_1.BarWidth;
for i=1:length(data_1(:, 1))
    row = data_1(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
error_1 = (high_erro_1 - low_erro_1)./2  % compute error half interval 
barwitherr(error_1, data_1)
plot_line = zeros(1,6) % line plot have to match demension
plot_line(:,:) = 0.1786 % upside Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.1786 % down side Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
set(gca,'xticklabel',xlabel);
ylim([-0.8 +0.8]);
set(gca, 'YTick', [-0.8:0.05:+0.8], 'YTickLabel',[-0.8:0.05:0.8])
ylabel('Correaltion Coefficent');
title('1989~1998 FULL region');

ax3 = subplot(1,3,3);
hb_2 = bar(ax3,data_2)
width = hb_2.BarWidth;
for i=1:length(data_2(:, 1))
    row = data_2(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
error_2 = (high_erro_2 - low_erro_2)./2  % compute error half interval 
barwitherr(error_2, data_2)
plot_line = zeros(1,6) % line plot have to match demension
plot_line(:,:) = 0.1337 % upside Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.1337 % down side Siginificant R critical value 
line((0:5),plot_line,'color','r') % plot Siginificant R critical value 
set(gca,'xticklabel',xlabel);
ylim([-0.8 +0.8]);
set(gca, 'YTick', [-0.8:0.05:+0.8], 'YTickLabel',[-0.8:0.05:0.8])
ylabel('Correaltion Coefficent');
title('1999~2016 FULL region');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CLIMATE index
clear all; clc; close all;
load ex_hadl_R_result_FULL_combine.mat
load ex_hadl_rho_interval_result_full
%%%%%% 0 period
    % R
AO_0 = R_AO_88(1:5);
ENSO_0 = R_ENSO_88(1:5);
PDO_0 = R_PDO_88(1:5);
data_0 = [AO_0; ENSO_0; PDO_0];
data_0 = round(data_0,2); % cut below 10^-2 

    % rho - low
low_AO_0 = rho_low_AO_88(1:5);
low_ENSO_0 = rho_low_ENSO_88(1:5);
low_PDO_0 = rho_low_PDO_88(1:5);
low_erro_0 = [low_AO_0; low_ENSO_0; low_PDO_0];
low_erro_0 = round(low_erro_0,2); % cut below 10^-2 

    % rho - high
high_AO_0 = rho_high_AO_88(1:5);
high_ENSO_0 = rho_high_ENSO_88(1:5);
high_PDO_0 = rho_high_PDO_88(1:5);
high_erro_0 = [high_AO_0; high_ENSO_0; high_PDO_0];
high_erro_0 = round(high_erro_0,2); % cut behigh 10^-2 

%%%%%% 1 period
AO_1 = R_AO_98(1:5);
ENSO_1 = R_ENSO_98(1:5);
PDO_1 = R_PDO_98(1:5);
data_1 = [AO_1; ENSO_1; PDO_1];
data_1 = round(data_1,2); % cut below 10^-2 
% 3 period
AO_2 = R_AO_16(1:5);
ENSO_2 = R_ENSO_16(1:5);
PDO_2 = R_PDO_16(1:5);
data_2 = [AO_2; ENSO_2; PDO_2];
data_2 = round(data_2,2); % cut below 10^-2 

% rho - low
low_AO_1 = rho_low_AO_98(1:5);
low_ENSO_1 = rho_low_ENSO_98(1:5);
low_PDO_1 = rho_low_PDO_98(1:5);
low_erro_1 = [low_AO_1; low_ENSO_1; low_PDO_1];
low_erro_1 = round(low_erro_1,2); % cut below 10^-2 

    % rho - high
high_AO_1 = rho_high_AO_98(1:5);
high_ENSO_1 = rho_high_ENSO_98(1:5);
high_PDO_1 = rho_high_PDO_98(1:5);
high_erro_1 = [high_AO_1; high_ENSO_1; high_PDO_1];
high_erro_1 = round(high_erro_1,2); % cut behigh 10^-2 

% 3 period
AO_2 = R_AO_16(1:5);
ENSO_2 = R_ENSO_16(1:5);
PDO_2 = R_PDO_16(1:5);
data_2 = [AO_2; ENSO_2; PDO_2];
data_2 = round(data_2,2); % cut below 10^-2 

% rho - low
low_AO_2 = rho_low_AO_16(1:5);
low_ENSO_2 = rho_low_ENSO_16(1:5);
low_PDO_2 = rho_low_PDO_16(1:5);
low_erro_2 = [low_AO_2; low_ENSO_2; low_PDO_2];
low_erro_2 = round(low_erro_2,2); % cut below 10^-2 

    % rho - high
high_AO_2 = rho_high_AO_16(1:5);
high_ENSO_2 = rho_high_ENSO_16(1:5);
high_PDO_2 = rho_high_PDO_16(1:5);
high_erro_2 = [high_AO_2; high_ENSO_2; high_PDO_2];
high_erro_2 = round(high_erro_2,2); % cut behigh 10^-2 

myaxis_0 = linspace(0,1,2);
myaxis_1 = linspace(3,1,5);
myaxis_2 = linspace(6,1,8);

myaxis_0 = linspace(0,1,2);
myaxis_1 = linspace(3,1,5);
myaxis_2 = linspace(6,1,8);

% % % hist(myaxis_0,AO_P1_0)
% % bar(AO_P1_0,'stack');
% % ENSO_1_0 = bar(ENSO_P1_0,'stack');
% % subplot(ENSO_1_0);


figure(1)
ax1 = subplot(1,3,1); 
hb_0 = bar(ax1,data_0)
width = hb_0.BarWidth;
for i=1:length(data_0(:, 1))
    row = data_0(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
% bar(ax1,high_erro_0)
% bar(ax1,low_erro_0)
error_0 = (high_erro_0 - low_erro_0)./2 % compute error half interval 
barwitherr(error_0, data_0)
plot_line = zeros(1,5) % line plot have to match demension
plot_line(:,:) = 0.2159 % upside Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.2159 % down side Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
% set(hb_0(1), 'FaceColor','r')
% set(hb_0(2), 'FaceColor','b')
% set(hb_0(3), 'FaceColor','g')
xlabel = {'AO','ENSO','PDO'};
set(gca,'xticklabel',xlabel);
ylim([-0.7 +0.7]);
set(gca, 'YTick', [-0.7:0.05:+0.7], 'YTickLabel',[-0.7:0.05:0.7])
ylabel('Correaltion Coefficent');
title('1982~1988 FULL region');
% title('1982~1991.07 P1 region')

ax2 = subplot(1,3,2); 
hb_1 = bar(ax2,data_1)
width = hb_1.BarWidth;
for i=1:length(data_1(:, 1))
    row = data_1(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
error_1 = (high_erro_1 - low_erro_1)./2  % compute error half interval 
barwitherr(error_1, data_1)
plot_line = zeros(1,5) % line plot have to match demension
plot_line(:,:) = 0.1786 % upside Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.1786 % down side Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
% set(hb_1(1), 'FaceColor','r')
% set(hb_1(2), 'FaceColor','b')
% set(hb_1(3), 'FaceColor','g')
set(gca,'xticklabel',xlabel);
ylim([-0.7 +0.7]);
set(gca, 'YTick', [-0.7:0.05:+0.7], 'YTickLabel',[-0.7:0.05:0.7])
ylabel('Correaltion Coefficent');
title('1989~1998 FULL region');

ax3 = subplot(1,3,3);
hb_2 = bar(ax3,data_2)
width = hb_2.BarWidth;
for i=1:length(data_2(:, 1))
    row = data_2(i, :);
    % 0.5 is approximate net width of white spacings per group
    offset = ((width + 2) / length(row)) / 2;
    x = linspace(i-offset, i+offset, length(row));
    text(x,row,num2str(row'),'vert','bottom','horiz','center');
end
error_2 = (high_erro_2 - low_erro_2)./2  % compute error half interval 
barwitherr(error_2, data_2);
plot_line = zeros(1,5) % line plot have to match demension
plot_line(:,:) = 0.1337 % upside Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
plot_line(:,:) = - 0.1337 % down side Siginificant R critical value 
line((0:4),plot_line,'color','r') % plot Siginificant R critical value 
% set(hb_2(1), 'FaceColor','r')
% set(hb_2(2), 'FaceColor','b')
% set(hb_2(3), 'FaceColor','g')
set(gca,'xticklabel',xlabel);
ylim([-0.7 +0.7]);
set(gca, 'YTick', [-0.7:0.05:+0.7], 'YTickLabel',[-0.7:0.05:0.7])
ylabel('Correaltion Coefficent');
title('1999~2016 FULL region');
% legend(hb_2(:),{'1mode','2mode','3mode','4mode','5mode'});

