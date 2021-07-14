
clear all; clc; close all;

CI = 0
for i = 0 : 0.0001 : 0.3
 
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

t_sample_88 = i * sqrt(N_88 - 2) / sqrt(1- i^2);
t_sample_98 = i * sqrt(N_98 - 2) / sqrt(1- i^2);
t_sample_16 = i * sqrt(N_16 - 2) / sqrt(1- i^2);


 
if( abs(t_sample_88) > t_c_88 )
    disp(i) % indicate R 
disp('biger than this R will be significant 88')
end
if( abs(t_sample_98) > t_c_98 ) 
    disp(i) % indicate R 
disp('biger than this R will be significant 98')
end
if( abs(t_sample_16) > t_c_16 )
    disp(i) % indicate R 
disp('biger than this R will be significant 16')
end
if i == i(end)
    diary('log_FULL_confidence') % save the display massage to file
end
end
