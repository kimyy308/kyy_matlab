clear; clc

addpath(genpath('/home/jhjung/Matlab/toolbox'));

for yyyy = 2014:2014
for mm = 4:12
    clearvars -except mm yyyy
    HYCOM_compile_monthly 
end
end
