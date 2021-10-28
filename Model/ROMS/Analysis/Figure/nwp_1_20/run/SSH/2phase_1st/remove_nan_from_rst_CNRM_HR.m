clc; close all; clear all;

addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

rootdir='D:\Data\Model\ROMS\nwp_1_20\input\test2120\';
scenname='historical';
modelname ='CNRM-CM6-1-HR';
expname ='nwp_1_20';
% filedir=[rootdir, scenname, '/', modelname, '/'];
filedir=[rootdir];

year = 2010;
% filename = [filedir, filesep, 'ocean_avg_0365.nc'];
% writefilename = [filedir, filesep, 'ocean_rst_from_avg_0365.nc'];
filename = [filedir, filesep, 'ocean_rst.nc'];
writefilename = [filedir, filesep, 'ocean_rst_fix_nan.nc'];


fileinfo=ncinfo(writefilename);


% vars= {'zeta','ubar','vbar','u','v','w','temp','salt', ...
%     'h', 'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'hc', 'lon_rho', 'lat_rho'};

clear vars
for i= 1:length(fileinfo.Variables(:))
    vars{i}=fileinfo.Variables(i).Name;
end

for i=1:length(vars)
    switch (vars{i})
        case {'Vtransform', 'Vstretching', 'theta_s', 'theta_b', 'hc'}
            vardata.(vars{i})= ncread(filename, vars{i});            
        case {'h', 'f', 'pm', 'pn',  ...
                'lon_rho', 'lat_rho',  'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', ...
                'mask_rho', 'mask_u', 'mask_v', 'mask_psi'}
            vardata.(vars{i})= ncread(filename, vars{i} , [770, 375], [31, 30]);            
        case {'zeta', 'ubar', 'vbar', 'Hsbl', 'Hbbl'}
            vardata.(vars{i})= ncread(filename, vars{i} , [770, 375, 1], [31, 30, 1]);            
        case {'u', 'v', 'temp', 'salt', 'rho'}
            vardata.(vars{i})= ncread(filename, vars{i} , [770, 375, 1, 1], [31, 30, 40, 1]);
        case {'w', 'AKv', 'AKt', 'AKs'} %avg file has 'w' but rst file doesn't have it
            vardata.(vars{i})= ncread(filename, vars{i} , [770, 375, 1, 1], [31, 30, 41, 1]);
    end
end
vardata.N=size(vardata.temp,3);



vardata.lon_rho_3d=repmat(vardata.lon_rho, [1 1 40]);
vardata.lat_rho_3d=repmat(vardata.lat_rho, [1 1 40]);

vardata.lon_rho_3dw=repmat(vardata.lon_rho, [1 1 41]);
vardata.lat_rho_3dw=repmat(vardata.lat_rho, [1 1 41]);

maskvars={'mask_rho', 'mask_u', 'mask_v', 'mask_psi'};
for i=1:length(maskvars)
    vardata.(maskvars{i})(vardata.(maskvars{i})==0)=1;
    ncwrite([writefilename], maskvars{i}, vardata.(maskvars{i}), [770,375]);
end

outvars= {'zeta', 'ubar', 'vbar', 'Hsbl', 'Hbbl'};
for i= 1: length(outvars)
    ismask=isnan(vardata.(outvars{i}));
    vardata.(outvars{i})(ismask)= ...
        griddata(vardata.lon_rho(~ismask), vardata.lat_rho(~ismask), vardata.(outvars{i})(~ismask), ...
        vardata.lon_rho(ismask), vardata.lat_rho(ismask), 'nearest');
        
    ncwrite([writefilename], outvars{i}, vardata.(outvars{i}), [770,375, 1]);
end

vardata.zr = zlevs(vardata.Vtransform, vardata.Vstretching, vardata.h, vardata.zeta, ...
    vardata.theta_s, vardata.theta_b, vardata.hc,vardata.N,'r');
vardata.zw = zlevs(vardata.Vtransform, vardata.Vstretching, vardata.h, vardata.zeta, ...
    vardata.theta_s, vardata.theta_b, vardata.hc,vardata.N,'w');
vardata.zr_p=permute(vardata.zr, [2 3 1]);
vardata.zw_p=permute(vardata.zw, [2 3 1]);

outvars= {'u', 'v', 'temp', 'salt', 'rho'};
for i= 1: length(outvars)
    ismask=isnan(vardata.(outvars{i}));
    vardata.(outvars{i})(ismask)= ...
        griddata(vardata.lon_rho_3d(~ismask), vardata.lat_rho_3d(~ismask), vardata.zr_p(~ismask), ...
        vardata.(outvars{i})(~ismask), vardata.lon_rho_3d(ismask), vardata.lat_rho_3d(ismask), vardata.zr_p(ismask), 'nearest');
    
    ncwrite([writefilename], outvars{i}, vardata.(outvars{i}), [770,375,1,1]);
end


outvars= {'AKv', 'AKt', 'AKs'}; %avg file has 'w' but rst file doesn't have it
for i= 1: length(outvars)
    ismask=isnan(vardata.(outvars{i}));
    vardata.(outvars{i})(ismask)= ...
        griddata(vardata.lon_rho_3dw(~ismask), vardata.lat_rho_3dw(~ismask), vardata.zw_p(~ismask), ...
        vardata.(outvars{i})(~ismask), vardata.lon_rho_3dw(ismask), vardata.lat_rho_3dw(ismask), vardata.zw_p(ismask), 'nearest');
    
    ncwrite([writefilename], outvars{i}, vardata.(outvars{i}), [770,375,1,1]);
end



