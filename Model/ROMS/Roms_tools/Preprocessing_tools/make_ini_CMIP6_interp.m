%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Build a ROMS initial file from CMIP6 CNRM-ESM2 data
%
%  Extrapole and interpole temperature, salinity, u, v and ssh from a
%  Climatology to get initial conditions for
%  ROMS (initial netcdf files) .
%
%  Data input format (netcdf):
%     temperature(T, Z, Y, X)
%     T : time [Months]
%     Z : st_ocean(depth) [m]
%     Y : yt_ocean [degree north]
%     X : xt_ocean [degree east]
%
%     Updated 22-Jul-2021 by Yong-Yub Kim
%     Updated 27-Jul-2021 by Yong-Yub Kim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
%%%%%%%%%%%%%%%%%%%%% USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%
%
% Common parameters
%
romstools_param



CMIP6_info.rootdir='E:\Data\Model\CMIP6';
% CMIP6_info.modelname='CNRM-ESM2-1';
% CMIP6_info.modelname='EC-Earth3-Veg';
% CMIP6_info.modelname='ACCESS-CM2';
% CMIP6_info.modelname='CNRM-CM6-1-HR';
CMIP6_info.modelname='CMCC-ESM2';

CMIP6_info.varname{1} = 'thetao';
CMIP6_info.varname{2} = 'so';
CMIP6_info.varname{3} = 'uo';
CMIP6_info.varname{4} = 'vo';
CMIP6_info.varname{5} = 'zos';
CMIP6_info.scenname = 'historical';
% CMIP6_info.ensname = 'r1i1p1f2';
CMIP6_info.ensname=Func_0015_GCM_CMIP6_ensname(CMIP6_info.modelname);
CMIP6_info.iniyear = '1985';

ROMS_info.varname{1} = 'temp';
ROMS_info.varname{2} = 'salt';
ROMS_info.varname{3} = 'u';
ROMS_info.varname{4} = 'v';
ROMS_info.varname{5} = 'zeta';



for vari=1:length(CMIP6_info.varname)
    CMIP6_info.filename{vari}=[CMIP6_info.rootdir, filesep, CMIP6_info.varname{vari}, filesep, ...
        CMIP6_info.scenname, filesep, 'interp', filesep, CMIP6_info.modelname, filesep, ...
        CMIP6_info.varname{vari}, '_interp_', CMIP6_info.modelname, '_', ...
        CMIP6_info.scenname, '_', CMIP6_info.ensname, '_', CMIP6_info.iniyear, '.nc'];
    switch (CMIP6_info.varname{vari})
        case {'thetao', 'so', 'zos'}
            ROMS_info.horgrid{vari}='r';
        case {'uo'}
            ROMS_info.horgrid{vari}='u';
        case {'vo'}
            ROMS_info.horgrid{vari}='v';
    end
end


title=['ROMS initial file from CMIP6 ', CMIP6_info.modelname, ' 1 month averaged result'];

%
%
%%%%%%%%%%%%%%%%%%% END USERS DEFINED VARIABLES %%%%%%%%%%%%%%%%%%%%%%%
%
% Title
%
disp(' ')
disp([' Making initial file: ',ininame])
disp(' ')
disp([' Title: ',title])
%
% Initial file
%
% create_inifile(ininame,grdname,title,...
%                theta_s,theta_b,hc,N,...
%                tini,'clobber');

create_inifile_Y(ininame,grdname,title,...
               Vtransform, Vstretching, theta_s,theta_b,hc,N,...
               tini,'clobber');
           
%
% Horizontal and vertical interp/extrapolations 
%
disp(' ')
disp(' Interpolations / extrapolations')


for vari=1:length(CMIP6_info.varname)
    disp(' ')
    disp(CMIP6_info.varname{vari})
    ext_datas_ini_CMIP6_interp(ininame,grdname,CMIP6_info.filename{vari}, ...
                CMIP6_info.varname{vari},ROMS_info.varname{vari}, ...
                ROMS_info.horgrid{vari},tini, Vtransform, Vstretching, CMIP6_info.modelname);  
end
            

u=ncread('D:\Data\Model\ROMS\nwp_1_20\input\test2199\roms_nwp_ini_test2199.nc', 'u');
pcolor(u(:,:,40)); shading flat; colorbar;            
%
% Make a few plots
%
disp(' ')
disp(' Make a few plots...')
test_clim(ininame,grdname,'temp',1,coastfileplot)
figure
test_clim(ininame,grdname,'salt',1,coastfileplot)
%
% End
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
