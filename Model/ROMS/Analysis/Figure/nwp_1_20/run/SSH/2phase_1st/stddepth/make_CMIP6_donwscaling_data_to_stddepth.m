%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% program  make_CMIP6_dwonscaling_data_to_stddepth
%
% converting ROMS output monthly data (sigma coordinate) to z-coordinate data
%
%  input:
%  data                   monthly data (temp, salt, zeta, u, v, w [lon, lat, s, t])
%
%  output:
%  data                   monthly data (temp, salt, zeta, u, v, w [lon, lat, s, t])
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    25-Jun-2021 by Yong-Yub Kim (now editing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all;

all_RCM.number={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
all_RCM.name ={'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC', 'RCM-CNRM-HR', 'RCM-CMCC'};
all_RCM.model = 'nwp_1_20';
all_RCM.root = '/data1/RCM/CMIP6/';
all_RCM.progress = 'run';

years = [1985 : 2014];

flags.make_std_nc = 0;
flags.make_std_mat = 1;


for ind_testname=1:length(all_RCM.name)
    testname = all_RCM.number{ind_testname};
    
% %     temppath=
    
    
    stddepth = [-10,-20,-30,-50,-75,-100, ...
        -150,-200,-300,-400,-500, ...
        -700,-1000,-1250,-1500, ...
        -1750,-2000,-2250,-2500, ...
        -3000,-3500,-4000,-4500,-5000];
    
end

for test = 1
    Vstretching = 4; Vtransform = 2;
    
    gridpath = ['/data2/ykang/',modelname,'/input/'];
    gridname = ['ES_grid.nc'];
    nc = netcdf([gridpath, gridname]);
    h = nc{'h'}(:);
    mask_rho = nc{'mask_rho'}(:);
    lon_rho = nc{'lon_rho'}(:);
    lat_rho = nc{'lat_rho'}(:);
    angle = nc{'angle'}(:);
    domaxis = [min(min(lon_rho)) max(max(lon_rho)) min(min(lat_rho)) max(max(lat_rho))];
    close(nc);
    
    gpath = gridpath;
    gname = ['ES_grid_info.nc'];
    nw = netcdf([gpath, gname]);
    hc = nc{'hc'}(:);
    theta_s = nc{'theta_s'}(:);
    theta_b = nc{'theta_b'}(:);
    close(nw);
    
    for year = syear : eyear
        yy = num2str(year);
        
        for month = 1 : 12
            mm = num2char(month,2);
            disp([testname])
            disp([yy,'.',mm])
            tstart = tic;
            
            filepath = ['/data2/ykang/',modelname,'/output/',testname,'/monthly/',yy,'/'];
            savepath = [filepath, 'stddepth/'];
            if exist(savepath) == 0
                mkdir(savepath)
            end
            
            filename = [filepath, 'ES_',testname,'_monthly_',yy,'_',mm,'.nc'];
            nc = netcdf([filename]);
            temp = nc{'temp'}(:); salt = nc{'salt'}(:); zeta = nc{'zeta'}(:);
            w = nc{'w'}(:);
            close(nc);
            [N,L,M] = size(temp);
            
            uname = [filepath, 'ES_',testname,'_monthly_u_',yy,'_',mm,'.nc'];
            nc = netcdf([uname]);
            u = nc{'u'}(:); close(nc);
            
            vname = [filepath, 'ES_',testname,'_monthly_v_',yy,'_',mm,'.nc'];
            nc = netcdf([vname]);
            v = nc{'v'}(:); close(nc);
            
            u = u2rho_3d(u); v = v2rho_3d(v);
            z = zlevs(h,zeta,theta_s,theta_b,hc,N,Vstretching,Vtransform,'temp');
            zw = zlevs(h,zeta,theta_s,theta_b,hc,N,Vstretching,Vtransform,'w');
            
            for i = 1 : L
                if mod(i,100) == 0
                    disp(['i = ',num2str(i),'/',num2str(L)])
                end
                for j = 1 : M
                    org_depth = z(:,i,j);
                    org_wdepth = zw(:,i,j);
                    
                    org_temp = temp(:,i,j); org_salt = salt(:,i,j);
                    org_u = u(:,i,j); org_v = v(:,i,j); org_w = w(:,i,j);
                    %org_rho = rho(:,i,j);
                    
                    re_temp(:,i,j) = interp1(org_depth, org_temp, stddepth);
                    re_salt(:,i,j) = interp1(org_depth, org_salt, stddepth);
                    re_u(:,i,j) = interp1(org_depth, org_u, stddepth);
                    re_v(:,i,j) = interp1(org_depth, org_v, stddepth);
                    re_w(:,i,j) = interp1(org_wdepth, org_w, stddepth);
                    %re_rho(:,i,j) = interp1(org_depth, org_rho, stddepth);
                end
            end
            
            if make_std_nc == 1
                savename = [savepath,testname,'_monthly_std_',yy,'_',mm,'.nc'];
                make_stddepth_ncfile('temp',stddepth,savename,filename);
                nw = netcdf([savename],'w');
                nw{'temp'}(:) = re_temp; nw{'salt'}(:) = re_salt;
                nw{'u'}(:) = re_u; nw{'v'}(:) = re_v;  nw{'w'}(:) = re_w;
                nw{'depth'}(:) = stddepth;
                close(nw);
            end
            
            if make_std_mat == 1
                savename = [savepath,testname,'_monthly_std_',yy,'_',mm,'.mat'];
                temp = re_temp; salt = re_salt; %rho = re_rho;
                u = re_u; v = re_v;  w = re_w; depth = stddepth;
                
                save(savename,'temp','salt','u','v','w','depth','zeta');
            end
            tend = toc(tstart);
            disp(['takes ',num2str(round(tend,2)),'sec for a month'])
        end
    end
end
