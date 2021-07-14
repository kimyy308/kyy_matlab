clc; clear all; close all;
addpath(genpath('/home/kisti_backup/'));
addpath(genpath('/home/ykang/function/'));
modelname = 'CMIP6';
headpath = ['/home/ykang/',modelname,'/'];
tt = [2106];

for test = 1 : length(tt)
    testname = ['test',num2str(tt(test))];
    progress = 'run';
    
    make_std_nc = 1;
    make_std_mat = 0;
    
    stdpath = headpath;
    stdname = 'stddepth.mat';
    load([stdpath,stdname])
    
    Vstretching = 4; Vtransform = 2;
    
    gridpath = [headpath,'input/'];
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
        
        filepath = [headpath,'output/',testname,'/daily/',yy,'/'];
        savepath = [filepath, 'stddepth/'];
        if exist(savepath) == 0
            mkdir(savepath)
        end
        
        if leapyear(year) == 1
            days = 366;
        else days = 365;
        end
        
        % days=365;
        
        for day = 1 : days
            tstart2 = tic;
            
            dd = num2char(day,4);
            disp(testname)
            disp([yy,'.',dd])
            
            filename = [filepath, 'ES_',testname,'_daily_',yy,'_',dd,'.nc'];
            nc = netcdf(filename);
            temp = nc{'temp'}(:); salt = nc{'salt'}(:);
            zeta = nc{'zeta'}(:);  w = nc{'w'}(:);
            close(nc);
            [N,L,M] = size(temp);
            
            uname = [filepath, 'ES_',testname,'_daily_u_',yy,'_',dd,'.nc'];
            nc = netcdf(uname);
            u = nc{'u'}(:); close(nc);
            
            vname = [filepath, 'ES_',testname,'_daily_v_',yy,'_',dd,'.nc'];
            nc = netcdf(vname);
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
                    
                    re_temp(:,i,j) = interp1(org_depth, org_temp, stddepth);
                    re_salt(:,i,j) = interp1(org_depth, org_salt, stddepth);
                    re_u(:,i,j) = interp1(org_depth, org_u, stddepth);
                    re_v(:,i,j) = interp1(org_depth, org_v, stddepth);
                    re_w(:,i,j) = interp1(org_wdepth, org_w, stddepth);
                end
            end
            
            if make_std_nc == 1
                savename = [savepath,testname,'_daily_std_',yy,'_',dd,'.nc'];
                make_stddepth_ncfile('temp',stddepth,savename,filename);
                nw = netcdf([savename],'w');
                nw{'temp'}(:) = re_temp; nw{'salt'}(:) = re_salt;  nw{'zeta'}(:) = zeta;
                nw{'u'}(:) = re_u; nw{'v'}(:) = re_v;  nw{'w'}(:) = re_w;
                nw{'depth'}(:) = stddepth;
                close(nw);
            end
            
            if make_std_mat == 1
                savename = [savepath,testname,'_daily_std_',yy,'_',dd,'.mat'];
                temp = re_temp; salt = re_salt;
                u = re_u; v = re_v;  w = re_w; depth = stddepth;
                save(savename,'temp','salt','u','v','w','depth','zeta');
            end
         tend2 = toc(tstart2);   
         disp(['takes ',num2str(round(tend2,3)),'sec for a year'])
        end
    end
end
