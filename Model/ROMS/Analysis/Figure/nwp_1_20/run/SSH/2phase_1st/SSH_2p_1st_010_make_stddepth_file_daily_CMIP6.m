addpath(genpath('/data1/ykang/function/'))

clc; clear all; close all;
modelname = 'CMIP6';
testnames = {'test2117','test2118','test2119','test2120','test2121'};
syear = 1985; ss = num2str(syear);

for i = 1 : length(testnames)
    testname = testnames{i};
    
%     switch testname
%         case {'test2101','test2102'}
%             filedir = 'output_hdd1';
%         case {'test2103','test2104'}
%             filedir = 'output_hdd2';
%         case {'test2105','test2106'}
%             filedir = 'output_hdd3';
%     end
    
    fpath = ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_daily/', testname, '/'];
    flist = dir(fpath);
    ee = flist(end).name; eyear = str2num(ee);
    
%     stdpath = '/data1/ykang/data/';
%     stdname = 'stddepth.mat';
%     load([stdpath, stdname])
    stddepth=[0, -10,-20,-30, -40, -50, ...
        -60, -70, -80, -90, -100, ...
        -110, -120, -130, -140, -150, ...
        -160, -170, -180, -190, -200, ...
        -220, -240, -260, -280, -300, ...
        -350, -400,-500, ...
        -700,-1000,-1250,-1500,...
        -1750,-2000,-2250,-2500,-3000,-3500];
    for year = syear : eyear
        yy = num2str(year);
        
        filepath = [fpath,yy,'/'];
        
        inputday = [1:91, 334:eomday(year,2)+337]; 

        for dayij = 1 : length(inputday)
%             mm = num2char(month,2);
            tempday=inputday(dayij);
            daystr=num2str(tempday, '%04i');
            disp([testname,' : ',yy,'.',daystr])
            
            filename = ['pck_ES_',testname,'_daily_',yy,'_',daystr,'.nc'];
            nc = netcdf([filepath, filename]);
            h = nc{'h'}(:);
            hc = nc{'hc'}(:);
            Vstretching = nc{'Vstretching'}(:);
            Vtransform = nc{'Vtransform'}(:);
            theta_s = nc{'theta_s'}(:);
            theta_b = nc{'theta_b'}(:);
            N = length(nc{'Cs_r'}(:));
            mask_rho = nc{'mask_rho'}(:);
            lon_rho = nc{'lon_rho'}(:);
            lat_rho = nc{'lat_rho'}(:);
            angle = nc{'angle'}(:);
            
            RCM_grid.lon_rho=lon_rho';
            RCM_grid.lat_rho=lat_rho';
            RCM_grid.depth=-stddepth;
            RCM_grid.size_depth = length(RCM_grid.depth);
            RCM_grid.size_lon_rho = size(RCM_grid.lon_rho, 1);
            RCM_grid.size_lat_rho = size(RCM_grid.lat_rho, 2);
            
            domaxis = [min(min(lon_rho)) max(max(lon_rho)) min(min(lat_rho)) max(max(lat_rho))];
            
            temp = nc{'temp'}(:); salt = nc{'salt'}(:); zeta = nc{'zeta'}(:);
            u = nc{'u'}(:); v = nc{'v'}(:); w = nc{'w'}(:);
            uwind=nc{'Uwind'}(:); vwind=nc{'Vwind'}(:);
            
            temp(temp>=10e10)=NaN;
            salt(salt>=10e10)=NaN;
            zeta(zeta>=10e10)=NaN;
            u(u>=10e10)=NaN;
            v(v>=10e10)=NaN;

            close(nc);
            [N,L,M] = size(temp);
            
            u = u2rho_3d(u); v = v2rho_3d(v);
            z = zlevs(h,zeta,theta_s,theta_b,hc,N,Vstretching,Vtransform,'temp');
            z(end+1,:,:)=10e5;
            zw = zlevs(h,zeta,theta_s,theta_b,hc,N,Vstretching,Vtransform,'w');
            zw(end+1,:,:)=10e5;
            
            for lonij=1:size(RCM_grid.lon_rho,1)
                for latij=1:size(RCM_grid.lat_rho,2)
                    if(isfinite(temp(end,latij,lonij)))
                        re_temp(lonij,latij,1:RCM_grid.size_depth)=interp1(z(:,latij, lonij), ...
                            [temp(:,latij,lonij); temp(end,latij,lonij)], stddepth);
                        re_salt(lonij,latij,1:RCM_grid.size_depth)=interp1(z(:,latij, lonij), ...
                            [salt(:,latij,lonij); salt(end, latij, lonij)], stddepth);
                        re_u(lonij,latij,1:RCM_grid.size_depth)=interp1(z(:,latij, lonij), ...
                            [u(:,latij,lonij); u(end, latij, lonij)], stddepth);
                        re_v(lonij,latij,1:RCM_grid.size_depth)=interp1(z(:,latij, lonij), ...
                            [v(:,latij,lonij); v(end, latij, lonij)], stddepth);
                        re_w(lonij,latij,1:RCM_grid.size_depth)=interp1(zw(:,latij, lonij), ...
                            [w(:,latij,lonij); w(end, latij, lonij)], stddepth);
                    else
                        re_temp(lonij,latij,1:RCM_grid.size_depth)=NaN;
                        re_salt(lonij,latij,1:RCM_grid.size_depth)=NaN;
                        re_u(lonij,latij,1:RCM_grid.size_depth)=NaN;
                        re_v(lonij,latij,1:RCM_grid.size_depth)=NaN;
                        re_w(lonij,latij,1:RCM_grid.size_depth)=NaN;
                    end
                end
%                 lonij
            end
            
%             pcolor(squeeze(re_temp(:,:,1))'); shading flat; colorbar;

            
%             tt = interp3(temp,lon_rho(1,:),lat_rho(:,1),stddepth);
%             ss = interp3(salt,lon_rho(1,:),lat_rho(:,1),stddepth);
%             uu = interp3(u,lon_rho(1,:),lat_rho(:,1),stddepth);
%             vv = interp3(v,lon_rho(1,:),lat_rho(:,1),stddepth);
%             ww = interp3(w,lon_rho(1,:),lat_rho(:,1),stddepth);
%             
%             re_temp = permute(tt,[2 1 3]); re_salt = permute(ss,[2 1 3]);
%             re_u = permute(uu,[2 1 3]); re_v = permute(vv,[2 1 3]);
%             re_w = permute(ww,[2 1 3]);
            
            
%             savepath = ['/data1/ykang/CMIP6/output/',testname,'/stddepth/monthly/',yy,'/'];
            savepath = ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_stddepth_mat_daily/', testname, '/'];

            if exist(savepath) == 0
                mkdir(savepath)
            end
            
            savename = [savepath,testname,'_monthly_std_',yy,'_',daystr,'.mat'];
            temp = re_temp; salt = re_salt;
            u = re_u; v = re_v;  w = re_w; depth = stddepth;
            
            
            save(savename,'temp','salt','u','v','w','zeta','uwind', 'vwind', 'RCM_grid');
            
%             RCM_time.ftime= datenum(year,month,15) - datenum(1900,12,31);
            RCM_time.ftime= datenum(year,1,1) +(tempday-1) - datenum(1900,12,31);

            % % % % % %         make ncfile
            ncsavepath = ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_stddepth_daily/', testname, '/', yy, '/'];

            savencname = [ncsavepath,testname,'_daily_std_',yy,'_',daystr,'.nc'];
            if exist(ncsavepath) == 0
                mkdir(ncsavepath)
            end
            
            

            ncid = netcdf.create(savencname,'NETCDF4');

            lon_dimid = netcdf.defDim(ncid, 'xi_rho', RCM_grid.size_lon_rho);
            lat_dimid = netcdf.defDim(ncid,'eta_rho',RCM_grid.size_lat_rho);
            depth_dimid = netcdf.defDim(ncid, 'depth', RCM_grid.size_depth);
            time_dimid = netcdf.defDim(ncid, 'time', 0);

            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'type', ['NWP 1/20 _ ', testname, 'model, ES stddepth interpolated file']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'title', [' daily SSH analysis (', num2str(min(year)), '-', num2str(max(year)) ,') ']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'source', [' ROMS NWP 1/20 data from _ ',testname ]);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'author', 'Created by Y.Y.Kim');
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'date', date);

            timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
            netcdf.putAtt(ncid,timevarid,'long_name','time');
            netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
            netcdf.putAtt(ncid,timevarid,'calendar','gregorian');           

            lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
            netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
            netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

            lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
            netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
            netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
            
            tempvarid=netcdf.defVar(ncid, 'temp', 'NC_FLOAT', [lon_dimid lat_dimid depth_dimid time_dimid]);
            netcdf.putAtt(ncid,tempvarid,'long_name','temp');
            netcdf.putAtt(ncid,tempvarid,'units','Celsius degree');  
            
            saltvarid=netcdf.defVar(ncid, 'salt', 'NC_FLOAT', [lon_dimid lat_dimid depth_dimid time_dimid]);
            netcdf.putAtt(ncid,saltvarid,'long_name','salt');
            netcdf.putAtt(ncid,saltvarid,'units',' ');
            
            uvarid=netcdf.defVar(ncid, 'u', 'NC_FLOAT', [lon_dimid lat_dimid depth_dimid time_dimid]);
            netcdf.putAtt(ncid,uvarid,'long_name','u');
            netcdf.putAtt(ncid,uvarid,'units','m/s');
            
            vvarid=netcdf.defVar(ncid, 'v', 'NC_FLOAT', [lon_dimid lat_dimid depth_dimid time_dimid]);
            netcdf.putAtt(ncid,vvarid,'long_name','v');
            netcdf.putAtt(ncid,vvarid,'units','m/s');
            
            wvarid=netcdf.defVar(ncid, 'w', 'NC_FLOAT', [lon_dimid lat_dimid depth_dimid time_dimid]);
            netcdf.putAtt(ncid,wvarid,'long_name','w');
            netcdf.putAtt(ncid,wvarid,'units','m/s');
            
            zetavarid=netcdf.defVar(ncid, 'zeta', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
            netcdf.putAtt(ncid,zetavarid,'long_name','zeta');
            netcdf.putAtt(ncid,zetavarid,'units','m');
            
            uwindvarid=netcdf.defVar(ncid, 'uwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
            netcdf.putAtt(ncid,uwindvarid,'long_name','Uwind');
            netcdf.putAtt(ncid,uwindvarid,'units','m/s');
            
            vwindvarid=netcdf.defVar(ncid, 'vwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
            netcdf.putAtt(ncid,vwindvarid,'long_name','Vwind');
            netcdf.putAtt(ncid,vwindvarid,'units','m/s');
            
            depthvarid=netcdf.defVar(ncid, 'depth', 'NC_FLOAT', [depth_dimid]);
            netcdf.putAtt(ncid,depthvarid,'long_name','depth');
            netcdf.putAtt(ncid,depthvarid,'units','m');

            netcdf.endDef(ncid);

            netcdf.putVar(ncid, timevarid, 0, length(RCM_time.ftime), RCM_time.ftime);
            netcdf.putVar(ncid, lon_rhovarid, [0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho], RCM_grid.lon_rho);
            netcdf.putVar(ncid, lat_rhovarid, [0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho], RCM_grid.lat_rho);
            netcdf.putVar(ncid, tempvarid, [0 0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho RCM_grid.size_depth length(RCM_time.ftime)], temp);
            netcdf.putVar(ncid, saltvarid, [0 0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho RCM_grid.size_depth length(RCM_time.ftime)], salt);
            netcdf.putVar(ncid, uvarid, [0 0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho RCM_grid.size_depth length(RCM_time.ftime)], u);
            netcdf.putVar(ncid, vvarid, [0 0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho RCM_grid.size_depth length(RCM_time.ftime)], v);
            netcdf.putVar(ncid, wvarid, [0 0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho RCM_grid.size_depth length(RCM_time.ftime)], w);
            netcdf.putVar(ncid, depthvarid, [0], [RCM_grid.size_depth], RCM_grid.depth);
            netcdf.putVar(ncid, zetavarid, [0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho length(RCM_time.ftime)], zeta');
            netcdf.putVar(ncid, uwindvarid, [0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho length(RCM_time.ftime)], uwind');
            netcdf.putVar(ncid, vwindvarid, [0 0 0], [RCM_grid.size_lon_rho RCM_grid.size_lat_rho length(RCM_time.ftime)], vwind');

            netcdf.close(ncid);

            
            
        end % month
    end % year
end % testname
