close all; clear all;  clc;

all_testname = {'test2117'};

% all_region ={'pollock_egg', 'pollock_egg2', 'ES'};
all_region ={'ES_KHOA'};
% all_region ={'ES_KHOA','YS_KHOA', 'SS_KHOA'};


for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname
        % % % 
        dropboxpath='/home/kimyy/Dropbox/';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/mat_tool']));
        addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path

        shadlev = [-2 2];
        % rms_shadlev = [0 4];
        trend_shadlev = [0 4];
        % bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/10;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1995:2014]; % % put year which you want to plot [year year ...]
%         inputyear = [2081:2100]; % % put year which you want to plot [year year ...]

%         inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        RCM_info.season='JF-';
        [inputmonth,error_status] = Func_0019_get_month_from_season(RCM_info.season);
%         inputmonth = [1,2]; % % put month which you want to plot [month month ...]
        
        varname ='zeta';
        variable='zeta';
        tmp.regionname=all_region{regionind};
        run('nwp_polygon_point.m');
        tmp.dropboxpath = '/home/kimyy/Dropbox';
        addpath(genpath([tmp.dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'function']));
        [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, filesep, 'source', filesep, 'matlab', filesep, 'Model' ...
            filesep, 'ROMS', filesep, 'Analysis', filesep, 'Figure', filesep, 'nwp_1_20', filesep ...
            'run', filesep, 'SSH', filesep, '2phase_2nd', filesep, 'subroutine']));
% % %         switch region
        [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);

        
%         switch testname
%             case {'test06'}
%                 drivename='D:/';
%         end
        param_script ='/home/kimyy/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_20/run/fig_param/fig_param2_kyy_ES_KHOA.m';
        filedir = strcat('/home/kimyy/Model/ROMS/nwp_1_20/output/', testname, '/run/'); % % where data files are
        ltransdir = ['/home/kimyy/Model/LTRANS/LTRANSv2b_auto/output/', testname, filesep];
        savedir = strcat('/home/kimyy/Model/LTRANS/LTRANSv2b_auto/output/', testname, filesep, 'analysis', filesep, testname, filesep);
        inputdir = ['/home/kimyy/Model/ROMS/nwp_1_20/output/', testname, '/run/' ];
        mkdir(savedir);
        
        
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='get pollock hatchery data';

            for flagi=1:10
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=2;
        end
        
% % %         get model data      
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            run(param_script);
            for yearij = 1:length(inputyear)
                tempyear = inputyear(yearij);
                tempyearstr = num2str(tempyear, '%04i');
%                 tairname = [inputdir, 'SBC/old/nwp_1_10_', num2str(tempyear, '%04i'), '_', 'Tair.nc'];
                for monthij = 1:length(inputmonth)
                    tempmonth = inputmonth(monthij);
                    ncname = [savedir,testname,'_',tmp.regionname,'_model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
%                     if (exist(ncname , 'file') ~= 2 || fig_flag==2)   
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        ind=1;
                        if tempmonth==12
                            [status, lastday_m_str]=system(['date -d "', num2str(tempyear+1), num2str(1, '%02i'), '01 1 days ago" +%d']); 
                        else
                            [status, lastday_m_str]=system(['date -d "', num2str(tempyear), num2str(tempmonth+1, '%02i'), '01 1 days ago" +%d']); 
                        end
                        lastday_m=str2num(lastday_m_str);
                        for dayij = 1:lastday_m
                            [status, filenum_str]=system(['date -d "', ...
                                num2str(tempyear), num2str(tempmonth, '%02i'), num2str(dayij, '%02i'), ' 0 days ago" +%j']);
                            filenum=str2num(filenum_str);
                            
                            filename = strcat(filedir, tempyearstr, '/', ...
                                    'pck_', testname, '_', tempyearstr, '_daily_avg_', num2str(filenum,'%04i'), '.nc');
                            ltransfilename = [ltransdir, num2str(tempyear, '%04i'), '/', '0030d_', ...
                                num2str(tempyear, '%04i'), num2str(tempmonth, '%02i'), num2str(dayij, '%02i'), '/', ...
                                'output_0030d_', num2str(tempyear, '%04i'), num2str(tempmonth, '%02i'), num2str(dayij, '%02i'), '.nc'];
                            
%                             ncinfo(ltransfilename)
                            
                            
                            
                            % read model data
                            if (exist('lon')==0)
                                modelinfo=ncinfo(filename);
%                                 lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
%                                 lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
                                
                                lon = ncread(filename, 'lon_rho');
                                lat = ncread(filename, 'lat_rho');
                                mask_rho2 = ncread(filename, 'mask_rho');

                                switch(tmp.regionname)
                                    case('NWP')
                                        mask_model2(1:size(lon,1),1:size(lon,2))=1;
                                    otherwise
                                        mask_model2 = double(inpolygon(lon,lat,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                                        mask_model2(mask_model2==0)=NaN;
                                end
                                mask_rho2 = mask_rho2 .* mask_model2;
                                size_lon = size(mask_rho2, 1);
                                size_lat = size(mask_rho2, 2);
% % %                                 for es grid (ykang)
%                                 lon_min=min(mod(find(mask_rho2(:,:)==1),size_lon));
%                                 lon_max=max(mod(find(mask_rho2(:,:)==1),size_lon));
%                                 lat_min=min(mod(find(mask_rho2(:,:)'==1),size_lat));
%                                 lat_max=max(mod(find(mask_rho2(:,:)'==1),size_lat));
                                
                                lon_min=min(mod(find(mask_model2(:,:)==1),size_lon));
                                lon_max=max(mod(find(mask_model2(:,:)==1),size_lon));
                                lat_min=min(mod(find(mask_model2(:,:)'==1),size_lat));
                                lat_max=max(mod(find(mask_model2(:,:)'==1),size_lat));
                                
                                
                                loncount= lon_max(1)-lon_min(1)+1;
                                latcount = lat_max(1)-lat_min(1)+1;
                                
                                lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [loncount latcount]);
                                lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [loncount latcount]);

                                switch(tmp.regionname)
                                    case('NWP') %% North western Pacific
                                        mask_model(1:size(lon,1),1:size(lon,2))=1;
                                    otherwise
                                        mask_model = double(inpolygon(lon,lat,RCM_grid.refpolygon(:,1),RCM_grid.refpolygon(:,2)));
                                        mask_model(mask_model==0)=NaN;
%                                         ref2polygon = [127 36; 127 41; 133 41; 133 36]; %%%% just for the EKB analysis
                                        ref2polygon = [0 -90; 0 -90; 360 90; 360 90]; %%%% just for analysis
                                        egg_mask_model = double(inpolygon(lon,lat,ref2polygon(:,1),ref2polygon(:,2)));
                                        egg_mask_model(egg_mask_model==0)=NaN;
                                end
                                max_slev=length(ncread(filename, 's_rho'));
                                mask_rho=ncread(filename, 'mask_rho', [lon_min(1) lat_min(1)], [loncount latcount]);
                                h = ncread(filename, 'h', [lon_min(1) lat_min(1)], [loncount latcount]);
                                pm = ncread(filename, 'pm', [lon_min(1) lat_min(1)], [loncount latcount]);
                                pn = ncread(filename, 'pn', [lon_min(1) lat_min(1)], [loncount latcount]);
                            end
%                             data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
                            
                            time_5=48 + 24* 5 +1;
                            time_10=48 + 24 * 10 +1;
                            time_15=48 + 24 * 15 +1;
                            time_30=48 + 24 * 30 +1;
                            
                            mask_5day=zeros(loncount, latcount);
                            mask_10day=zeros(loncount, latcount);
                            mask_15day=zeros(loncount, latcount);
                            mask_30day=zeros(loncount, latcount);
                            xdist = 1./pm;
                            ydist = 1./pn; 
                            
                            dob=ncread(ltransfilename, 'dob');

                            if isempty(dob)
%                                 lon_par_5=NaN;
%                                 lat_par_5=NaN;
%                                 lon_par_10=NaN;
%                                 lat_par_10=NaN;
                                lon_par_15=NaN;
                                lat_par_15=NaN;
                                lon_par_30=NaN;
                                lat_par_30=NaN;
                            else
%                                 lon_par_5=ncread(ltransfilename, 'lon', [1, time_5], [inf 1]);
%                                 lat_par_5=ncread(ltransfilename, 'lat', [1, time_5], [inf 1]);
%                                 lon_par_10=ncread(ltransfilename, 'lon', [1, time_10], [inf 1]);
%                                 lat_par_10=ncread(ltransfilename, 'lat', [1, time_10], [inf 1]);
                                lon_par_15=ncread(ltransfilename, 'lon', [1, time_15], [inf 1]);
                                lat_par_15=ncread(ltransfilename, 'lat', [1, time_15], [inf 1]);
                                lon_par_30=ncread(ltransfilename, 'lon', [1, time_30], [inf 1]);
                                lat_par_30=ncread(ltransfilename, 'lat', [1, time_30], [inf 1]);
                                
                                for pari=1:length(lon_par_15)

%                                     if (lon_par_5(pari)> RCM_grid.domain(1) & lon_par_5(pari) < RCM_grid.domain(2) & lat_par_5(pari) > RCM_grid.domain(3) & lat_par_5(pari) < RCM_grid.domain(4))
%                                         for loni= 1:loncount
%                                             dist_lon5(loni)=abs(lon(loni,1)-lon_par_5(pari));
%                                         end
%                                         ind_lon5(pari)=find(dist_lon5(:)==min(dist_lon5(:)));
%                                         for lati = 1:latcount
%                                             dist_lat5(lati)=abs(lat(1,lati)-lat_par_5(pari));
%                                         end
%                                         ind_lat5(pari)=find(dist_lat5(:)==min(dist_lat5(:)));
%                                         mask_5day(ind_lon5(pari), ind_lat5(pari))=mask_5day(ind_lon5(pari), ind_lat5(pari))+1;
%                                     end
% 
%                                     if (lon_par_10(pari)> RCM_grid.domain(1) & lon_par_10(pari) < RCM_grid.domain(2) & lat_par_10(pari) > RCM_grid.domain(3) & lat_par_10(pari) < RCM_grid.domain(4))
%                                         for loni= 1:loncount
%                                             dist_lon10(loni)=abs(lon(loni,1)-lon_par_10(pari));
%                                         end
%                                         ind_lon10(pari)=find(dist_lon10(:)==min(dist_lon10(:)));
%                                         for lati = 1:latcount
%                                             dist_lat10(lati)=abs(lat(1,lati)-lat_par_10(pari));
%                                         end
%                                         ind_lat10(pari)=find(dist_lat10(:)==min(dist_lat10(:)));
%                                         mask_10day(ind_lon10(pari), ind_lat10(pari))=mask_10day(ind_lon10(pari), ind_lat10(pari))+1;
%                                     end

                                    if (lon_par_15(pari)> RCM_grid.domain(1) & lon_par_15(pari) < RCM_grid.domain(2) & lat_par_15(pari) > RCM_grid.domain(3) & lat_par_15(pari) < RCM_grid.domain(4))
                                        for loni= 1:loncount
                                            dist_lon15(loni)=abs(lon(loni,1)-lon_par_15(pari));
                                        end
                                        ind_lon15(pari)=find(dist_lon15(:)==min(dist_lon15(:)));
                                        for lati = 1:latcount
                                            dist_lat15(lati)=abs(lat(1,lati)-lat_par_15(pari));
                                        end
                                        ind_lat15(pari)=find(dist_lat15(:)==min(dist_lat15(:)));
                                        mask_15day(ind_lon15(pari), ind_lat15(pari))=mask_15day(ind_lon15(pari), ind_lat15(pari))+1;
                                    end

                                    if (lon_par_30(pari)> RCM_grid.domain(1) & lon_par_30(pari) < RCM_grid.domain(2) & lat_par_30(pari) > RCM_grid.domain(3) & lat_par_30(pari) < RCM_grid.domain(4))
                                        for loni= 1:loncount
                                            dist_lon30(loni)=abs(lon(loni,1)-lon_par_30(pari));
                                        end
                                        ind_lon30(pari)=find(dist_lon30(:)==min(dist_lon30(:)));
                                        for lati = 1:latcount
                                            dist_lat30(lati)=abs(lat(1,lati)-lat_par_30(pari));
                                        end
                                        ind_lat30(pari)=find(dist_lat30(:)==min(dist_lat30(:)));
                                        mask_30day(ind_lon30(pari), ind_lat30(pari))=mask_30day(ind_lon30(pari), ind_lat30(pari))+1;
                                    end
                                end
                            end

                            temp = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [loncount latcount max_slev 1]);
                            
                            lonlat_ind=[lon_min(1), lon_max(1), lat_min(1), lat_max(1)];
                            temp_surf= temp(:,:,end);
                            temp_50 = get_hslice_pollock(filename, filename, temp, -50, 'temp', lonlat_ind)';
                            
                            egg_mask=zeros(size(mask_rho)).*egg_mask_model;
                            idx_egg = find(-h >= -500 & temp_50 >= 2 & temp_50 <= 5);
                            egg_mask(idx_egg)=1;

                            uwind = ncread(filename,'Uwind',[lon_min(1) lat_min(1) 1], [loncount latcount  1]);
                            vwind = ncread(filename,'Vwind',[lon_min(1) lat_min(1) 1], [loncount latcount  1]);
                            shflux = ncread(filename, 'shflux', [lon_min(1) lat_min(1) 1], [loncount latcount 1]);
                            u_surf = ncread(filename,'u',[lon_min(1) lat_min(1) max_slev 1], [loncount-1 latcount  1 1]);
                            v_surf = ncread(filename,'v',[lon_min(1) lat_min(1) max_slev 1], [loncount latcount-1  1 1]);
%                             tair = ncread(tairname, 'Tair', [lon_min(1) lat_min(1) filenum], [loncount latcount 1]);
                            pair = ncread(filename, 'Pair', [lon_min(1) lat_min(1) 1], [loncount latcount 1]);
                            zeta = ncread(filename, 'zeta', [lon_min(1) lat_min(1) 1], [loncount latcount 1]);
                            
                            [uwstr, vwstr]=ra_windstr(uwind,vwind);
                            
                            u_rho=NaN(loncount,latcount);
                            u_rho(2:loncount-1,:) = 0.5 *(u_surf(1:loncount-2,:) + u_surf(2:loncount-1,:));
                            u_rho(1,:)=u_surf(1,:);
                            u_rho(loncount,:)=u_surf(loncount-1,:);
                            v_rho=NaN(loncount,latcount);
                            v_rho(:,2:latcount-1) = 0.5 *(v_surf(:,1:latcount-2) + v_surf(:,2:latcount-1));
                            v_rho(:,1)=v_surf(:,1);
                            v_rho(:,latcount)=v_surf(:,latcount-1);
                            
                            for i=1:loncount-2
                                for j=1:latcount-2
                                    vort(i+1,j+1) = (v_rho(i+2,j+1)-v_rho(i,j+1)) ...
                                        / (xdist(i+1, j+1)*2) + ...
                                        (u_rho(i+1,j+2)-u_rho(i+1,j)) ...
                                        / (xdist(i+1, j+1)*2);
                                    wind_curl(i+1,j+1) = (vwind(i+2,j+1)-vwind(i,j+1)) ...
                                        / (xdist(i+1, j+1)*2) + ...
                                        (uwind(i+1,j+2)-uwind(i+1,j)) ...
                                        / (xdist(i+1, j+1)*2);
                                    wstr_curl(i+1,j+1) = (vwstr(i+2,j+1)-vwstr(i,j+1)) ...
                                        / (xdist(i+1, j+1)*2) + ...
                                        (uwstr(i+1,j+2)-uwstr(i+1,j)) ...
                                        / (xdist(i+1, j+1)*2);
                                end
                            end
                            vort(1,2:latcount-2)=vort(2,2:latcount-2);
                            vort(loncount, 2:latcount-2) =vort(loncount-1,2:latcount-2);
                            vort(2:loncount-2,1)=vort(2:loncount-2,2);
                            vort(2:loncount-2,latcount)=vort(2:loncount-2,latcount-1);
                            vort(1,1)=vort(2,2);
                            vort(1,latcount)=vort(2,latcount-1);
                            vort(loncount,1)=vort(loncount-1,2);
                            vort(loncount,latcount)=vort(loncount-1,latcount-1);
                            wind_curl(1,2:latcount-2)=wind_curl(2,2:latcount-2);
                            wind_curl(loncount, 2:latcount-2) =wind_curl(loncount-1,2:latcount-2);
                            wind_curl(2:loncount-2,1)=wind_curl(2:loncount-2,2);
                            wind_curl(2:loncount-2,latcount)=wind_curl(2:loncount-2,latcount-1);
                            wind_curl(1,1)=wind_curl(2,2);
                            wind_curl(1,latcount)=wind_curl(2,latcount-1);
                            wind_curl(loncount,1)=wind_curl(loncount-1,2);
                            wind_curl(loncount,latcount)=wind_curl(loncount-1,latcount-1);
                            wstr_curl(1,2:latcount-2)=wstr_curl(2,2:latcount-2);
                            wstr_curl(loncount, 2:latcount-2) =wstr_curl(loncount-1,2:latcount-2);
                            wstr_curl(2:loncount-2,1)=wstr_curl(2:loncount-2,2);
                            wstr_curl(2:loncount-2,latcount)=wstr_curl(2:loncount-2,latcount-1);
                            wstr_curl(1,1)=wstr_curl(2,2);
                            wstr_curl(1,latcount)=wstr_curl(2,latcount-1);
                            wstr_curl(loncount,1)=wstr_curl(loncount-1,2);
                            wstr_curl(loncount,latcount)=wstr_curl(loncount-1,latcount-1);
                               
%                             if (exist('comb_spatial_meanmodel')==0)
%                                 comb_spatial_meanmodel=(zeros([len_lon_model,len_lat_model,12]));
%                                 del=(zeros([len_lon_model,len_lat_model,12]));
%                             end

                            comb_temp_surf(:,:,ind) = single(temp_surf);
                            comb_temp_50(:,:,ind) = single(temp_50);
                            comb_uwind(:,:,ind) = single(uwind);
                            comb_vwind(:,:,ind) = single(vwind);
                            comb_uwstr(:,:,ind) = single(uwstr);
                            comb_vwstr(:,:,ind) = single(vwstr);
                            comb_wind_curl(:,:,ind) = single(wind_curl);
                            comb_wstr_curl(:,:,ind) = single(wstr_curl);
                            comb_shflux(:,:,ind) = single(shflux);
                            comb_u_rho(:,:,ind) = single(u_rho);
                            comb_v_rho(:,:,ind) = single(v_rho);
                            comb_vort(:,:,ind) = single(vort);
                            comb_egg_mask(:,:,ind) = single(egg_mask);
%                             comb_tair(:,:,ind) = single(tair);
                            comb_pair(:,:,ind) = single(pair);
                            comb_zeta(:,:,ind) = single(zeta);
                            comb_mask_5day(:,:,ind) = single(mask_5day);
                            comb_mask_10day(:,:,ind) = single(mask_10day);
                            comb_mask_15day(:,:,ind) = single(mask_15day);
                            comb_mask_30day(:,:,ind) = single(mask_30day);

                            comb_time(ind) = datenum(tempyear, tempmonth, dayij) -datenum(1900,12,31);
                            ind = ind + 1;
                        end
                        
                        % % %         make ncfile
                        ncid = netcdf.create(ncname,'NETCDF4');

                        lon_dimid = netcdf.defDim(ncid, 'lon', loncount);
                        lat_dimid = netcdf.defDim(ncid,'lat',latcount);
                        time_dimid = netcdf.defDim(ncid, 'time', 0);

                        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                            'type', ['NWP 1/20 _ ', testname, 'model, pollock data analysis file']);
                        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                            'title', [' daily pollock analysis (', num2str(tempyear), '-', num2str(tempmonth,'%02i') ,') ']);
                        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                            'source', [' ROMS NWP 1/10 DA data from _ ',testname ]);
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
                        netcdf.defVarChunking(ncid, lon_rhovarid, 'CHUNKED', [loncount/10 latcount/10]);
                        netcdf.defVarDeflate(ncid, lon_rhovarid, true, true, 1);
                        
                        lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                        netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
                        netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
%                         netcdf.defVarChunking(ncid, lat_rhovarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarChunking(ncid, lat_rhovarid, 'CHUNKED', [loncount/10 latcount/10]);
                        netcdf.defVarDeflate(ncid, lat_rhovarid, true, true, 1);

                        temp_surfvarid=netcdf.defVar(ncid, 'temp_surf', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,temp_surfvarid,'long_name','temp_surf');
                        netcdf.putAtt(ncid,temp_surfvarid,'units','Celcius degree');
                        netcdf.defVarChunking(ncid, temp_surfvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, temp_surfvarid, true, true, 1);
                        
                        temp_50varid=netcdf.defVar(ncid, 'temp_50', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,temp_50varid,'long_name','temp_50');
                        netcdf.putAtt(ncid,temp_50varid,'units','Celcius degree');
                        netcdf.defVarChunking(ncid, temp_50varid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, temp_50varid, true, true, 1);
                        
                        uwindvarid=netcdf.defVar(ncid, 'uwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,uwindvarid,'long_name','uwind');
                        netcdf.putAtt(ncid,uwindvarid,'units','m/s');
                        netcdf.defVarChunking(ncid, uwindvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, uwindvarid, true, true, 1);
                        
                        vwindvarid=netcdf.defVar(ncid, 'vwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,vwindvarid,'long_name','vwind');
                        netcdf.putAtt(ncid,vwindvarid,'units','m/s');
                        netcdf.defVarChunking(ncid, vwindvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, vwindvarid, true, true, 1);
                        
                        uwstrvarid=netcdf.defVar(ncid, 'uwstr', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,uwstrvarid,'long_name','uwind stress');
                        netcdf.putAtt(ncid,uwstrvarid,'units','N/m^2');
                        netcdf.defVarChunking(ncid, uwstrvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, uwstrvarid, true, true, 1);
                        
                        vwstrvarid=netcdf.defVar(ncid, 'vwstr', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,vwstrvarid,'long_name','vwstr');
                        netcdf.putAtt(ncid,vwstrvarid,'units','N/m^2');
                        netcdf.defVarChunking(ncid, vwstrvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, vwstrvarid, true, true, 1);
                        
                        wind_curlvarid=netcdf.defVar(ncid, 'wind_curl', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,wind_curlvarid,'long_name','wind_curl');
                        netcdf.putAtt(ncid,wind_curlvarid,'units',' /s');
                        netcdf.defVarChunking(ncid, wind_curlvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, wind_curlvarid, true, true, 1);
                        
                        wstr_curlvarid=netcdf.defVar(ncid, 'wstr_curl', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,wstr_curlvarid,'long_name','wind stress curl');
                        netcdf.putAtt(ncid,wstr_curlvarid,'units',' /s');
                        netcdf.defVarChunking(ncid, wstr_curlvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, wstr_curlvarid, true, true, 1);
                        
                        shfluxvarid=netcdf.defVar(ncid, 'shflux', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,shfluxvarid,'long_name','shflux');
                        netcdf.putAtt(ncid,shfluxvarid,'units','N/m^2');
                        netcdf.defVarChunking(ncid, shfluxvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, shfluxvarid, true, true, 1);
                        
%                         tairvarid=netcdf.defVar(ncid, 'tair', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
%                         netcdf.putAtt(ncid,tairvarid,'long_name','tair');
%                         netcdf.putAtt(ncid,tairvarid,'units','^oC ');
                        
                        pairvarid=netcdf.defVar(ncid, 'pair', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,pairvarid,'long_name','pair');
                        netcdf.putAtt(ncid,pairvarid,'units','millibar');
                        netcdf.defVarChunking(ncid, pairvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, pairvarid, true, true, 1);
                        
                        zetavarid=netcdf.defVar(ncid, 'zeta', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,zetavarid,'long_name','zeta');
                        netcdf.putAtt(ncid,zetavarid,'units','m');
                        netcdf.defVarChunking(ncid, zetavarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, zetavarid, true, true, 1);
                        
                        u_rhovarid=netcdf.defVar(ncid, 'u_rho', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,u_rhovarid,'long_name','u_rho');
                        netcdf.putAtt(ncid,u_rhovarid,'units','m/s');
                        netcdf.defVarChunking(ncid, u_rhovarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, u_rhovarid, true, true, 1);
                        
                        v_rhovarid=netcdf.defVar(ncid, 'v_rho', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,v_rhovarid,'long_name','v_rho');
                        netcdf.putAtt(ncid,v_rhovarid,'units','m.s');
                        netcdf.defVarChunking(ncid, v_rhovarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, v_rhovarid, true, true, 1);
                        
                        vortvarid=netcdf.defVar(ncid, 'vort', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,vortvarid,'long_name','vort');
                        netcdf.putAtt(ncid,vortvarid,'units',' /s');
                        netcdf.defVarChunking(ncid, vortvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, vortvarid, true, true, 1);
                        
                        egg_maskvarid=netcdf.defVar(ncid, 'egg_mask', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,egg_maskvarid,'long_name','egg_mask');
                        netcdf.putAtt(ncid,egg_maskvarid,'units','m.s');
                        netcdf.defVarChunking(ncid, egg_maskvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, egg_maskvarid, true, true, 1);
                        
%                         mask_5dayvarid=netcdf.defVar(ncid, 'mask_5day', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
%                         netcdf.putAtt(ncid,mask_5dayvarid,'long_name','mask_5day');
%                         netcdf.putAtt(ncid,mask_5dayvarid,'units','m.s');
%                         netcdf.defVarChunking(ncid, varid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
%                         netcdf.defVarDeflate(ncid, varid, true, true, 1);
%                         
%                         mask_10dayvarid=netcdf.defVar(ncid, 'mask_10day', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
%                         netcdf.putAtt(ncid,mask_10dayvarid,'long_name','mask_10day');
%                         netcdf.putAtt(ncid,mask_10dayvarid,'units','m.s');
%                         netcdf.defVarChunking(ncid, varid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
%                         netcdf.defVarDeflate(ncid, varid, true, true, 1);
                        
                        mask_15dayvarid=netcdf.defVar(ncid, 'mask_15day', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,mask_15dayvarid,'long_name','mask_15day');
                        netcdf.putAtt(ncid,mask_15dayvarid,'units','m.s');
                        netcdf.defVarChunking(ncid, mask_15dayvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, mask_15dayvarid, true, true, 1);
                        
                        mask_30dayvarid=netcdf.defVar(ncid, 'mask_30day', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                        netcdf.putAtt(ncid,mask_30dayvarid,'long_name','mask_30day');
                        netcdf.putAtt(ncid,mask_30dayvarid,'units','m.s');
                        netcdf.defVarChunking(ncid, mask_30dayvarid, 'CHUNKED', [loncount/10 latcount/10 length(comb_time)/10]);
                        netcdf.defVarDeflate(ncid, mask_30dayvarid, true, true, 1);
                        
                        netcdf.endDef(ncid);

                        netcdf.putVar(ncid, timevarid, 0, length(comb_time), comb_time);
                        netcdf.putVar(ncid, lon_rhovarid, [0 0], [loncount latcount], lon);
                        netcdf.putVar(ncid, lat_rhovarid, [0 0], [loncount latcount], lat);
                        netcdf.putVar(ncid, temp_surfvarid, [0 0 0], [loncount latcount length(comb_time)], comb_temp_surf);
                        netcdf.putVar(ncid, temp_50varid, [0 0 0], [loncount latcount length(comb_time)], comb_temp_50);
                        netcdf.putVar(ncid, uwindvarid, [0 0 0], [loncount latcount length(comb_time)], comb_uwind);
                        netcdf.putVar(ncid, vwindvarid, [0 0 0], [loncount latcount length(comb_time)], comb_vwind);
                        netcdf.putVar(ncid, uwstrvarid, [0 0 0], [loncount latcount length(comb_time)], comb_uwstr);
                        netcdf.putVar(ncid, vwstrvarid, [0 0 0], [loncount latcount length(comb_time)], comb_vwstr);
                        netcdf.putVar(ncid, wind_curlvarid, [0 0 0], [loncount latcount length(comb_time)], comb_wind_curl);
                        netcdf.putVar(ncid, wstr_curlvarid, [0 0 0], [loncount latcount length(comb_time)], comb_wstr_curl);
                        netcdf.putVar(ncid, shfluxvarid, [0 0 0], [loncount latcount length(comb_time)], comb_shflux);
%                         netcdf.putVar(ncid, tairvarid, [0 0 0], [loncount latcount length(comb_time)], comb_tair);
                        netcdf.putVar(ncid, pairvarid, [0 0 0], [loncount latcount length(comb_time)], comb_pair);
                        netcdf.putVar(ncid, zetavarid, [0 0 0], [loncount latcount length(comb_time)], comb_zeta);
                        netcdf.putVar(ncid, u_rhovarid, [0 0 0], [loncount latcount length(comb_time)], comb_u_rho);
                        netcdf.putVar(ncid, v_rhovarid, [0 0 0], [loncount latcount length(comb_time)], comb_v_rho);
                        netcdf.putVar(ncid, vortvarid, [0 0 0], [loncount latcount length(comb_time)], comb_vort);
                        netcdf.putVar(ncid, egg_maskvarid, [0 0 0], [loncount latcount length(comb_time)], comb_egg_mask);
%                         netcdf.putVar(ncid, mask_5dayvarid, [0 0 0], [loncount latcount length(comb_time)], comb_mask_5day);
%                         netcdf.putVar(ncid, mask_10dayvarid, [0 0 0], [loncount latcount length(comb_time)], comb_mask_10day);
                        netcdf.putVar(ncid, mask_15dayvarid, [0 0 0], [loncount latcount length(comb_time)], comb_mask_15day);
                        netcdf.putVar(ncid, mask_30dayvarid, [0 0 0], [loncount latcount length(comb_time)], comb_mask_30day);
                        
                        netcdf.close(ncid);
                        
                        clear comb_temp_surf comb_temp_50 comb_uwind comb_vwind comb_shflux comb_u_rho comb_v_rho comb_egg_mask comb_time comb_tair
                        clear comb_mask_5day comb_mask_10day comb_mask_15day comb_mask_30day comb_wind_curl comb_vort comb_uwstr comb_vwstr comb_wstr_curl comb_pair comb_zeta
                        toc;
%                     else
%                         disp('abc')
%                     end
                end
            end
            fig_flag=0;
        end
        
    % % %         time set
        for folding=1:1
            tind=1;
            for yearij = 1:length(inputyear)
                for month=1:12 
                    tempyear = inputyear(yearij);
                    ftime(tind) = datenum(tempyear,month,15) - datenum(1900,12,31);
                    tind=tind+1;
                end
            end
            for month=1:12 
                    tempyear = inputyear(yearij);
                    climtime(month) = datenum(1950,month,15) - datenum(1900,12,31);
            end

            for i =1:length(inputyear) 
                tempyear=inputyear(i);
                for month=1:12
                    xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
                end
            end

            trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
            trendtime_yearly=[inputyear(1) : inputyear(end)]-inputyear(1);
        end     

        fig_flag=0;
    end
end

