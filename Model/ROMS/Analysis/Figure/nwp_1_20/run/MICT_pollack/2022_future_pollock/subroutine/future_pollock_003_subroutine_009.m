% %  Updated    28-Apr-2023 by Yong-Yub Kim   % make

if testnameind==1
    
    for checkti=1:length(RCM_info.checktime)
        tmp.checktime=RCM_info.checktime(checkti);

        %% historical plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
                    num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        tmp.testname_prefix=tmp.testname_his(1:4); 
        tmp.testlen=length(RCM_info.name);
        tmp.tlen_his=0;

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            for sub_testnameind=1:length(RCM_info.name)
                tmp_sub.testname_ssp=RCM_info.name{sub_testnameind};
                [tmp_sub.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp_sub.testname_ssp);
                dirs_sub.filedir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_his, '/pollock/'); % % where data files are
                dirs_sub.savedir_his = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_his, '/pollock/');

                %%    initialization
                if sub_testnameind==1
                    for yearij = 1:length(RCM_info.years_his)
                        tmp.tempyear = RCM_info.years_his(yearij);
                        for monthij = 1:length(RCM_info.months)
                            tmp.tempmonth = RCM_info.months(monthij);
                            tmp.ncname = [dirs_sub.savedir_his,tmp_sub.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                            tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
                        end
                    end
                    tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
                    tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
                    tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his.*tmp.testlen);
%                     tmp.mov_dist_lon_mean=NaN(1, tmp.tlen_his.*tmp.testlen); % [total releasing days]
%                     tmp.mov_dist_lat_mean=NaN(1, tmp.tlen_his.*tmp.testlen);
                    tmp.mov_dist_lon_km_mean=NaN(1, tmp.tlen_his.*tmp.testlen);
                    tmp.mov_dist_lat_km_mean=NaN(1, tmp.tlen_his.*tmp.testlen);
                    tmp.NoE = NaN(1,tmp.tlen_his.*tmp.testlen); % Number of released eggs per day
                end

            %% read grids, eggs    
                for yearij = 1:length(RCM_info.years_his)
                    tmp.tempyear = RCM_info.years_his(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs_sub.savedir_his,tmp_sub.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        
                        tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                        tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                    
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
                        tmp.ini_egg_mask=squeeze(ncread(tmp.ncname, 'egg_mask', [1 1 1], [inf inf inf]));
%                         tmp.md_lon=squeeze(ncread(tmp.ncname, 'mov_dist_lon_mean', [1 tmp.ind_checktime], [inf 1]));
%                         tmp.md_lat=squeeze(ncread(tmp.ncname, 'mov_dist_lat_mean', [1 tmp.ind_checktime], [inf 1]));
                        tmp.md_lon_km=squeeze(ncread(tmp.ncname, 'mov_dist_lon_km_mean', [1 tmp.ind_checktime], [inf 1]));
                        tmp.md_lat_km=squeeze(ncread(tmp.ncname, 'mov_dist_lat_km_mean', [1 tmp.ind_checktime], [inf 1]));
                        tmp.lastday_m=size(tmp.egg_mask,3);
                        if yearij==1 && monthij==1 && sub_testnameind==1
                            tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                            tmp.NoE(1:tmp.lastday_m) = sum(sum(tmp.ini_egg_mask,1),2);
%                             tmp.mov_dist_lon_mean(1:tmp.lastday_m) = tmp.md_lon;
%                             tmp.mov_dist_lat_mean(1:tmp.lastday_m) = tmp.md_lat;
                            tmp.mov_dist_lon_km_mean(1:tmp.lastday_m) = tmp.md_lon_km;
                            tmp.mov_dist_lat_km_mean(1:tmp.lastday_m) = tmp.md_lat_km;
                            tmp.endij=tmp.lastday_m;
                        else
                            tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.egg_mask;
                            tmp.NoE(tmp.endij+1:tmp.endij+tmp.lastday_m) = sum(sum(tmp.ini_egg_mask,1),2);
%                             tmp.mov_dist_lon_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lon;
%                             tmp.mov_dist_lat_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lat;
                            tmp.mov_dist_lon_km_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lon_km;
                            tmp.mov_dist_lat_km_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lat_km;
                            tmp.endij=tmp.endij+tmp.lastday_m;
                        end
                    end
                end
                RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
                RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');


            end
            tmp.comb_egg_mask_testsep=reshape(tmp.comb_egg_mask, ...
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), tmp.testlen, size(tmp.comb_egg_mask,3)/tmp.testlen]);
            tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,4));

            tmp.mean_md_lon_km = sum(tmp.mov_dist_lon_km_mean(isfinite(tmp.mov_dist_lon_km_mean)) ...
                .* tmp.NoE(isfinite(tmp.mov_dist_lon_km_mean))) ...
                ./ sum(tmp.NoE(isfinite(tmp.mov_dist_lon_km_mean))); % particle number weighted mean of travel distance
            tmp.mean_md_lat_km = sum(tmp.mov_dist_lat_km_mean(isfinite(tmp.mov_dist_lat_km_mean)) ...
                .* tmp.NoE(isfinite(tmp.mov_dist_lat_km_mean))) ...
                ./ sum(tmp.NoE(isfinite(tmp.mov_dist_lat_km_mean))); % particle number weighted mean of travel distance
            tmp.mean_data = sum(tmp.comb_egg_mask,3)./tmp.testlen;
            tmp.mean_data_his=tmp.mean_data;
            tmp.mean_md_lon_km_his(checkti) = tmp.mean_md_lon_km
            tmp.mean_md_lat_km_his(checkti) = tmp.mean_md_lat_km         
            tmp.std_data = 2*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.mean_data(tmp.mean_data==0)=NaN;
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;

        end



        %% future(ssp) plot
        tmp.tifname=strcat(dirs.figdir_allmean, 'all_mean','_',tmp.regionname, '_loc_num_', num2str(tmp.checktime, '%02i'),'days', ...
                    num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                    RCM_info.season, '.tif'); 

        tmp.testname_prefix=tmp.testname_ssp(1:4); 
        tmp.testlen=length(RCM_info.name);
        tmp.tlen_ssp=0;

        if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)

            run(tmp.param_script);

            for sub_testnameind=1:length(RCM_info.name)
                tmp_sub.testname_ssp=RCM_info.name{sub_testnameind};
                dirs_sub.filedir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_ssp, '/pollock/'); % % where data files are
                dirs_sub.savedir_ssp = strcat('/Volumes/kyy_raid/Data/Model/ROMS/nwp_1_20/', tmp_sub.testname_ssp, '/pollock/');

                %%    initialization
                if sub_testnameind==1
                    for yearij = 1:length(RCM_info.years_ssp)
                        tmp.tempyear = RCM_info.years_ssp(yearij);
                        for monthij = 1:length(RCM_info.months)
                            tmp.tempmonth = RCM_info.months(monthij);
                            tmp.ncname = [dirs_sub.savedir_ssp,tmp_sub.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                            tmp.tlen_ssp=tmp.tlen_ssp + length(ncread(tmp.ncname, 'time'));
                        end
                    end
                    tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
                    tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
                    tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp.*tmp.testlen);
                    tmp.mov_dist_lon_km_mean=NaN(1, tmp.tlen_ssp.*tmp.testlen);
                    tmp.mov_dist_lat_km_mean=NaN(1, tmp.tlen_ssp.*tmp.testlen);
                    tmp.NoE = NaN(1,tmp.tlen_ssp.*tmp.testlen); % Number of released eggs per day
                end

            %% read grids, eggs    
                for yearij = 1:length(RCM_info.years_ssp)
                    tmp.tempyear = RCM_info.years_ssp(yearij);
                    for monthij = 1:length(RCM_info.months)
                        tmp.tempmonth = RCM_info.months(monthij);
                        tmp.ncname = [dirs_sub.savedir_ssp,tmp_sub.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        
                        tmp.all_checktime=ncread(tmp.ncname, 'checktime');
                        tmp.ind_checktime=find(tmp.all_checktime==tmp.checktime);
                    
                        tmp.egg_mask=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
                        tmp.ini_egg_mask=squeeze(ncread(tmp.ncname, 'egg_mask', [1 1 1], [inf inf inf]));
%                         tmp.md_lon=squeeze(ncread(tmp.ncname, 'mov_dist_lon_mean', [1 tmp.ind_checktime], [inf 1]));
%                         tmp.md_lat=squeeze(ncread(tmp.ncname, 'mov_dist_lat_mean', [1 tmp.ind_checktime], [inf 1]));
                        tmp.md_lon_km=squeeze(ncread(tmp.ncname, 'mov_dist_lon_km_mean', [1 tmp.ind_checktime], [inf 1]));
                        tmp.md_lat_km=squeeze(ncread(tmp.ncname, 'mov_dist_lat_km_mean', [1 tmp.ind_checktime], [inf 1]));

                        tmp.lastday_m=size(tmp.egg_mask,3);
                        if yearij==1 && monthij==1 && sub_testnameind==1
                            tmp.comb_egg_mask(:,:,1:tmp.lastday_m)=tmp.egg_mask;
                            tmp.NoE(1:tmp.lastday_m) = sum(sum(tmp.ini_egg_mask,1),2);
%                             tmp.mov_dist_lon_mean(1:tmp.lastday_m) = tmp.md_lon;
%                             tmp.mov_dist_lat_mean(1:tmp.lastday_m) = tmp.md_lat;
                            tmp.mov_dist_lon_km_mean(1:tmp.lastday_m) = tmp.md_lon_km;
                            tmp.mov_dist_lat_km_mean(1:tmp.lastday_m) = tmp.md_lat_km;
                            tmp.endij=tmp.lastday_m;
                        else
                            tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.egg_mask;
                            tmp.NoE(tmp.endij+1:tmp.endij+tmp.lastday_m) = sum(sum(tmp.ini_egg_mask,1),2);
%                             tmp.mov_dist_lon_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lon;
%                             tmp.mov_dist_lat_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lat;
                            tmp.mov_dist_lon_km_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lon_km;
                            tmp.mov_dist_lat_km_mean(tmp.endij+1:tmp.endij+tmp.lastday_m) = tmp.md_lat_km;
                            tmp.endij=tmp.endij+tmp.lastday_m;
                        end
                    end
                end
                RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
                RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');


            end
            tmp.mean_md_lon_km = sum(tmp.mov_dist_lon_km_mean(isfinite(tmp.mov_dist_lon_km_mean)) ...
                .* tmp.NoE(isfinite(tmp.mov_dist_lon_km_mean))) ...
                ./ sum(tmp.NoE(isfinite(tmp.mov_dist_lon_km_mean))); % particle number weighted mean of travel distance
            tmp.mean_md_lat_km = sum(tmp.mov_dist_lat_km_mean(isfinite(tmp.mov_dist_lat_km_mean)) ...
                .* tmp.NoE(isfinite(tmp.mov_dist_lat_km_mean))) ...
                ./ sum(tmp.NoE(isfinite(tmp.mov_dist_lat_km_mean))); % particle number weighted mean of travel distance
            tmp.mean_data = sum(tmp.comb_egg_mask,3)./tmp.testlen;
            tmp.mean_data_ssp=tmp.mean_data;
            tmp.mean_md_lon_km_ssp(checkti) = tmp.mean_md_lon_km
            tmp.mean_md_lat_km_ssp(checkti) = tmp.mean_md_lat_km     
            
            tmp.mean_md_lon_km_diff(checkti) = tmp.mean_md_lon_km_ssp(checkti) - tmp.mean_md_lon_km_his(checkti)
            tmp.mean_md_lat_km_diff(checkti) = tmp.mean_md_lat_km_ssp(checkti) - tmp.mean_md_lat_km_his(checkti)

            tmp.mean_data(tmp.mean_data==0)=NaN;
            
            tmp.explen=length(RCM_info.name);
            tmp.comb_egg_mask_testsep=reshape(tmp.comb_egg_mask, ...
                [size(tmp.comb_egg_mask,1) size(tmp.comb_egg_mask,2), tmp.explen, size(tmp.comb_egg_mask,3)/tmp.explen]);
            tmp.comb_egg_mask_testsep=squeeze(sum(tmp.comb_egg_mask_testsep,4));
            tmp.std_data = 2*std(tmp.comb_egg_mask_testsep,0,3);
            tmp.agreement_data=NaN(size(tmp.mean_data));
            tmp.agreement_data(tmp.mean_data>tmp.std_data)=1;
            
        end
        
    end

end