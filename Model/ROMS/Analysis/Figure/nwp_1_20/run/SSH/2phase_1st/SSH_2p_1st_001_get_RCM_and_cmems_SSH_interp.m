close all; clear all;  clc;
% %  Updated 05-Jul-2021 by Yong-Yub Kim, structure

% % configuration of RCM
RCM_info.name={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};
RCM_info.abbs = {'RCM-CNRM', 'RCM-EC-Veg', 'RCM-ACC', 'RCM-CNRM-HR', 'RCM-CMCC'};
RCM_info.model = 'nwp_1_20';
% RCM_info.dataroot = '/data1/RCM/CMIP6/';
RCM_info.dataroot = ['E:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'output', filesep];
% RCM_info.saveroot = '/data1/RCM/CMIP6/';
RCM_info.saveroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'ROMS', filesep, 'nwp_1_20', filesep, 'output', filesep];
RCM_info.progress = 'run';
RCM_info.region = {'AKP4'};
RCM_info.years = 1993:2014;
RCM_info.months = 1:12;
RCM_grid.dl = 1/20;

% % configuration of GCM
GCM_info.name={'CNRM-ESM2-1', 'EC-Earth3-Veg', 'ACCESS-CM2', 'CNRM-CM6-1-HR', 'CMCC-CM2-HR4'};
GCM_info.abbs = {'GCM-CNRM', 'GCM-EC-Veg', 'GCM-ACC', 'GCM-CNRM-HR', 'GCM-CMCC'};
GCM_info.model = GCM_info.name;
GCM_info.dataroot = ['D:', filesep, 'Data', filesep, 'Model', filesep, ...
    'CMIP6', filesep, 'zos', filesep];
GCM_info.saveroot = GCM_info.dataroot;
GCM_info.progress = RCM_info.progress;
GCM_info.region = RCM_info.region;
GCM_info.years = RCM_info.years;
GCM_info.months = RCM_info.months;
GCM_grid.dl = 1/2;

% % configuration of flags (make file)
flags.make_std_mat = 1;
flags.make_std_nc = 1;

% % configuration of figure levels
param.fig_lev_shad = [-2 2];
param.fig_lev_shad_trend = [0 4];
param.fig_lev_shad_bias = [-4 4];
param.fig_lev_con = 0:5:35;
param.fig_lev_rms = [0 4];

% %  working
for testnameind=1:length(RCM_info.name)
    for regionind=1:length(RCM_info.region)
        clearvars '*' -except regionind testnameind RCM_info flags flg_lev RCM_grid param ...
            GCM_info GCM_grid
        
        tmp.fs = filesep; % file separator win = '\', linux = '/'
        
% %     set temporary variables (testname, regionname, filesep, ...)
        RCM_info.testname = RCM_info.name{testnameind};
        RCM_info.regionname = RCM_info.region{regionind};
        RCM_info.abb = RCM_info.abbs{testnameind};
        [RCM_info.scenario, tmp.error_status] = Func_0013_RCM_CMIP6_scenname(RCM_info.testname);
        
        GCM_info.testname = GCM_info.name{testnameind};
        GCM_info.regionname = RCM_info.regionname;
        GCM_info.abb = GCM_info.abbs{testnameind};
        GCM_info.scenario = RCM_info.scenario;
        
        variable = 'zeta';
        % %     set dropbox path
        if (strcmp(computer,'PCWIN64'))
            tmp.dropboxpath = 'C:\Users\User\Dropbox';
        else
            tmp.dropboxpath = '/home/kimyy/Dropbox';
        end
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
        
        [tmp.error_status, tmp.dropboxpath] = Func_0008_set_dropbox_path(computer);
        addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
            tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
            'run', tmp.fs, 'SSH', tmp.fs, '2phase_1st', tmp.fs, 'subroutine']));
        
        [tmp.error_status, RCM_grid.refpolygon, RCM_grid.domain] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

        tmp.param_script = [tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', ...
            tmp.fs, 'Model', tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', ...
            tmp.fs, 'nwp_1_20', tmp.fs, 'run', tmp.fs, 'fig_param', tmp.fs, ...
            'fig_param2_kyy_', RCM_info.regionname, '.m'];
        RCM_info.filedir = [RCM_info.dataroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            'packed_monthly', tmp.fs];
        RCM_info.savedir = [RCM_info.saveroot, RCM_info.testname, tmp.fs, 'run', tmp.fs, ...
            'packed_monthly', tmp.fs];
        CMEMS_info.filedir = 'Z:\내 드라이브\Data\Observation\CMEMS\';
        
% % %         flag configuration (process)
        for folding=1:1
            flags.fig_name{1}='get model data and cmemsstructed data';
            flags.fig_name{2}='sea level analysis';
            flags.fig_name{3}='cmems sea level trend analysis';
            flags.fig_name{4}='interped sea level trend analysis';
            flags.fig_name{5}='interped sea level correlation analysis';
            flags.fig_name{6}='low pass filtered interped sea level trend analysis';
            flags.fig_name{7}='low pass filtered interped sea level correlation analysis';
            flags.fig_name{8}='detrended sea level correlation analysis';
            flags.fig_name{9}='moving averaged interped sea level trend analysis';
            flags.fig_name{10}='moving averaged interped sea level correlation analysis';

            for flagi=1:length(flags.fig_name)
                flags.fig_switch(flagi)=0;
            end
            flags.fig_switch(1)=1;
            flags.fig_switch(2)=1;
            flags.fig_switch(3)=1;
            flags.fig_switch(4)=1;
            flags.fig_switch(5)=1;
            flags.fig_switch(6)=0;
            flags.fig_switch(7)=0;
            flags.fig_switch(8)=0;
            flags.fig_switch(9)=1;
            flags.fig_switch(10)=1;
        end
        
        if flags.fig_switch(1) > 0
            flags.fig_tmp = flags.fig_switch(1);
            SSH_2p_1st_001_sub_001_get_SSH;
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
            trendtime_yearly=inputyear(1) : inputyear(end);
        end     
    
% % %         sea level analysis
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
            % % %         trend
                trend(1:len_lon_model,1:len_lat_model)=NaN;
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        p=polyfit(trendtime,squeeze(comb_data(i,j,:))',1);
                        trend(i,j)=p(1);
                    end
                end
                trend = trend * 1000.0; %% m/y -> mm/y
               disp('trend complete') 

            % % %         trend (seasonal filtered)
                for t=1:length(inputyear)
                    comb_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(comb_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel);
                end

                trend_filtered(1:len_lon_model,1:len_lat_model)=NaN;
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        p=polyfit(trendtime,squeeze(comb_data_filtered(i,j,:))',1);
                        trend_filtered(i,j)=p(1);
                    end
                end
                trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y
                disp('trend_filtered complete') 

            % % %         climatological trend 
                comb_spatial_data=reshape(comb_data, [len_lon_model len_lat_model 12 length(inputyear)]);
                climtrendtime=inputyear(1):inputyear(end);
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        for k=1:12  % month
                            p=polyfit(climtrendtime,squeeze(comb_spatial_data(i,j,k,:))',1);
                            trend_clim(i,j,k)=p(1);
                        end
                    end
                end
                disp('climatological trend complete') 

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
                lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
                netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

                lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
                netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

                raw_sshvarid=netcdf.defVar(ncid, 'raw_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,raw_sshvarid,'long_name','raw_ssh');
                netcdf.putAtt(ncid,raw_sshvarid,'units','m');

                ssh_filteredvarid=netcdf.defVar(ncid, 'ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,ssh_filteredvarid,'long_name','ssh_filtered');
                netcdf.putAtt(ncid,ssh_filteredvarid,'units','m');

                trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,trendvarid,'long_name','trend');
                netcdf.putAtt(ncid,trendvarid,'units','mm/year');

                trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
                netcdf.putAtt(ncid,trend_filteredvarid,'units','mm/year');

                clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_sshvarid,'long_name','clim_ssh');
                netcdf.putAtt(ncid,clim_sshvarid,'units','m');

                clim_ssh_trendvarid=netcdf.defVar(ncid, 'clim_ssh_trend', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','clim_ssh_trend');
                netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
                netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
                netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
                netcdf.putVar(ncid, ssh_filteredvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data_filtered);
                netcdf.putVar(ncid, trendvarid, [0 0], [len_lon_model, len_lat_model], trend);
                netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon_model, len_lat_model], trend_filtered);
                netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [len_lon_model, len_lat_model length(climtime)], comb_spatial_meanmodel);
                netcdf.putVar(ncid, clim_ssh_trendvarid, [0 0 0], [len_lon_model, len_lat_model length(climtime)], trend_clim);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     cmems sea level trend analysis
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   

                len_lon_cmems=cmems_lonsize_cut;
                len_lat_cmems=cmems_latsize_cut;

                cmems_sla=comb_cmems_data;
                cmems_sla_divided=reshape(cmems_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                cmems_sla_yearly=squeeze(mean(cmems_sla_divided,3,'omitnan')).*100;
                clim_cmems_sla=mean(cmems_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    cmems_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(cmems_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_cmems_sla);
                end

            % % %         cmems trend 
                comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

                cmems_trend(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_cmems_data(i,j,:))',1);
                        cmems_trend(i,j)=p(1) * 1000.0 ;
                    end
                end
               disp('cmems trend complete') 
            
           % % %         cmems yearly trend 
            cmems_trend_yearly(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    p=polyfit(trendtime_yearly,squeeze(cmems_sla_yearly(i,j,:))',1);
                    cmems_trend_yearly(i,j)=p(1) * 10.0 ;
                end
            end
           disp('cmems yearly trend complete') 
           
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    for k=1:size(cmems_sla_yearly,3)
                        cmems_sla_yearly_detrended(i,j,k)=cmems_sla_yearly(i,j,k)-(cmems_trend_yearly(i,j)/10.0)*(k-1);
                    end
                end
            end
%             plot(squeeze(mean(mean(cmems_sla_yearly_detrended,1,'omitnan'),2,'omitnan')))
%             hold on
%             plot(squeeze(mean(mean(cmems_sla_yearly,1,'omitnan'),2,'omitnan')))
%             hold off

             % % %         interped exponential fitting yearly
                corr_for_exp=abs(2*(min(cmems_sla_yearly(:))));
                tic
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        if isfinite(cmems_sla_yearly(i,j,1))
                            [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:))+corr_for_exp,'exp1');
                            cmems_sla_yearly_exp_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
                            cmems_sla_yearly_exp_fit_rsquare(i,j)=gof_exp.rsquare;
                            cmems_sla_yearly_exp_fit_rmse(i,j)=gof_exp.rmse;
                        else
                            cmems_sla_yearly_exp_fit(i,j,1:length(trendtime_yearly))=NaN;
                            cmems_sla_yearly_exp_fit_rsquare(i,j)=NaN;
                            cmems_sla_yearly_exp_fit_rmse(i,j)=NaN;
                        end
%                         p=polyfit(trendtime_yearly,squeeze(cmems_sla_yearly(i,j,:))',1);
%                         cmems_trend_yearly(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('cmems exponential fitting yearly complete') 
                toc;
                
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:size(cmems_sla_yearly,3)
                            cmems_sla_yearly_exp_fit_detrended(i,j,k)=cmems_sla_yearly(i,j,k)-cmems_sla_yearly_exp_fit(i,j,k);
                        end
                    end
                end
                
                % % %         cmems poly1 fitting yearly
                tic;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        if isfinite(cmems_sla_yearly(i,j,1))
                            [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(cmems_sla_yearly(i,j,:)),'poly1');
                            cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
                            cmems_sla_yearly_poly1_fit_rsquare(i,j)=gof_exp.rsquare;
                            cmems_sla_yearly_poly1_fit_rmse(i,j)=gof_exp.rmse;
                        else
                            cmems_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=NaN;
                            cmems_sla_yearly_poly1_fit_rsquare(i,j)=NaN;
                            cmems_sla_yearly_poly1_fit_rmse(i,j)=NaN;
                        end
%                         p=polyfit(trendtime_yearly,squeeze(cmems_sla_yearly(i,j,:))',1);
%                         cmems_trend_yearly(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('cmems poly1 fitting yearly complete') 
                toc;
                
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:size(cmems_sla_yearly,3)
                            cmems_sla_yearly_poly1_fit_detrended(i,j,k)=cmems_sla_yearly(i,j,k)-cmems_sla_yearly_poly1_fit(i,j,k);
                        end
                    end
                end


            % % %         cmems trend (seasonal filtered)
                for t=1:length(inputyear)
                    comb_cmems_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_cmems_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meancmems;
                end

                cmems_trend_filtered(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_cmems_data_filtered(i,j,:))',1);
                        cmems_trend_filtered(i,j)=p(1) * 1000.0 ;
                    end
                end
               disp('cmems trend_filtered complete') 

            % % %         climatological cmems trend 

                clim_cmems_trend_divided(1:cmems_lonsize_cut,1:cmems_latsize_cut,1:12)=NaN;
                for k=1:12
                    clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            p=polyfit(clim_trendtime(:,k)',squeeze(comb_cmems_clim_divided(i,j,k,:))',1);
                            clim_cmems_trend_divided(i,j,k)=p(1) * 1000.0 ;
                        end
                    end
                end
               disp('cmems climatological trend complete') 

                mean_cmems_trend=mean(mean(cmems_trend,'omitnan'),'omitnan');
                mean_cmems_trend_filtered=mean(mean(cmems_trend_filtered,'omitnan'),'omitnan');
                mean_clim_cmems_trend_divided=mean(mean(clim_cmems_trend_divided,'omitnan'),'omitnan');

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
                yearly_time_dimid = netcdf.defDim(ncid, 'yearly_time', size(cmems_sla_yearly,3));

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');
                
                yearly_timevarid=netcdf.defVar(ncid, 'yearly_time', 'NC_DOUBLE', yearly_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','yearly_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                cmems_slavarid=netcdf.defVar(ncid, 'cmems_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,cmems_slavarid,'long_name','cmems_sla');
                netcdf.putAtt(ncid,cmems_slavarid,'units','m ');
                
                cmems_adtvarid=netcdf.defVar(ncid, 'cmems_adt', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,cmems_adtvarid,'long_name','cmems_adt');
                netcdf.putAtt(ncid,cmems_adtvarid,'units','m ');
                
                cmems_sla_yearlyvarid=netcdf.defVar(ncid, 'cmems_sla_yearly', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearlyvarid,'long_name','cmems_sla_yearly');
                netcdf.putAtt(ncid,cmems_sla_yearlyvarid,'units','m ');
                
                cmems_sla_yearly_detrendedvarid=netcdf.defVar(ncid, 'cmems_sla_yearly_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_detrendedvarid,'long_name','cmems_sla_yearly_detrended');
                netcdf.putAtt(ncid,cmems_sla_yearly_detrendedvarid,'units','m ');
                
                cmems_sla_yearly_exp_fitvarid=netcdf.defVar(ncid, 'cmems_sla_yearly_exp_fit', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fitvarid,'long_name','cmems_sla_yearly_exp_fit');
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fitvarid,'units','cm ');
                
                cmems_sla_yearly_poly1_fitvarid=netcdf.defVar(ncid, 'cmems_sla_yearly_poly1_fit', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fitvarid,'long_name','cmems_sla_yearly_poly1_fit');
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fitvarid,'units','cm ');
                
                cmems_sla_yearly_exp_fit_detrendedvarid=netcdf.defVar(ncid, 'cmems_sla_yearly_exp_fit_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_detrendedvarid,'long_name','cmems_sla_yearly_exp_fit_detrended');
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_detrendedvarid,'units','cm ');
                 
                cmems_sla_yearly_poly1_fit_detrendedvarid=netcdf.defVar(ncid, 'cmems_sla_yearly_poly1_fit_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_detrendedvarid,'long_name','cmems_sla_yearly_poly1_fit_detrended');
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_detrendedvarid,'units','cm ');
                
                cmems_sla_filteredvarid=netcdf.defVar(ncid, 'cmems_sla_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_filteredvarid,'long_name','cmems_sla_filtered');
                netcdf.putAtt(ncid,cmems_sla_filteredvarid,'units','m ');

                clim_cmemsvarid=netcdf.defVar(ncid, 'clim_cmems_ssh', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_cmemsvarid,'long_name','clim_cmems_ssh');
                netcdf.putAtt(ncid,clim_cmemsvarid,'units','m ');

                cmems_trendvarid=netcdf.defVar(ncid, 'cmems_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_trendvarid,'long_name','cmems_trend');
                netcdf.putAtt(ncid,cmems_trendvarid,'units','mm /year');
                                
                cmems_trend_yearlyvarid=netcdf.defVar(ncid, 'cmems_trend_yearly', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_trend_yearlyvarid,'long_name','cmems_trend_yearly');
                netcdf.putAtt(ncid,cmems_trend_yearlyvarid,'units','mm /year');

                cmems_sla_yearly_exp_fit_rsquarevarid=netcdf.defVar(ncid, 'cmems_sla_yearly_exp_fit_rsquare', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_rsquarevarid,'long_name','cmems_sla_yearly_exp_fit_rsquare');
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_rsquarevarid,'units',' ');
                
                cmems_sla_yearly_exp_fit_rmsevarid=netcdf.defVar(ncid, 'cmems_sla_yearly_exp_fit_rmse', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_rmsevarid,'long_name','cmems_sla_yearly_exp_fit_rmse');
                netcdf.putAtt(ncid,cmems_sla_yearly_exp_fit_rmsevarid,'units','cm');
                
                cmems_sla_yearly_poly1_fit_rsquarevarid=netcdf.defVar(ncid, 'cmems_sla_yearly_poly1_fit_rsquare', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_rsquarevarid,'long_name','cmems_sla_yearly_poly1_fit_rsquare');
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_rsquarevarid,'units',' ');
                
                cmems_sla_yearly_poly1_fit_rmsevarid=netcdf.defVar(ncid, 'cmems_sla_yearly_poly1_fit_rmse', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_rmsevarid,'long_name','cmems_sla_yearly_poly1_fit_rmse');
                netcdf.putAtt(ncid,cmems_sla_yearly_poly1_fit_rmsevarid,'units','cm');
                
                cmems_trend_filteredvarid=netcdf.defVar(ncid, 'cmems_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_trend_filteredvarid,'long_name','cmems_trend_filtered');
                netcdf.putAtt(ncid,cmems_trend_filteredvarid,'units','mm /year');

                clim_cmems_trend_dividedvarid=netcdf.defVar(ncid, 'clim_cmems_trend_divided', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'long_name','clim_cmems_trend_divided');
                netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'units','mm /year');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, yearly_timevarid, 0, length(trendtime_yearly), trendtime_yearly);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                netcdf.putVar(ncid, cmems_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data);
                netcdf.putVar(ncid, cmems_adtvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_adt);
                netcdf.putVar(ncid, cmems_sla_yearlyvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly);
                netcdf.putVar(ncid, cmems_sla_yearly_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly_detrended);
                netcdf.putVar(ncid, cmems_sla_yearly_exp_fitvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly_exp_fit);                
                netcdf.putVar(ncid, cmems_sla_yearly_poly1_fitvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly_poly1_fit);                
                netcdf.putVar(ncid, cmems_sla_yearly_exp_fit_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly_exp_fit_detrended);                
                netcdf.putVar(ncid, cmems_sla_yearly_poly1_fit_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], cmems_sla_yearly_poly1_fit_detrended);                
                netcdf.putVar(ncid, cmems_sla_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data_filtered);
                netcdf.putVar(ncid, clim_cmemsvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], comb_spatial_meancmems);
                netcdf.putVar(ncid, cmems_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend);
                netcdf.putVar(ncid, cmems_trend_yearlyvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend_yearly);
                netcdf.putVar(ncid, cmems_sla_yearly_exp_fit_rsquarevarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_sla_yearly_exp_fit_rsquare);                
                netcdf.putVar(ncid, cmems_sla_yearly_exp_fit_rmsevarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_sla_yearly_exp_fit_rmse);                
                netcdf.putVar(ncid, cmems_sla_yearly_poly1_fit_rsquarevarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_sla_yearly_poly1_fit_rsquare);                
                netcdf.putVar(ncid, cmems_sla_yearly_poly1_fit_rmsevarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_sla_yearly_poly1_fit_rmse); 
                netcdf.putVar(ncid, cmems_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend_filtered);
                netcdf.putVar(ncid, clim_cmems_trend_dividedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], clim_cmems_trend_divided);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     interped sea level trend analysis
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
               cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;
                             
                cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               cmems_sla_mean = mean(cmems_sla,3);
               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                interped_sla_yearly=squeeze(mean(interped_sla_divided,3,'omitnan')).*100;
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % % %         interped trend 
                comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

                interped_trend(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                        interped_trend(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('interped trend complete') 
                
                % % %         interped trend yearly
                interped_trend_yearly(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime_yearly,squeeze(interped_sla_yearly(i,j,:))',1);
                        interped_trend_yearly(i,j)=p(1) * 10.0 ;
                    end
                end
                disp('interped trend yearly complete') 

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:size(cmems_sla_yearly,3)
                            interped_sla_yearly_detrended(i,j,k)=interped_sla_yearly(i,j,k)-(interped_trend_yearly(i,j)/10.0)*(k-1);
                        end
                    end
                end
                
                % % %         interped exponential fitting yearly
                corr_for_exp=abs(2*(min(interped_sla_yearly(:))));
                tic
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        if isfinite(interped_sla_yearly(i,j,1))
                            [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(interped_sla_yearly(i,j,:))+corr_for_exp,'exp1');
                            interped_sla_yearly_exp_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
                            interped_sla_yearly_exp_fit_rsquare(i,j)=gof_exp.rsquare;
                            interped_sla_yearly_exp_fit_rmse(i,j)=gof_exp.rmse;
                        else
                            interped_sla_yearly_exp_fit(i,j,1:length(trendtime_yearly))=NaN;
                            interped_sla_yearly_exp_fit_rsquare(i,j)=NaN;
                            interped_sla_yearly_exp_fit_rmse(i,j)=NaN;
                        end
%                         p=polyfit(trendtime_yearly,squeeze(interped_sla_yearly(i,j,:))',1);
%                         interped_trend_yearly(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('interped exponential fitting yearly complete') 
                toc;
                
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:size(cmems_sla_yearly,3)
                            interped_sla_yearly_exp_fit_detrended(i,j,k)=interped_sla_yearly(i,j,k)-interped_sla_yearly_exp_fit(i,j,k);
                        end
                    end
                end
                
                % % %         interped poly1 fitting yearly
                tic;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        if isfinite(interped_sla_yearly(i,j,1))
                            [f_exp, gof_exp] = fit(trendtime_yearly',squeeze(interped_sla_yearly(i,j,:)),'poly1');
                            interped_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=f_exp(trendtime_yearly);
                            interped_sla_yearly_poly1_fit_rsquare(i,j)=gof_exp.rsquare;
                            interped_sla_yearly_poly1_fit_rmse(i,j)=gof_exp.rmse;
                        else
                            interped_sla_yearly_poly1_fit(i,j,1:length(trendtime_yearly))=NaN;
                            interped_sla_yearly_poly1_fit_rsquare(i,j)=NaN;
                            interped_sla_yearly_poly1_fit_rmse(i,j)=NaN;
                        end
%                         p=polyfit(trendtime_yearly,squeeze(interped_sla_yearly(i,j,:))',1);
%                         interped_trend_yearly(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('interped poly1 fitting yearly complete') 
                toc;
                
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:size(cmems_sla_yearly,3)
                            interped_sla_yearly_poly1_fit_detrended(i,j,k)=interped_sla_yearly(i,j,k)-interped_sla_yearly_poly1_fit(i,j,k);
                        end
                    end
                end
                
                
        % % %         interped trend (seasonal filtered)
                for t=1:length(inputyear)
                    comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meaninterped;
                end

                interped_trend_filtered(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                        interped_trend_filtered(i,j)=p(1) * 1000.0 ;
                    end
                end
                disp('interped trend_filtered complete') 
               mean_trend_filtered=mean(interped_trend_filtered(:), 'omitnan');
                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:size(comb_data,3)
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
                yearly_time_dimid = netcdf.defDim(ncid, 'yearly_time', size(interped_sla_yearly,3));

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');
                
                yearly_timevarid=netcdf.defVar(ncid, 'yearly_time', 'NC_DOUBLE', yearly_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','yearly_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');
                
                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                interped_sshvarid=netcdf.defVar(ncid, 'interped_ssh', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_sshvarid,'long_name','interped_ssh');
                netcdf.putAtt(ncid,interped_sshvarid,'units','m ');

                interped_trendvarid=netcdf.defVar(ncid, 'interped_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_trendvarid,'long_name','interped_trend');
                netcdf.putAtt(ncid,interped_trendvarid,'units','mm /year');
                
                interped_trend_yearlyvarid=netcdf.defVar(ncid, 'interped_trend_yearly', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_trend_yearlyvarid,'long_name','interped_trend_yearly');
                netcdf.putAtt(ncid,interped_trend_yearlyvarid,'units','mm /year');
                
                interped_sla_yearly_exp_fit_rsquarevarid=netcdf.defVar(ncid, 'interped_sla_yearly_exp_fit_rsquare', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_rsquarevarid,'long_name','interped_sla_yearly_exp_fit_rsquare');
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_rsquarevarid,'units',' ');
                
                interped_sla_yearly_exp_fit_rmsevarid=netcdf.defVar(ncid, 'interped_sla_yearly_exp_fit_rmse', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_rmsevarid,'long_name','interped_sla_yearly_exp_fit_rmse');
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_rmsevarid,'units','cm');
                
                interped_sla_yearly_poly1_fit_rsquarevarid=netcdf.defVar(ncid, 'interped_sla_yearly_poly1_fit_rsquare', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_rsquarevarid,'long_name','interped_sla_yearly_poly1_fit_rsquare');
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_rsquarevarid,'units',' ');
                
                interped_sla_yearly_poly1_fit_rmsevarid=netcdf.defVar(ncid, 'interped_sla_yearly_poly1_fit_rmse', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_rmsevarid,'long_name','interped_sla_yearly_poly1_fit_rmse');
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_rmsevarid,'units','cm');
   
                interped_trend_filteredvarid=netcdf.defVar(ncid, 'interped_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_trend_filteredvarid,'long_name','interped_trend_filtered_trend');
                netcdf.putAtt(ncid,interped_trend_filteredvarid,'units','mm /year');

                corrected_interped_sshvarid=netcdf.defVar(ncid, 'corrected_interped_ssh', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,corrected_interped_sshvarid,'long_name','corrected_interped_ssh');
                netcdf.putAtt(ncid,corrected_interped_sshvarid,'units','m ');

                interped_slavarid=netcdf.defVar(ncid, 'interped_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_slavarid,'long_name','interped_sla');
                netcdf.putAtt(ncid,interped_slavarid,'units','m ');
                
                interped_sla_yearlyvarid=netcdf.defVar(ncid, 'interped_sla_yearly', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearlyvarid,'long_name','interped_sla_yearly');
                netcdf.putAtt(ncid,interped_sla_yearlyvarid,'units','m ');
                
                interped_sla_yearly_detrendedvarid=netcdf.defVar(ncid, 'interped_sla_yearly_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_detrendedvarid,'long_name','interped_sla_yearly_detrended');
                netcdf.putAtt(ncid,interped_sla_yearly_detrendedvarid,'units','m ');
                
                corrected_interped_slavarid=netcdf.defVar(ncid, 'corrected_interped_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,corrected_interped_slavarid,'long_name','corrected_interped_sla');
                netcdf.putAtt(ncid,corrected_interped_slavarid,'units','m ');
                
                interped_sla_yearly_exp_fitvarid=netcdf.defVar(ncid, 'interped_sla_yearly_exp_fit', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fitvarid,'long_name','interped_sla_yearly_exp_fit');
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fitvarid,'units','cm ');
                
                interped_sla_yearly_poly1_fitvarid=netcdf.defVar(ncid, 'interped_sla_yearly_poly1_fit', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fitvarid,'long_name','interped_sla_yearly_poly1_fit');
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fitvarid,'units','cm ');
                
                interped_sla_yearly_exp_fit_detrendedvarid=netcdf.defVar(ncid, 'interped_sla_yearly_exp_fit_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_detrendedvarid,'long_name','interped_sla_yearly_exp_fit_detrended');
                netcdf.putAtt(ncid,interped_sla_yearly_exp_fit_detrendedvarid,'units','cm ');
                 
                interped_sla_yearly_poly1_fit_detrendedvarid=netcdf.defVar(ncid, 'interped_sla_yearly_poly1_fit_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid yearly_time_dimid]);
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_detrendedvarid,'long_name','interped_sla_yearly_poly1_fit_detrended');
                netcdf.putAtt(ncid,interped_sla_yearly_poly1_fit_detrendedvarid,'units','cm ');
                
                interped_sla_filteredvarid=netcdf.defVar(ncid, 'interped_sla_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_sla_filteredvarid,'long_name','interped_sla_filtered');
                netcdf.putAtt(ncid,interped_sla_filteredvarid,'units','m ');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, yearly_timevarid, 0, length(trendtime_yearly), trendtime_yearly);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                netcdf.putVar(ncid, interped_sshvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_interped_data);
                netcdf.putVar(ncid, corrected_interped_sshvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_interped_data_corrected);
                netcdf.putVar(ncid, interped_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], interped_sla);
                netcdf.putVar(ncid, interped_sla_yearlyvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly);
                netcdf.putVar(ncid, interped_sla_yearly_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly_detrended);
                netcdf.putVar(ncid, interped_sla_yearly_exp_fitvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly_exp_fit);                
                netcdf.putVar(ncid, interped_sla_yearly_poly1_fitvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly_poly1_fit);                
                netcdf.putVar(ncid, interped_sla_yearly_exp_fit_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly_exp_fit_detrended);                
                netcdf.putVar(ncid, interped_sla_yearly_poly1_fit_detrendedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(trendtime_yearly)], interped_sla_yearly_poly1_fit_detrended);                

                netcdf.putVar(ncid, corrected_interped_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], corrected_interped_sla);
                netcdf.putVar(ncid, interped_sla_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], interped_sla_filtered);
                netcdf.putVar(ncid, interped_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend);
                netcdf.putVar(ncid, interped_trend_yearlyvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend_yearly);
                netcdf.putVar(ncid, interped_sla_yearly_exp_fit_rsquarevarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sla_yearly_exp_fit_rsquare);                
                netcdf.putVar(ncid, interped_sla_yearly_exp_fit_rmsevarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sla_yearly_exp_fit_rmse);                
                netcdf.putVar(ncid, interped_sla_yearly_poly1_fit_rsquarevarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sla_yearly_poly1_fit_rsquare);                
                netcdf.putVar(ncid, interped_sla_yearly_poly1_fit_rmsevarid, [0 0], [len_lon_cmems len_lat_cmems], interped_sla_yearly_poly1_fit_rmse); 
                netcdf.putVar(ncid, interped_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend_filtered);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     interped sea level correlation analysis
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;
                
               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:size(comb_data,3)
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

               % % %   correlation coefficient between model and cmems
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corrcoef(squeeze(comb_interped_data(i,j,:))',squeeze(comb_cmems_data(i,j,:))');
                    corr_interped(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_filtered complete') 

    % % %   correlation coefficient between model and cmems (corrected)
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corrcoef(squeeze(comb_interped_data_corrected(i,j,:))',squeeze(comb_cmems_data(i,j,:))');
                    corr_corrected_interped(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_corrected complete') 

    % % %   correlation coefficient between model_filtered and cmems_filtered 
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corrcoef(squeeze(comb_interped_data_filtered(i,j,:))',squeeze(comb_cmems_data_filtered(i,j,:))');
                    corr_interped_filtered(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_filtered complete') 

    % % %   correlation coefficient between model_climatology and cmems_climatology 
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corrcoef(squeeze(comb_spatial_meaninterped(i,j,:))',squeeze(comb_spatial_meancmems(i,j,:))');
                    corr_spatial_mean(i,j)=temp_corr(1,2);
                end
            end
            disp('corr coef_spatial_mean complete')         

            % % %   correlation coefficient between climatological ssh and climatological cmems ssh 
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    for k=1:12
                        temp_corr=corrcoef(squeeze(comb_interped_clim_divided(i,j,k,:))',squeeze(comb_cmems_clim_divided(i,j,k,:))');
                        corr_clim(i,j,k)=temp_corr(1,2);
                    end
                end
            end
            disp('corr coef_clim complete') 


            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');   

                corr_interpedvarid=netcdf.defVar(ncid, 'corr_interped', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_interpedvarid,'long_name','corr_interped');
                netcdf.putAtt(ncid,corr_interpedvarid,'units',' ');

                corr_interped_filteredvarid=netcdf.defVar(ncid, 'corr_interped_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_interped_filteredvarid,'long_name','corr_interped_filtered');
                netcdf.putAtt(ncid,corr_interped_filteredvarid,'units',' ');

                corr_corrected_interpedvarid=netcdf.defVar(ncid, 'corr_corrected_interped', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_corrected_interpedvarid,'long_name','corr_corrected_interped');
                netcdf.putAtt(ncid,corr_corrected_interpedvarid,'units',' ');

                corr_spatial_meanvarid=netcdf.defVar(ncid, 'corr_spatial_mean', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_spatial_meanvarid,'long_name','corr_spatial_mean');
                netcdf.putAtt(ncid,corr_spatial_meanvarid,'units',' ');

                corr_climvarid=netcdf.defVar(ncid, 'corr_clim', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,corr_climvarid,'long_name','corr_clim');
                netcdf.putAtt(ncid,corr_climvarid,'units',' ');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                netcdf.putVar(ncid, corr_interpedvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped);
                netcdf.putVar(ncid, corr_interped_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_filtered);
                netcdf.putVar(ncid, corr_corrected_interpedvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_corrected_interped);
                netcdf.putVar(ncid, corr_spatial_meanvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_spatial_mean);
                netcdf.putVar(ncid, corr_climvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], corr_clim);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     low pass filtered interped sea level trend analysis
        fig_flag=fig_flags{6,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;
                cmems_lonsize_cut=len_lon_cmems;
                cmems_latsize_cut=len_lat_cmems;
                
               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:length(inputyear)*12
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

               % %        2, 3, 5, 10year lowpass filter 
                sample_freq=1;   %sampling frequency
                filt_order=5;  %filtering order
                nq_freq=sample_freq/2;  %Nyquist frequency
                ftype='low';  %filter type 

                nyears=[1,2,3,4,5];
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    cutoff_freq=1/(12*nyear);
                    [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                         eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_interped_data(i,j,:)-mean(comb_interped_data(i,j,:)))+mean(comb_interped_data(i,j,:));'])
                            eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_cmems_data(i,j,:)-mean(comb_cmems_data(i,j,:)))+mean(comb_cmems_data(i,j,:));'])
%                             eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
%                             eval(['comb_interped_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
%                             eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
%                             eval(['comb_cmems_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end
                
%                 for l=1:288
%                     tempmsl=interped_sla(:,:,l);
%                     msl(l)=mean(tempmsl(:),'omitnan');
%         %             tempmsl_lp=comb_interped_sla_2y_lowpass(:,:,l);
%                     tempmsl_lp=comb_interped_5y_lowpass(:,:,l);
%                     msl_lp(l)=mean(tempmsl_lp(:),'omitnan');
%                 end
%                 plot(msl,'k');
%                 hold on;
%                 plot(msl_lp-mean(msl_lp(:), 'omitnan'),'r');
%                 hold off;
        
                % % %          sea level anomaly low pass filter
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    cutoff_freq=1/(12*nyear);
                    [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,interped_sla(i,j,:)-mean(interped_sla(i,j,:)))+mean(interped_sla(i,j,:));'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,cmems_sla(i,j,:)-mean(cmems_sla(i,j,:)))+mean(cmems_sla(i,j,:));'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end

                % % %          corrected sea level anomaly low pass filter
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    cutoff_freq=1/(12*nyear);
                    [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,corrected_interped_sla(i,j,:)-mean(corrected_interped_sla(i,j,:)))+mean(corrected_interped_sla(i,j,:));'])
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end       

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=[nc_varname_prefix, num2str(nyear), 'y_lowpass'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid time_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', 'm', '''', ');']);
                    end
                end

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_', nc_varname, ');']);
                    end
                end

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     low pass filtered interped sea level correlation analysis
        fig_flag=fig_flags{7,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:size(comb_data,3)
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

            % % %   correlation coefficient between model_low_passed and cmems_low_passed
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_',num2str(nyear),'y_lowpass(i,j,:)));']);
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');']);
                            eval(['numnan_interped=sum(~isnan(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['numnan_cmems=sum(~isnan(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['tlen=length(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:));']);
                            if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                                eval(['corr_interped_',num2str(nyear),'y_lowpass(i,j)=NaN;'])
                            else
                                eval(['corr_interped_',num2str(nyear),'y_lowpass(i,j)=temp_corr(1,2);'])
                            end
                        end
                    end
        %             eval(['m_corr_',num2str(nyear),'y_lowpass=mean(corr_interped_',num2str(nyear),'y_lowpass(:),','''','omitnan','''',');']);
                end
                disp('corr coef_low_passed complete')    

                % % %   correlation coefficient between model_low_passed and cmems_low_passed (corrected)
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_',num2str(nyear),'y_lowpass(i,j,:)));']);
                            eval(['temp_corr=corrcoef(squeeze(comb_corrected_interped_sla_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');']);
                            eval(['numnan_interped=sum(~isnan(comb_corrected_interped_sla_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['numnan_cmems=sum(~isnan(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['tlen=length(comb_cmems_', num2str(nyear), 'y_lowpass(i,j,:));']);
                            if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                                eval(['corr_corrected_interped_sla_',num2str(nyear),'y_lowpass(i,j)=NaN;'])
                            else
                                eval(['corr_corrected_interped_sla_',num2str(nyear),'y_lowpass(i,j)=temp_corr(1,2);'])
                            end
                        end
                    end
        %             eval(['m_corr_',num2str(nyear),'y_lowpass=mean(corr_interped_',num2str(nyear),'y_lowpass(:),','''','omitnan','''',');']);
                end
                disp('corr coef_low_passed_corrected complete')     

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=['corr_interped_',num2str(nyear),'y_lowpass'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
                end

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=['corr_corrected_interped_sla_',num2str(nyear),'y_lowpass'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
                end

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));


                nc_varname_prefixes={'corr_interped_', 'corr_corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix,num2str(nyear),'y_lowpass'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0], [len_lon_cmems len_lat_cmems],', nc_varname, ');']);
                    end
                end
                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     detrended sea level correlation analysis
        fig_flag=fig_flags{8,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'detrended_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
               cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        comb_interped_detrended(i,j,:)=comb_interped_data_filtered(i,j,:)-mean(comb_interped_data_filtered(i,j,:));
                        ppp=polyfit(xData,squeeze(comb_interped_detrended(i,j,:))',1);
                        comb_interped_linear(i,j,:)=xData*ppp(1)+ppp(2);
                        comb_interped_detrended(i,j,:)=comb_interped_detrended(i,j,:)-comb_interped_linear(i,j,:);
                        comb_cmems_detrended(i,j,:)=comb_cmems_data_filtered(i,j,:)-mean(comb_cmems_data_filtered(i,j,:));
                        ppp=polyfit(xData,squeeze(comb_cmems_detrended(i,j,:))',1);
                        comb_cmems_linear(i,j,:)=xData*ppp(1)+ppp(2);
                        comb_cmems_detrended(i,j,:)=comb_cmems_detrended(i,j,:)-comb_cmems_linear(i,j,:);
                    end
                end

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        temp_corr=corrcoef(squeeze(comb_interped_detrended(i,j,:))',squeeze(comb_cmems_detrended(i,j,:))');
                        corr_interped_detrended(i,j)=temp_corr(1,2);
                    end
                end
                disp('corr coef_detrended complete') 

                % %        2, 3, 5, 10year detrended data lowpass filter 

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    cutoff_freq=1/(12*nyear);
                    [filt_coef_b,filt_coef_a]=butter(filt_order,cutoff_freq/nq_freq, ftype); %butterworth filter
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_interped_detrended(i,j,:));'])
                            eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,:)=filter(filt_coef_b,filt_coef_a,comb_cmems_detrended(i,j,:));'])
                            eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_detrended_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                            eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_cmems_detrended_',num2str(nyear),'y_lowpass(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end

                % % %   correlation coefficient between detrended model_low_passed and detrended cmems_low_passed
        %         nyears=[2,3,5];
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));'])
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),squeeze(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
                            eval(['numnan_interped=sum(~isnan(comb_interped_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['numnan_cmems=sum(~isnan(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:)));']);
                            eval(['tlen=length(comb_cmems_detrended_', num2str(nyear), 'y_lowpass(i,j,:));']);
                            if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                                eval(['corr_interped_detrended_',num2str(nyear),'y_lowpass(i,j)=NaN;'])
                            else
                                eval(['corr_interped_detrended_',num2str(nyear),'y_lowpass(i,j)=temp_corr(1,2);'])
                            end
                        end
                    end 
        %             eval(['m_corr_',num2str(nyear),'y_lowpass=mean(corr_interped_',num2str(nyear),'y_lowpass(:),','''','omitnan','''',');']);
                end
                disp('corr detrended coef_low_passed complete') 


            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                nc_varname_prefixes={'cmems_detrended_', 'interped_detrended_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=[nc_varname_prefix, num2str(nyear), 'y_lowpass'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid time_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', 'm', '''', ');']);
                    end
                end

                corr_interped_detrendedvarid=netcdf.defVar(ncid, 'corr_interped_detrended', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_interped_detrendedvarid,'long_name','corr_interped_detrended');
                netcdf.putAtt(ncid,corr_interped_detrendedvarid,'units',' ');

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=['corr_interped_detrended_',num2str(nyear),'y_lowpass'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
                end

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                nc_varname_prefixes={'cmems_detrended_', 'interped_detrended_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_lowpass'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_', nc_varname, ');']);
                    end
                end
                nc_varname_prefixes={'corr_interped_detrended_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix,num2str(nyear),'y_lowpass'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0], [len_lon_cmems len_lat_cmems],', nc_varname, ');']);
                    end
                end

                netcdf.putVar(ncid, corr_interped_detrendedvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_detrended);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end   

% % %     moving averaged interped sea level trend analysis
        fig_flag=fig_flags{9,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;
                cmems_lonsize_cut=len_lon_cmems;
                cmems_latsize_cut=len_lat_cmems;
                
               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:length(inputyear)*12
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

               % %        2, 3, 5, 10year moving average
                sample_freq=1;   %sampling frequency
                filt_order=5;  %filtering order
                nq_freq=sample_freq/2;  %Nyquist frequency
                ftype='low';  %filter type 

                nyears=[1,2,3,4,5];
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                         eval(['comb_interped_',num2str(nyear),'y_movmean(i,j,:)=movmean(comb_interped_data(i,j,:)-mean(comb_interped_data(i,j,:)),',num2str(12*nyear), ')+mean(comb_interped_data(i,j,:));'])
                            eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,:)=movmean(comb_cmems_data(i,j,:)-mean(comb_cmems_data(i,j,:)),', num2str(12*nyear), ')+mean(comb_cmems_data(i,j,:));'])
                            eval(['comb_interped_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
                            eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end
                
% %                 for l=1:288
% %                     tempmsl=interped_sla(:,:,l);
% %                     msl(l)=mean(tempmsl(:),'omitnan');
% %         %             tempmsl_lp=comb_interped_sla_2y_lowpass(:,:,l);
% %                     tempmsl_lp=comb_interped_2y_movmean(:,:,l);
% %                     msl_lp(l)=mean(tempmsl_lp(:),'omitnan');
% %                 end
% %                 plot(msl,'k');
% %                 hold on;
% %                 plot(msl_lp-mean(msl_lp(:), 'omitnan'),'r');
% %                 hold off;
        
                % % %          sea level anomaly moving average
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,:)=movmean(interped_sla(i,j,:)-mean(interped_sla(i,j,:)),',num2str(12*nyear), ')+mean(interped_sla(i,j,:));'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,:)=movmean(cmems_sla(i,j,:)-mean(cmems_sla(i,j,:)),',num2str(12*nyear), ')+mean(cmems_sla(i,j,:));'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end

                % % %          corrected sea level anomaly moving average
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_movmean(i,j,:)=movmean(corrected_interped_sla(i,j,:)-mean(corrected_interped_sla(i,j,:)),',num2str(12*nyear), ')+mean(corrected_interped_sla(i,j,:));'])
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_corrected_interped_sla_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
                        end
                    end
                end       

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=[nc_varname_prefix, num2str(nyear), 'y_movmean'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid time_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', 'm', '''', ');']);
                    end
                end

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                nc_varname_prefixes={'cmems_', 'interped_', 'interped_sla_', 'corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix, num2str(nyear),'y_movmean'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_', nc_varname, ');']);
                    end
                end

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     moving averaged interped sea level correlation analysis
        fig_flag=fig_flags{10,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, 'moving_averaged_cmems_interped_ssh_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(savedir, testname,'_',regionname, 'cmems_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(savedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_ssh');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_ssh=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected ssh
               mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
               diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
               for t=1:size(comb_data,3)
                   comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
               end

               for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
                        corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
               end

            % % %   correlation coefficient between model_low_passed and cmems_low_passed
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_movmean(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_',num2str(nyear),'y_movmean(i,j,:)));']);
                            eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_movmean(i,j,:)),','''', 'rows','''',',','''','complete','''',');']);
                            eval(['numnan_interped=sum(~isnan(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)));']);
                            eval(['numnan_cmems=sum(~isnan(comb_cmems_', num2str(nyear), 'y_movmean(i,j,:)));']);
                            eval(['tlen=length(comb_cmems_', num2str(nyear), 'y_movmean(i,j,:));']);
                            if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                                eval(['corr_interped_',num2str(nyear),'y_movmean(i,j)=NaN;'])
                            else
                                eval(['corr_interped_',num2str(nyear),'y_movmean(i,j)=temp_corr(1,2);'])
                            end
                        end
                    end
        %             eval(['m_corr_',num2str(nyear),'y_movmean=mean(corr_interped_',num2str(nyear),'y_movmean(:),','''','omitnan','''',');']);
                end
                disp('corr coef_low_passed complete')    

                % % %   correlation coefficient between model_low_passed and cmems_low_passed (corrected)
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_movmean(i,j,:)),','''', 'rows','''',',','''','complete','''',');'])
        %                     eval(['temp_corr=corrcoef(squeeze(comb_interped_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_',num2str(nyear),'y_movmean(i,j,:)));']);
                            eval(['temp_corr=corrcoef(squeeze(comb_corrected_interped_sla_', num2str(nyear), 'y_movmean(i,j,:)),squeeze(comb_cmems_', num2str(nyear),'y_movmean(i,j,:)),','''', 'rows','''',',','''','complete','''',');']);
                            eval(['numnan_interped=sum(~isnan(comb_corrected_interped_sla_', num2str(nyear), 'y_movmean(i,j,:)));']);
                            eval(['numnan_cmems=sum(~isnan(comb_cmems_', num2str(nyear), 'y_movmean(i,j,:)));']);
                            eval(['tlen=length(comb_cmems_', num2str(nyear), 'y_movmean(i,j,:));']);
                            if (numnan_interped < tlen-nyear*12 || numnan_cmems < tlen-nyear*12 )
                                eval(['corr_corrected_interped_sla_',num2str(nyear),'y_movmean(i,j)=NaN;'])
                            else
                                eval(['corr_corrected_interped_sla_',num2str(nyear),'y_movmean(i,j)=temp_corr(1,2);'])
                            end
                        end
                    end
        %             eval(['m_corr_',num2str(nyear),'y_movmean=mean(corr_interped_',num2str(nyear),'y_movmean(:),','''','omitnan','''',');']);
                end
                disp('corr coef_low_passed_corrected complete')     

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                clim_timevarid=netcdf.defVar(ncid, 'clim_time', 'NC_DOUBLE', clim_time_dimid);
                netcdf.putAtt(ncid,clim_timevarid,'long_name','clim_time');
                netcdf.putAtt(ncid,clim_timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,clim_timevarid,'calendar','gregorian');

                lon_cmemsvarid=netcdf.defVar(ncid, 'lon_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lon_cmemsvarid,'long_name','longitude');
                netcdf.putAtt(ncid,lon_cmemsvarid,'units','degree_east');

                lat_cmemsvarid=netcdf.defVar(ncid, 'lat_cmems', 'NC_DOUBLE', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,lat_cmemsvarid,'long_name','latitude');
                netcdf.putAtt(ncid,lat_cmemsvarid,'units','degree_north');

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=['corr_interped_',num2str(nyear),'y_movmean'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
                end

                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    nc_varname=['corr_corrected_interped_sla_',num2str(nyear),'y_movmean'];
                    eval([nc_varname,'varid=netcdf.defVar(ncid,','''', ...
                        nc_varname,'''',',', '''', 'NC_FLOAT','''',',', '[lon_cmems_dimid lat_cmems_dimid]);'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'long_name','''', ',', '''', nc_varname,'''',');'])
                    eval(['netcdf.putAtt(ncid,', nc_varname,'varid,', '''', 'units'    ,'''', ',', '''', ' ', '''', ');']);
                end

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));


                nc_varname_prefixes={'corr_interped_', 'corr_corrected_interped_sla_'};
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    for nyeari=1:length(nyears)
                        nyear=nyears(nyeari);
                        nc_varname=[nc_varname_prefix,num2str(nyear),'y_movmean'];
                        eval(['netcdf.putVar(ncid,', nc_varname,'varid, [0 0], [len_lon_cmems len_lat_cmems],', nc_varname, ');']);
                    end
                end
                netcdf.close(ncid);
            end
            fig_flag=0;
        end

    end
end
