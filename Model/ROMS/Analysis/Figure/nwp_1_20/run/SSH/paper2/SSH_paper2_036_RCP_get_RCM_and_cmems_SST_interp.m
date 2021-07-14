close all; clear all;  clc;
% %  get cmemsstructed sst. compare model and resst. save. 

all_region ={'AKP4'}
% all_region ={'YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}

% all_region ={'AKP4'}
% all_region ={'NWP', 'AKP2', 'ES', 'YS', 'NES', 'SES'}
% all_testname = {'test11', 'test12'};
% all_testname = {'test57', 'test58', 'test59', 'test60', 'ens07'};
% all_testname = {'test65', 'test66', 'test67', 'test68'};
% all_testname = {'test61', 'test62', 'test63', 'test64', 'ens09'};
% all_testname = {'test65', 'test66', 'test67', 'test68', 'ens10'};
% all_testname = {'test61', 'test62', 'test63', 'test64', 'ens09', 'test65', 'test66', 'test67', 'test68', 'ens10'};

% all_testname = {'test57', 'test58', 'test59', 'test60', 'test65', 'test66', 'test67', 'test68'};
% all_testname = {'test65', 'test66', 'test67', 'test68'};
% all_testname = {'test53', 'test54', 'test55', 'test56'};
all_testname = {'ens08', 'ens09', 'ens10'};

% all_region ={'AKP4'};
for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname
        % % % 

        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\User\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

        shadlev = [-2 2];
        % rms_shadlev = [0 4];
        trend_shadlev = [0 4];
        % bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/20;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [2006:2100]; % % put year which you want to plot [year year ...];
%         inputyear = [1976:2005]; % % put year which you want to plot [year year ...];
        
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
%         inputmonth = [1]; % % put month which you want to plot [month month ...]

        regionname=all_region{regionind};
        run('nwp_polygon_point.m');
        
% % %         switch region
        regionname=all_region{regionind};
        [error_status, refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);

        switch testname
            case {'test53', 'test54', 'test55', 'test56'}
                drivename='J:\';
            case {'test61', 'test62', 'test63', 'test64'}
                drivename='H:\';
            case {'test57', 'test58', 'test59', 'test60'}
                drivename='I:\';
            case {'test65', 'test66', 'test67', 'test68'}
                drivename='G:\';
            case {'ens08', 'ens09', 'ens10'}
                drivename='F:\';
        end
            
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
%             figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat(drivename, 'Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='get model data and cmemsstructed data';
            fig_flags{2,1}='sea surface temperature analysis';
            fig_flags{3,1}='cmems sea surface temperature trend analysis';
            fig_flags{4,1}='interped sea surface temperature trend analysis';
            fig_flags{5,1}='interped sea surface temperature correlation analysis';
            fig_flags{6,1}='low pass filtered interped sea surface temperature trend analysis';
            fig_flags{7,1}='low pass filtered interped sea surface temperature correlation analysis';
            fig_flags{8,1}='detrended sea surface temperature correlation analysis';
            fig_flags{9,1}='moving averaged interped sea surface temperature trend analysis';
            fig_flags{10,1}='moving averaged interped sea surface temperature correlation analysis';

            for flagi=1:10
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=1;
            fig_flags{2,2}=0;
            fig_flags{3,2}=0;
            fig_flags{4,2}=1;
            fig_flags{5,2}=0;
            fig_flags{6,2}=0;
            fig_flags{7,2}=0;
            fig_flags{8,2}=0;
            fig_flags{9,2}=0;
            fig_flags{10,2}=0;
        end
        
        
        variable = 'SST';
% % %         get model data and cmems data      
        fig_flag=fig_flags{1,2};
        while (fig_flag)
%             run(param_script);
            ind=1;
            matname = [filedir,testname,'_',regionname,'model_cmems_sst_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];

            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                        % read model data
                        if (exist('lon')==0)
                            modelinfo=ncinfo(filename);
                            lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                            lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

                            lon_west = abs(lon - (lonlat(1)-1));
                            min_lon_west=min(lon_west);
                            lon_east = abs(lon - (lonlat(2)+1));
                            min_lon_east=min(lon_east);
                            lat_south = abs(lat - (lonlat(3)-1));
                            min_lat_south=min(lat_south);
                            lat_north = abs(lat - (lonlat(4)+1));
                            min_lat_north=min(lat_north);

                            lon_min = find(lon_west == min_lon_west);
                            lon_max = find(lon_east == min_lon_east);
                            lat_min = find(lat_south == min_lat_south);
                            lat_max = find(lat_north == min_lat_north);

                            lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_model(1:size(lon,1),1:size(lon,2))=1;
                                otherwise
                                    mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model(mask_model==0)=NaN;
                            end
                        end


                        data_info = ncinfo(filename, 'temp');  %% [lon lat depth time] -> [1601 1201 33 1]
                        
                        varname ='temp';
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) modelinfo.Dimensions(2).Length 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
                        data=data.*mask_model;
                        
                        len_lon_model = size(data,1);
                        len_lat_model = size(data,2);

                        if (exist('comb_spatial_meanmodel')==0)
                            comb_spatial_mean_data=(zeros([len_lon_model,len_lat_model,12]));
                            del=(zeros([len_lon_model,len_lat_model,12]));
                        end

                        comb_data(:,:,ind) = single(data);
                        comb_spatial_mean_data(:,:,monthij)=comb_spatial_mean_data(:,:,monthij)+data/double(length(inputyear));

        % % % %                 % read Reconstructed DATA
                        cmemsfilename = strcat(cmemsdir,'cmems_nwp_ssh_', num2str(2006,'%04i'),'.nc');
                        if (exist('cmems_lon')==0)
                            cmemsinfo=ncinfo(cmemsfilename);
                            cmems_lonname='longitude';
                            cmems_latname='latitude';
                            cmems_varname='ssha';
                            cmemsinfo_lon=ncinfo(cmemsfilename,cmems_lonname);
                            cmemsinfo_lat=ncinfo(cmemsfilename,cmems_latname);
                            cmems_lon = ncread(cmemsfilename,cmems_lonname,1,cmemsinfo_lon.Dimensions.Length);
                            cmems_lat = ncread(cmemsfilename,cmems_latname,1,cmemsinfo_lat.Dimensions.Length);
                            [cmems_lat2, cmems_lon2]= meshgrid(cmems_lat, cmems_lon);
                            [cmems_lon_min, cmems_lon_max, cmems_lat_min, cmems_lat_max] = findind_Y(1/20, lonlat(1:4), cmems_lon2, cmems_lat2);
                            cmems_lon2 = cmems_lon2(cmems_lon_min(1):cmems_lon_max(1), cmems_lat_min(1):cmems_lat_max(1));
                            cmems_lat2 = cmems_lat2(cmems_lon_min(1):cmems_lon_max(1), cmems_lat_min(1):cmems_lat_max(1));

                            cmems_lonsize_cut=cmems_lon_max(1)-cmems_lon_min(1)+1;
                            cmems_latsize_cut=cmems_lat_max(1)-cmems_lat_min(1)+1;
                            len_lon=cmems_lonsize_cut;
                            len_lat=cmems_latsize_cut;
        %                     comb_spatial_meanrms=(zeros([length(cmems_lon),length(cmems_lat),12]));
        %                     comb_spatial_meanbias=(zeros([length(cmems_lon),length(cmems_lat),12]));
                            comb_spatial_meancmems=(zeros([cmems_lonsize_cut,cmems_latsize_cut,12]));
                            comb_spatial_meaninterped_data=(zeros([cmems_lonsize_cut,cmems_latsize_cut,12]));
                            
        %                     comb_spatial_meanmodel=(zeros([length(cmems_lon),length(cmems_lat),12]));
        
                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_cmems(1:cmems_lonsize_cut,1:cmems_latsize_cut)=1;
                                otherwise
                                    mask_cmems = double(inpolygon(cmems_lon2,cmems_lat2,refpolygon(:,1),refpolygon(:,2)));
                                    mask_cmems(mask_cmems==0)=NaN;
                            end
                        end

                        interped_data = griddata(double(lon), double(lat), data, double(cmems_lon2),double(cmems_lat2));   
                        
                        comb_interped_data(:,:,ind) = interped_data;
                        comb_spatial_meaninterped_data(:,:,monthij)=comb_spatial_meaninterped_data(:,:,monthij)+interped_data/double(length(inputyear));

                        ind = ind + 1;
                        toc;
                    end
                end
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_data', ...
                    'cmems_lonsize_cut', 'cmems_latsize_cut', 'comb_spatial_mean_data', ...
                    'mask_cmems', 'comb_interped_data', 'comb_spatial_meaninterped_data', ...
                    'cmems_lon2', 'cmems_lat2', 'lon', 'lat', '-v7.3');
            else
                load(matname);
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
    end     
    
% % %         sea surface temperature analysis
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, '_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                raw_sstvarid=netcdf.defVar(ncid, 'raw_sst', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,raw_sstvarid,'long_name','raw_sst');
                netcdf.putAtt(ncid,raw_sstvarid,'units','m');

                sst_filteredvarid=netcdf.defVar(ncid, 'sst_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,sst_filteredvarid,'long_name','sst_filtered');
                netcdf.putAtt(ncid,sst_filteredvarid,'units','m');

                trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,trendvarid,'long_name','trend');
                netcdf.putAtt(ncid,trendvarid,'units','mm/year');

                trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
                netcdf.putAtt(ncid,trend_filteredvarid,'units','mm/year');

                clim_sstvarid=netcdf.defVar(ncid, 'clim_sst', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_sstvarid,'long_name','clim_sst');
                netcdf.putAtt(ncid,clim_sstvarid,'units','m');

                clim_sst_trendvarid=netcdf.defVar(ncid, 'clim_sst_trend', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_sst_trendvarid,'long_name','clim_sst_trend');
                netcdf.putAtt(ncid,clim_sst_trendvarid,'units','m');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
                netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
                netcdf.putVar(ncid, raw_sstvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
                netcdf.putVar(ncid, sst_filteredvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data_filtered);
                netcdf.putVar(ncid, trendvarid, [0 0], [len_lon_model, len_lat_model], trend);
                netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon_model, len_lat_model], trend_filtered);
                netcdf.putVar(ncid, clim_sstvarid, [0 0 0], [len_lon_model, len_lat_model length(climtime)], comb_spatial_meanmodel);
                netcdf.putVar(ncid, clim_sst_trendvarid, [0 0 0], [len_lon_model, len_lat_model length(climtime)], trend_clim);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     cmems sea surface temperature trend analysis
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   

                len_lon_cmems=cmems_lonsize_cut;
                len_lat_cmems=cmems_latsize_cut;

                cmems_sla=comb_cmems_data;
                cmems_sla_divided=reshape(cmems_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
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

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                cmems_slavarid=netcdf.defVar(ncid, 'cmems_sla', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,cmems_slavarid,'long_name','cmems_sla');
                netcdf.putAtt(ncid,cmems_slavarid,'units','m ');

                cmems_sla_filteredvarid=netcdf.defVar(ncid, 'cmems_sla_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,cmems_sla_filteredvarid,'long_name','cmems_sla_filtered');
                netcdf.putAtt(ncid,cmems_sla_filteredvarid,'units','m ');

                clim_cmemsvarid=netcdf.defVar(ncid, 'clim_cmems_sst', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_cmemsvarid,'long_name','clim_cmems_sst');
                netcdf.putAtt(ncid,clim_cmemsvarid,'units','m ');

                cmems_trendvarid=netcdf.defVar(ncid, 'cmems_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_trendvarid,'long_name','cmems_trend');
                netcdf.putAtt(ncid,cmems_trendvarid,'units','mm /year');

                cmems_trend_filteredvarid=netcdf.defVar(ncid, 'cmems_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,cmems_trend_filteredvarid,'long_name','cmems_trend_filtered');
                netcdf.putAtt(ncid,cmems_trend_filteredvarid,'units','mm /year');

                clim_cmems_trend_dividedvarid=netcdf.defVar(ncid, 'clim_cmems_trend_divided', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'long_name','clim_cmems_trend_divided');
                netcdf.putAtt(ncid,clim_cmems_trend_dividedvarid,'units','mm /year');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                netcdf.putVar(ncid, cmems_slavarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data);
                netcdf.putVar(ncid, cmems_sla_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_cmems_data_filtered);
                netcdf.putVar(ncid, clim_cmemsvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], comb_spatial_meancmems);
                netcdf.putVar(ncid, cmems_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend);
                netcdf.putVar(ncid, cmems_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_trend_filtered);
                netcdf.putVar(ncid, clim_cmems_trend_dividedvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(climtime)], clim_cmems_trend_divided);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     interped sea surface temperature trend analysis
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir,  testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                len_lon_cmems=size(cmems_lon2,1);
                len_lat_cmems=size(cmems_lon2,2);
                             
               cmems_sla_mean = zeros(size(cmems_lon2));

               interped_data_divided=reshape(comb_interped_data,[size(cmems_lon2,1), size(cmems_lon2,2), 12, length(inputyear)]);
                clim_interped_data=mean(interped_data_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_data);
                end

                % % %         interped trend 
                comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

                interped_trend(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                        interped_trend(i,j)=p(1);
                    end
                end
                disp('interped trend complete') 

                
                comb_spatial_meaninterped_data=mean(comb_interped_clim_divided, 4, 'omitnan');
        % % %         interped trend (seasonal filtered)
                for t=1:length(inputyear)
                    comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meaninterped_data;
                end

                interped_trend_filtered(1:cmems_lonsize_cut,1:cmems_latsize_cut)=NaN;

                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                        interped_trend_filtered(i,j)=p(1);
                    end
                end
                disp('interped trend_filtered complete') 
               mean_trend_filtered=mean(interped_trend_filtered(:), 'omitnan');
               
               % % %         climatological trend 
                comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);
                climtrendtime=inputyear(1):inputyear(end);
                for i=1:cmems_lonsize_cut
                    for j=1:cmems_latsize_cut
                        for k=1:12  % month
                            p=polyfit(climtrendtime,squeeze(comb_interped_clim_divided(i,j,k,:))',1);
                            interped_trend_clim(i,j,k)=p(1);
                        end
                    end
                end
                disp('climatological interped trend complete')       
                
% %                 % %        get trend corrected sst
% %                mean_cmems_trend_filtered=mean(cmems_trend_filtered(:), 'omitnan');
% %                diff_trend_correction = (mean_cmems_trend_filtered - mean_trend_filtered)/1000.0;
% %                for t=1:size(comb_data,3)
% %                    comb_interped_data_corrected(:,:,t)=comb_interped_data(:,:,t)+(t-1)*diff_trend_correction/12.0;
% %                end
% % 
% %                for sla_i=1:size(cmems_sla,1)
% %                     for sla_j=1:size(cmems_sla,2)
% %                         corrected_interped_sla_mean(sla_i,sla_j)=mean(comb_interped_data_corrected(sla_i,sla_j,:),'omitnan');
% %                         corrected_interped_sla(sla_i,sla_j,:)=comb_interped_data_corrected(sla_i,sla_j,:)-(corrected_interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
% %                     end
% %                end

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_cmems_dimid = netcdf.defDim(ncid, 'lon_cmems', len_lon_cmems);
                lat_cmems_dimid = netcdf.defDim(ncid,'lat_cmems', len_lat_cmems);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
                clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                interped_sstvarid=netcdf.defVar(ncid, 'interped_sst', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_sstvarid,'long_name','interped_sst');
                netcdf.putAtt(ncid,interped_sstvarid,'units','m/s ');
                
                interped_sst_filteredvarid=netcdf.defVar(ncid, 'interped_sst_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_sst_filteredvarid,'long_name','interped_sst_filtered');
                netcdf.putAtt(ncid,interped_sst_filteredvarid,'units','m/s ');

                interped_trendvarid=netcdf.defVar(ncid, 'interped_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_trendvarid,'long_name','interped_trend');
                netcdf.putAtt(ncid,interped_trendvarid,'units','m/s /year');

                interped_trend_filteredvarid=netcdf.defVar(ncid, 'interped_trend_filtered', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,interped_trend_filteredvarid,'long_name','interped_trend_filtered_trend');
                netcdf.putAtt(ncid,interped_trend_filteredvarid,'units','m/s /year');
                
                clim_interped_sst_trendvarid=netcdf.defVar(ncid, 'clim_interped_sst_trend', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid clim_time_dimid]);
                netcdf.putAtt(ncid,clim_interped_sst_trendvarid,'long_name','clim_interped_sst_trend');
                netcdf.putAtt(ncid,clim_interped_sst_trendvarid,'units','m/s');

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lon2(:));
                netcdf.putVar(ncid, lat_cmemsvarid, [0 0], [len_lon_cmems len_lat_cmems], cmems_lat2(:));

                netcdf.putVar(ncid, interped_sstvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], comb_interped_data);
                netcdf.putVar(ncid, interped_sst_filteredvarid, [0 0 0], [len_lon_cmems len_lat_cmems length(ftime)], interped_data_filtered);
                netcdf.putVar(ncid, interped_trendvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend);
                netcdf.putVar(ncid, interped_trend_filteredvarid, [0 0], [len_lon_cmems len_lat_cmems], interped_trend_filtered);
                netcdf.putVar(ncid, clim_interped_sst_trendvarid, [0 0 0], [len_lon_cmems len_lat_cmems 12], interped_trend_clim);

                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %     interped sea surface temperature correlation analysis
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;
                
               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected sst
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

            % % %   correlation coefficient between climatological sst and climatological cmems sst 
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

% % %     low pass filtered interped sea surface temperature trend analysis
        fig_flag=fig_flags{6,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
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

               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected sst
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
        
                % % %          sea surface temperature anomaly low pass filter
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

                % % %          corrected sea surface temperature anomaly low pass filter
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

% % %     low pass filtered interped sea surface temperature correlation analysis
        fig_flag=fig_flags{7,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'low_pass_filtered_cmems_interped_sst_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected sst
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

% % %     detrended sea surface temperature correlation analysis
        fig_flag=fig_flags{8,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'detrended_cmems_interped_sst_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
               cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

% % %     moving averaged interped sea surface temperature trend analysis
        fig_flag=fig_flags{9,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                %                cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%                 lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=size(cmems_lon2,1);
                len_lat_cmems=size(cmems_lon2,2);

               cmems_sla_mean = zeros(size(cmems_lon2));
               
                cmems_lonsize_cut=len_lon_cmems;
                cmems_latsize_cut=len_lat_cmems;


               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_lon2,1)
                    for sla_j=1:size(cmems_lon2,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_lon2,1), size(cmems_lon2,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
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
%                             eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,:)=movmean(comb_cmems_data(i,j,:)-mean(comb_cmems_data(i,j,:)),', num2str(12*nyear), ')+mean(comb_cmems_data(i,j,:));'])
                            eval(['comb_interped_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
%                             eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
%                             eval(['comb_cmems_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
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
        
                % % %          sea surface temperature anomaly moving average
                for nyeari=1:length(nyears)
                    nyear=nyears(nyeari);
                    for i=1:cmems_lonsize_cut
                        for j=1:cmems_latsize_cut
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,:)=movmean(interped_sla(i,j,:)-mean(interped_sla(i,j,:)),',num2str(12*nyear), ')+mean(interped_sla(i,j,:));'])
%                             eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,:)=movmean(cmems_sla(i,j,:)-mean(cmems_sla(i,j,:)),',num2str(12*nyear), ')+mean(cmems_sla(i,j,:));'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
                            eval(['comb_interped_sla_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
%                             eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,1:6*nyear)=NaN;'])
%                             eval(['comb_cmems_sla_',num2str(nyear),'y_movmean(i,j,end-6*nyear+1:end)=NaN;'])
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                nc_varname_prefixes={'interped_', 'interped_sla_'};
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

                nc_varname_prefixes={'interped_', 'interped_sla_'};
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

% % %     moving averaged interped sea surface temperature correlation analysis
        fig_flag=fig_flags{10,2};
        while (fig_flag)
            ncoutfilename = strcat(filedir, testname,'_',regionname, 'moving_averaged_cmems_interped_sst_corr_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                cmemsfilename = strcat(filedir, testname,'_',regionname, 'cmems_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                lon_info=ncinfo(cmemsfilename, 'lon_cmems');
                len_lon_cmems=lon_info.Dimensions(1).Length;
                len_lat_cmems=lon_info.Dimensions(2).Length;

               cmems_sla=ncread(cmemsfilename, 'cmems_sla');
               cmems_sla_mean = mean(cmems_sla,3);
               cmems_trend_filtered=ncread(cmemsfilename, 'cmems_trend_filtered');
               comb_cmems_data=ncread(cmemsfilename, 'cmems_sla');
               comb_cmems_data_filtered=ncread(cmemsfilename, 'cmems_sla_filtered');
               comb_cmems_clim_divided=reshape(comb_cmems_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_sst_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
               interped_trend_filtered= ncread(interpedfilename, 'interped_trend_filtered');
               mean_trend_filtered= mean(interped_trend_filtered(:), 'omitnan');
               comb_interped_data_filtered = ncread(interpedfilename, 'interped_sla_filtered');
               comb_interped_data = ncread(interpedfilename, 'interped_sst');
               comb_interped_clim_divided=reshape(comb_interped_data,[cmems_lonsize_cut, cmems_latsize_cut, 12, length(inputyear)]);

               interped_sst=comb_interped_data;
                for sla_i=1:size(cmems_sla,1)
                    for sla_j=1:size(cmems_sla,2)
                        interped_sla_mean(sla_i,sla_j)=mean(interped_sst(sla_i,sla_j,:),'omitnan');
                        interped_sla(sla_i,sla_j,:)=interped_sst(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
                    end
                end
                interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
                clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
                for t=1:length(inputyear)
                    interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
                end

                % %        get trend corrected sst
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
                    'type', ['NWP 1/20 _ ', testname, 'model, cmems monthly sst analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' monthly sst analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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
        
%         %%% get sea surface temperature anomaly
%         cmems_sla = comb_cmems_data;
%         for sla_i=1:size(cmems_sla,1)
%             for sla_j=1:size(cmems_sla,2)
%                 cmems_sla_mean(sla_i,sla_j)=mean(cmems_sla(sla_i,sla_j,:),'omitnan');
%             end
%         end
        
%         for l=1:288
%             tempmsl=interped_sla(:,:,l);
%             msl(l)=mean(tempmsl(:),'omitnan');
% %             tempmsl_lp=comb_interped_sla_2y_lowpass(:,:,l);
%             tempmsl_lp=comb_interped_2y_lowpass(:,:,l);
%             msl_lp(l)=mean(tempmsl_lp(:),'omitnan');
%         end
%         plot(msl,'k');
%         hold on;
%         plot(msl_lp,'r');
%         hold off;

% 
% % % %         mean trend (raw, seasonal filtered)
%         mean_trend=mean(trend(:),'omitnan');
%         mean_trend_filtered=mean(trend_filtered(:),'omitnan');
        
        
% % %   correlation coefficient between model_sst_filtered_detrended and cmems_filtered_detrended
        
        

        
% %         figure;
% %         pcolor(corr_clim(:,:,2)')
% %         shading flat
% %         colorbar
% %         disp('corr coef_detrended complete') 
             
%         for l=1:300
%             tempm=squeeze(comb_interped_data_filtered(:,:,l));
%             m_interped(l)=mean(tempm(:),'omitnan');
%             tempm=squeeze(comb_cmems_data_filtered(:,:,l));
%             m_cmems(l)=mean(tempm(:),'omitnan');
%         end
%         corrcoef(m_interped, m_cmems)
%         pcolor(corr_interped')
%         shading flat
%         colorbar;

%         save([filedir,testname,'_',regionname,'sst_trend_corr_coef_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'], '-v7.3');

    end
end
