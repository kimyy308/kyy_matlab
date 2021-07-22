close all; clear all;  clc;
% %  get cmemsstructed SSH. compare model and reSSH. save. 

% all_region ={'AKP4'}
% all_region ={'YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}

% all_region ={'NWP', 'AKP4'}
all_region ={'AKP4'}
% all_testname = {'test11', 'test12'};
% all_testname = {'test2102'};
all_testname = {'test2102', 'test2103', 'test2104', 'test2105', 'test2106'};

% all_testname = {'ens03'};

% all_region ={'AKP4'};
for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname
        % % % 
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\USER\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([dropboxpath '\source\matlab\Common\seawater_ver3_2']));
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
        dl=1/10;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1993:2014]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
        variable = 'zeta';
        varname ='zeta';
        regionname=all_region{regionind}
        run('nwp_polygon_point.m');
        
% % %         switch region
        for folding=1:1
            switch(regionname)
                case('NWP') %% North western Pacific
                    lonlat = [115, 164, 15, 52];  %% whole data area
                    refpolygon(1,1)=lonlat(1);
                    refpolygon(2,1)=lonlat(2);
                    refpolygon(1,2)=lonlat(3);
                    refpolygon(2,2)=lonlat(4);
                case('NWP2') %% North western Pacific
                    lonlat = [115, 145, 25, 52];  %% whole data area
                    refpolygon(1,1)=lonlat(1);
                    refpolygon(2,1)=lonlat(2);
                    refpolygon(1,2)=lonlat(3);
                    refpolygon(2,2)=lonlat(4);
                case('ES') %% East Sea
                    refpolygon=espolygon;
                case('NES') %% Northern East Sea
                    refpolygon=nespolygon;
                case('SES') %% Southern East Sea
                    refpolygon=sespolygon;
                case('SS') %% South Sea
                    refpolygon=sspolygon;
                case('YS') %% Yellow Sea
                    refpolygon=yspolygon;
                case('ECS') %% East China Sea
                    refpolygon=ecspolygon;
                case('ECS2') %% East China Sea2
                    refpolygon=ecs2polygon;
                case('YSECS') %% YS & East China Sea
                    refpolygon=ysecspolygon;
                case('AKP') %% Around Korea Peninsula
                    refpolygon=akppolygon;
                case('AKP2') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('AKP3') %% Around Korea Peninsula
                    refpolygon=akp3polygon;
                case('AKP4') %% Around Korea Peninsula
                    refpolygon=akp4polygon;
                case('CA') %% Coastal Area around korea peninsula
                    refpolygon=capolygon;
                case('EKB') %% Coastal Area around korea peninsula
                    refpolygon=ekbpolygon;
                case('BOH') %% Coastal Area around korea peninsula
                    refpolygon=bohpolygon;
                case('TEST') %% for debugging
                    refpolygon=testpolygon;
                otherwise
                    ('?')
            end
                lonlat(1)=min(refpolygon(:,1));
                lonlat(2)=max(refpolygon(:,1));
                lonlat(3)=min(refpolygon(:,2));
                lonlat(4)=max(refpolygon(:,2));
        end
        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
%             figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
%             filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\'); % % where data files are
            monfiledir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\');
            griddir = strcat('E:\Data\Model\ROMS\nwp_1_20\output\', testname, '\run\packed_monthly\'); % % where data files are
            savedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='get model data and cmemsstructed data';
            fig_flags{2,1}='sea level analysis';
            fig_flags{3,1}='cmems sea level trend analysis';
            fig_flags{4,1}='interped sea level trend analysis';
            fig_flags{5,1}='interped sea level correlation analysis';
            fig_flags{6,1}='low pass filtered interped sea level trend analysis';
            fig_flags{7,1}='low pass filtered interped sea level correlation analysis';
            fig_flags{8,1}='detrended sea level correlation analysis';
            fig_flags{9,1}='moving averaged interped sea level trend analysis';
            fig_flags{10,1}='moving averaged interped sea level correlation analysis';
            fig_flags{11,1}='steric sea level analysis';
            fig_flags{12,1}='steric sea level data nc save';
            fig_flags{13,1}='RCM and GCM steric sea level interped data nc save';
            fig_flags{14,1}='thermosteric sea level analysis';
            fig_flags{15,1}='steric sea level data nc save';
            fig_flags{16,1}='halosteric sea level analysis';
            
            for flagi=1:20
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=1;
            fig_flags{2,2}=1;
            fig_flags{3,2}=1;
            fig_flags{4,2}=1;
            fig_flags{5,2}=1;
            fig_flags{6,2}=0;
            fig_flags{7,2}=0;
            fig_flags{8,2}=0;
            fig_flags{9,2}=0;
            fig_flags{10,2}=0;
            fig_flags{11,2}=1;
            fig_flags{12,2}=1;
            fig_flags{13,2}=0;
            fig_flags{14,2}=1;
            fig_flags{15,2}=1;
            fig_flags{16,2}=1;
            fig_flags{17,2}=1;
        end
        
% % %         get model data and cmemsstructed data      
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            run(param_script);
            ind=1;
            matname = [savedir,testname,'_',regionname,'model_cmems_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];

            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
%                         filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
%                                 'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                        filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_zeta_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                        gridname = strcat(griddir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                        % read model data
                        if (exist('lon')==0)
                            modelinfo=ncinfo(gridname);
                            lon = ncread(gridname,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                            lat = ncread(gridname,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

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

                            lon = ncread(gridname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            lat = ncread(gridname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_model(1:size(lon,1),1:size(lon,2))=1;
                                otherwise
                                    mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model(mask_model==0)=NaN;
                            end
                        end


                        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        data=data.*mask_model;

                        len_lon_model = size(data,1);
                        len_lat_model = size(data,2);

                        if (exist('comb_spatial_meanmodel')==0)
                            comb_spatial_meanmodel=(zeros([len_lon_model,len_lat_model,12]));
                            del=(zeros([len_lon_model,len_lat_model,12]));
                        end

                        comb_data(:,:,ind) = single(data);

                        comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+data/double(length(inputyear));



        % % % %                 % read Reconstructed DATA
                        cmemsfilename = strcat(cmemsdir,'cmems_nwp_ssh_', num2str(tempyear,'%04i'),'.nc');
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
                            comb_spatial_meaninterped=(zeros([cmems_lonsize_cut,cmems_latsize_cut,12]));

        %                     comb_spatial_meanmodel=(zeros([length(cmems_lon),length(cmems_lat),12]));

                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_cmems(1:cmems_lonsize_cut,1:cmems_latsize_cut)=1;
                                otherwise
                                    mask_cmems = double(inpolygon(cmems_lon2,cmems_lat2,refpolygon(:,1),refpolygon(:,2)));
                                    mask_cmems(mask_cmems==0)=NaN;
                            end
                        end

                        if tempmonth==1
                            cmems_daily_adt=ncread(cmemsfilename,'adt',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 31]);
                            cmems_daily_data=ncread(cmemsfilename,'sla',[cmems_lon_min(1) cmems_lat_min(1) 1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 31]);
                        else
                            cmems_daily_adt=ncread(cmemsfilename,'adt',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 eomday(tempyear,tempmonth)]);
                            cmems_daily_data=ncread(cmemsfilename,'sla',[cmems_lon_min(1) cmems_lat_min(1) sum(eomday(tempyear,1:tempmonth-1))+1], [cmems_lon_max(1)-cmems_lon_min(1)+1 cmems_lat_max(1)-cmems_lat_min(1)+1 eomday(tempyear,tempmonth)]);
                        end
                        cmems_adt = mean(cmems_daily_adt,3,'omitnan');
                        cmems_data = mean(cmems_daily_data,3,'omitnan');
            %             cmems_data = ncread(cmemsfilename,varname,[cmems_lon_min(1) cmems_lat_min(1) tempmonth], [cmems_lon_max(1)-cmems_lon_min(1) cmems_lat_max(1)-cmems_lat_min(1) 1]);
                        cmems_data(cmems_data<-1000)=NaN;
                        cmems_data(cmems_data>1000)=NaN;
                        cmems_data=cmems_data.*mask_cmems;

                        comb_cmems_adt(:,:,ind) = cmems_adt;
                        comb_cmems_data(:,:,ind) = cmems_data;
                        comb_spatial_meancmems(:,:,monthij)=comb_spatial_meancmems(:,:,monthij)+cmems_adt/double(length(inputyear));        %                 comb_spatial_meancmems(:,:,monthij)=comb_spatial_meancmems(:,:,monthij)+cmems_adt/double(length(inputyear));

                        interped_data = griddata(double(lon), double(lat), data,double(cmems_lon2),double(cmems_lat2));   

                        comb_interped_data(:,:,ind) = interped_data;
                        comb_spatial_meaninterped(:,:,monthij)=comb_spatial_meaninterped(:,:,monthij)+interped_data/double(length(inputyear));

                        ind = ind + 1;
                        toc;
                    end
                end
                if (exist(strcat(savedir) , 'dir') ~= 7)
                    mkdir(strcat(savedir));
                end 
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_spatial_meanmodel', 'comb_data', ...
                    'cmems_lonsize_cut', 'cmems_latsize_cut', 'comb_spatial_meancmems', 'comb_spatial_meaninterped', ...
                    'mask_cmems', 'comb_cmems_data', 'comb_cmems_adt', 'comb_interped_data', 'comb_spatial_meaninterped', ...
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
    
    % % %   correlation coefficient between model and cmems (Spearman)
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corr(squeeze(comb_interped_data(i,j,:)),squeeze(comb_cmems_data(i,j,:)), 'Type', 'Spearman');
                    corr_interped_spearman(i,j)=temp_corr;
                end
            end
            disp('corr coef_spearman complete') 
            
    % % %   correlation coefficient between model and cmems (Kendall)
            for i=1:cmems_lonsize_cut
                for j=1:cmems_latsize_cut
                    temp_corr=corr(squeeze(comb_interped_data(i,j,:)),squeeze(comb_cmems_data(i,j,:)), 'Type', 'Kendall');
                    corr_interped_kendall(i,j)=temp_corr;
                end
            end
            disp('corr coef_kendall complete') 

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
                
                corr_interped_spearmanvarid=netcdf.defVar(ncid, 'corr_interped_spearman', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_interped_spearmanvarid,'long_name','corr_interped_spearman');
                netcdf.putAtt(ncid,corr_interped_spearmanvarid,'units',' ');
                
                corr_interped_kendallvarid=netcdf.defVar(ncid, 'corr_interped_kendall', 'NC_FLOAT', [lon_cmems_dimid lat_cmems_dimid]);
                netcdf.putAtt(ncid,corr_interped_kendallvarid,'long_name','corr_interped_kendall');
                netcdf.putAtt(ncid,corr_interped_kendallvarid,'units',' ');

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
                netcdf.putVar(ncid, corr_interped_spearmanvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_spearman);
                netcdf.putVar(ncid, corr_interped_kendallvarid, [0 0], [len_lon_cmems len_lat_cmems], corr_interped_kendall);                
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
        
% % %         get ts & calculate zosto
        fig_flag=fig_flags{11,2};
        while (fig_flag)
            run(param_script);
            ind=1;
            matname = [savedir,testname,'_',regionname,'model_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    if (exist('rho_0_va')==0)
                        for monthij = 1:length(inputmonth)
                            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                            tic;
                            tempyear = inputyear(yearij);
                            tempmonth = inputmonth(monthij);
                            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
%                             filename = strcat(filedir, num2str(tempyear,'%04i'), '\', ...
%                                     testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                            filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
%                             gridname = strcat(griddir, num2str(tempyear,'%04i'), '\', ...
%                                 'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                                
                            % read model data
                            if (exist('lon_min')==0)
                                modelinfo=ncinfo(filename);
                                lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                                lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
                                
                                [lon_min, lon_max, lat_min, lat_max] = ...
                                    Func_0012_findind_Y(1, lonlat, lon, lat);
%                                 
%                                 lon_west = abs(lon - (lonlat(1)-1));
%                                 [min_lon_west, lon_min]=min(lon_west);
%                                 lon_east = abs(lon - (lonlat(2)+1));
%                                 [min_lon_east, lon_max]=min(lon_east);
%                                 lat_south = abs(lat - (lonlat(3)-1));
%                                 [min_lat_south, lat_min]=min(lat_south);
%                                 lat_north = abs(lat - (lonlat(4)+1));
%                                 [min_lat_north, lat_max]=min(lat_north);

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

                            data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                            if (exist('h')==0)
                                h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                dA=1./pm .* 1./pn;
                                Vtransform = ncread(filename, 'Vtransform');
                                Vstretching = ncread(filename, 'Vstretching');
                                theta_s = ncread(filename, 'theta_s');
                                theta_b = ncread(filename, 'theta_b');
                                hc = ncread(filename, 'hc');
                                N = length(ncread(filename, 's_rho'));
                            end
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                            PT_src(:,:,:,monthij) = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(:,:,:,monthij) = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(S_src<=0)=0;

                            vtype ='w';
                            zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                            dz = diff(zw);
                            dz = permute(dz, [2,3,1]);
    %                         vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
    %                         area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
                            vtype ='r';
                            P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                            P = -permute(P, [2,3,1]);
                            comb_P(:,:,:,monthij)=P;
                        end
                            PT_src = mean(PT_src, 4 );
                            S_src = mean(S_src, 4 );
                            P = mean(comb_P, 4);
                            T_src = sw_temp(S_src, PT_src, P, zeros(size(P)) );
                            rho_src = sw_dens(S_src, T_src, P );
                            dep_0 = h + data;
                            rho_0_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_0;       % vertical(local) average, 2-D
                    end
                    
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
%                         filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
%                                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                        filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                        if (exist('h')==0)
                            h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            dA=1./pm .* 1./pn;
                            Vtransform = ncread(filename, 'Vtransform');
                            Vstretching = ncread(filename, 'Vstretching');
                            theta_s = ncread(filename, 'theta_s');
                            theta_b = ncread(filename, 'theta_b');
                            hc = ncread(filename, 'hc');
                            N = length(ncread(filename, 's_rho'));
                        end
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        PT_src = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                        S_src = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                        S_src(S_src<=0)=0;

                        vtype ='w';
                        zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                        dz = diff(zw);
                        dz = permute(dz, [2,3,1]);

                        vtype ='r';
                        P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                        P = -permute(P, [2,3,1]);
                        T_src = sw_temp (S_src, PT_src, P, zeros(size(P)) );
                        rho_src = sw_dens (S_src, T_src, P );

                        dep_n = h + data;
                        rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
%                         rho_n_va2 = sum( rho_src.*dz, 3 , 'omitnan') ./ (h+data);

%                         rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
                        zosto = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
%                         zosto_2 = dep_n       .* ( 1 - rho_n_va ./ rho_0_va );
                        comb_zosto(:,:,ind) = zosto;

                        len_lon_model = size(data,1);
                        len_lat_model = size(data,2);

                        GLO_zostoname='D:\Data\Model\CMIP5\zosto\historical\zosto_interp_IPSL-CM5A-LR_historical_1976-2005.nc';

                        if (exist('lon_glo')==0)
                            GLO_modelinfo=ncinfo(GLO_zostoname);
                            lon_glo = ncread(GLO_zostoname,'lon',[1 1],[inf,1]);
                            lat_glo = ncread(GLO_zostoname,'lat',[1 1],[1,inf]);

                            lon_west_glo = abs(lon_glo - (lonlat(1)-1));
                            [min_lon_west_glo, lon_min_glo]=min(lon_west_glo);
                            lon_east_glo = abs(lon_glo - (lonlat(2)+1));
                            [min_lon_east_glo, lon_max_glo]=min(lon_east_glo);
                            lat_south_glo = abs(lat_glo - (lonlat(3)-1));
                            [min_lat_south_glo, lat_min_glo]=min(lat_south_glo);
                            lat_north_glo = abs(lat_glo - (lonlat(4)+1));
                            [min_lat_north_glo, lat_max_glo]=min(lat_north_glo);

                            lon_glo = ncread(GLO_zostoname,'lon', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
                            lat_glo = ncread(GLO_zostoname,'lat', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);

                            switch(regionname)
                                otherwise
                                    mask_model_glo = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model_glo(mask_model_glo==0)=NaN;
                            end
                        end

                        interped_zosto = griddata(double(lon), double(lat), zosto,double(lon_glo),double(lat_glo));   

                        comb_interped_zosto(:,:,ind) = interped_zosto;

                        ind = ind + 1;
                        toc;
                    end
                end
                
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_zosto',  ...
                    'mask_model_glo', 'comb_interped_zosto', ...
                    'lon_glo', 'lat_glo', 'lon', 'lat', '-v7.3');
                
                
                
            else
                load(matname);
% % % % % % % % % % % % % % % % % % %                     landname = [savedir,testname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     filename = strcat('D:\Data\Model\ROMS\nwp_1_20\test53\run\2005\ocean_rst2.nc');
% % % % % % % % % % % % % % % % % % %                     modelinfo=ncinfo(filename);
% % % % % % % % % % % % % % % % % % %                     lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
% % % % % % % % % % % % % % % % % % %                     lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     lon_west = abs(lon - (lonlat(1)-1));
% % % % % % % % % % % % % % % % % % %                     [min_lon_west, lon_min]=min(lon_west);
% % % % % % % % % % % % % % % % % % %                     lon_east = abs(lon - (lonlat(2)+1));
% % % % % % % % % % % % % % % % % % %                     [min_lon_east, lon_max]=min(lon_east);
% % % % % % % % % % % % % % % % % % %                     lat_south = abs(lat - (lonlat(3)-1));
% % % % % % % % % % % % % % % % % % %                     [min_lat_south, lat_min]=min(lat_south);
% % % % % % % % % % % % % % % % % % %                     lat_north = abs(lat - (lonlat(4)+1));
% % % % % % % % % % % % % % % % % % %                     [min_lat_north, lat_max]=min(lat_north);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
% % % % % % % % % % % % % % % % % % %                     lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
% % % % % % % % % % % % % % % % % % %                     mask_rho = ncread(filename,'mask_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
% % % % % % % % % % % % % % % % % % %                     interped_mask_rho = griddata(double(lon), double(lat), mask_rho,double(lon_glo),double(lat_glo));   
% % % % % % % % % % % % % % % % % % %                     land_model_interped=interped_mask_rho;
% % % % % % % % % % % % % % % % % % %                     land_model_interped(interped_mask_rho==0)=1;
% % % % % % % % % % % % % % % % % % %                     interped_mask_rho(interped_mask_rho==0)=NaN;
% % % % % % % % % % % % % % % % % % %                     land_model_interped(isfinite(interped_mask_rho))=NaN;
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     [scenname, error_status] = Func_0013_RCM_CMIP6_scenname(testname);
% % % % % % % % % % % % % % % % % % %                     gcmname = Func_0004_get_GCMname_from_RCM(testname);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     GLO_zostoname = ['D:\Data\Model\CMIP5\zosto\', scenname, '\', ...
% % % % % % % % % % % % % % % % % % %                         'zosto_interp_', gcmname, '_', scenname, '_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.nc'];
% % % % % % % % % % % % % % % % % % %                     lon_glo = ncread(GLO_zostoname,'lon',[1 1],[inf,1]);
% % % % % % % % % % % % % % % % % % %                     lat_glo = ncread(GLO_zostoname,'lat',[1 1],[1,inf]);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     lon_west_glo = abs(lon_glo - (lonlat(1)-1));
% % % % % % % % % % % % % % % % % % %                     [min_lon_west_glo, lon_min_glo]=min(lon_west_glo);
% % % % % % % % % % % % % % % % % % %                     lon_east_glo = abs(lon_glo - (lonlat(2)+1));
% % % % % % % % % % % % % % % % % % %                     [min_lon_east_glo, lon_max_glo]=min(lon_east_glo);
% % % % % % % % % % % % % % % % % % %                     lat_south_glo = abs(lat_glo - (lonlat(3)-1));
% % % % % % % % % % % % % % % % % % %                     [min_lat_south_glo, lat_min_glo]=min(lat_south_glo);
% % % % % % % % % % % % % % % % % % %                     lat_north_glo = abs(lat_glo - (lonlat(4)+1));
% % % % % % % % % % % % % % % % % % %                     [min_lat_north_glo, lat_max_glo]=min(lat_north_glo);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     lon_glo = ncread(GLO_zostoname,'lon', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
% % % % % % % % % % % % % % % % % % %                     lat_glo = ncread(GLO_zostoname,'lat', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
% % % % % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % %                     switch(regionname)
% % % % % % % % % % % % % % % % % % %                         otherwise
% % % % % % % % % % % % % % % % % % %                             mask_model_glo = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
% % % % % % % % % % % % % % % % % % %                             mask_model_glo(mask_model_glo==0)=NaN;
% % % % % % % % % % % % % % % % % % %                     end
% % % % % % % % % % % % % % % % % % %                     land_model_glo2=ncread(GLO_zostoname, 'zosto',[lon_min_glo(1) lat_min_glo(1), 1], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1, 1]);
% % % % % % % % % % % % % % % % % % %                     land_model_glo=land_model_glo2;
% % % % % % % % % % % % % % % % % % %                     land_model_glo(isnan(land_model_glo2))=1;
% % % % % % % % % % % % % % % % % % %                     land_model_glo(isfinite(land_model_glo2))=NaN;
% % % % % % % % % % % % % % % % % % %                     save(landname, 'land_model_interped', 'land_model_glo', '-v7.3');
%                 end
            end
            fig_flag=0;
        end

% % %         steric sea level analysis
        fig_flag=fig_flags{12,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, '_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
            
                comb_zosto_correct = comb_zosto;
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        raw_zosto=squeeze(comb_zosto(i,j,:));
                        std5_zosto = 5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                for i=1:size(comb_interped_zosto,1)
                    for j=1:size(comb_interped_zosto,2)
                        raw_zosto=squeeze(comb_interped_zosto(i,j,:));
                        std5_zosto =5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_interped_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                rawsshfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                mask_rho = ncread(rawsshfilename, 'trend');
                mask_rho(isfinite(mask_rho))=1;
%                 for tt= 1:size(ftime)
                    comb_zosto_correct = comb_zosto_correct .* mask_rho;
%                 end
%                 plot(squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan')))

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
                lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
%                 clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' steric SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                steric_sshvarid=netcdf.defVar(ncid, 'steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,steric_sshvarid,'long_name','steric_ssh');
                netcdf.putAtt(ncid,steric_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,steric_sshvarid,true,true,1);

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
%                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
                netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
                netcdf.putVar(ncid, steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_zosto_correct);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end
        
% % %         steric sea level analysis (interp)
        fig_flag=fig_flags{13,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, '_interped_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                comb_interped_zosto_correct = comb_interped_zosto;
                
                rawsshfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

                raw_ssh = ncread(rawsshfilename, 'raw_ssh');
                clear comb_interped_ssh_correct
%                 oldncoutfilename = strcat(savedir, 'old_', testname,'_',regionname, '_interped_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
%                 comb_interped_ssh_correct = ncread(oldncoutfilename, 'interped_ssh');
                for tt=1:size(ftime,2)
                    comb_interped_ssh_correct(:,:,tt) = griddata(double(lon), double(lat), double(squeeze(raw_ssh(:,:,tt))),double(lon_glo),double(lat_glo));   
                end

                mean_interped_zosto_correct = squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan'));
%                 diff_mean_interped_zosto_correct = diff( mean_interped_zosto_correct);
                slope_mzos = polyfit(ftime, mean_interped_zosto_correct,1);
                fit_mzos = polyval(slope_mzos, ftime);
                noise_ind = find(mean_interped_zosto_correct > fit_mzos'+std(mean_interped_zosto_correct));
                comb_interped_zosto_correct(:,:,noise_ind) = NaN;
                
                for i=1:size(comb_interped_zosto,1)
                    for j=1:size(comb_interped_zosto,2)
                        raw_zosto = squeeze(comb_interped_zosto_correct(i,j,:));
                        nanx = isnan(raw_zosto);
                        if sum(nanx) ~= length(ftime)
                            raw_zosto(nanx) = interp1(ftime(~nanx), raw_zosto(~nanx), ftime(nanx));
                            comb_interped_zosto_correct(i,j,:) = raw_zosto;
                            raw_ssh = squeeze(comb_interped_ssh_correct(i,j,:));
                            raw_ssh(nanx) = interp1(ftime(~nanx), raw_ssh(~nanx), ftime(nanx));
                            comb_interped_ssh_correct(i,j,:) = raw_ssh;
                        end
                    end
                end
                
                mask_rho = ncread(filename,'mask_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                interped_mask_rho = griddata(double(lon), double(lat), mask_rho,double(lon_glo),double(lat_glo));   
                land_model_interped=interped_mask_rho;
                land_model_interped(interped_mask_rho==0)=1;
                interped_mask_rho(interped_mask_rho==0)=NaN;
                land_model_interped(isfinite(interped_mask_rho))=NaN;
                
                ocean_mask_model = land_model_interped;
                ocean_mask_model(isnan(land_model_interped))=1;
                ocean_mask_model(isfinite(land_model_interped))=NaN;
                
                comb_interped_zosto_correct= comb_interped_zosto_correct .* ocean_mask_model .* mask_model_glo;
                
                
                
%                 mask_rho = squeeze(std(permute(comb_interped_zosto_correct, [3 1 2])));
%                 mask_rho (mask_rho < mask_rho(2,2).*5) = NaN;
%                 pcolor(mask_rho'); shading flat; colorbar

                plot(squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan')))
                       
                switch testname
                    case {'test53', 'test54', 'test55', 'test56'}
                        scenname='historical';
                    case {'test61', 'test62', 'test63', 'test64'}
                        scenname='rcp26';
                    case {'test57', 'test58', 'test59', 'test60'}
                        scenname='rcp45';
                    case {'test65', 'test66', 'test67', 'test68'}
                        scenname='rcp85';
                end
                switch testname
                    case {'test53', 'test61', 'test57', 'test65'}
                        gcmname='IPSL-CM5A-LR';
                    case {'test54', 'test62', 'test58', 'test66'}
                        gcmname='IPSL-CM5A-MR';
                    case {'test55', 'test63', 'test59', 'test67'}
                        gcmname='NorESM1-M';
                    case {'test56', 'test64', 'test60', 'test68'}
                        gcmname='MPI-ESM-LR';
                end
                
                GLO_zostoname = ['D:\Data\Model\CMIP5\zosto\', scenname, '\', ...
                    'zosto_interp_', gcmname, '_', scenname, '_', num2str(min(inputyear),'%04i'), '-', num2str(max(inputyear),'%04i'), '.nc'];

                if (exist('lon_min_glo')==0)
                    GLO_modelinfo=ncinfo(GLO_zostoname);
                    lon_glo = ncread(GLO_zostoname,'lon',[1 1],[inf,1]);
                    lat_glo = ncread(GLO_zostoname,'lat',[1 1],[1,inf]);

                    lon_west_glo = abs(lon_glo - (lonlat(1)-1));
                    [min_lon_west_glo, lon_min_glo]=min(lon_west_glo);
                    lon_east_glo = abs(lon_glo - (lonlat(2)+1));
                    [min_lon_east_glo, lon_max_glo]=min(lon_east_glo);
                    lat_south_glo = abs(lat_glo - (lonlat(3)-1));
                    [min_lat_south_glo, lat_min_glo]=min(lat_south_glo);
                    lat_north_glo = abs(lat_glo - (lonlat(4)+1));
                    [min_lat_north_glo, lat_max_glo]=min(lat_north_glo);

                    lon_glo = ncread(GLO_zostoname,'lon', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
                    lat_glo = ncread(GLO_zostoname,'lat', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);

                    switch(regionname)
                        otherwise
                            mask_model_glo = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
                            mask_model_glo(mask_model_glo==0)=NaN;
                    end
                end
                comb_gcm_zosto = ncread(GLO_zostoname,'zosto', [lon_min_glo(1) lat_min_glo(1) 1 ], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1 inf]);
                
                comb_gcm_zosto = comb_gcm_zosto .* mask_model_glo;
                
                len_lon_model = size( lon_glo,1);
                len_lat_model = size( lat_glo,2);
                % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
                lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
%                 clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'GCM and RCM monthly interped steric SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' steric SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'source', [' ROMS NWP 1/20 data and from _ ',testname, ' and ', gcmname]);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'author', 'Created by Y.Y.Kim');
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'date', date);

                timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
                netcdf.putAtt(ncid,timevarid,'long_name','time');
                netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
                netcdf.putAtt(ncid,timevarid,'calendar','gregorian');

                lon_glovarid=netcdf.defVar(ncid, 'lon_glo', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,lon_glovarid,'long_name','lon_model');
                netcdf.putAtt(ncid,lon_glovarid,'units','degree_east');

                lat_glovarid=netcdf.defVar(ncid, 'lat_glo', 'NC_DOUBLE', [lon_dimid lat_dimid]);
                netcdf.putAtt(ncid,lat_glovarid,'long_name','lat_model');
                netcdf.putAtt(ncid,lat_glovarid,'units','degree_north');

                interped_steric_sshvarid=netcdf.defVar(ncid, 'interped_steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_steric_sshvarid,'long_name','interped_steric_ssh');
                netcdf.putAtt(ncid,interped_steric_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,interped_steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,interped_steric_sshvarid,true,true,1);
                
                interped_sshvarid=netcdf.defVar(ncid, 'interped_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,interped_sshvarid,'long_name','interped_ssh');
                netcdf.putAtt(ncid,interped_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,interped_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,interped_sshvarid,true,true,1);

                gcm_steric_sshvarid=netcdf.defVar(ncid, 'gcm_steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,gcm_steric_sshvarid,'long_name',['gcm_steric_ssh from ', gcmname]);
                netcdf.putAtt(ncid,gcm_steric_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,gcm_steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,gcm_steric_sshvarid,true,true,1);
                
                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
%                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_glovarid, [0 0], [len_lon_model len_lat_model], lon_glo);
                netcdf.putVar(ncid, lat_glovarid, [0 0], [len_lon_model len_lat_model], lat_glo);
                netcdf.putVar(ncid, gcm_steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_gcm_zosto);
                netcdf.putVar(ncid, interped_steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_interped_zosto_correct);
                netcdf.putVar(ncid, interped_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_interped_ssh_correct);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end

% % %         get ts & calculate zosto_thermo
        fig_flag=fig_flags{14,2};
        while (fig_flag)
            run(param_script);
            ind=1;
            matname = [savedir,testname,'_',regionname,'model_thermosteric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
            clear rho_0_va
            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    if (exist('rho_0_va')==0)
                        for monthij = 1:length(inputmonth)
                            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                            tic;
                            tempyear = inputyear(yearij);
                            tempmonth = inputmonth(monthij);
                            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                            filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                                
                            % read model data
                            if (exist('lon_min')==0)
                                modelinfo=ncinfo(filename);
                                lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                                lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
                                
                                [lon_min, lon_max, lat_min, lat_max] = ...
                                    Func_0012_findind_Y(1, lonlat, lon, lat);

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

                            data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                            if (exist('h')==0)
                                h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                dA=1./pm .* 1./pn;
                                Vtransform = ncread(filename, 'Vtransform');
                                Vstretching = ncread(filename, 'Vstretching');
                                theta_s = ncread(filename, 'theta_s');
                                theta_b = ncread(filename, 'theta_b');
                                hc = ncread(filename, 'hc');
                                N = length(ncread(filename, 's_rho'));
                            end
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                            PT_src(:,:,:,monthij) = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(:,:,:,monthij) = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(S_src<=0)=0;

                            vtype ='w';
                            zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                            dz = diff(zw);
                            dz = permute(dz, [2,3,1]);
    %                         vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
    %                         area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
                            vtype ='r';
                            P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                            P = -permute(P, [2,3,1]);
                            comb_P(:,:,:,monthij)=P;
                        end
                            PT_src = mean(PT_src, 4 );
                            S_src = mean(S_src, 4 );
                            P = mean(comb_P, 4);
                            T_src = sw_temp(S_src, PT_src, P, zeros(size(P)) );
                            rho_src = sw_dens(S_src, T_src, P );
                            dep_0 = h + data;
                            rho_0_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_0;       % vertical(local) average, 2-D
                    end
                    
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                        if (exist('h')==0)
                            h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            dA=1./pm .* 1./pn;
                            Vtransform = ncread(filename, 'Vtransform');
                            Vstretching = ncread(filename, 'Vstretching');
                            theta_s = ncread(filename, 'theta_s');
                            theta_b = ncread(filename, 'theta_b');
                            hc = ncread(filename, 'hc');
                            N = length(ncread(filename, 's_rho'));
                        end
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % % %                        S should be used as climate value
                        PT_src = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
%                         S_src = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                        S_src(S_src<=0)=0;

                        vtype ='w';
                        zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                        dz = diff(zw);
                        dz = permute(dz, [2,3,1]);

                        vtype ='r';
                        P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                        P = -permute(P, [2,3,1]);
                        T_src = sw_temp (S_src, PT_src, P, zeros(size(P)) );
                        rho_src = sw_dens (S_src, T_src, P );

                        dep_n = h + data;
                        rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
%                         rho_n_va2 = sum( rho_src.*dz, 3 , 'omitnan') ./ (h+data);

%                         rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
                        zosto = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
%                         zosto_2 = dep_n       .* ( 1 - rho_n_va ./ rho_0_va );
                        comb_zosto_thermo(:,:,ind) = zosto;

                        len_lon_model = size(data,1);
                        len_lat_model = size(data,2);

                        GLO_zostoname='D:\Data\Model\CMIP5\zosto\historical\zosto_interp_IPSL-CM5A-LR_historical_1976-2005.nc';

                        if (exist('lon_glo')==0)
                            GLO_modelinfo=ncinfo(GLO_zostoname);
                            lon_glo = ncread(GLO_zostoname,'lon',[1 1],[inf,1]);
                            lat_glo = ncread(GLO_zostoname,'lat',[1 1],[1,inf]);

                            lon_west_glo = abs(lon_glo - (lonlat(1)-1));
                            [min_lon_west_glo, lon_min_glo]=min(lon_west_glo);
                            lon_east_glo = abs(lon_glo - (lonlat(2)+1));
                            [min_lon_east_glo, lon_max_glo]=min(lon_east_glo);
                            lat_south_glo = abs(lat_glo - (lonlat(3)-1));
                            [min_lat_south_glo, lat_min_glo]=min(lat_south_glo);
                            lat_north_glo = abs(lat_glo - (lonlat(4)+1));
                            [min_lat_north_glo, lat_max_glo]=min(lat_north_glo);

                            lon_glo = ncread(GLO_zostoname,'lon', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
                            lat_glo = ncread(GLO_zostoname,'lat', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);

                            switch(regionname)
                                otherwise
                                    mask_model_glo = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model_glo(mask_model_glo==0)=NaN;
                            end
                        end

                        interped_zosto = griddata(double(lon), double(lat), zosto,double(lon_glo),double(lat_glo));   

                        comb_interped_zosto_thermo(:,:,ind) = interped_zosto;

                        ind = ind + 1;
                        toc;
                    end
                end
                
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_zosto_thermo',  ...
                    'mask_model_glo', 'comb_interped_zosto_thermo', ...
                    'lon_glo', 'lat_glo', 'lon', 'lat', '-v7.3');
            else
                load(matname);
            end
            fig_flag=0;
        end
        
% % %         thermosteric sea level analysis (ncsave)
        fig_flag=fig_flags{15,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, '_thermosteric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                comb_zosto_correct = comb_zosto_thermo;
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        raw_zosto=squeeze(comb_zosto_thermo(i,j,:));
                        std5_zosto = 5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                for i=1:size(comb_interped_zosto,1)
                    for j=1:size(comb_interped_zosto,2)
                        raw_zosto=squeeze(comb_interped_zosto(i,j,:));
                        std5_zosto =5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_interped_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                rawsshfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                mask_rho = ncread(rawsshfilename, 'trend');
                mask_rho(isfinite(mask_rho))=1;
%                 for tt= 1:size(ftime)
                    comb_zosto_correct = comb_zosto_correct .* mask_rho;
%                 end
%                 plot(squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan')))

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
                lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
%                 clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' steric SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                steric_sshvarid=netcdf.defVar(ncid, 'steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,steric_sshvarid,'long_name','steric_ssh');
                netcdf.putAtt(ncid,steric_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,steric_sshvarid,true,true,1);

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
%                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
                netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
                netcdf.putVar(ncid, steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_zosto_correct);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end
        
% % %         get ts & calculate zosto_halo
        fig_flag=fig_flags{16,2};
        while (fig_flag)
            run(param_script);
            ind=1;
            matname = [savedir,testname,'_',regionname,'model_halosteric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];
            clear rho_0_va
            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    if (exist('rho_0_va')==0)
                        for monthij = 1:length(inputmonth)
                            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                            tic;
                            tempyear = inputyear(yearij);
                            tempmonth = inputmonth(monthij);
                            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                            filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                                
                            % read model data
                            if (exist('lon_min')==0)
                                modelinfo=ncinfo(filename);
                                lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                                lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);
                                
                                [lon_min, lon_max, lat_min, lat_max] = ...
                                    Func_0012_findind_Y(1, lonlat, lon, lat);

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

                            data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                            if (exist('h')==0)
                                h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                                dA=1./pm .* 1./pn;
                                Vtransform = ncread(filename, 'Vtransform');
                                Vstretching = ncread(filename, 'Vstretching');
                                theta_s = ncread(filename, 'theta_s');
                                theta_b = ncread(filename, 'theta_b');
                                hc = ncread(filename, 'hc');
                                N = length(ncread(filename, 's_rho'));
                            end
                            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                            PT_src(:,:,:,monthij) = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(:,:,:,monthij) = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                            S_src(S_src<=0)=0;

                            vtype ='w';
                            zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                            dz = diff(zw);
                            dz = permute(dz, [2,3,1]);
    %                         vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
    %                         area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
                            vtype ='r';
                            P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                            P = -permute(P, [2,3,1]);
                            comb_P(:,:,:,monthij)=P;
                        end
                            PT_src = mean(PT_src, 4 );
                            S_src = mean(S_src, 4 );
                            P = mean(comb_P, 4);
                            T_src = sw_temp(S_src, PT_src, P, zeros(size(P)) );
                            rho_src = sw_dens(S_src, T_src, P );
                            dep_0 = h + data;
                            rho_0_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_0;       % vertical(local) average, 2-D
                    end
                    
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                                'pck_', testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

                        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

                        if (exist('h')==0)
                            h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            dA=1./pm .* 1./pn;
                            Vtransform = ncread(filename, 'Vtransform');
                            Vstretching = ncread(filename, 'Vstretching');
                            theta_s = ncread(filename, 'theta_s');
                            theta_b = ncread(filename, 'theta_b');
                            hc = ncread(filename, 'hc');
                            N = length(ncread(filename, 's_rho'));
                        end
                        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
% % % % %                        PT should be used as climate value
%                         PT_src = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                        S_src = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                        S_src(S_src<=0)=0;

                        vtype ='w';
                        zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

                        dz = diff(zw);
                        dz = permute(dz, [2,3,1]);

                        vtype ='r';
                        P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
                        P = -permute(P, [2,3,1]);
                        T_src = sw_temp (S_src, PT_src, P, zeros(size(P)) );
                        rho_src = sw_dens (S_src, T_src, P );

                        dep_n = h + data;
                        rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
%                         rho_n_va2 = sum( rho_src.*dz, 3 , 'omitnan') ./ (h+data);

%                         rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
                        zosto = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
%                         zosto_2 = dep_n       .* ( 1 - rho_n_va ./ rho_0_va );
                        comb_zosto_halo(:,:,ind) = zosto;

                        len_lon_model = size(data,1);
                        len_lat_model = size(data,2);

                        GLO_zostoname='D:\Data\Model\CMIP5\zosto\historical\zosto_interp_IPSL-CM5A-LR_historical_1976-2005.nc';

                        if (exist('lon_glo')==0)
                            GLO_modelinfo=ncinfo(GLO_zostoname);
                            lon_glo = ncread(GLO_zostoname,'lon',[1 1],[inf,1]);
                            lat_glo = ncread(GLO_zostoname,'lat',[1 1],[1,inf]);

                            lon_west_glo = abs(lon_glo - (lonlat(1)-1));
                            [min_lon_west_glo, lon_min_glo]=min(lon_west_glo);
                            lon_east_glo = abs(lon_glo - (lonlat(2)+1));
                            [min_lon_east_glo, lon_max_glo]=min(lon_east_glo);
                            lat_south_glo = abs(lat_glo - (lonlat(3)-1));
                            [min_lat_south_glo, lat_min_glo]=min(lat_south_glo);
                            lat_north_glo = abs(lat_glo - (lonlat(4)+1));
                            [min_lat_north_glo, lat_max_glo]=min(lat_north_glo);

                            lon_glo = ncread(GLO_zostoname,'lon', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);
                            lat_glo = ncread(GLO_zostoname,'lat', [lon_min_glo(1) lat_min_glo(1)], [lon_max_glo(1)-lon_min_glo(1)+1 lat_max_glo(1)-lat_min_glo(1)+1]);

                            switch(regionname)
                                otherwise
                                    mask_model_glo = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model_glo(mask_model_glo==0)=NaN;
                            end
                        end

                        interped_zosto = griddata(double(lon), double(lat), zosto,double(lon_glo),double(lat_glo));   

                        comb_interped_zosto_halo(:,:,ind) = interped_zosto;

                        ind = ind + 1;
                        toc;
                    end
                end
                
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_zosto_halo',  ...
                    'mask_model_glo', 'comb_interped_zosto_halo', ...
                    'lon_glo', 'lat_glo', 'lon', 'lat', '-v7.3');
            else
                load(matname);
            end
            fig_flag=0;
        end

% % %         halosteric sea level analysis (ncsave)
        fig_flag=fig_flags{17,2};
        while (fig_flag)
            ncoutfilename = strcat(savedir, testname,'_',regionname, '_halosteric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
                comb_zosto_correct = comb_zosto_halo;
                for i=1:len_lon_model
                    for j=1:len_lat_model
                        raw_zosto=squeeze(comb_zosto_halo(i,j,:));
                        std5_zosto = 5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                for i=1:size(comb_interped_zosto,1)
                    for j=1:size(comb_interped_zosto,2)
                        raw_zosto=squeeze(comb_interped_zosto(i,j,:));
                        std5_zosto =5.* std(raw_zosto);
                        diff_zosto = diff(raw_zosto);
                        if diff_zosto(end)> std5_zosto;
                            raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
                            diff_zosto(end)=diff_zosto(end-1);
                        end
                        raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
                        continuity_error=diff(find(diff_zosto>std5_zosto));
                        continuity_error(continuity_error==1)=55555;
                        if continuity_error==55555
                            'error'
                            break
                        end
                        comb_interped_zosto_correct(i,j,:)=raw_zosto;
                    end
                end
                
                rawsshfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
                mask_rho = ncread(rawsshfilename, 'trend');
                mask_rho(isfinite(mask_rho))=1;
%                 for tt= 1:size(ftime)
                    comb_zosto_correct = comb_zosto_correct .* mask_rho;
%                 end
%                 plot(squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan')))

            % % %         make ncfile
                ncid = netcdf.create(ncoutfilename,'NETCDF4');

                lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
                lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
                time_dimid = netcdf.defDim(ncid, 'time', 0);
%                 clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'type', ['NWP 1/20 _ ', testname, 'model monthly SSH analysis file']);
                netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                    'title', [' steric SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
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

                steric_sshvarid=netcdf.defVar(ncid, 'steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
                netcdf.putAtt(ncid,steric_sshvarid,'long_name','steric_ssh');
                netcdf.putAtt(ncid,steric_sshvarid,'units','m');
                netcdf.defVarChunking(ncid,steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
                netcdf.defVarDeflate(ncid,steric_sshvarid,true,true,1);

                netcdf.endDef(ncid);

                netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
%                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
                netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
                netcdf.putVar(ncid, steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_zosto_correct);
                netcdf.close(ncid);
            end
            fig_flag=0;
        end
        
    end
end
