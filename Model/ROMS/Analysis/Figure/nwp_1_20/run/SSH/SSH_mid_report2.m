close all; clear all;  clc;
% %  get satellite SSH data and compared data
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'AKP'}
all_testname = {'test55', 'test56'};


for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname

        % % % 
        % % % Read Model SST
        % % % interp
        % % % get RMS
        % % % get BIAS
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\KYY\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
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
        inputyear = [1993:2008]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
        regionname=all_region{regionind}
        run('nwp_polygon_point.m');
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
            case('SS') %% South Sea
                refpolygon=sspolygon;
            case('YS') %% Yellow Sea
                refpolygon=yspolygon;
            case('ECS') %% East China Sea
                refpolygon=ecspolygon;
            case('AKP') %% Around Korea Peninsula
                refpolygon=akppolygon;
            otherwise
                ('?')
        end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));

        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            avisodir='E:\Data\Observation\AVISO\monthly\';
        elseif (strcmp(system_name,'GLNXA64'))
        end

        run(param_script);
        ind=1;
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


                data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

        %         data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
                data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                data=data.*mask_model;

                % read OISST DATA
                avisofilename = strcat(avisodir,'dt_global_allsat_msla_h_y', num2str(tempyear,'%04i'), ...
                        '_m',num2str(tempmonth,'%02i'), '.nc');
                if (exist('aviso_lon')==0)
                    avisoinfo=ncinfo(avisofilename);
                    aviso_lon = ncread(avisofilename,'lon');
                    aviso_lat = ncread(avisofilename,'lat');

                    aviso_lon_west = abs(aviso_lon - (lonlat(1)));
                    min_aviso_lon_west=min(aviso_lon_west);
                    aviso_lon_east = abs(aviso_lon - (lonlat(2)));
                    min_aviso_lon_east=min(aviso_lon_east);
                    aviso_lat_south = abs(aviso_lat - (lonlat(3)));
                    min_aviso_lat_south=min(aviso_lat_south);
                    aviso_lat_north = abs(aviso_lat - (lonlat(4)));
                    min_aviso_lat_north=min(aviso_lat_north);

                    aviso_lon_min = find(aviso_lon_west == min_aviso_lon_west);
                    aviso_lon_max = find(aviso_lon_east == min_aviso_lon_east);
                    aviso_lat_min = find(aviso_lat_south == min_aviso_lat_south);
                    aviso_lat_max = find(aviso_lat_north == min_aviso_lat_north);

            %         ncinfo('E:\Data\Observation\OISST\monthly\aviso_monthly1983_11.nc');

                    aviso_lon = ncread(avisofilename,'lon', [aviso_lon_min(1)], [aviso_lon_max(1)-aviso_lon_min(1)]);
                    aviso_lat = ncread(avisofilename,'lat', [aviso_lat_min(1)], [aviso_lat_max(1)-aviso_lat_min(1)]);

        %             comb_spatial_meanrms=(zeros([length(aviso_lon),length(aviso_lat),12]));
        %             comb_spatial_meanbias=(zeros([length(aviso_lon),length(aviso_lat),12]));
                    comb_spatial_meanaviso=(zeros([length(aviso_lon),length(aviso_lat),12]));
                    comb_spatial_meanmodel=(zeros([length(aviso_lon),length(aviso_lat),12]));

                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_aviso(1:size(aviso_lon,1),1:size(aviso_lon,2))=1;
                        otherwise
                            [aviso_lat2 aviso_lon2]=meshgrid(aviso_lat, aviso_lon);
                            mask_aviso = double(inpolygon(aviso_lon2,aviso_lat2,refpolygon(:,1),refpolygon(:,2)));
                            mask_aviso(mask_aviso==0)=NaN;
                    end
                end
                len_lon = length(aviso_lon);
                len_lat = length(aviso_lat);
                len_lon_model = size(data,1);
                len_lat_model = size(data,2);

                aviso_data = ncread(avisofilename,'sla',[aviso_lon_min(1) aviso_lat_min(1) 1], [aviso_lon_max(1)-aviso_lon_min(1) aviso_lat_max(1)-aviso_lat_min(1) 1]);
                aviso_data=aviso_data.*mask_aviso;

                interped_data = griddata(double(lon), double(lat), data,double(aviso_lon),double(aviso_lat)')';   
        %         bias = aviso_data-interped_data;  
        %         rms = sqrt((aviso_data-interped_data).^2);  
        %         meanbias = mean(mean(bias,'omitnan'),'omitnan');
        %         meanrms = mean(mean(rms,'omitnan'),'omitnan');
                comb_data(:,:,ind) = data;
                comb_interped_data(:,:,ind) = interped_data;
                comb_aviso_data(:,:,ind) = aviso_data;
        %         comb_rms_data(:,:,ind) = rms;
        %         comb_bias_data(:,:,ind) = bias;

        %         comb_meanrms(yearij,monthij)=meanrms;
        %         comb_meanbias(yearij,monthij)=meanbias;
        %         comb_spatial_meanrms(:,:,monthij)=comb_spatial_meanrms(:,:,monthij)+rms/double(length(inputyear));
        %         comb_spatial_meanbias(:,:,monthij)=comb_spatial_meanbias(:,:,monthij)+bias/double(length(inputyear));
                comb_spatial_meanaviso(:,:,monthij)=comb_spatial_meanaviso(:,:,monthij)+aviso_data/double(length(inputyear));
                comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+interped_data/double(length(inputyear));
                ind = ind + 1;
                toc;
            end
        end

        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data(i,j,:))',1);
                trend(i,j)=p(1);
            end
        end
        trend = trend * 1000.0; %% m/y -> mm/y

        aviso_trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_aviso_data(i,j,:))',1);
                aviso_trend(i,j)=p(1);
            end
        end
        aviso_trend = aviso_trend * 1000.0; %% m/y -> mm/y


        for t=1:length(inputyear)
            comb_interped_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
            comb_aviso_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_aviso_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanaviso;
        end


        trend_filtered(1:len_lon,1:len_lat)=NaN;
        aviso_trend_filtered(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered(i,j,:))',1);
                trend_filtered(i,j)=p(1);
            end
        end
        trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y

        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_aviso_data_filtered(i,j,:))',1);
                aviso_trend_filtered(i,j)=p(1);
            end
        end
        aviso_trend_filtered = aviso_trend_filtered * 1000.0; %% m/y -> mm/y

        mean_trend=mean(mean(trend,'omitnan'),'omitnan');
        mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
        mean_aviso_trend=mean(mean(aviso_trend,'omitnan'),'omitnan');
        mean_aviso_trend_filtered=mean(mean(aviso_trend_filtered,'omitnan'),'omitnan');

        % % pcolor(trend_filtered')
        % % shading interp
        % % mean(mean(trend,'omitnan'),'omitnan')
        % % colorbar;



        figdir=[figrawdir,'CLIM\'];
        outfile = strcat(figdir,regionname);
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        % 
        % rmsplot=plot(mean(comb_meanrms,1),'k')
        % jpgname=strcat(outfile, '_', testname, '_climrms', '.jpg'); %% ~_year_month.jpg
        % xlabel('month')
        % ylabel('rms(^o)')
        % title(['rms, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
        % ylim([0 4])
        % set(rmsplot,'LineWidth',2);
        % grid on
        % saveas(gcf,jpgname,'jpg');
        % grid off
        % 
        % biasplot=plot(mean(comb_meanbias,1) ,'k')
        % jpgname=strcat(outfile, '_', testname, '_climbias', '.jpg'); %% ~_year_month.jpg
        % xlabel('month')
        % ylabel('bias(^o)')
        % title(['bias, climatology(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),')'])
        % ylim([-4 4])
        % set(biasplot,'LineWidth',2);
        % grid on
        % saveas(gcf,jpgname,'jpg');
        % grid off

        save([filedir,testname,'_',regionname,'ssh_aviso_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        ncid = netcdf.create(strcat(filedir, testname,'_',regionname, '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
        xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
        etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['NWP 1/20 _ ', testname, 'model, aviso monthly SSH analysis file']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'title', ' monthly SSH analysis (1993-2009) ');
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

        lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', lon_dimid);
        netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
        netcdf.putAtt(ncid,lonvarid,'units','degree_east');

        latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', lat_dimid);
        netcdf.putAtt(ncid,latvarid,'long_name','latitude');
        netcdf.putAtt(ncid,latvarid,'units','degree_north');

        lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
        netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

        lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xidimid etadimid]);
        netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
        netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

        raw_sshvarid=netcdf.defVar(ncid, 'raw_ssh', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_sshvarid,'long_name','raw_ssh');
        netcdf.putAtt(ncid,raw_sshvarid,'units','m');

        sshvarid=netcdf.defVar(ncid, 'ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sshvarid,'long_name','ssh');
        netcdf.putAtt(ncid,sshvarid,'units','m');

        ssh_filteredvarid=netcdf.defVar(ncid, 'ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,ssh_filteredvarid,'long_name','ssh_filtered');
        netcdf.putAtt(ncid,ssh_filteredvarid,'units','m');

        aviso_sshvarid=netcdf.defVar(ncid, 'aviso_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,aviso_sshvarid,'long_name','aviso_ssh');
        netcdf.putAtt(ncid,aviso_sshvarid,'units','m');

        aviso_ssh_filteredvarid=netcdf.defVar(ncid, 'aviso_ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,aviso_ssh_filteredvarid,'long_name','aviso_ssh_filtered');
        netcdf.putAtt(ncid,aviso_ssh_filteredvarid,'units','m');

        trendvarid=netcdf.defVar(ncid, 'trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trendvarid,'long_name','trend');
        netcdf.putAtt(ncid,trendvarid,'units','mm/year');

        trend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,trend_filteredvarid,'long_name','trend_filtered');
        netcdf.putAtt(ncid,trend_filteredvarid,'units','mm/year');

        aviso_trendvarid=netcdf.defVar(ncid, 'aviso_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,aviso_trendvarid,'long_name','aviso_trend');
        netcdf.putAtt(ncid,aviso_trendvarid,'units','mm/year');

        aviso_trend_filteredvarid=netcdf.defVar(ncid, 'aviso_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,aviso_trend_filteredvarid,'long_name','aviso_trend_filtered');
        netcdf.putAtt(ncid,aviso_trend_filteredvarid,'units','mm/year');

        mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
        netcdf.putAtt(ncid,mean_trendvarid,'units','mm/year');

        mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','mm/year');

        mean_aviso_trendvarid=netcdf.defVar(ncid, 'mean_aviso_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_aviso_trendvarid,'long_name','mean_aviso_trend');
        netcdf.putAtt(ncid,mean_aviso_trendvarid,'units','mm/year');

        mean_aviso_trend_filteredvarid=netcdf.defVar(ncid, 'mean_aviso_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_aviso_trend_filteredvarid,'long_name','mean_aviso_trend_filtered');
        netcdf.putAtt(ncid,mean_aviso_trend_filteredvarid,'units','mm/year');

        clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_sshvarid,'long_name','clim_ssh');
        netcdf.putAtt(ncid,clim_sshvarid,'units','m');

        clim_avisovarid=netcdf.defVar(ncid, 'clim_aviso_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_avisovarid,'long_name','clim_aviso_ssh');
        netcdf.putAtt(ncid,clim_avisovarid,'units','m');


        netcdf.endDef(ncid);

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

        netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
        netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
        netcdf.putVar(ncid, lonvarid, 0, len_lon, aviso_lon);
        netcdf.putVar(ncid, latvarid, 0, len_lat, aviso_lat);
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
        netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
        netcdf.putVar(ncid, sshvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data);
        netcdf.putVar(ncid, ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered);
        netcdf.putVar(ncid, aviso_ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_aviso_data_filtered);
        netcdf.putVar(ncid, aviso_sshvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_aviso_data);
        netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
        netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
        netcdf.putVar(ncid, aviso_trendvarid, [0 0], [len_lon len_lat], aviso_trend);
        netcdf.putVar(ncid, aviso_trend_filteredvarid, [0 0], [len_lon len_lat], aviso_trend_filtered);
        netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
        netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
        netcdf.putVar(ncid, mean_aviso_trendvarid, [0], [1], mean_aviso_trend);
        netcdf.putVar(ncid, mean_aviso_trend_filteredvarid, [0], [1], mean_aviso_trend_filtered);
        netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
        netcdf.putVar(ncid, clim_avisovarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanaviso);

        netcdf.close(ncid);
    end
end