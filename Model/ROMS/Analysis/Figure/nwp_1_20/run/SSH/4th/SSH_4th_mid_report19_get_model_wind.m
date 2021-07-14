close all; clear all;  clc;
% %  get wind, interpolate it to reSSH grid. save. 

% all_region ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region ={'AKP2', 'NWP', 'ES', 'YS'}
% all_testname = {'ens03','test53','test54','test55','test56'};
all_testname = {'test52'};

all_region ={'AKP2'};
for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind testnameind all_region all_testname

        % % % 
        % % % Read Model SSH
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
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROM` S/Grid_kyy']));
        end

        shadlev = [-2 2];
        % rms_shadlev = [0 4];
        trend_shadlev = [0 4];
        % bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/20;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1980:2008]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

%         varname ='zeta'
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
            case('AKP2') %% Around Korea Peninsula
                refpolygon=akp2polygon;
            case('CA') %% Coastal Area around korea peninsula
                refpolygon=capolygon;
            case('EKB') %% Coastal Area around korea peninsula
                refpolygon=ekbpolygon;
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
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('E:\Data\Model\ROMS\nwp_1_20\input\', testname, '\'); % % where data files are
            recondir='E:\Data\Reanalysis\CSEOF_reconstructed_SSHA\';
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
                filename = strcat(filedir,  ...
                        'nwp_1_20_', num2str(tempyear,'%04i'), '_Uwind.nc');
                filename2 = strcat(filedir,  ...
                        'nwp_1_20_', num2str(tempyear,'%04i'), '_Vwind.nc');
                gridname = strcat(filedir,  ...
                        'roms_grid_nwp_1_20_',testname,'.nc');
                % read model data
                if (exist('lon')==0)
                    modelinfo=ncinfo(gridname);
                    lon = ncread(gridname,'lon_rho',[1 1],[modelinfo.Variables(23).Dimensions(1).Length,1]);
                    lat = ncread(gridname,'lat_rho',[1 1],[1,modelinfo.Variables(27).Dimensions(2).Length]);

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


                data_info_u = ncinfo(filename, 'Uwind');  %% [lon lat depth time] -> [1601 1201 33 1]
                data_info_v = ncinfo(filename2, 'Vwind');  %% [lon lat depth time] -> [1601 1201 33 1]

        %         data = ncread(filename,varname,[lon_min(1) lat_min(1) 40 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1 1]);
                if (tempmonth==1)
                    sind=1;
                    eind=31;
                    cind=31;
                else
                    sind=sum(eomday(tempyear,1:tempmonth-1)) + 1;
                    eind=sind + eomday(tempyear,tempmonth) - 2;
                    cind=eomday(tempyear,tempmonth);
                end
                clear Udata_daily Vdata_daily
                Udata_daily = ncread(filename,'Uwind',[lon_min(1) lat_min(1) sind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 cind]);
                Vdata_daily = ncread(filename2,'Vwind',[lon_min(1) lat_min(1) sind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 cind]);
                Udata=squeeze(mean(Udata_daily,3));
                Vdata=squeeze(mean(Vdata_daily,3));
                
                Udata=Udata.*mask_model;
                Vdata=Vdata.*mask_model;
                SWdata=Udata*cos(pi/4) + Vdata*cos(pi/4);
                
                % read reconstructed SSH DATA
                reconfilename = strcat(recondir,'CCAR_recon_sea_level_monthly_merged', '.nc');
                if (exist('recon_lon')==0)
                    reconinfo=ncinfo(reconfilename);
                    recon_lon = ncread(reconfilename,'lon');
                    recon_lat = ncread(reconfilename,'lat');

                    recon_lon_west = abs(recon_lon - (lonlat(1)));
                    min_recon_lon_west=min(recon_lon_west);
                    recon_lon_east = abs(recon_lon - (lonlat(2)));
                    min_recon_lon_east=min(recon_lon_east);
                    recon_lat_south = abs(recon_lat - (lonlat(3)));
                    min_recon_lat_south=min(recon_lat_south);
                    recon_lat_north = abs(recon_lat - (lonlat(4)));
                    min_recon_lat_north=min(recon_lat_north);

                    recon_lon_min = find(recon_lon_west == min_recon_lon_west);
                    recon_lon_max = find(recon_lon_east == min_recon_lon_east);
                    recon_lat_min = find(recon_lat_south == min_recon_lat_south);
                    recon_lat_max = find(recon_lat_north == min_recon_lat_north);

                    recon_lon = ncread(reconfilename,'lon', [recon_lon_min(1)], [recon_lon_max(1)-recon_lon_min(1)]);
                    recon_lat = ncread(reconfilename,'lat', [recon_lat_min(1)], [recon_lat_max(1)-recon_lat_min(1)]);

                    comb_spatial_mean_U=(zeros([length(recon_lon),length(recon_lat),12]));
                    comb_spatial_mean_V=(zeros([length(recon_lon),length(recon_lat),12]));
                    comb_spatial_mean_SW=(zeros([length(recon_lon),length(recon_lat),12]));

                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_recon(1:size(recon_lon,1),1:size(recon_lon,2))=1;
                        otherwise
                            [recon_lat2 recon_lon2]=meshgrid(recon_lat, recon_lon);
                            mask_recon = double(inpolygon(recon_lon2,recon_lat2,refpolygon(:,1),refpolygon(:,2)));
                            mask_recon(mask_recon==0)=NaN;
                    end
                end
                len_lon = length(recon_lon);
                len_lat = length(recon_lat);
                len_lon_model = size(Udata,1);
                len_lat_model = size(Udata,2);

                interped_Udata = griddata(double(lon), double(lat), double(Udata),double(recon_lon),double(recon_lat)')';   
                interped_Vdata = griddata(double(lon), double(lat), double(Vdata),double(recon_lon),double(recon_lat)')';   
                interped_SWdata = griddata(double(lon), double(lat), double(SWdata),double(recon_lon),double(recon_lat)')';  
                
                comb_Udata(:,:,ind) = Udata;
                comb_Vdata(:,:,ind) = Vdata;
                comb_SWdata(:,:,ind) = SWdata;
                
                comb_interped_Udata(:,:,ind) = interped_Udata;
                comb_interped_Vdata(:,:,ind) = interped_Vdata;
                comb_interped_SWdata(:,:,ind) = interped_SWdata;
                
                comb_spatial_mean_U(:,:,monthij)=comb_spatial_mean_U(:,:,monthij)+interped_Udata/double(length(inputyear));
                comb_spatial_mean_V(:,:,monthij)=comb_spatial_mean_V(:,:,monthij)+interped_Vdata/double(length(inputyear));
                comb_spatial_mean_SW(:,:,monthij)=comb_spatial_mean_SW(:,:,monthij)+interped_SWdata/double(length(inputyear));
                
                ind = ind + 1;
                toc;
            end
        end
        comb_interped_clim_divided_U=reshape(comb_interped_Udata,[len_lon, len_lat, 12, length(inputyear)]);
        comb_interped_clim_divided_V=reshape(comb_interped_Vdata,[len_lon, len_lat, 12, length(inputyear)]);
        comb_interped_clim_divided_SW=reshape(comb_interped_SWdata,[len_lon, len_lat, 12, length(inputyear)]);
        
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend_U(1:len_lon,1:len_lat)=NaN;
        trend_V(1:len_lon,1:len_lat)=NaN;
        trend_SW(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_Udata(i,j,:))',1);
                trend_U(i,j)=p(1);
                p=polyfit(trendtime,squeeze(comb_interped_Vdata(i,j,:))',1);
                trend_V(i,j)=p(1);
                p=polyfit(trendtime,squeeze(comb_interped_SWdata(i,j,:))',1);
                trend_SW(i,j)=p(1);
            end
        end

        for t=1:length(inputyear)
            comb_interped_data_filtered_U(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_Udata(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_mean_U;
            comb_interped_data_filtered_V(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_Vdata(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_mean_V;
            comb_interped_data_filtered_SW(:,:,(t-1)*12+1:(t-1)*12+12)=comb_interped_SWdata(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_mean_SW;
        end

        trend_filtered_U(1:len_lon,1:len_lat)=NaN;
        trend_filtered_V(1:len_lon,1:len_lat)=NaN;
        trend_filtered_SW(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered_U(i,j,:))',1);
                trend_filtered_U(i,j)=p(1);
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered_V(i,j,:))',1);
                trend_filtered_V(i,j)=p(1);
                p=polyfit(trendtime,squeeze(comb_interped_data_filtered_SW(i,j,:))',1);
                trend_filtered_SW(i,j)=p(1);
            end
        end
        
        clim_interped_trend_divided_U(1:len_lon,1:len_lat,1:12)=NaN;
        clim_interped_trend_divided_V(1:len_lon,1:len_lat,1:12)=NaN;
        clim_interped_trend_divided_SW(1:len_lon,1:len_lat,1:12)=NaN;
        for k=1:12
            clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
            for i=1:len_lon
                for j=1:len_lat
                    p=polyfit(clim_trendtime(:,k)',squeeze(comb_interped_clim_divided_U(i,j,k,:))',1);
                    clim_interped_trend_divided_U(i,j,k)=p(1);
                    p=polyfit(clim_trendtime(:,k)',squeeze(comb_interped_clim_divided_V(i,j,k,:))',1);
                    clim_interped_trend_divided_V(i,j,k)=p(1);
                    p=polyfit(clim_trendtime(:,k)',squeeze(comb_interped_clim_divided_SW(i,j,k,:))',1);
                    clim_interped_trend_divided_SW(i,j,k)=p(1);
                end
            end
        end
   
        
        mean_trend_U=mean(mean(trend_U,'omitnan'),'omitnan');
        mean_trend_V=mean(mean(trend_V,'omitnan'),'omitnan');
        mean_trend_SW=mean(mean(trend_SW,'omitnan'),'omitnan');
        mean_trend_filtered_U=mean(mean(trend_filtered_U,'omitnan'),'omitnan');
        mean_trend_filtered_V=mean(mean(trend_filtered_V,'omitnan'),'omitnan');
        mean_trend_filtered_SW=mean(mean(trend_filtered_SW,'omitnan'),'omitnan');

        % % pcolor(trend_filtered')
        % % shading interp
        % % mean(mean(trend,'omitnan'),'omitnan')
        % % colorbar;



%         figdir=[figrawdir,'CLIM\'];
%         outfile = strcat(figdir,regionname);
%         if (exist(strcat(figdir) , 'dir') ~= 7)
%             mkdir(strcat(figdir));
%         end 
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

        save([filedir,testname,'_',regionname,'wind_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        ncid = netcdf.create(strcat(filedir, testname,'_',regionname, '_wind_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
        xidimid = netcdf.defDim(ncid, 'xi_rho', len_lon_model);
        etadimid = netcdf.defDim(ncid,'eta_rho',len_lat_model);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['NWP 1/20 _ ', testname, 'model, recon monthly SSH analysis file']);
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

        raw_uvarid=netcdf.defVar(ncid, 'raw_Uwind', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_uvarid,'long_name','raw_Uwind');
        netcdf.putAtt(ncid,raw_uvarid,'units','m/s');
        
        raw_vvarid=netcdf.defVar(ncid, 'raw_Vwind', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_vvarid,'long_name','raw_Vwind');
        netcdf.putAtt(ncid,raw_vvarid,'units','m/s');
        
        raw_swvarid=netcdf.defVar(ncid, 'raw_SWwind', 'NC_FLOAT', [xidimid etadimid time_dimid]);
        netcdf.putAtt(ncid,raw_swvarid,'long_name','raw_SWwind');
        netcdf.putAtt(ncid,raw_swvarid,'units','m/s');

        uvarid=netcdf.defVar(ncid, 'Uwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,uvarid,'long_name','Uwind');
        netcdf.putAtt(ncid,uvarid,'units','m/s');
        
        vvarid=netcdf.defVar(ncid, 'Vwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,vvarid,'long_name','Vwind');
        netcdf.putAtt(ncid,vvarid,'units','m/s');
        
        swvarid=netcdf.defVar(ncid, 'SWwind', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,swvarid,'long_name','SWwind');
        netcdf.putAtt(ncid,swvarid,'units','m/s');

        u_filteredvarid=netcdf.defVar(ncid, 'Uwind_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,u_filteredvarid,'long_name','Uwind_filtered');
        netcdf.putAtt(ncid,u_filteredvarid,'units','m/s');
        
        v_filteredvarid=netcdf.defVar(ncid, 'Vwind_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,v_filteredvarid,'long_name','Vwind_filtered');
        netcdf.putAtt(ncid,v_filteredvarid,'units','m/s');
        
        sw_filteredvarid=netcdf.defVar(ncid, 'SWwind_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,sw_filteredvarid,'long_name','SWwind_filtered');
        netcdf.putAtt(ncid,sw_filteredvarid,'units','m/s');
        
        utrendvarid=netcdf.defVar(ncid, 'trend_Uwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,utrendvarid,'long_name','trend_Uwind');
        netcdf.putAtt(ncid,utrendvarid,'units','m/s/year');
        
        vtrendvarid=netcdf.defVar(ncid, 'trend_Vwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,vtrendvarid,'long_name','trend_Vwind');
        netcdf.putAtt(ncid,vtrendvarid,'units','m/s/year');
        
        swtrendvarid=netcdf.defVar(ncid, 'trend_SWwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,swtrendvarid,'long_name','trend_SWwind');
        netcdf.putAtt(ncid,swtrendvarid,'units','m/s/year');

        utrend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered_Uwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,utrend_filteredvarid,'long_name','trend_filtered_Uwind');
        netcdf.putAtt(ncid,utrend_filteredvarid,'units','m/s/year');
        
        vtrend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered_Vwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,vtrend_filteredvarid,'long_name','trend_filtered_Vwind');
        netcdf.putAtt(ncid,vtrend_filteredvarid,'units','m/s/year');
        
        swtrend_filteredvarid=netcdf.defVar(ncid, 'trend_filtered_SWwind', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,swtrend_filteredvarid,'long_name','trend_filtered_SWwind');
        netcdf.putAtt(ncid,swtrend_filteredvarid,'units','m/s/year');
        
        clim_model_utrend_dividedvarid=netcdf.defVar(ncid, 'clim_model_trend_divided_Uwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_model_utrend_dividedvarid,'long_name','clim_model_trend_divided_Uwind');
        netcdf.putAtt(ncid,clim_model_utrend_dividedvarid,'units','m/s/year');
        
        clim_model_vtrend_dividedvarid=netcdf.defVar(ncid, 'clim_model_trend_divided_Vwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_model_vtrend_dividedvarid,'long_name','clim_model_trend_divided_Vwind');
        netcdf.putAtt(ncid,clim_model_vtrend_dividedvarid,'units','m/s/year');
        
        clim_model_swtrend_dividedvarid=netcdf.defVar(ncid, 'clim_model_trend_divided_SWwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_model_swtrend_dividedvarid,'long_name','clim_model_trend_divided_SWwind');
        netcdf.putAtt(ncid,clim_model_swtrend_dividedvarid,'units','m/s/year');
       
        
        mean_utrendvarid=netcdf.defVar(ncid, 'mean_trend_Uwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_utrendvarid,'long_name','mean_trend_Uwind');
        netcdf.putAtt(ncid,mean_utrendvarid,'units','m/s/year');
        
        mean_vtrendvarid=netcdf.defVar(ncid, 'mean_trend_Vwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_vtrendvarid,'long_name','mean_trend_Vwind');
        netcdf.putAtt(ncid,mean_vtrendvarid,'units','m/s/year');
        
        mean_swtrendvarid=netcdf.defVar(ncid, 'mean_trend_SWwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_swtrendvarid,'long_name','mean_trend_SWwind');
        netcdf.putAtt(ncid,mean_swtrendvarid,'units','m/s/year');

        mean_utrend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered_Uwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_utrend_filteredvarid,'long_name','mean_trend_filtered_Uwind');
        netcdf.putAtt(ncid,mean_utrend_filteredvarid,'units','m/s/year');
        
        mean_vtrend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered_Vwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_vtrend_filteredvarid,'long_name','mean_trend_filtered_Vwind');
        netcdf.putAtt(ncid,mean_vtrend_filteredvarid,'units','m/s/year');
        
        mean_swtrend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered_SWwind', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_swtrend_filteredvarid,'long_name','mean_trend_filtered_SWwind');
        netcdf.putAtt(ncid,mean_swtrend_filteredvarid,'units','m/s/year')

        clim_uvarid=netcdf.defVar(ncid, 'clim_Uwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_uvarid,'long_name','clim_Uwind');
        netcdf.putAtt(ncid,clim_uvarid,'units','m/s');
        
        clim_vvarid=netcdf.defVar(ncid, 'clim_Vwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_vvarid,'long_name','clim_Vwind');
        netcdf.putAtt(ncid,clim_vvarid,'units','m/s');
        
        clim_swvarid=netcdf.defVar(ncid, 'clim_SWwind', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_swvarid,'long_name','clim_SWwind');
        netcdf.putAtt(ncid,clim_swvarid,'units','m/s');

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
        netcdf.putVar(ncid, lonvarid, 0, len_lon, recon_lon);
        netcdf.putVar(ncid, latvarid, 0, len_lat, recon_lat);
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
        netcdf.putVar(ncid, raw_uvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_Udata);
        netcdf.putVar(ncid, raw_vvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_Vdata);
        netcdf.putVar(ncid, raw_swvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_SWdata);
        netcdf.putVar(ncid, uvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_Udata);
        netcdf.putVar(ncid, vvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_Vdata);
        netcdf.putVar(ncid, swvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_SWdata);
        netcdf.putVar(ncid, u_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered_U);
        netcdf.putVar(ncid, v_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered_V);
        netcdf.putVar(ncid, sw_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_interped_data_filtered_SW);
        netcdf.putVar(ncid, utrendvarid, [0 0], [len_lon len_lat], trend_U);
        netcdf.putVar(ncid, vtrendvarid, [0 0], [len_lon len_lat], trend_V);
        netcdf.putVar(ncid, swtrendvarid, [0 0], [len_lon len_lat], trend_SW);
        netcdf.putVar(ncid, utrend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered_U);
        netcdf.putVar(ncid, vtrend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered_V);
        netcdf.putVar(ncid, swtrend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered_SW);
        netcdf.putVar(ncid, clim_model_utrend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_interped_trend_divided_U);
        netcdf.putVar(ncid, clim_model_vtrend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_interped_trend_divided_V);
        netcdf.putVar(ncid, clim_model_swtrend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_interped_trend_divided_SW);
        netcdf.putVar(ncid, mean_utrendvarid, [0], [1], mean_trend_U);
        netcdf.putVar(ncid, mean_vtrendvarid, [0], [1], mean_trend_V);
        netcdf.putVar(ncid, mean_swtrendvarid, [0], [1], mean_trend_SW);
        netcdf.putVar(ncid, mean_utrend_filteredvarid, [0], [1], mean_trend_filtered_U);
        netcdf.putVar(ncid, mean_vtrend_filteredvarid, [0], [1], mean_trend_filtered_V);
        netcdf.putVar(ncid, mean_swtrend_filteredvarid, [0], [1], mean_trend_filtered_SW);
        netcdf.putVar(ncid, clim_uvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_mean_U);
        netcdf.putVar(ncid, clim_vvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_mean_V);
        netcdf.putVar(ncid, clim_swvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_mean_SW);
        
        netcdf.close(ncid);
    end
end
% SSH_4th_mid_report7