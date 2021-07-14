close all; clear all;  clc;   
% all_region ={'NWP','ES', 'SS', 'YS', 'ECS'}
all_region ={'NWP'}
% all_testname = {'IPSL-CM6A-LR'};
% all_testname = {'MIROC6', 'MPI-ESM1-2-LR', 'NorESM2-LM'};
% all_testname = {'CNRM-CM6-1-HR','MPI-ESM1-2-HR'};
all_testname = {'MPI-ESM1-2-HR'};

for testnameind=1:length(all_testname)
    for regionind=1:length(all_region)
        clearvars '*' -except regionind all_region testnameind all_testname

        
        % % % 
        % % % Read Model ssh
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


        shadlev = [0 35];
        rms_shadlev = [0 4];
        bias_shadlev = [-4 4];
        conlev  = 0:5:35;
        dl=1/20;
        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]



        varname ='ssh';  %% reference variable -> temperature
        run('nwp_polygon_point.m');
        regionname=all_region{regionind}
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
            case('EKB') %% North western Pacific
                lonlat = [127, 131, 37, 42];  %% East Korea Bay
                refpolygon(1,1)=lonlat(1);
                refpolygon(2,1)=lonlat(2);
                refpolygon(1,2)=lonlat(3);
                refpolygon(2,2)=lonlat(4);
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

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('E:\Data\Observation\SSH\monthly\'); % % where data files are
            CMIP6dir='E:\Data\Model\CMIP6\';
        end

        scenname='historical';
        varname='zos';
        variable='SSH';

        run(param_script);
        ind=1;
        for yearij = 1:length(inputyear)
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tic;
                tempyear = inputyear(yearij);
                tempmonth = inputmonth(monthij);
                % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc

                % read SSH DATA
                filedir = strcat(CMIP6dir, varname, '\', scenname, '\Omon\', testname, '\'); % % where data files are
                flag_file_in = false;
                list = dir( [ filedir, '\', varname, '*' ]); 
                for kk = 1 : length( list )
                    fname_in    = list(kk).name;
                    fname_split = strsplit( fname_in, {'_','.'} );
                    fyear_str   = strsplit( fname_split{end-1}, '-' );
                    fyear_start = str2num( fyear_str{1}(1:4) );
                    fyear_end   = str2num( fyear_str{2}(1:4) );
                    if( tempyear >= fyear_start && tempyear <= fyear_end &&     ...                 
                            strcmp( fname_split{2}, 'Omon' ) &&         ...
                            strcmp( fname_split{3}, testname ) &&      ...                 
                            strcmp( fname_split{4}, scenname ) )
                        flag_file_in = true;            break;
                    end         
                end         
                if( ~flag_file_in )
                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                end
                filename=[filedir, '\', fname_in];
                tind=(tempyear-fyear_start)*12+tempmonth;

                CMIP6filename = filename;
                if (exist('CMIP6_lon')==0)
                    CMIP6info=ncinfo(CMIP6filename);
                    if (exist('CMIP6_lon' , 'var') ~= 1)
                        switch(testname)
                            case('IPSL-CM6A-LR') 
                                lon_name='nav_lon';
                                lat_name='nav_lat';
                            case('CNRM-CM6-1-HR')
                                lon_name='lon';
                                lat_name='lat';
                            otherwise
                                lon_name='longitude';
                                lat_name='latitude';
                        end
                        CMIP6_lon = ncread(CMIP6filename,lon_name,[1, 1],[inf, inf]);
                        CMIP6_lat = ncread(CMIP6filename,lat_name,[1, 1],[inf, inf]);
                        if (mean(CMIP6_lon(:))>180)
                            lon_inv_ind=CMIP6_lon>180;
                            CMIP6_lon(lon_inv_ind)=180-CMIP6_lon(lon_inv_ind);
                        end
                        if (sum(diff(CMIP6_lon(:,1)))<-50)
                            ind_reverse_lon=1;
                            for loni=1:size(CMIP6_lon,1)
                                CMIP6_lon2(loni,:)=CMIP6_lon(size(CMIP6_lon,1)-loni+1,:);
                                CMIP6_lat2(loni,:)=CMIP6_lat(size(CMIP6_lon,1)-loni+1,:);
                            end
                            CMIP6_lon=CMIP6_lon2;
                            CMIP6_lat=CMIP6_lat2;
                            clear CMIP6_lon2 CMIP6_lat2
                        else
                            ind_reverse_lon=0;
                        end
                        
                        if (sum(diff(CMIP6_lat(1,:)))<0)
                            ind_reverse_lat=1;
                            for lati=1:size(CMIP6_lat,2)
                                CMIP6_lon2(:,lati)=CMIP6_lon(:,size(CMIP6_lat,2)-lati+1);
                                CMIP6_lat2(:,lati)=CMIP6_lat(:,size(CMIP6_lat,2)-lati+1);
                            end
                            CMIP6_lon=CMIP6_lon2;
                            CMIP6_lat=CMIP6_lat2;
                            clear CMIP6_lon2 CMIP6_lat2
                        else
                            ind_reverse_lat=0;
                        end
                        [CMIP6_lon_min, CMIP6_lon_max, CMIP6_lat_min, CMIP6_lat_max] = findind_Y(1/20, lonlat(1:4), CMIP6_lon, CMIP6_lat, 1);
                    end

%                     CMIP6_lon = ncread(CMIP6filename,lon_name, [CMIP6_lon_min(1), CMIP6_lat_min(1)], [CMIP6_lon_max(1)-CMIP6_lon_min(1), CMIP6_lat_max(1)-CMIP6_lat_min(1)]);
%                     CMIP6_lat = ncread(CMIP6filename,lat_name, [CMIP6_lon_min(1), CMIP6_lat_min(1)], [CMIP6_lon_max(1)-CMIP6_lon_min(1), CMIP6_lat_max(1)-CMIP6_lat_min(1)]);
                    CMIP6_lon = CMIP6_lon(CMIP6_lon_min(1):CMIP6_lon_max(1), CMIP6_lat_min(1):CMIP6_lat_max(1));
                    CMIP6_lat = CMIP6_lat(CMIP6_lon_min(1):CMIP6_lon_max(1), CMIP6_lat_min(1):CMIP6_lat_max(1));

                    comb_spatial_meanrms=(zeros([size(CMIP6_lon,1),size(CMIP6_lat,2),12]));
                    comb_spatial_meanbias=(zeros([size(CMIP6_lon,1),size(CMIP6_lat,2),12]));
                    comb_spatial_meanCMIP6=(zeros([size(CMIP6_lon,1),size(CMIP6_lat,2),12]));
                    comb_spatial_meanmodel=(zeros([size(CMIP6_lon,1),size(CMIP6_lat,2),12]));

                    switch(regionname)
                        case('NWP') %% North western Pacific
                            mask_CMIP6(1:size(CMIP6_lon,1),1:size(CMIP6_lon,2))=1;
                        otherwise
                            mask_CMIP6 = double(inpolygon(CMIP6_lon,CMIP6_lat,refpolygon(:,1),refpolygon(:,2)));
                            mask_CMIP6(mask_CMIP6==0)=NaN;
                    end
                end
                len_lon = size(CMIP6_lon,1);
                len_lat = size(CMIP6_lat,2);
                
                if(ind_reverse_lon==0 && ind_reverse_lat==0)
                    CMIP6_zos = ncread(CMIP6filename,varname,[CMIP6_lon_min(1) CMIP6_lat_min(1) tind], [CMIP6_lon_max(1)-CMIP6_lon_min(1)+1 CMIP6_lat_max(1)-CMIP6_lat_min(1)+1 1]);
                else
                    CMIP6_zos2 = ncread(CMIP6filename, varname, [1 1 tind], [inf inf 1]);
                    if (ind_reverse_lon==1)
                        for loni=1:size(CMIP6_lon,1)
                            CMIP6_zos3(loni,:,1)=CMIP6_zos2(size(CMIP6_lon,1)-loni+1,:,1) ;
                        end
                        CMIP6_zos=CMIP6_zos3(:,CMIP6_lat_min(1):CMIP6_lat_max(1),1);
                    end
                    if (ind_reverse_lat==1)
                        for lati=1:size(CMIP6_lon,2)
%                             CMIP6_zos3(:,lati,1)=CMIP6_zos2(:,CMIP6_lat_max(1)-lati+1,1) ;
%                             CMIP6_zos3(:,lati,1)=CMIP6_zos2(:,lati,1) ;
                            CMIP6_zos3(:,lati,1)=CMIP6_zos2(:,size(CMIP6_zos2,2)- CMIP6_lat_min(1)-lati+1,1) ;
                        end
                        CMIP6_zos=CMIP6_zos3(CMIP6_lon_min(1):CMIP6_lon_max(1),:,1);
                    end
                    clear CMIP6_zos2 CMIP6_zos3
                end
                CMIP6_zos(CMIP6_zos<-1000)=NaN;
                CMIP6_zos(CMIP6_zos>1000)=NaN;
                CMIP6_zos=CMIP6_zos.*mask_CMIP6;
                
                % read zostoga DATA
                ztfiledir = strcat(CMIP6dir, 'zostoga', '\', scenname, '\Omon\', testname, '\'); % % where data files are
                ztflag_file_in = false;
                list = dir( [ ztfiledir, '\', 'zostoga', '*' ]); 
                for kk = 1 : length( list )
                    ztfname_in    = list(kk).name;
                    ztfname_split = strsplit( ztfname_in, {'_','.'} );
                    ztfyear_str   = strsplit( ztfname_split{end-1}, '-' );
                    ztfyear_start = str2num( ztfyear_str{1}(1:4) );
                    ztfyear_end   = str2num( ztfyear_str{2}(1:4) );
                    if( tempyear >= ztfyear_start && tempyear <= ztfyear_end &&     ...                 
                            strcmp( ztfname_split{2}, 'Omon' ) &&         ...
                            strcmp( ztfname_split{3}, testname ) &&      ...                 
                            strcmp( ztfname_split{4}, scenname ) )
                        ztflag_file_in = true;            break;
                    end         
                end         
                if( ~flag_file_in )
                    fprintf('Source File for %04i does not Exist. Continue...\n',tempyear);   continue;
                end
                ztfilename=[ztfiledir, '\', ztfname_in];
                tind=(tempyear-fyear_start)*12+tempmonth;
                
                switch(testname)
                    case('IPSL-CM6A-LR') %% North western Pacific
                        CMIP6_zostoga = ncread(ztfilename,'zostoga',[1 1 tind], [1 1 1]);
                    otherwise
                        CMIP6_zostoga = ncread(ztfilename,'zostoga',[tind], [1]);
                end

    %             CMIP6_data = ncread(CMIP6filename,varname,[CMIP6_lon_min(1) CMIP6_lat_min(1) tind], [CMIP6_lon_max(1)-CMIP6_lon_min(1) CMIP6_lat_max(1)-CMIP6_lat_min(1) 1]);
                CMIP6_data = CMIP6_zos+CMIP6_zostoga;

                CMIP6_data(CMIP6_data<-1000)=NaN;
                CMIP6_data(CMIP6_data>1000)=NaN;
                CMIP6_data=CMIP6_data.*mask_CMIP6;

                comb_CMIP6_data(:,:,ind) = CMIP6_data;
                comb_CMIP6_zos(:,:,ind) = CMIP6_zos;
                comb_CMIP6_zostoga(ind) = CMIP6_zostoga;
                comb_spatial_meanCMIP6(:,:,monthij)=comb_spatial_meanCMIP6(:,:,monthij)+CMIP6_data/double(length(inputyear));

                ind = ind + 1;
                toc;
            end
        end
    %     plot(squeeze(mean(mean(comb_CMIP6_data- mean(comb_CMIP6_zostoga(1)),1,'omitnan'),2,'omitnan')))
    %     hold on
    %     plot(squeeze(mean(mean(comb_CMIP6_zos,1,'omitnan'),2,'omitnan')))
    %     hold off
        comb_CMIP6_clim_divided=reshape(comb_CMIP6_data,[len_lon, len_lat, 12, length(inputyear)]);

        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        CMIP6_trend(1:len_lon,1:len_lat)=NaN;
        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_CMIP6_data(i,j,:))',1);
                CMIP6_trend(i,j)=p(1) * 1000.0 ;
            end
        end

        for t=1:length(inputyear)
            comb_CMIP6_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_CMIP6_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanCMIP6;
        end

        CMIP6_trend_filtered(1:len_lon,1:len_lat)=NaN;

        for i=1:len_lon
            for j=1:len_lat
                p=polyfit(trendtime,squeeze(comb_CMIP6_data_filtered(i,j,:))',1);
                CMIP6_trend_filtered(i,j)=p(1) * 1000.0 ;
            end
        end

        clim_CMIP6_trend_filtered(1:len_lon,1:len_lat,1:12)=NaN;
        for k=1:12
            clim_trendtime(:,k)=inputyear(1)+(k-1)/12 : 1 : inputyear(end)+(k-1)/12;
            for i=1:len_lon
                for j=1:len_lat
                    p=polyfit(clim_trendtime(:,k)',squeeze(comb_CMIP6_clim_divided(i,j,k,:))',1);
                    clim_CMIP6_trend_divided(i,j,k)=p(1) * 1000.0 ;
                end
            end
        end

        mean_CMIP6_trend=mean(mean(CMIP6_trend,'omitnan'),'omitnan');
        mean_CMIP6_trend_filtered=mean(mean(CMIP6_trend_filtered,'omitnan'),'omitnan');
        mean_clim_CMIP6_trend_divided=mean(mean(clim_CMIP6_trend_divided,'omitnan'),'omitnan');
    %     figdir=[figrawdir,'CLIM\'];
    %     outfile = strcat(figdir,regionname);
    %     if (exist(strcat(figdir) , 'dir') ~= 7)
    %         mkdir(strcat(figdir));
    %     end 

        save([filedir,regionname,'CMIP6_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

        ncid = netcdf.create(strcat(filedir,regionname,'CMIP6_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

        onedimid = netcdf.defDim(ncid,'one', 1);
        lon_dimid = netcdf.defDim(ncid, 'lon', len_lon);
        lat_dimid = netcdf.defDim(ncid,'lat',len_lat);
        time_dimid = netcdf.defDim(ncid, 'time', 0);
        clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);

        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'type', ['CMIP6(ssh) _ ', 'monthly ssh trend file']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'title', [' CMIP6 ssh trend (',num2str(inputyear(1)),'-',num2str(inputyear(end)),') ']);
        netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
            'source', [' CMIP6(SSH) satellite data ' ]);
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

        lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
        netcdf.putAtt(ncid,lonvarid,'units','degree_east');

        latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,latvarid,'long_name','latitude');
        netcdf.putAtt(ncid,latvarid,'units','degree_north');

        CMIP6_sshvarid=netcdf.defVar(ncid, 'CMIP6_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,CMIP6_sshvarid,'long_name','CMIP6_ssh');
        netcdf.putAtt(ncid,CMIP6_sshvarid,'units','m ');

        CMIP6_zosvarid=netcdf.defVar(ncid, 'CMIP6_zos', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,CMIP6_zosvarid,'long_name','CMIP6_zos');
        netcdf.putAtt(ncid,CMIP6_zosvarid,'units','m ');

        CMIP6_zostogavarid=netcdf.defVar(ncid, 'CMIP6_zostoga', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,CMIP6_zostogavarid,'long_name','CMIP6_toga');
        netcdf.putAtt(ncid,CMIP6_zostogavarid,'units','m ');

        CMIP6_ssh_filteredvarid=netcdf.defVar(ncid, 'CMIP6_ssh_filtered', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
        netcdf.putAtt(ncid,CMIP6_ssh_filteredvarid,'long_name','CMIP6_ssh_filtered');
        netcdf.putAtt(ncid,CMIP6_ssh_filteredvarid,'units','mm ');

        CMIP6_trendvarid=netcdf.defVar(ncid, 'CMIP6_trend', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,CMIP6_trendvarid,'long_name','CMIP6_trend');
        netcdf.putAtt(ncid,CMIP6_trendvarid,'units','mm /year');

        CMIP6_trend_filteredvarid=netcdf.defVar(ncid, 'CMIP6_trend_filtered', 'NC_FLOAT', [lon_dimid lat_dimid]);
        netcdf.putAtt(ncid,CMIP6_trend_filteredvarid,'long_name','CMIP6_trend_filtered');
        netcdf.putAtt(ncid,CMIP6_trend_filteredvarid,'units','mm /year');

        clim_CMIP6_trend_dividedvarid=netcdf.defVar(ncid, 'clim_CMIP6_trend_divided', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_CMIP6_trend_dividedvarid,'long_name','clim_CMIP6_trend_divided');
        netcdf.putAtt(ncid,clim_CMIP6_trend_dividedvarid,'units','mm /year');

        mean_CMIP6_trendvarid=netcdf.defVar(ncid, 'mean_CMIP6_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_CMIP6_trendvarid,'long_name','mean_CMIP6_trend');
        netcdf.putAtt(ncid,mean_CMIP6_trendvarid,'units','mm /year');

        mean_CMIP6_trend_filteredvarid=netcdf.defVar(ncid, 'mean_CMIP6_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_CMIP6_trend_filteredvarid,'long_name','mean_CMIP6_trend_filtered');
        netcdf.putAtt(ncid,mean_CMIP6_trend_filteredvarid,'units','mm /year');

        mean_clim_CMIP6_trend_dividedvarid=netcdf.defVar(ncid, 'mean_clim_CMIP6_trend_divided', 'NC_FLOAT', clim_time_dimid);
        netcdf.putAtt(ncid,mean_clim_CMIP6_trend_dividedvarid,'long_name','mean_clim_CMIP6_trend_divided');
        netcdf.putAtt(ncid,mean_clim_CMIP6_trend_dividedvarid,'units','mm /year');

        clim_CMIP6varid=netcdf.defVar(ncid, 'clim_CMIP6_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_CMIP6varid,'long_name','clim_CMIP6_ssh');
        netcdf.putAtt(ncid,clim_CMIP6varid,'units','mm ');

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
        netcdf.putVar(ncid, lonvarid, [0 0], [len_lon len_lat], CMIP6_lon);
        netcdf.putVar(ncid, latvarid, [0 0], [len_lon len_lat], CMIP6_lat);
        netcdf.putVar(ncid, CMIP6_sshvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_CMIP6_data);
        netcdf.putVar(ncid, CMIP6_zosvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_CMIP6_zos);
        netcdf.putVar(ncid, CMIP6_zostogavarid, [0], [length(ftime)], comb_CMIP6_zostoga);
        netcdf.putVar(ncid, CMIP6_ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_CMIP6_data_filtered);
        netcdf.putVar(ncid, clim_CMIP6varid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanCMIP6);
        netcdf.putVar(ncid, CMIP6_trendvarid, [0 0], [len_lon len_lat], CMIP6_trend);
        netcdf.putVar(ncid, CMIP6_trend_filteredvarid, [0 0], [len_lon len_lat], CMIP6_trend_filtered);
        netcdf.putVar(ncid, clim_CMIP6_trend_dividedvarid, [0 0 0], [len_lon len_lat length(climtime)], clim_CMIP6_trend_divided);
        netcdf.putVar(ncid, mean_CMIP6_trendvarid, [0], [1], mean_CMIP6_trend);
        netcdf.putVar(ncid, mean_CMIP6_trend_filteredvarid, [0], [1], mean_CMIP6_trend_filtered);
        netcdf.putVar(ncid, mean_clim_CMIP6_trend_dividedvarid, [0], [length(climtime)], mean_clim_CMIP6_trend_divided);

        netcdf.close(ncid);
    end
end