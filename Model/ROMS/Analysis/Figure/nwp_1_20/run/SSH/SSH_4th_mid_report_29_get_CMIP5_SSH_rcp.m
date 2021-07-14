close all; clear all;  clc;
% %  get reconstructed SSH. compare model and reSSH. save. 

all_region ={'NWP', 'YS', 'AKP2'}
% all_region ={'AKP2', 'ES'}
all_testname = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
% all_testname = {'NorESM1-M'};

% all_region ={'AKP2'};
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
            dropboxpath='C:/Users/KYY/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
            addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path
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
        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zos'
        regionname=all_region{regionind}
        run('nwp_polygon_point.m');
        switch(regionname)
            case('NWP') %% North western Pacific
                lonlat = [115, 164, 15, 52];  %% whole data area
%                 refpolygon(1,1)=lonlat(1);
%                 refpolygon(2,1)=lonlat(2);
%                 refpolygon(1,2)=lonlat(3);
%                 refpolygon(2,2)=lonlat(4);
                refpolygon=nwppolygon;
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
            
        elseif (strcmp(system_name,'GLNXA64'))
            figrawdir =strcat('/home/kimyy/SSH_project_source/fig/',testname,'/',regionname,'/'); % % where figure files will be saved
            param_script ='/home/kimyy/SSH_project_source/fig_param_kyy_EKB_RMS.m'
            filedir = strcat('/home/auto/ext_hde/CMIP5/', varname, '/rcp45/Omon/', testname, '/'); % % where data files are
            stfiledir = strcat('/home/auto/ext_hde/CMIP5/', 'zostoga', '/rcp45/Omon/', testname, '/'); % % where data files are
            outputdir = '/home/kimyy/SSH_project_source/data/';
        end

        run(param_script);
        
        
        glaciername='/home/auto/IPCC/glacier_data/model_mean.txt';
        fid = fopen(glaciername);
        glacier_all = textscan(fid, '%f%f');
        glacier_time = glacier_all{1}'+.5;
        glacier_sle = glacier_all{2}'*.001;
        
        
% %         next to the last value -> NaN to extrapolated value (If NaN)
% %         [val, idx] = func
        [~,idx_end_yr] = min(abs( glacier_time-inputyear(end)-.5 )); 
        if( isnan(glacier_sle(idx_end_yr+1)) )
            glacier_sle(idx_end_yr+1) = glacier_sle(idx_end_yr) + ... 
                (glacier_sle(idx_end_yr)   -glacier_sle(idx_end_yr-1)) / ... 
                (glacier_time(idx_end_yr)  -glacier_time(idx_end_yr-1)) * ... 
                (glacier_time(idx_end_yr+1)-glacier_time(idx_end_yr));
        end
        
        %---------- land water correction
        landwater_time = [1995.5; 2100.5; 2101.5];
        landwater_sle = [0; 0.05; 0.05*(1+1/105)];
    %---------- ice sheets correction
        icesheets_time = [1995.5; 2100.5; 2101.5];
        icesheets_sle = [0; 0.14; 0.14*(1+1/105)];     % Greenland(0.09) + Antarctic(0.05)
        
        % %     calculate correction values
        interptime = linspace(inputyear(1), inputyear(end), length(inputyear)*12);
        glacier   = interp1( glacier_time,   glacier_sle,   interptime );
        landwater = interp1( landwater_time, landwater_sle, interptime );
        icesheets = interp1( icesheets_time, icesheets_sle, interptime );

        ind=1;
        for yearij = 1:length(inputyear)
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tic;
                tempyear = inputyear(yearij);
                tempmonth = inputmonth(monthij);
                % ex : rootdir/test37/data/2001/test37_monthly_2001_01.nc
                yearstr=num2str(tempyear, '%04i');
                monthstr=num2str(tempmonth, '%04i');
                switch testname
                    case 'IPSL-CM5A-LR'
                        filename = strcat(filedir, varname, '_Omon_', testname, '_rcp45_r1i1p1_', num2str(200601), '-', ...
                            num2str(230012,'%04i'), '.nc');
                        stfilename = strcat(stfiledir, 'zostoga', '_Omon_', testname, '_rcp45_r1i1p1_', num2str(200601), '-', ...
                            num2str(230012,'%04i'), '.nc');
                    otherwise
                        filename = strcat(filedir, varname, '_Omon_', testname, '_rcp45_r1i1p1_', num2str(200601), '-', ...
                            num2str(210012,'%04i'), '.nc');
                        stfilename = strcat(stfiledir, 'zostoga', '_Omon_', testname, '_rcp45_r1i1p1_', num2str(200601), '-', ...
                            num2str(210012,'%04i'), '.nc');
                end
                 

                
                % read model data
                if (exist('lon')==0)
                    modelinfo=ncinfo(filename);
%                     lon = ncread(filename,'lon',[1 1],[inf,1]);
%                     lat = ncread(filename,'lat',[1 1],[1,inf]);
                    lon = ncread(filename,'lon');
                    lat = ncread(filename,'lat');
                    
                    [lon_min, lon_max, lat_min, lat_max] = findind_Y(dl, lonlat, lon, lat, 1);
                    
                    lon = ncread(filename,'lon', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                    lat = ncread(filename,'lat', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                    
%                     switch(regionname)
%                         case('NWP') %% North western Pacific
%                             mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
%                             mask_model(mask_model==0)=NaN;
%                         otherwise
                            mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                            mask_model(mask_model==0)=NaN;
%                     end
                end


                data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]
                    
                data = ncread(filename,varname,[lon_min(1) lat_min(1) ind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                data=data.*mask_model;
                glo_data = ncread(filename, varname, [1 1 ind], [inf inf 1]);
                glo_mean = mean(mean(glo_data, 1, 'omitnan'),2, 'omitnan');
                stdata = ncread(stfilename,'zostoga',ind, 1);
               
                comb_rawdata(:,:,ind) = data;
%                 switch testname
%                     case 'NorESM1-M'
%                         data=data + stdata + glacier(ind) + icesheets(ind) + landwater(ind) -glo_mean;
%                     otherwise
                        data=data + stdata + glacier(ind) + icesheets(ind) + landwater(ind) -glo_mean;
%                 end
                
                len_lon_model = size(data,1);
                len_lat_model = size(data,2);
                len_lon=len_lon_model;
                len_lat=len_lat_model;
                if (exist('comb_spatial_meanmodel')==0)
                    comb_spatial_meanmodel=(zeros([len_lon_model,len_lat_model,12]));
                end
                
                comb_data(:,:,ind) = data;
                comb_stdata(ind) = stdata;
                comb_glo_mean(ind) = glo_mean;
                comb_spatial_meanmodel(:,:,monthij)=comb_spatial_meanmodel(:,:,monthij)+data/double(length(inputyear));
                ind = ind + 1;
                toc;
            end
        end
        
        mlength=(max(inputyear)-min(inputyear)+1)*12;
        gla_trend=polyfit(1:mlength,glacier,1)*12000;  % monthly to yearly trend -> *12 , m to mm -> *1000
        icesheets_trend=polyfit(1:mlength,icesheets,1)*12000;
        landwater_trend=polyfit(1:mlength,landwater,1)*12000;
        steric_trend=polyfit(1:mlength,comb_stdata,1)*12000;
        other_trend=gla_trend(1) + icesheets_trend(1) + landwater_trend(1)
        
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        rawtrend(1:len_lon_model,1:len_lat_model)=NaN;
        for i=1:len_lon_model
            for j=1:len_lat_model
                p=polyfit(trendtime,squeeze(comb_rawdata(i,j,:))',1);
                rawtrend(i,j)=p(1);
            end
        end
        rawtrend = rawtrend * 1000.0; %% m/y -> mm/y
        
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:len_lon_model,1:len_lat_model)=NaN;
        for i=1:len_lon_model
            for j=1:len_lat_model
                p=polyfit(trendtime,squeeze(comb_data(i,j,:))',1);
                trend(i,j)=p(1);
            end
        end
        trend = trend * 1000.0; %% m/y -> mm/y
        
        for t=1:length(inputyear)
            comb_data_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=comb_data(:,:,(t-1)*12+1:(t-1)*12+12)-comb_spatial_meanmodel;
        end
        
        trend_filtered(1:len_lon_model,1:len_lat_model)=NaN;
        for i=1:len_lon_model
            for j=1:len_lat_model
                p=polyfit(trendtime,squeeze(comb_data_filtered(i,j,:))',1);
                trend_filtered(i,j)=p(1);
            end
        end
        trend_filtered = trend_filtered * 1000.0; %% m/y -> mm/y

        mean_trend=mean(mean(trend,'omitnan'),'omitnan');
        mean_trend_filtered=mean(mean(trend_filtered,'omitnan'),'omitnan');
        
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
        trend_clim = trend_clim * 1000.0;
        
        disp('stop');
        
        % % pcolor(trend_filtered')
        % % shading interp
        % % mean(mean(trend,'omitnan'),'omitnan')
        % % colorbar;

%         figdir=[figrawdir,'CLIM/'];
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

        save([outputdir,testname,'_',regionname,'ssh_trend_rcp45_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'], '-v7.3', '-nocompression');

        ncid = netcdf.create(strcat(outputdir, testname,'_',regionname, '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc'),'NETCDF4');

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

        raw_sshvarid=netcdf.defVar(ncid, 'raw_ssh', 'NC_FLOAT', [xidimid etadimid time_dimid]);
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

        mean_trendvarid=netcdf.defVar(ncid, 'mean_trend', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trendvarid,'long_name','mean_trend');
        netcdf.putAtt(ncid,mean_trendvarid,'units','mm/year');

        mean_trend_filteredvarid=netcdf.defVar(ncid, 'mean_trend_filtered', 'NC_FLOAT', onedimid);
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'long_name','mean_trend_filtered');
        netcdf.putAtt(ncid,mean_trend_filteredvarid,'units','mm/year');

        clim_sshvarid=netcdf.defVar(ncid, 'clim_ssh', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_sshvarid,'long_name','clim_ssh');
        netcdf.putAtt(ncid,clim_sshvarid,'units','m');
        
        clim_ssh_trendvarid=netcdf.defVar(ncid, 'clim_ssh_trend', 'NC_FLOAT', [lon_dimid lat_dimid clim_time_dimid]);
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','clim_ssh_trend');
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');
        
        glaciervarid=netcdf.defVar(ncid, 'glacier', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','glacier');
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');
        
        icesheetsvarid=netcdf.defVar(ncid, 'icesheets', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','icesheets');
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');
        
        landwatervarid=netcdf.defVar(ncid, 'landwater', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'long_name','landwater');
        netcdf.putAtt(ncid,clim_ssh_trendvarid,'units','m');
        
        stericvarid=netcdf.defVar(ncid, 'steric', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,stericvarid,'long_name','steric');
        netcdf.putAtt(ncid,stericvarid,'units','m');
        
        glo_meanvarid=netcdf.defVar(ncid, 'glo_mean', 'NC_FLOAT', [time_dimid]);
        netcdf.putAtt(ncid,glo_meanvarid,'long_name','glo_mean');
        netcdf.putAtt(ncid,glo_meanvarid,'units','m');
        
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
        netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
        netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
        netcdf.putVar(ncid, raw_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_data);
        netcdf.putVar(ncid, ssh_filteredvarid, [0 0 0], [len_lon len_lat length(ftime)], comb_data_filtered);
        netcdf.putVar(ncid, trendvarid, [0 0], [len_lon len_lat], trend);
        netcdf.putVar(ncid, trend_filteredvarid, [0 0], [len_lon len_lat], trend_filtered);
        netcdf.putVar(ncid, mean_trendvarid, [0], [1], mean_trend);
        netcdf.putVar(ncid, mean_trend_filteredvarid, [0], [1], mean_trend_filtered);
        netcdf.putVar(ncid, clim_sshvarid, [0 0 0], [len_lon len_lat length(climtime)], comb_spatial_meanmodel);
        netcdf.putVar(ncid, clim_ssh_trendvarid, [0 0 0], [len_lon len_lat length(climtime)], trend_clim);
        netcdf.putVar(ncid, glaciervarid, 0, length(ftime), glacier);
        netcdf.putVar(ncid, icesheetsvarid, 0, length(ftime), icesheets);
        netcdf.putVar(ncid, landwatervarid, 0, length(ftime), landwater);
        netcdf.putVar(ncid, stericvarid, 0, length(ftime), comb_stdata);
        netcdf.putVar(ncid, stericvarid, 0, length(ftime), comb_glo_mean);
        netcdf.close(ncid);
    end
end