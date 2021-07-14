close all; clear all;  clc;
% %  get cmemsstructed SSH. compare model and reSSH. save. 

% all_region ={'AKP4'}
% all_region ={'YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}

% all_region ={'NWP', 'AKP4'}
all_region ={'ES_deep', 'NWES', 'NEES', 'SES_deep'};
% all_region ={'ES_deep', 'SES_deep', 'NES_deep', 'CES_deep', 'NES'};

% all_testname = {'test11', 'test12'};
all_testname = {'v52'};
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
        dl=1/30;

        % for snu_desktop
        testname=all_testname{testnameind} 
        inputyear = [1994:2012]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...];

        varname ='wind';
        regionname=all_region{regionind}
        run('es_polygon_point.m');
        
% % %         switch region
        for folding=1:1
            switch(regionname)
                case('ES') %% East Sea
                    refpolygon=espolygon;
                case('ES_deep') %% East Sea
                    refpolygon=es_deeppolygon;
                case('ES_KHOA') %% East Sea
                    refpolygon=es_khoapolygon;
                case('NES') %% Northern East Sea
                    refpolygon=nespolygon;
                case('NWES') %% North western East Sea
                    refpolygon=nwespolygon;
                case('NEES') %% North eastern East Sea
                    refpolygon=neespolygon;
                case('CES_deep') %% Northern East Sea (Center)
                    refpolygon=ces_deeppolygon;
                case('NES_deep') %% Northern East Sea
                    refpolygon=nes_deeppolygon;
                case('SES_deep') %% Southern East Sea
                    refpolygon=ses_deeppolygon;
                case('SES') %% Southern East Sea
                    refpolygon=sespolygon;
                case('CA') %% Coastal Area around korea peninsula
                    refpolygon=capolygon;
                case('EKB') %% Coastal Area around korea peninsula
                    refpolygon=ekbpolygon;
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
            param_script =['C:\Users\USER\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\es_1_30\run\fig_param\fig_param_kyy_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\es_1_30\', testname, '\output\run\short_monthly\'); % % where data files are
            savedir = strcat('D:\Data\Model\ROMS\es_1_30\', testname, '\output\run\short_monthly\'); % % where data files are
            cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='get model data and cmemsstructed data';
            fig_flags{2,1}='analysis of atm variable and mean energies';

            for flagi=1:10
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=1;
            fig_flags{2,2}=1;
            fig_flags{3,2}=0;
            fig_flags{4,2}=0;
            fig_flags{5,2}=0;

        end
        
        variable='wind';
% % %         get model data and cmemsstructed data      
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            run(param_script);
            ind=1;
            matname = [savedir,testname,'_',regionname,'model_atm_variable_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat'];

            if (exist(matname , 'file') ~= 2 || fig_flag==2)   
                for yearij = 1:length(inputyear)
                    for monthij = 1:length(inputmonth)
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                        tic;
                        tempyear = inputyear(yearij);
                        tempmonth = inputmonth(monthij);
                        % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
                        filename = strcat(filedir, num2str(tempyear,'%04i'), '\short_', ...
                                testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');
                        gridname = [filedir, 'grid_deep_eastsea_', testname, '.nc'];
                        % read model data
                        if (exist('lon')==0)
                            modelinfo=ncinfo(gridname);
                            lon = ncread(gridname,'lon_rho');
                            lat = ncread(gridname,'lat_rho');
                            mask_rho = ncread(gridname,'mask_rho');
                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_model2(1:size(lon,1),1:size(lon,2))=1;
                                otherwise
                                    mask_model2 = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model2(mask_model2==0)=NaN;
                            end
                            mask_rho2=mask_rho.*mask_model2;

                            size_lon=size(mask_rho2,1);
                            size_lat=size(mask_rho2,2);
                            lon_min=min(mod(find(mask_rho2(:,:)==1),size_lon));
                            lon_max=max(mod(find(mask_rho2(:,:)==1),size_lon));
                            lat_min=min(mod(find(mask_rho2(:,:)'==1),size_lat));
                            lat_max=max(mod(find(mask_rho2(:,:)'==1),size_lat));
                            
                            lon = ncread(gridname,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            lat = ncread(gridname,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            mask_rho = ncread(gridname,'mask_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            mask_u = ncread(gridname,'mask_u', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            mask_v = ncread(gridname,'mask_v', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                            
                            pcolor(lon,lat,mask_rho);
                            shading flat;
                            switch(regionname)
                                case('NWP') %% North western Pacific
                                    mask_model(1:size(lon,1),1:size(lon,2))=1;
                                otherwise
                                    mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                    mask_model(mask_model==0)=NaN;
                            end
                            mask_model(mask_rho==0)=NaN;
                            mask_model(mask_u==0)=NaN;
                            mask_model(mask_v==0)=NaN;
                        end


                        uwind_info = ncinfo(filename, 'Uwind');  %% [lon lat depth time] -> [1601 1201 33 1]

                        uwind = ncread(filename,'Uwind',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        uwind=uwind.*mask_model;
                        vwind = ncread(filename,'Vwind',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        vwind=vwind.*mask_model;
                        shflux = ncread(filename,'shflux',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        shflux=shflux.*mask_model;
                        latent = ncread(filename,'latent',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        latent=latent.*mask_model;
                        sensible = ncread(filename,'sensible',[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                        sensible=sensible.*mask_model;
                        
                        speed=sqrt(uwind.^2 + vwind.^2);
                        
                        pcolor(lon,lat,uwind);
                        shading flat;

                        len_lon_model = size(uwind,1);
                        len_lat_model = size(uwind,2);

                        if (exist('comb_spatial_meanmodel_uwind')==0)
                            comb_spatial_meanmodel_uwind=(zeros([len_lon_model,len_lat_model,12]));
                            comb_spatial_meanmodel_vwind=(zeros([len_lon_model,len_lat_model,12]));
                            comb_spatial_meanmodel_shflux=(zeros([len_lon_model,len_lat_model,12]));
                            comb_spatial_meanmodel_latent=(zeros([len_lon_model,len_lat_model,12]));
                            comb_spatial_meanmodel_sensible=(zeros([len_lon_model,len_lat_model,12]));
                            comb_spatial_meanmodel_speed=(zeros([len_lon_model,len_lat_model,12]));
%                             del=(zeros([len_lon_model,len_lat_model,12]));
                        end

                        comb_uwind(:,:,ind) = single(uwind);
                        comb_vwind(:,:,ind) = single(vwind);
                        comb_speed(:,:,ind) = single(speed);
                        comb_shflux(:,:,ind) = single(shflux);
                        comb_latent(:,:,ind) = single(latent);
                        comb_sensible(:,:,ind) = single(sensible);
                        
                        mean_uwind(ind) = mean(uwind(:), 'omitnan');
                        mean_vwind(ind) = mean(vwind(:), 'omitnan');
                        mean_shflux(ind) = mean(shflux(:), 'omitnan');
                        mean_latent(ind) = mean(latent(:), 'omitnan');
                        mean_sensible(ind) = mean(sensible(:), 'omitnan');
                        mean_speed(ind) = mean(speed(:), 'omitnan');
                        
                        comb_spatial_meanmodel_uwind(:,:,monthij)=comb_spatial_meanmodel_uwind(:,:,monthij)+uwind/double(length(inputyear));
                        comb_spatial_meanmodel_vwind(:,:,monthij)=comb_spatial_meanmodel_vwind(:,:,monthij)+vwind/double(length(inputyear));
                        comb_spatial_meanmodel_shflux(:,:,monthij)=comb_spatial_meanmodel_shflux(:,:,monthij)+shflux/double(length(inputyear));
                        comb_spatial_meanmodel_latent(:,:,monthij)=comb_spatial_meanmodel_latent(:,:,monthij)+latent/double(length(inputyear));
                        comb_spatial_meanmodel_sensible(:,:,monthij)=comb_spatial_meanmodel_sensible(:,:,monthij)+sensible/double(length(inputyear));
                        comb_spatial_meanmodel_speed(:,:,monthij)=comb_spatial_meanmodel_speed(:,:,monthij)+speed/double(length(inputyear));


%                         interped_data = griddata(double(lon), double(lat), data,double(cmems_lon2),double(cmems_lat2));   

%                         comb_interped_data(:,:,ind) = interped_data;
%                         comb_spatial_meaninterped(:,:,monthij)=comb_spatial_meaninterped(:,:,monthij)+interped_data/double(length(inputyear));

                        ind = ind + 1;
                        toc;
                    end
                end
                save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', ...
                    'comb_spatial_meanmodel_uwind', 'comb_spatial_meanmodel_vwind', 'comb_spatial_meanmodel_shflux', ...
                    'comb_spatial_meanmodel_latent', 'comb_spatial_meanmodel_sensible', 'comb_spatial_meanmodel_speed',...
                    'comb_uwind', 'comb_vwind', 'comb_shflux', 'comb_latent', 'comb_sensible', 'comb_speed', ...
                    'mean_uwind', 'mean_vwind', 'mean_shflux', 'mean_latent', 'mean_sensible', 'mean_speed', ...
                     'lon', 'lat', 'mask_rho', 'mask_model', '-v7.3');
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
    
% % %     analysis of atm var and mean energies
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            corrdir=[savedir,'\corr'];
            if (exist(strcat(corrdir) , 'dir') ~= 7)
                mkdir(strcat(corrdir));
            end 
            load(matname);
%             load([filedir, '\..\MKE\monthly_MKE_199401.mat']);  %MKE
            stddepth=[-1,-10,-20,-30,-50,-75,-100,-125,-150,-200,-300,-400,-500,-1000,-1500,-2000,-2500,-3000];
%             ana_stddepth=-[1, 50, 100, 150, 200, 300, 400, 500, 1000, 1500, 2000, 2500, 3000];
%             ana_stddepth=-[1, 500, 1000, 2000];
            ana_stddepth=stddepth;
            for depthi=1:length(ana_stddepth)
                tempdepth=ana_stddepth(depthi);
                stddepthi=find(stddepth==tempdepth);
                tempdepthstr=num2str(-tempdepth,'%04i');
%                 variable_energy={'MKE', 'MPE', 'EKE', 'EPE'};
                variable_energy={'MKE', 'MPE', 'EKE', 'EPE'};
                variable_atm={'speed', 'shflux'};
                for var_atmi=1:length(variable_atm)
                    var_atmname=variable_atm{var_atmi};
                    for var_ei=1:length(variable_energy)
                        ind=1;
                        var_ename=variable_energy{var_ei};

                        ncoutfilename = strcat(corrdir, '\', testname,'_',regionname, '_', var_atmname, ...
                            '_corr_',var_ename, '_', tempdepthstr, 'm_', num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc')
                        if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)
                            
                            gridname = [filedir, 'grid_deep_eastsea_', testname, '.nc'];
                            % read model data
%                             if (exist('lon_min')==0)
                                modelinfo=ncinfo(gridname);
                                lon = ncread(gridname,'lon_rho');
                                lat = ncread(gridname,'lat_rho');
                                mask_rho = ncread(gridname,'mask_rho');
                                switch(regionname)
                                    case('NWP') %% North western Pacific
                                        mask_model2(1:size(lon,1),1:size(lon,2))=1;
                                    otherwise
                                        mask_model2 = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                        mask_model2(mask_model2==0)=NaN;
                                end
                                mask_rho2=mask_rho.*mask_model2;

                                size_lon=size(mask_rho2,1);
                                size_lat=size(mask_rho2,2);
                                lon_min=min(mod(find(mask_rho2(:,:)==1),size_lon));
                                lon_max=max(mod(find(mask_rho2(:,:)==1),size_lon));
                                lat_min=min(mod(find(mask_rho2(:,:)'==1),size_lat));
                                lat_max=max(mod(find(mask_rho2(:,:)'==1),size_lat));
            %                     mask_model_region=mask_model(lon_min(1):lon_max(1), lat_min(1):lat_max(1));

%                             end
                                cut_lon=lon(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                cut_lat=lat(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                                
                                len_xi=size(cut_lon,1);
                                len_eta=size(cut_lat,2);
                            
                            if (exist('comb_spatial_meanmodel_e')==0)
                                comb_spatial_meanmodel_e=(zeros([len_lon_model,len_lat_model,12]));
                            end
                            for yearij = 1:length(inputyear)
                                for monthij = 1:length(inputmonth)
                                   yearstr=num2str(inputyear(yearij), '%04i');
                                   monthstr=num2str(inputmonth(monthij), '%02i');
                                   var_efilename=[filedir, '\..\', var_ename, '\monthly_', var_ename, '_', yearstr, monthstr, '.mat'];

                                   load(var_efilename);
                                   temp_var_e=eval(['squeeze(', var_ename, '(stddepthi, lat_min(1):lat_max(1), lon_min(1):lon_max(1)))'';']);
                                   fillval=round(temp_var_e(1,1),10);
                                   temp_var_e(temp_var_e==fillval)=NaN;
                                   temp_var_e=temp_var_e.*mask_model;
    %                                eval(['comb_',var_ename, '(:,:,ind)=temp_var_e;']);
                                   comb_e(:,:,ind)=temp_var_e;
                                   comb_spatial_meanmodel_e(:,:,monthij)=comb_spatial_meanmodel_e(:,:,monthij)+temp_var_e/double(length(inputyear));
                                   ind=ind+1;
                                end
                            end

                           % % %   correlation coefficient between atm variable and energies
                            for i=1:len_xi
                                for j=1:len_eta
        %                             eval(['temp_corr=corrcoef(squeeze(comb_curl(i,j,:))'',','squeeze(comb_', var_ename, '(i,j,:))'');']);
                                    eval(['temp_corr=corrcoef(squeeze(comb_', var_atmname, '(i,j,:))'',','squeeze(comb_e(i,j,:))'');']);

%                                     temp_corr=corrcoef(squeeze(comb_curl(i,j,:))',squeeze(comb_e(i,j,:))');
                                    corr_var(i,j)=temp_corr(1,2);
                                end
                            end
                            disp(['corr coef_', var_atmname, ' complete']) 
        % pcolor(cut_lon', cut_lat',corr_curl'); shading flat; colorbar; colormap(bwrmap); caxis([-1 1]);

                    % % %   correlation coefficient between model_climatology and cmems_climatology 
                            for i=1:len_xi
                                for j=1:len_eta
%                                     temp_corr=corrcoef(squeeze(comb_spatial_meanmodel_curl(i,j,:))',squeeze(comb_spatial_meanmodel_e(i,j,:))');
                                    eval(['temp_corr=corrcoef(squeeze(comb_spatial_meanmodel_', var_atmname, '(i,j,:))'',squeeze(comb_spatial_meanmodel_e(i,j,:))'');']);
                                    corr_spatial_mean(i,j)=temp_corr(1,2);
                                end
                            end
                            disp('corr coef_spatial_mean complete')         
        % pcolor(cut_lon', cut_lat',corr_spatial_mean'); shading flat; colorbar; colormap(bwrmap); caxis([-1 1]);
                            
                            eval(['comb_clim_divided_', var_atmname, '= reshape(comb_', var_atmname, ', [len_xi, len_eta, 12, length(inputyear)]);']);
                            comb_clim_divided_e= reshape(comb_e, [len_xi, len_eta, 12, length(inputyear)]);
                            % % %   correlation coefficient between climatological ssh and climatological cmems ssh 
                            for i=1:len_xi
                                for j=1:len_eta
                                    for k=1:12
%                                         temp_corr=corrcoef(squeeze(comb_clim_divided_curl(i,j,k,:))',squeeze(comb_clim_divided_e(i,j,k,:))');
                                        eval(['temp_corr=corrcoef(squeeze(comb_clim_divided_', var_atmname, '(i,j,k,:))'',squeeze(comb_clim_divided_e(i,j,k,:))'');']);
                                        corr_clim(i,j,k)=temp_corr(1,2);
                                    end
                                end
                            end
                            disp('corr coef_clim complete') 
        % pcolor(cut_lon', cut_lat',corr_clim(:,:,2)'); shading flat; colorbar; colormap(bwrmap); caxis([-1 1]);

                            % % %   lag correlation coefficient between curl and energies
                            lag_month=[-3:-1, 1:3];
                            max_lag_month=max(abs(lag_month));
                            for i=1:len_xi
                                for j=1:len_eta
                                    for k=1:length(lag_month)
                                        var_lagtime_range=1+max_lag_month+lag_month(k):size(comb_e,3)-max_lag_month+lag_month(k);
                                        e_lagtime_range=1+max_lag_month:size(comb_e,3)-max_lag_month;
                                        eval(['temp_corr=corrcoef(squeeze(comb_', var_atmname, '(i,j,var_lagtime_range))'',squeeze(comb_e(i,j,e_lagtime_range))'');']);
%                                         temp_corr=corrcoef(squeeze(comb_curl(i,j,var_lagtime_range))',squeeze(comb_e(i,j,e_lagtime_range))');
                                        corr_lag(i,j,k)=temp_corr(1,2);
                                    end
                                end
                            end
                            disp('corr coef_lag complete') 


                        % % %         make ncfile
                            ncid = netcdf.create(ncoutfilename,'NETCDF4');

                            xi_dimid = netcdf.defDim(ncid, 'xi', len_xi);
                            eta_dimid = netcdf.defDim(ncid,'eta', len_eta);
                            time_dimid = netcdf.defDim(ncid, 'time', 0);
                            clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
                            lag_time_dimid = netcdf.defDim(ncid, 'lag_time', length(lag_month));

                            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                                'type', ['ES 1/30 _ ', testname, 'model, monthly ', var_atmname, ' correlation analysis file']);
                            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                                'title', [' monthly analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
                            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                                'source', [' ROMS ES 1/30 data from _ ',testname ]);
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

                            lag_timevarid=netcdf.defVar(ncid, 'lag_time', 'NC_DOUBLE', lag_time_dimid);
                            netcdf.putAtt(ncid,lag_timevarid,'long_name','lag_time');
                            netcdf.putAtt(ncid,lag_timevarid,'units','days since 1900-12-31 00:00:00');
                            netcdf.putAtt(ncid,lag_timevarid,'calendar','gregorian');
                            netcdf.putAtt(ncid,lag_timevarid,'sign','atm variable is criterion, -1 means atm variable is 1 month earlier than e');

                            lonvarid=netcdf.defVar(ncid, 'lon', 'NC_DOUBLE', [xi_dimid eta_dimid]);
                            netcdf.putAtt(ncid,lonvarid,'long_name','longitude');
                            netcdf.putAtt(ncid,lonvarid,'units','degree_east');

                            latvarid=netcdf.defVar(ncid, 'lat', 'NC_DOUBLE', [xi_dimid eta_dimid]);
                            netcdf.putAtt(ncid,latvarid,'long_name','latitude');
                            netcdf.putAtt(ncid,latvarid,'units','degree_north');   

                            corr_varvarid=netcdf.defVar(ncid, 'corr_var', 'NC_FLOAT', [xi_dimid eta_dimid]);
                            netcdf.putAtt(ncid,corr_varvarid,'long_name','corr_interped');
                            netcdf.putAtt(ncid,corr_varvarid,'units',' ');

                            corr_spatial_meanvarid=netcdf.defVar(ncid, 'corr_spatial_mean', 'NC_FLOAT', [xi_dimid eta_dimid]);
                            netcdf.putAtt(ncid,corr_spatial_meanvarid,'long_name','corr_spatial_mean');
                            netcdf.putAtt(ncid,corr_spatial_meanvarid,'units',' ');

                            corr_climvarid=netcdf.defVar(ncid, 'corr_clim', 'NC_FLOAT', [xi_dimid eta_dimid clim_time_dimid]);
                            netcdf.putAtt(ncid,corr_climvarid,'long_name','corr_clim');
                            netcdf.putAtt(ncid,corr_climvarid,'units',' ');

                            corr_lagvarid=netcdf.defVar(ncid, 'corr_lag', 'NC_FLOAT', [xi_dimid eta_dimid lag_time_dimid]);
                            netcdf.putAtt(ncid,corr_lagvarid,'long_name','corr_lag');
                            netcdf.putAtt(ncid,corr_lagvarid,'units',' ');

                            var_evarid=netcdf.defVar(ncid, 'var_e', 'NC_FLOAT', [xi_dimid eta_dimid time_dimid]);
                            netcdf.putAtt(ncid,var_evarid,'long_name',var_ename);

                            netcdf.endDef(ncid);

                            netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
                            netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
                            netcdf.putVar(ncid, lag_timevarid, 0, length(lag_month), lag_month);
                            netcdf.putVar(ncid, lonvarid, [0 0], [len_xi len_eta], cut_lon);
                            netcdf.putVar(ncid, latvarid, [0 0], [len_xi len_eta], cut_lat);

                            netcdf.putVar(ncid, corr_varvarid, [0 0], [len_xi len_eta], corr_var);
                            netcdf.putVar(ncid, corr_spatial_meanvarid, [0 0], [len_xi len_eta], corr_spatial_mean);
                            netcdf.putVar(ncid, corr_climvarid, [0 0 0], [len_xi len_eta length(climtime)], corr_clim);
                            netcdf.putVar(ncid, corr_lagvarid, [0 0 0], [len_xi len_eta length(lag_month)], corr_lag);
                            netcdf.putVar(ncid, var_evarid, [0 0 0], [len_xi len_eta length(ftime)], comb_e);

                            netcdf.close(ncid);
                        end
                    end
                end
            end

            fig_flag=0;
        end
    

    end
end
