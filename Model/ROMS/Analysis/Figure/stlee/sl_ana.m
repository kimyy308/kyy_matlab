close all; clc; clear all;

%% set path
[error_status, tmp.hostname] = system('hostname');
tmp.hostname=tmp.hostname(1:end-1);
switch tmp.hostname
    case 'Yong-Yubs-iMac-Pro.local'
        tmp.dropboxpath = '/Volumes/kyy_raid/kimyy/Dropbox';
    case {'da1', 'da2', 'da3', 'da4'}
        tmp.dropboxpath = '/mnt/lustre/proj/kimyy/Dropbox';
end
tmp.fs=filesep;
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'mca']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'order']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Common', tmp.fs, 'seawater_ver3_2']));
            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);


monfiledir='/mnt/lustre/proj/kimyy/Model/ROMS/drifter_ROMS/monthly/';
savedir='/mnt/lustre/proj/kimyy/Model/ROMS/drifter_ROMS/mat/';

for i=1:11
    for j=1:2
        fig_flags{i,j}=1;
    end
end
inputyear=1982:2019;
inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
regionname='stlee_drifter';


%% polygon set
run('nwp_polygon_point.m');
switch(regionname)
    case('stlee_drifter') %% for debugging
        refpolygon=stlee_drifter_polygon;
    otherwise
        ('?')
end
lonlat(1)=min(refpolygon(:,1));
lonlat(2)=max(refpolygon(:,1));
lonlat(3)=min(refpolygon(:,2));
lonlat(4)=max(refpolygon(:,2));

%% get ts & calculate zosto
fig_flag=fig_flags{11,2};
while (fig_flag)
%     run(param_script);
    ind=1;
    matname = [savedir,regionname,'_model_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), ...
        '_', num2str(min(inputmonth),'%02i'),'_',num2str(max(inputmonth),'%02i'),'.mat'];
    if (exist(matname , 'file') ~= 2 || fig_flag==2) 
        for yearij = 1:length(inputyear)
            tempyear = inputyear(yearij);
            %% file open
            filename = strcat(monfiledir, num2str(tempyear,'%04i'), '/', ...
                       'roms_monthly_avg.nc');
            ncid=netcdf.open(filename, 'NOWRITE');

            %% get initial rho information to caluculate steric sea-leval in relative to initial rho
            if (exist('rho_0_va')==0) %  rho at initial time, yearly mean or seasonal mean
                for monthij = 1:length(inputmonth)
                    disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                    tic;
                    tempmonth = inputmonth(monthij);
                    %% read model data and masking with the polygon information
                    if (exist('lon_min')==0)
                        modelinfo=ncinfo(filename);
                        lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(2).Length,1]);
                        lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(3).Length]);
                        
                        [lon_min, lon_max, lat_min, lat_max] = ...
                            Func_0012_findind_Y(1, lonlat, lon, lat);

                        lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                        lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                        %% masking
                        switch(regionname)
                            case('NWP') %% North western Pacific
                                mask_model(1:size(lon,1),1:size(lon,2))=1;
                            otherwise
                                mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                                mask_model(mask_model==0)=NaN;
                        end
                    end

                    data_info = ncinfo(filename, 'zeta');  %% [lon lat depth time] -> [470 460 40 12]
                    
                    %% read s-coord information
                    if (exist('h')==0)
                        h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                        pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                        pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                        dA=1./pm .* 1./pn;
%                                 Vtransform = ncread(filename, 'Vtransform');
%                                 Vstretching = ncread(filename, 'Vstretching');
%                                 theta_s = ncread(filename, 'theta_s');
%                                 theta_b = ncread(filename, 'theta_b');
%                                 hc = ncread(filename, 'hc');
%                                 N = length(ncread(filename, 's_rho'));
                        Vtransform= 2;
                        Vstretching = 4;
                        theta_s=10;
                        theta_b=1;
                        hc=250;
                        N=40;
                    end

                    %% read zeta
                    varid.zeta=netcdf.inqVarID(ncid, 'zeta');
                    zeta = netcdf.getVar(ncid, varid.zeta, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%                     data = ncread(filename,varname, [lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                    
                    %% read potential temperature
                    varid.temp=netcdf.inqVarID(ncid, 'temp');
%                     PT_src(:,:,:,monthij) = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                    PT_src(:,:,:,monthij) = netcdf.getVar(ncid, varid.temp, [lon_min(1)-1 lat_min(1)-1 0 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 40 1]);
%                     S_src(:,:,:,monthij) = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
                    
                    %% read salinity
                    varid.salt=netcdf.inqVarID(ncid, 'salt');
                    S_src(:,:,:,monthij) = netcdf.getVar(ncid, varid.salt, [lon_min(1)-1 lat_min(1)-1 0 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 40 1]);
                    S_src(S_src<=0)=0;
                    
                    %% get zlev information
                    vtype ='w';
                    zw=zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);

                    dz = diff(zw);
                    dz = permute(dz, [2,3,1]);
%                         vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
%                         area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
                    vtype ='r';

                    %% pressure calculation
                    P = zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);
                    P = -permute(P, [2,3,1]);
                    comb_P(:,:,:,monthij)=P;
                end
                    %% yearly or seasonal mean
                    PT_src = mean(PT_src, 4 );
                    PT_src_0=PT_src;
                    S_src = mean(S_src, 4 );
                    S_src_0=S_src;
                    P = mean(comb_P, 4);
                    T_src = sw_temp(S_src, PT_src, P, zeros(size(P)) );
                    T_src_0 =T_src;
                    rho_src = sw_dens(S_src, T_src, P );
                    dep_0 = h + zeta;
                    rho_0_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_0;       % vertical(local) average, 2-D
            end
            
            %% steric sea-level calculation
            for monthij = 1:length(inputmonth)
                disp([num2str(yearij), 'y_',num2str(monthij),'m'])
                tic;
                tempmonth = inputmonth(monthij);
                
                %% read zeta
                varid.zeta=netcdf.inqVarID(ncid, 'zeta');
                zeta = netcdf.getVar(ncid, varid.zeta, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);

                %% read sustr
                varid.sustr=netcdf.inqVarID(ncid, 'sustr');
                sustr = netcdf.getVar(ncid, varid.sustr, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);

                %% read svstr
                varid.svstr=netcdf.inqVarID(ncid, 'svstr');
                svstr = netcdf.getVar(ncid, varid.svstr, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                
                %% read Uwind
                varid.Uwind=netcdf.inqVarID(ncid, 'Uwind');
                Uwind = netcdf.getVar(ncid, varid.Uwind, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);

                %% read Vwind
                varid.Vwind=netcdf.inqVarID(ncid, 'Vwind');
                Vwind = netcdf.getVar(ncid, varid.Vwind, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);

                %% read ubar
                varid.ubar=netcdf.inqVarID(ncid, 'ubar');
                ubar = netcdf.getVar(ncid, varid.ubar, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);

                %% read vbar
                varid.vbar=netcdf.inqVarID(ncid, 'vbar');
                vbar = netcdf.getVar(ncid, varid.vbar, [lon_min(1)-1 lat_min(1)-1 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
                
                %% read potential temperature
                varid.temp=netcdf.inqVarID(ncid, 'temp');
                PT_src = netcdf.getVar(ncid, varid.temp, [lon_min(1)-1 lat_min(1)-1 0 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 40 1]);
                
                %% read salinity
                varid.salt=netcdf.inqVarID(ncid, 'salt');
                S_src = netcdf.getVar(ncid, varid.salt, [lon_min(1)-1 lat_min(1)-1 0 tempmonth-1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 40 1]);
                S_src(S_src<=0)=0;

                vtype ='w';
                zw=zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);

                dz = diff(zw);
                dz = permute(dz, [2,3,1]);

                %% get steric sea-lev anomaly
                vtype ='r';
                P = zlevs(Vtransform, Vstretching, h,zeta,theta_s,theta_b,hc,N,vtype);
                P = -permute(P, [2,3,1]);
                T_src = sw_temp (S_src, PT_src, P, zeros(size(P)) );
                rho_src = sw_dens (S_src, T_src, P );

                dep_n = h + zeta;
                rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
%                         rho_n_va2 = sum( rho_src.*dz, 3 , 'omitnan') ./ (h+data);
%                         rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
                zosto = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
%                         zosto_2 = dep_n       .* ( 1 - rho_n_va ./ rho_0_va );

                %% get thermosteric sla
                T_src = sw_temp (S_src_0, PT_src, P, zeros(size(P)) );
                rho_src = sw_dens (S_src_0, T_src, P );
                dep_n = h + zeta;
                rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
                zosto_thermo = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );

                %% get halosteric sla
                rho_src = sw_dens (S_src, T_src_0, P );
                dep_n = h + zeta;
                rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
                zosto_halo = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );

                comb_data.zosto(:,:,ind) = zosto;
                comb_data.zosto_thermo(:,:,ind) = zosto_thermo;
                comb_data.zosto_halo(:,:,ind) = zosto_halo;
                comb_data.zeta(:,:,ind) = zeta;
                comb_data.sustr(:,:,ind) = sustr;
                comb_data.svstr(:,:,ind) = svstr;
                comb_data.Uwind(:,:,ind) = Uwind;
                comb_data.Vwind(:,:,ind) = Vwind;
                comb_data.ubar(:,:,ind) = ubar;
                comb_data.vbar(:,:,ind) = vbar;

                len_lon_model = size(zeta,1);
                len_lat_model = size(zeta,2);

                ind = ind + 1;
                toc;
            end
            netcdf.close(ncid);
        end
        
        save(matname, 'mask_model', 'len_lon_model', 'len_lat_model', 'comb_data',  ...
            'lon', 'lat', '-v7.3');
    else
        load(matname);
    end
    fig_flag=0;
end



% % % % % % % %         steric sea level analysis
% % % % %         fig_flag=fig_flags{12,2};
% % % % %         while (fig_flag)
% % % % %             ncoutfilename = strcat(savedir, testname,'_',regionname, '_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
% % % % %             if (exist(ncoutfilename , 'file') ~= 2 || fig_flag==2)   
% % % % %             
% % % % %                 comb_zosto_correct = comb_zosto;
% % % % %                 for i=1:len_lon_model
% % % % %                     for j=1:len_lat_model
% % % % %                         raw_zosto=squeeze(comb_zosto(i,j,:));
% % % % %                         std5_zosto = 5.* std(raw_zosto);
% % % % %                         diff_zosto = diff(raw_zosto);
% % % % %                         if diff_zosto(end)> std5_zosto;
% % % % %                             raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
% % % % %                             diff_zosto(end)=diff_zosto(end-1);
% % % % %                         end
% % % % %                         raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
% % % % %                         continuity_error=diff(find(diff_zosto>std5_zosto));
% % % % %                         continuity_error(continuity_error==1)=55555;
% % % % %                         if continuity_error==55555
% % % % %                             'error'
% % % % %                             break
% % % % %                         end
% % % % %                         comb_zosto_correct(i,j,:)=raw_zosto;
% % % % %                     end
% % % % %                 end
% % % % %                 
% % % % %                 for i=1:size(comb_interped_zosto,1)
% % % % %                     for j=1:size(comb_interped_zosto,2)
% % % % %                         raw_zosto=squeeze(comb_interped_zosto(i,j,:));
% % % % %                         std5_zosto =5.* std(raw_zosto);
% % % % %                         diff_zosto = diff(raw_zosto);
% % % % %                         if diff_zosto(end)> std5_zosto;
% % % % %                             raw_zosto(end) = raw_zosto(end-1) + diff_zosto(end-1);
% % % % %                             diff_zosto(end)=diff_zosto(end-1);
% % % % %                         end
% % % % %                         raw_zosto(find(diff_zosto>std5_zosto)+1)= (raw_zosto(find(diff_zosto>std5_zosto)) + raw_zosto(find(diff_zosto>std5_zosto)+2))./2.0;
% % % % %                         continuity_error=diff(find(diff_zosto>std5_zosto));
% % % % %                         continuity_error(continuity_error==1)=55555;
% % % % %                         if continuity_error==55555
% % % % %                             'error'
% % % % %                             break
% % % % %                         end
% % % % %                         comb_interped_zosto_correct(i,j,:)=raw_zosto;
% % % % %                     end
% % % % %                 end
% % % % %                 
% % % % %                 rawsshfilename = strcat(savedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
% % % % %                 mask_rho = ncread(rawsshfilename, 'trend');
% % % % %                 mask_rho(isfinite(mask_rho))=1;
% % % % % %                 for tt= 1:size(ftime)
% % % % %                     comb_zosto_correct = comb_zosto_correct .* mask_rho;
% % % % % %                 end
% % % % % %                 plot(squeeze(mean(mean(comb_interped_zosto_correct,1,'omitnan'),2,'omitnan')))
% % % % % 
% % % % %             % % %         make ncfile
% % % % %                 ncid = netcdf.create(ncoutfilename,'NETCDF4');
% % % % % 
% % % % %                 lon_dimid = netcdf.defDim(ncid, 'lon', len_lon_model);
% % % % %                 lat_dimid = netcdf.defDim(ncid,'lat',len_lat_model);
% % % % %                 time_dimid = netcdf.defDim(ncid, 'time', 0);
% % % % % %                 clim_time_dimid = netcdf.defDim(ncid, 'clim_time', 12);
% % % % % 
% % % % %                 netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % % % %                     'type', ['NWP 1/20 _ ', testname, 'model monthly SSH analysis file']);
% % % % %                 netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % % % %                     'title', [' steric SSH analysis (', num2str(min(inputyear)), '-', num2str(max(inputyear)) ,') ']);
% % % % %                 netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % % % %                     'source', [' ROMS NWP 1/20 data from _ ',testname ]);
% % % % %                 netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % % % %                     'author', 'Created by Y.Y.Kim');
% % % % %                 netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
% % % % %                     'date', date);
% % % % % 
% % % % %                 timevarid=netcdf.defVar(ncid, 'time', 'NC_DOUBLE', time_dimid);
% % % % %                 netcdf.putAtt(ncid,timevarid,'long_name','time');
% % % % %                 netcdf.putAtt(ncid,timevarid,'units','days since 1900-12-31 00:00:00');
% % % % %                 netcdf.putAtt(ncid,timevarid,'calendar','gregorian');
% % % % % 
% % % % %                 lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % % % %                 netcdf.putAtt(ncid,lon_rhovarid,'long_name','lon_model');
% % % % %                 netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');
% % % % % 
% % % % %                 lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [lon_dimid lat_dimid]);
% % % % %                 netcdf.putAtt(ncid,lat_rhovarid,'long_name','lat_model');
% % % % %                 netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
% % % % % 
% % % % %                 steric_sshvarid=netcdf.defVar(ncid, 'steric_ssh', 'NC_FLOAT', [lon_dimid lat_dimid time_dimid]);
% % % % %                 netcdf.putAtt(ncid,steric_sshvarid,'long_name','steric_ssh');
% % % % %                 netcdf.putAtt(ncid,steric_sshvarid,'units','m');
% % % % %                 netcdf.defVarChunking(ncid,steric_sshvarid,'CHUNKED',[len_lon_model/10, len_lat_model/10, length(ftime)/10]);
% % % % %                 netcdf.defVarDeflate(ncid,steric_sshvarid,true,true,1);
% % % % % 
% % % % %                 netcdf.endDef(ncid);
% % % % % 
% % % % %                 netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
% % % % % %                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);
% % % % %                 netcdf.putVar(ncid, lon_rhovarid, [0 0], [len_lon_model len_lat_model], lon);
% % % % %                 netcdf.putVar(ncid, lat_rhovarid, [0 0], [len_lon_model len_lat_model], lat);
% % % % %                 netcdf.putVar(ncid, steric_sshvarid, [0 0 0], [len_lon_model len_lat_model length(ftime)], comb_zosto_correct);
% % % % %                 netcdf.close(ncid);
% % % % %             end
% % % % %             fig_flag=0;
% % % % %         end