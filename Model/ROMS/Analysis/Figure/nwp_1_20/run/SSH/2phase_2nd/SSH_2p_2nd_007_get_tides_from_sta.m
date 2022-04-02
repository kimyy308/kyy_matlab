close all; clear all; clc;

dropboxpath='/home/kimyy/Dropbox';
addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
addpath(genpath([dropboxpath '/source/matlab/Common/t_tide_v1.3beta']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Roms_tools/Preprocessing_tools']));
addpath(genpath([matlabroot,'/toolbox/matlab/imagesci/'])); %% add new netcdf path

% all_region2 ={'YSECS'}
all_region2 ={'JSB'}
% all_var2 ={'zeta', 'ubar_eastward', 'vbar_northward'};
all_var2 ={'zeta'};
% % set info of station
for folding=1:1
station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

num_sta=length(station.name);

station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

end
% % get sum of 4 major tidal constitute (Observation)

all_testname2 = {'test2107'};
% inputyear=[1985:2014];
inputyear=[2033:2050];

% inputyear=[2100];

for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        for yearind2=1:length(inputyear)
            for varind2=1:length(all_var2)
                
            datetime
            close all;
            clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 inputyear yearind2 varind2 all_var2 varind2
            
            varname = all_var2{varind2};
            % % %         flag configuration
            for folding=1:1
                fig_flags{1,1}='Mokpo M2 time series';

            end
            for flagi=1:8
                fig_flags{flagi,2}=0;
            end

            fig_flags{1,2}=2;
            fig_flags{2,2}=2;
            fig_flags{3,2}=2;
            fig_flags{4,2}=2;

            testname=all_testname2{testnameind2};
            tempyear=inputyear(yearind2);
            switch testname
                case {'test61', 'test62', 'test63', 'test64'}
                    drivename='E:/';
                case {'test57', 'test58', 'test59', 'test60'}
                    drivename='I:/';
                case {'test65', 'test66', 'test67', 'test68'}
                    drivename='G:/';
            end
            outputdir=['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/JSB/', testname, '/',num2str(tempyear)];
            harmonic_outputdir=['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/JSB/', testname, '/harmonic_analysis'];
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
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
                case('SS') %% South Sea
                    refpolygon=sspolygon;
                case('YS') %% Yellow Sea
                    refpolygon=yspolygon;
                case('YSECS') %% Yellow Sea
                    refpolygon=ysecspolygon;
                case('ECS') %% East China Sea
                    refpolygon=ecspolygon;
                case('AKP') %% Around Korea Peninsula
                    refpolygon=akppolygon;
                case('AKP2') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('AKP4') %% Around Korea Peninsula
                    refpolygon=AKP4polygon;
                case('CA') %% Around Korea Peninsula
                    refpolygon=capolygon;
                case('EKB') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('JSB')
                    refpolygon(1,1)=121.5;
                    refpolygon(2,1)=122.5;
                    refpolygon(1,2)=31.6781;
                    refpolygon(2,2)=32.5254;
                otherwise
                    ('?')
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));

            lonfilename=['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/JSB/', testname, '/',num2str(tempyear),'/pck_his_lon_rho_JSB', '.nc'];
            latfilename=['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/JSB/', testname, '/',num2str(tempyear),'/pck_his_lat_rho_JSB', '.nc'];

             if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(lonfilename, 'lon_rho');
                lat_rho=ncread(latfilename, 'lat_rho');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                Model.tcon=cell(lon_max(1)-lon_min(1)+1,lat_max(1)-lat_min(1)+1);
                Model.amp_4con(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
                Model.amp_M2(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
                Model.amp_S2(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
                Model.amp_K1(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
                Model.amp_O1(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
             end
            end
             for lon_loop=lon_min(1) : lon_max(1)
                for lat_loop=lat_min(1) : lat_max(1)
                    loni=lon_loop-lon_min(1)+1;
                    lati=lat_loop-lat_min(1)+1;
                    Model.tcon{loni,lati}=NaN(59,4);
                end
             end
%             if (i == 2)
    %             mx = 28;
%                 tempyear=inputyear(1);
                fullday=365;
                if (mod(tempyear,4)==0)
    %                 mx =29;
                    fullday=366;
                    if (mod(tempyear,100)==0)
    %                     mx = 28;
                    fullday=365;
                    end
                end
    %         elseif i==4 || 6 || 9 || 11;  mx=30;
    %         elseif i==1 || 3 || 5 || 7 || 8 ||10||12; mx=31;
%             end
            tide_info.name{1}='M2  ';
            tide_info.name{2}='S2  ';
            tide_info.name{3}='K1  ';
            tide_info.name{4}='O1  ';

            tempvar(1:fullday*24)=NaN;
            tic;
%             Model.tcon{1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1}=NaN(lon_max(1)-lon_min(1)+1,lat_max(1)-lat_min(1)+1);
            
%             Model.amp_4con(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
%             Model.amp_M2(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
%             Model.amp_S2(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
%             Model.amp_K1(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
%             Model.amp_O1(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
            tic;
            for lon_loop=lon_min(1) : lon_max(1)
                for lat_loop=lat_min(1) : lat_max(1)
                    loni=lon_loop-lon_min(1)+1;
                    lati=lat_loop-lat_min(1)+1;
                    filename=[outputdir,'/pck_', testname, '_', num2str(tempyear), ...
                                '_daily_his_', varname, '_', regionname, '_', num2str(1, '%04i'), '.nc'];
                    varflag=ncread(filename,varname, [lon_loop, lat_loop, 1], [1, 1, 1]);
                    if(~isnan(varflag))
                        for dayi=1:fullday
%                             tic;
                            filename=[outputdir,'/pck_', testname, '_', num2str(tempyear), ...
                                '_daily_his_', varname, '_', regionname, '_', num2str(dayi, '%04i'), '.nc'];
                            tempvar((dayi-1)*24+1:dayi*24)=ncread(filename,varname, [lon_loop, lat_loop, 1], [1, 1, 24]);
%                             ncread(filename,'zeta', [loni, lati, 1], [1, 1, 24]);
%                             toc;
%                             dayi
                        end
                        [tname,tfreq,tcon,tout]=t_tide_noprint(tempvar(:),...
                               'interval',1, ...                     % hourly data
                               'start',datenum(tempyear,1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
                               'latitude',lat_rho(lon_loop,lat_loop),...               % Latitude of Model
                               'rayleigh',1, 'error','wboot', 'output', 'none');   

                        Model.tcon{loni,lati}=tcon;

                        tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio

                        num_tide_all=size(tname,1);
                        num_tide_tgt=length(tide_info.name);
                        for coni=1:num_tide_all
                            for tide_namei=1:num_tide_tgt
                                if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                                    tide_info.index(tide_namei)=coni;
                                end
                            end
                        end
                        Model.amp_4con(loni,lati)= sum(tcon(tide_info.index(:),1));
                        Model.amp_M2(loni,lati)=tcon(tide_info.index(1),1);
                        Model.amp_S2(loni,lati)=tcon(tide_info.index(2),1);
                        Model.amp_K1(loni,lati)=tcon(tide_info.index(3),1);
                        Model.amp_O1(loni,lati)=tcon(tide_info.index(4),1);
                    else
                        Model.tcon{loni,lati}=NaN;
                        Model.amp_4con(loni,lati)=NaN;
                        Model.amp_M2(loni,lati)=NaN;
                        Model.amp_S2(loni,lati)=NaN;
                        Model.amp_K1(loni,lati)=NaN;
                        Model.amp_O1(loni,lati)=NaN;
                    end
                end
                lon_loop/lon_max(1)
            end
            toc;
            save([harmonic_outputdir,'/', testname, '_', 'harmonic_analysis_',varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.mat'])

            ncoutfilename=[harmonic_outputdir,'/', testname, '_', 'harmonic_analysis_', varname, '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];

            % % %         make ncfile
            ncid = netcdf.create(ncoutfilename,'NETCDF4');

            cut_len_lon=size(cut_lon_rho,1);
            cut_len_lat=size(cut_lat_rho,2);
            xi_rho_dimid = netcdf.defDim(ncid, 'lon_rho', cut_len_lon);
            eta_rho_dimid = netcdf.defDim(ncid,'lat_rho', cut_len_lat);
            tide_dimid = netcdf.defDim(ncid, 'tide', size(tname,1));
            tcon_dimid = netcdf.defDim(ncid, 'tcon_kind', size(tname,2));

            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'type', ['NWP 1/20 _ ', testname, 'model, harmonic analysis file']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'title', [' monthly ', varname, ' analysis (', num2str(tempyear) ,') ']);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'source', [' ROMS NWP 1/20 data from _ ',testname ]);
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'author', 'Created by Y.Y.Kim');
            netcdf.putAtt(ncid,netcdf.getConstant('NC_GLOBAL'), ...
                'date', date);

            tidevarid=netcdf.defVar(ncid, 'tname', 'NC_CHAR', [tcon_dimid tide_dimid ]);
            netcdf.putAtt(ncid,tidevarid,'long_name','tide_name');

            lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
            netcdf.putAtt(ncid,lon_rhovarid,'long_name','longitude');
            netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

            lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
            netcdf.putAtt(ncid,lat_rhovarid,'long_name','latitude');
            netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');

            m2_ampvarid=netcdf.defVar(ncid, 'm2_amp', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
            netcdf.putAtt(ncid,m2_ampvarid,'long_name','m2_amplitude');
            netcdf.putAtt(ncid,m2_ampvarid,'units','m');

            tconvarid=netcdf.defVar(ncid, 'tcon', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid tide_dimid tcon_dimid]);
            netcdf.putAtt(ncid,tconvarid,'long_name','tcon(amp, amp_err, phase, phase_err)');

%                 cmems_trendvarid=netcdf.defVar(ncid, 'cmems_trend', 'NC_FLOAT', [lon_rho_dimid lat_rho_dimid]);
%                 netcdf.putAtt(ncid,cmems_trendvarid,'long_name','cmems_trend');
%                 netcdf.putAtt(ncid,cmems_trendvarid,'units','mm /year');

            netcdf.endDef(ncid);

%                 netcdf.putVar(ncid, timevarid, 0, length(ftime), ftime);
%                 netcdf.putVar(ncid, clim_timevarid, 0, length(climtime), climtime);

            netcdf.putVar(ncid, lon_rhovarid, [0 0], [cut_len_lon cut_len_lat], cut_lon_rho);
            netcdf.putVar(ncid, lat_rhovarid, [0 0], [cut_len_lon cut_len_lat], cut_lat_rho);
            netcdf.putVar(ncid, m2_ampvarid, [0 0], [cut_len_lon cut_len_lat], Model.amp_M2);

            netcdf.putVar(ncid, tidevarid, [0 0], [4 size(tname,1)], tname');
            for loni=1:cut_len_lon
                for lati=1:cut_len_lat
                    if isfinite(Model.tcon{loni, lati})
                        comb_tcon(loni,lati,1:size(tname,1), 1:4) = Model.tcon{loni, lati};
                    else
                        comb_tcon(loni,lati,1:size(tname,1), 1:4) = NaN;
                    end
                end
            end

            netcdf.putVar(ncid, tconvarid, [0 0 0 0], [cut_len_lon cut_len_lat size(tname,1) size(tname,2)], comb_tcon);

            netcdf.close(ncid);
            end
        end
    end
end
