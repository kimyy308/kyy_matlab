close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path


all_region2 ={'YSECS'}
% % set info of station

all_testname2 = {'test66'};
inputyear=[2006, 2100];

for testnameind2=1:length(all_testname2)
    for regionind2=1:length(all_region2)
        for yearind2=1:length(inputyear)
%             for varind2=1:length(all_var2)
                
            datetime
            close all;
            clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 inputyear yearind2 varind2 all_var2 varind2
            
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
            outputdir=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear)];
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
                    refpolygon=akp4polygon;
                case('CA') %% Around Korea Peninsula
                    refpolygon=capolygon;
                case('EKB') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                otherwise
                    ('?')
            end
            lonlat(1)=min(refpolygon(:,1));
            lonlat(2)=max(refpolygon(:,1));
            lonlat(3)=min(refpolygon(:,2));
            lonlat(4)=max(refpolygon(:,2));

            lonfilename=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lon_rho_AKP4', '.nc'];
            latfilename=['G:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(tempyear),'\ocean_his_lat_rho_AKP4', '.nc'];

             if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(lonfilename, 'lon_rho');
                lat_rho=ncread(latfilename, 'lat_rho');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
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
            
            all_var2 ={'zeta', 'ubar_eastward', 'vbar_northward'};
            for varind2= 1:length(all_var2)
                varname=all_var2{varind2};
                eval(['temp_', varname, '(1:fullday*24)=NaN;']);
            end
%             Model.tcon{1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1}=NaN(lon_max(1)-lon_min(1)+1,lat_max(1)-lat_min(1)+1);
            
            
            all_tidedir = {'flood', 'ebb'};
            all_var_analysis = {'speed', 'u', 'v', 'zeta', 'time'}
% % % % %             for tidedir_ind2=1:length(all_tidedir)
% % % % %                 for var_analysis_ind2=1:length(all_var_analysis)
% % % % % %                 Model.max_flood_speed(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_flood_u(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_flood_v(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_flood_zeta(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_flood_time(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_ebb_speed(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_ebb_u(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_ebb_v(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_ebb_zeta(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % % %                 Model.max_ebb_time(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;
% % % % %                     tidedir=all_tidedir{tidedir_ind2};
% % % % %                     var_analysis=all_var_analysis{var_analysis_ind2};
% % % % %                     eval(['Model.max_', tidedir,'_', var_analysis,'(1:lon_max(1)-lon_min(1)+1,1:lat_max(1)-lat_min(1)+1)=NaN;']);
% % % % %                 end
% % % % %             end
% % % % %             
            tic;
% % % % %             for lon_loop=lon_min(1) : lon_max(1)
% % % % %                 for lat_loop=lat_min(1) : lat_max(1)
% % % % %                     loni=lon_loop-lon_min(1)+1;
% % % % %                     lati=lat_loop-lat_min(1)+1;
% % % % %                     filename=[outputdir,'\ocean_his_', 'zeta', '_AKP4_', num2str(1, '%04i'), '.nc'];
% % % % %                     varflag=ncread(filename,'zeta', [lon_loop, lat_loop, 1], [1, 1, 1]);
% % % % %                     if(~isnan(varflag))
% % % % %                         for dayi=1:fullday
% % % % %                             all_var2 ={'zeta', 'ubar_eastward', 'vbar_northward'};
% % % % %                             for varind2= 1:length(all_var2)
% % % % %                                 varname=all_var2{varind2};
% % % % %                                 filename=[outputdir,'\ocean_his_', varname, '_AKP4_', num2str(dayi, '%04i'), '.nc'];
% % % % %                                 eval(['temp_', varname, '((dayi-1)*24+1:dayi*24)=ncread(filename,varname, [lon_loop, lat_loop, 1], [1, 1, 24]);']);
% % % % %                             end
% % % % % %                             dayi
% % % % %                         end
% % % % %                         diff_temp_zeta(2:fullday*24)=diff(temp_zeta);
% % % % %                         temp_speed=sqrt(temp_ubar_eastward.^2+temp_vbar_northward.^2);
% % % % %                         temp_ebb_speed=temp_speed;
% % % % %                         temp_ebb_speed((diff_temp_zeta>0))=NaN;
% % % % %                         temp_flood_speed=temp_speed;
% % % % %                         temp_flood_speed((diff_temp_zeta<0))=NaN;
% % % % %                         
% % % % %                         Model.max_ebb_speed(loni,lati)=max(temp_ebb_speed);
% % % % %                         Model.max_flood_speed(loni,lati)=max(temp_flood_speed);
% % % % %                         max_ebb_ind=find(temp_ebb_speed==Model.max_ebb_speed(loni,lati),1);
% % % % %                         max_flood_ind=find(temp_flood_speed==Model.max_flood_speed(loni,lati),1);
% % % % %                         Model.max_ebb_time(loni,lati)=max_ebb_ind / 24.0;
% % % % %                         Model.max_ebb_ubar_eastward(loni,lati)=temp_ubar_eastward(max_ebb_ind);
% % % % %                         Model.max_ebb_vbar_northward(loni,lati)=temp_vbar_northward(max_ebb_ind);
% % % % %                         Model.max_ebb_zeta(loni,lati)=temp_zeta(max_ebb_ind);
% % % % %                         Model.max_flood_time(loni,lati)=max_flood_ind / 24.0;
% % % % %                         Model.max_flood_ubar_eastward(loni,lati)=temp_ubar_eastward(max_flood_ind);
% % % % %                         Model.max_flood_vbar_northward(loni,lati)=temp_vbar_northward(max_flood_ind);
% % % % %                         Model.max_flood_zeta(loni,lati)=temp_zeta(max_flood_ind);
% % % % % 
% % % % %                     else
% % % % %                         for tidedir_ind2=1:length(all_tidedir)
% % % % %                             for var_analysis_ind2=1:length(all_var_analysis)
% % % % %                                 tidedir=all_tidedir{tidedir_ind2};
% % % % %                                 var_analysis=all_var_analysis{var_analysis_ind2};
% % % % %                                 eval(['Model.max_', tidedir,'_', var_analysis,'(loni,lati)=NaN;']);
% % % % %                             end
% % % % %                         end
% % % % %                     end
% % % % %                     
% % % % %                 end
% % % % %                 lon_loop/lon_max(1)
% % % % %             end
            toc;
% % % % %             save([outputdir,'\tide_vel_analysis', '_', regionname, '_', num2str(tempyear, '%04i'), '.mat'])
            

            load([outputdir,'\tide_vel_analysis', '_', regionname, '_', num2str(tempyear, '%04i'), '.mat'])
            Model.max_flood_u=Model.max_flood_ubar_eastward;
            Model.max_flood_u(Model.max_flood_u==0)=NaN;
            Model.max_flood_v=Model.max_flood_vbar_northward;
            Model.max_flood_v(Model.max_flood_v==0)=NaN;
            Model.max_ebb_u=Model.max_ebb_ubar_eastward;
            Model.max_ebb_u(Model.max_ebb_u==0)=NaN;
            Model.max_ebb_v=Model.max_ebb_vbar_northward;
            Model.max_ebb_v(Model.max_ebb_v==0)=NaN;
            all_var_analysis = {'speed', 'ubar_eastward', 'vbar_northward', 'zeta', 'time'};
            Model.max_flood_ubar_eastward(Model.max_flood_ubar_eastward==0)=NaN;
            Model.max_flood_vbar_northward(Model.max_flood_vbar_northward==0)=NaN;
            Model.max_ebb_ubar_eastward(Model.max_ebb_ubar_eastward==0)=NaN;
            Model.max_ebb_vbar_northward(Model.max_ebb_vbar_northward==0)=NaN;


                        
            ncoutfilename=[outputdir,'\tide_vel_analysis', '_', regionname, '_', num2str(tempyear, '%04i'), '.nc'];
            
            % % %         make ncfile
            ncid = netcdf.create(ncoutfilename,'NETCDF4');

            cut_len_lon=size(cut_lon_rho,1);
            cut_len_lat=size(cut_lat_rho,2);
            xi_rho_dimid = netcdf.defDim(ncid, 'lon_rho', cut_len_lon);
            eta_rho_dimid = netcdf.defDim(ncid,'lat_rho', cut_len_lat);

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

            lon_rhovarid=netcdf.defVar(ncid, 'lon_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
            netcdf.putAtt(ncid,lon_rhovarid,'long_name','longitude');
            netcdf.putAtt(ncid,lon_rhovarid,'units','degree_east');

            lat_rhovarid=netcdf.defVar(ncid, 'lat_rho', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
            netcdf.putAtt(ncid,lat_rhovarid,'long_name','latitude');
            netcdf.putAtt(ncid,lat_rhovarid,'units','degree_north');
            
            for tidedir_ind2=1:length(all_tidedir)
                for var_analysis_ind2=1:length(all_var_analysis)
                    tidedir=all_tidedir{tidedir_ind2};
                    var_analysis=all_var_analysis{var_analysis_ind2};
                    varvarname=['max_', tidedir,'_', var_analysis];
                    eval([varvarname, 'varid=netcdf.defVar(ncid,', '''', varvarname, '''',',', '''', 'NC_DOUBLE', '''', ',' '[xi_rho_dimid eta_rho_dimid]);']);
                end
            end
% 
%             tconvarid=netcdf.defVar(ncid, 'tcon', 'NC_DOUBLE', [xi_rho_dimid eta_rho_dimid]);
%             netcdf.putAtt(ncid,tconvarid,'long_name','max_ ebb');

            netcdf.endDef(ncid);
            netcdf.putVar(ncid, lon_rhovarid, [0 0], [cut_len_lon cut_len_lat], cut_lon_rho);
            netcdf.putVar(ncid, lat_rhovarid, [0 0], [cut_len_lon cut_len_lat], cut_lat_rho);
            for tidedir_ind2=1:length(all_tidedir)
                for var_analysis_ind2=1:length(all_var_analysis)
                    tidedir=all_tidedir{tidedir_ind2};
                    var_analysis=all_var_analysis{var_analysis_ind2};
                    varvarname=['max_', tidedir,'_', var_analysis];
                    eval(['netcdf.putVar(ncid, ', varvarname,'varid, [0 0], [cut_len_lon cut_len_lat], ', 'Model.',varvarname, ');']);
                end
            end
            netcdf.close(ncid);
%             end
        end
    end
end

% quiver(Model.max_flood_ubar_eastward(1:5:end,1:5:end)'*3, Model.max_flood_vbar_northward(1:5:end,1:5:end)'*3,  'AutoScale','off')