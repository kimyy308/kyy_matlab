clc; close all; clear all;

% /data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_stddepth/test2117/1985

% RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131', 'ens2201'};
RCM_info.name={'test2127', 'test2128', 'test2129', 'test2130', 'test2131'};

dirs.dataroot= ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_stddepth/'];
dirs.saveroot= ['/data2/RCM/CMIP6/Model/ROMS/nwp_1_20/cut_ES_stddepth_transport/'];
% RCM_info.region ={'pollock_egg3'};
RCM_info.region ={'ES'};
RCM_info.years = [2081:2100];
RCM_info.years_his = [1995:2014];
% seasons_group={'JF-'};
seasons_group={'all'};

RCM_grid.dl = 1/20;
RCM_grid.gridname = {'lon_rho', 'lat_rho', 'lon_u', 'lat_u', 'lon_v', 'lat_v', 'lon_psi', 'lat_psi', 'pm', 'pn', 'f'};
% RCM_grid.gridname = {'lon_rho', 'lat_rho', 'depth' 'pm', 'pn', 'f'};


for seasons_groupi=1:length(seasons_group)
    for testnameind2=1:length(RCM_info.name)
        for regionind2=1:length(RCM_info.region)
            close all;
            clearvars '*' -except RCM_info RCM_grid testnameind2 regionind2 years_groupi years_group ...
                seasons_group seasons_groupi season dirs cut_area cut_xdist cut_ydist all_data
            tmp.fs=filesep;  
            
            %%     set dropbox path
            tmp.dropboxpath = '/home/kimyy/Dropbox';
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'function']));
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab_git', tmp.fs, 'function']));

            [tmp.dropboxpath, tmp.error_status] = Func_0008_set_dropbox_path(computer);
            addpath(genpath([tmp.dropboxpath, tmp.fs, 'source', tmp.fs, 'matlab', tmp.fs, 'Model' ...
                tmp.fs, 'ROMS', tmp.fs, 'Analysis', tmp.fs, 'Figure', tmp.fs, 'nwp_1_20', tmp.fs ...
                'run', tmp.fs, 'SSH', tmp.fs, '2phase_2nd', tmp.fs, 'subroutine']));
            
            %% get season and month
            RCM_info.season=seasons_group{seasons_groupi};
            [RCM_info.months,tmp.error_status] = Func_0019_get_month_from_season(RCM_info.season);
            
            %% get testname and region polygon
            tmp.testname_ssp=RCM_info.name{testnameind2};   % % need to change
            
            [tmp.testname_his, tmp.error_status] = Func_0023_RCM_CMIP6_testname_his(tmp.testname_ssp);

            [tmp.abb, tmp.error_status] = Func_0018_abbname(tmp.testname_ssp);
            tmp.regionname=RCM_info.region{regionind2};
            [RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(tmp.regionname);
            
            dirs.filedir_ssp = [dirs.dataroot, tmp.testname_ssp, tmp.fs];
            dirs.filedir_his = [dirs.dataroot, tmp.testname_his, tmp.fs];            
            dirs.matdir = [dirs.saveroot, tmp.testname_ssp, tmp.fs];
            dirs.griddir = [dirs.dataroot, '../cut_ES/', tmp.testname_ssp, tmp.fs];
            
            dirs.filedir_shflux_ssp = [dirs.dataroot, '../cut_ES/', tmp.testname_ssp, tmp.fs];
            dirs.filedir_shflux_his = [dirs.dataroot, '../cut_ES/', tmp.testname_his, tmp.fs];
            
            if (exist(strcat(dirs.matdir) , 'dir') ~= 7)
                mkdir(strcat(dirs.matdir));
            end 
    
            tmp.matname = [dirs.matdir, tmp.testname_ssp, '_', tmp.testname_his, '_', tmp.regionname,...
                '_lat_transport_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
            '_', RCM_info.season, '.mat'];
            
            tmp.gridname = [dirs.griddir, '/2081/pck_ES_test2127_monthly_2081_01.nc'];
        
            if (exist(tmp.matname , 'file') ~= 2 || fig_flag==2)
                for yearij=1:length(RCM_info.years)  %% yearly loop
                    tmp.tempyear=RCM_info.years(yearij);
                    tmp.yearstr=num2str(tmp.tempyear, '%04i');
                    tmp.tempyear_his=RCM_info.years_his(yearij);
                    tmp.yearstr_his=num2str(tmp.tempyear_his, '%04i');
                    
                    disp(num2str(tmp.testname_ssp));
                    disp(num2str(yearij));
            
                    for monthij=1:length(RCM_info.months)  %% monthly loop
                        tmp.tempmonth=RCM_info.months(monthij);
                        tmp.monthstr=num2str(tmp.tempmonth, '%02i');

                        tmp.filename_ssp = [dirs.filedir_ssp, tmp.fs, tmp.yearstr, tmp.fs, ...
                            tmp.testname_ssp, '_monthly_std_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.filename_his = [dirs.filedir_his, tmp.fs, tmp.yearstr_his, tmp.fs, ...
                            tmp.testname_his, '_monthly_std_', tmp.yearstr_his, '_', tmp.monthstr, '.nc'];
                        
                        tmp.filename_shflux_ssp = [dirs.filedir_shflux_ssp, tmp.fs, tmp.yearstr, tmp.fs, ...
                            'pck_ES_', tmp.testname_ssp, '_monthly_', tmp.yearstr, '_', tmp.monthstr, '.nc'];
                        tmp.filename_shflux_his = [dirs.filedir_shflux_his, tmp.fs, tmp.yearstr_his, tmp.fs, ...
                            'pck_ES_', tmp.testname_his, '_monthly_', tmp.yearstr_his, '_', tmp.monthstr, '.nc'];

                        if (isfield(RCM_grid, 'lon_rho') ~= 1)

                            for gridi=1:length(RCM_grid.gridname)
                                RCM_grid.(RCM_grid.gridname{gridi})=ncread(tmp.gridname, RCM_grid.gridname{gridi});
                            end
                            RCM_grid.depth = ncread(tmp.filename_ssp, 'depth');

                            [RCM_grid.lon_min, RCM_grid.lon_max, RCM_grid.lat_min, RCM_grid.lat_max] = ...
                                findind_Y(1/20, RCM_grid.domain(1:4), RCM_grid.lon_rho, RCM_grid.lat_rho);
                            RCM_grid.cut_lon_rho = ...
                                RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                            RCM_grid.cut_lat_rho = ...
                                RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
%                             cut_lon_psi = ...
%                                 RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1)-1, RCM_grid.lat_min(1):RCM_grid.lat_max(1)-1);
%                             cut_lat_psi = ...
%                                 RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1)-1, RCM_grid.lat_min(1):RCM_grid.lat_max(1)-1);
                            RCM_grid.xdist=1./RCM_grid.pm;
                            RCM_grid.ydist=1./RCM_grid.pn;
                            cut_xdist= ...
                                RCM_grid.xdist(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                            cut_ydist= ...
                                RCM_grid.ydist(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
                            cut_area= cut_xdist.*cut_ydist;
                            
                        end
                        
                        tmp.temp_ssp = ncread(tmp.filename_ssp, 'temp',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
                        tmp.u_ssp = ncread(tmp.filename_ssp, 'u',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
                        tmp.v_ssp = ncread(tmp.filename_ssp, 'v',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
%                         tmp.w_ssp = ncread(tmp.filename_ssp, 'w',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 2 1], ...
%                             [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);  % 2nd layer
                        tmp.shflux_ssp = ncread(tmp.filename_shflux_ssp, 'shflux',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]); 
                       
                        
                        tmp.temp_his = ncread(tmp.filename_his, 'temp',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
                        tmp.u_his = ncread(tmp.filename_his, 'u',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
                        tmp.v_his = ncread(tmp.filename_his, 'v',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 inf 1]);
%                         tmp.w_his = ncread(tmp.filename_his, 'w',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 2 1], ...
%                             [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1 1]);  % 2nd layer
                        tmp.shflux_his = ncread(tmp.filename_shflux_his, 'shflux',[RCM_grid.lon_min(1) RCM_grid.lat_min(1) 1], ...
                            [RCM_grid.lon_max(1)-RCM_grid.lon_min(1)+1 RCM_grid.lat_max(1)-RCM_grid.lat_min(1)+1 1]);  

                        cut_mask=NaN(size(tmp.shflux_his));
                        cut_mask(isfinite(tmp.shflux_his))=1;
                        cut_area=cut_area.*cut_mask;
%                         tmp.c=4.1*10^6; % c
                        tmp.c=4.184; %heat capacity, J/g / celcius
                        tmp.c2=tmp.c * 1e3 * 1e3; % heat capacity, J/m^3 / celcius
%                         tmp.monthly_sec= 60*60*24*30;
%                         tmp.unit_watt=1e+15; % param unit for watt
                        RCM_grid.thickness=diff(RCM_grid.depth);
                        RCM_grid.thickness(end+1)=RCM_grid.thickness(end);
                        
                        RCM_grid.thickness_3d=repmat(RCM_grid.thickness,[1, size(cut_area,1), size(cut_area,2)]);
                        RCM_grid.thickness_3d=permute(RCM_grid.thickness_3d, [2 3 1]);
%                         RCM_grid.thickness_3d(:,:,3)
                        
                            %%ssp
                            tmp.heat_trans_u_ssp_3d= tmp.temp_ssp.*tmp.u_ssp.*RCM_grid.thickness_3d .*(tmp.c2).*cut_ydist; % Watt
                            tmp.heat_trans_v_ssp_3d= tmp.temp_ssp.*tmp.v_ssp.*RCM_grid.thickness_3d .*(tmp.c2).*cut_xdist; % Watt
                            tmp.heat_trans_u_ssp = sum(tmp.heat_trans_u_ssp_3d,3,'omitnan');
                            tmp.heat_trans_v_ssp = sum(tmp.heat_trans_v_ssp_3d,3,'omitnan');
    %                         tmp.heat_trans_w_ssp= tmp.temp_ssp.*tmp.w_ssp.*(tmp.c2).*cut_xdist.*cut_ydist; % Watt
                            tmp.heat_trans_shflux_ssp= tmp.shflux_ssp.*cut_area; % Watt

                            tmp.heat_trans_u_m2_ssp= tmp.temp_ssp.*tmp.u_ssp.*tmp.c2; % Watt / m^2
                            tmp.heat_trans_v_m2_ssp= tmp.temp_ssp.*tmp.v_ssp.*tmp.c2; % Watt / m^2
    %                         tmp.heat_trans_w_m2_ssp= tmp.temp_ssp.*tmp.w_ssp.*tmp.c2; % Watt / m^2

    %                         pcolor(tmp.heat_trans_w_m2'); shading flat; colorbar; caxis([-2000 2000]);
                            
                            %% historical
                            tmp.heat_trans_u_his_3d= tmp.temp_his.*tmp.u_his.*RCM_grid.thickness_3d .*(tmp.c2).*cut_ydist; % Watt
                            tmp.heat_trans_v_his_3d= tmp.temp_his.*tmp.v_his.*RCM_grid.thickness_3d .*(tmp.c2).*cut_xdist; % Watt
                            tmp.heat_trans_u_his = sum(tmp.heat_trans_u_his_3d,3,'omitnan');
                            tmp.heat_trans_v_his = sum(tmp.heat_trans_v_his_3d,3,'omitnan');
    %                         tmp.heat_trans_w_his= tmp.temp_his.*tmp.w_his.*(tmp.c2).*cut_xdist.*cut_ydist; % Watt
                            tmp.heat_trans_shflux_his= tmp.shflux_his.*cut_area; % Watt

                            tmp.heat_trans_u_m2_his= tmp.temp_his.*tmp.u_his.*tmp.c2; % Watt / m^2
                            tmp.heat_trans_v_m2_his= tmp.temp_his.*tmp.v_his.*tmp.c2; % Watt / m^2
    %                         tmp.heat_trans_w_m2_his= tmp.temp_his.*tmp.w_his.*tmp.c2; % Watt / m^2

    %                         pcolor(tmp.heat_trans_w_m2'); shading flat; colorbar; caxis([-2000 2000]);
    
                    switch tmp.regionname
                        case('pollock_egg3')
                            

                            tmp.sum_seg_ssp(1)=sum(tmp.heat_trans_u_ssp(1,:), 'omitnan'); % west
                            tmp.sum_seg_ssp(2)=-sum(tmp.heat_trans_u_ssp(end,:), 'omitnan'); % east
                            tmp.sum_seg_ssp(3)=sum(tmp.heat_trans_v_ssp(:,1), 'omitnan'); % south
                            tmp.sum_seg_ssp(4)=-sum(tmp.heat_trans_v_ssp(:,end), 'omitnan'); %north
                            tmp.sum_seg_ssp(5)=0; %vertical
                            tmp.sum_seg_ssp(6)=sum(sum(tmp.heat_trans_shflux_ssp, 'omitnan'),'omitnan'); %surface


                            %% historical

                            tmp.sum_seg_his(1)=sum(tmp.heat_trans_u_his(1,:), 'omitnan'); % west
                            tmp.sum_seg_his(2)=-sum(tmp.heat_trans_u_his(end,:), 'omitnan'); % east
                            tmp.sum_seg_his(3)=sum(tmp.heat_trans_v_his(:,1), 'omitnan'); % south
                            tmp.sum_seg_his(4)=-sum(tmp.heat_trans_v_his(:,end), 'omitnan'); %north
                            tmp.sum_seg_his(5)=0; %vertical
                            tmp.sum_seg_his(6)=sum(sum(tmp.heat_trans_shflux_his, 'omitnan'),'omitnan'); %surface
                        case ('ES')
                            tmp.sum_seg_ssp(1)=sum(tmp.heat_trans_v_ssp(35:122, 51), 'omitnan'); % Korea strait
                            tmp.sum_seg_ssp(2)=-sum(tmp.heat_trans_u_ssp(266, 181:214), 'omitnan'); % Tsugaru strait
                            tmp.sum_seg_ssp(3)=-sum(tmp.heat_trans_u_ssp(304, 310:342), 'omitnan'); % Soya strait
                            tmp.sum_seg_ssp(4)=0; 
                            tmp.sum_seg_ssp(5)=0;
                            [RCM_grid.ES_KHOApolygon, RCM_grid.ES_KHOAdomain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('ES_KHOA');
                            RCM_grid.mask_ES_KHOA = double(inpolygon(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, ...
                                RCM_grid.ES_KHOApolygon(:,1), RCM_grid.ES_KHOApolygon(:,2)));
                            RCM_grid.mask_ES_KHOA(RCM_grid.mask_ES_KHOA==0)=NaN;
                            tmp.heat_trans_shflux_ssp=tmp.heat_trans_shflux_ssp.*RCM_grid.mask_ES_KHOA;
                            tmp.sum_seg_ssp(6)=sum(sum(tmp.heat_trans_shflux_ssp, 'omitnan'),'omitnan'); %surface
                            %% historical
                            tmp.sum_seg_his(1)=sum(tmp.heat_trans_v_his(35:122, 51), 'omitnan'); % Korea strait
                            tmp.sum_seg_his(2)=-sum(tmp.heat_trans_u_his(266, 181:214), 'omitnan'); % Tsugaru strait
                            tmp.sum_seg_his(3)=-sum(tmp.heat_trans_u_his(304, 310:342), 'omitnan'); % Soya strait
                            tmp.sum_seg_his(4)=0; 
                            tmp.sum_seg_his(5)=0;
                            [RCM_grid.ES_KHOApolygon, RCM_grid.ES_KHOAdomain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname('ES_KHOA');
                            RCM_grid.mask_ES_KHOA = double(inpolygon(RCM_grid.cut_lon_rho, RCM_grid.cut_lat_rho, ...
                                RCM_grid.ES_KHOApolygon(:,1), RCM_grid.ES_KHOApolygon(:,2)));
                            RCM_grid.mask_ES_KHOA(RCM_grid.mask_ES_KHOA==0)=NaN;
                            tmp.heat_trans_shflux_his=tmp.heat_trans_shflux_his.*RCM_grid.mask_ES_KHOA;
                            tmp.sum_seg_his(6)=sum(sum(tmp.heat_trans_shflux_his, 'omitnan'),'omitnan'); %surface
                    end
                        
                        
                        
                        tmp.sum_seg_diff=tmp.sum_seg_ssp - tmp.sum_seg_his;
%                         sum(tmp.sum_seg_ssp) sum(tmp.sum_seg_his)
                        
%                         sum(sum(cut_area.*20, 'omitnan'), 'omitnan')

                        RCM_data.sum_seg_ssp(:,yearij,monthij)=tmp.sum_seg_ssp;
                        RCM_data.sum_seg_his(:,yearij,monthij)=tmp.sum_seg_his;
                        RCM_data.sum_seg_diff(:,yearij,monthij)=tmp.sum_seg_diff;
                    end
                    RCM_data.yearly_sum_seg_ssp(:,yearij)=mean(RCM_data.sum_seg_ssp(:,yearij,:),3);
                    RCM_data.yearly_sum_seg_his(:,yearij)=mean(RCM_data.sum_seg_his(:,yearij,:),3);
                    RCM_data.yearly_sum_seg_diff(:,yearij)=mean(RCM_data.sum_seg_diff(:,yearij,:),3);
                end
            end
            
            tmp.m_sum_seg_his=squeeze(mean(mean(RCM_data.sum_seg_his,2),3))
            tmp.m_sum_seg_ssp=squeeze(mean(mean(RCM_data.sum_seg_ssp,2),3))
            tmp.m_sum_seg_diff=squeeze(mean(mean(RCM_data.sum_seg_diff,2),3))
            sum(tmp.m_sum_seg_his(1:4)) + tmp.m_sum_seg_his(6)
            sum(tmp.m_sum_seg_ssp(1:4)) + tmp.m_sum_seg_ssp(6)
            sum(tmp.m_sum_seg_diff(1:4)) + tmp.m_sum_seg_diff(6)
            sum(tmp.m_sum_seg_ssp)
            
        end
        all_data.all_sum_seg_his(:, testnameind2)=tmp.m_sum_seg_his
        all_data.all_sum_seg_ssp(:, testnameind2)=tmp.m_sum_seg_ssp
        all_data.all_sum_seg_diff(:, testnameind2)=tmp.m_sum_seg_diff
        
        all_data.raw_all_sum_seg_ssp(:,:,:,testnameind2)=RCM_data.sum_seg_ssp;
        all_data.raw_all_sum_seg_his(:,:,:,testnameind2)=RCM_data.sum_seg_his;
        all_data.raw_all_sum_seg_diff(:,:,:,testnameind2)=RCM_data.sum_seg_diff;
        
        all_data.raw_yearly_all_sum_seg_ssp(:,:,testnameind2)=RCM_data.yearly_sum_seg_ssp;
        all_data.raw_yearly_all_sum_seg_his(:,:,testnameind2)=RCM_data.yearly_sum_seg_his;
        all_data.raw_yearly_all_sum_seg_diff(:,:,testnameind2)=RCM_data.yearly_sum_seg_diff;
        
    end
end

sum(all_data.all_sum_seg_diff(1:4,:),1)
all_data.all_sum_seg_diff(6,:)

tmp.data=mean(all_data.raw_yearly_all_sum_seg_ssp(6,:,1,:),4)
std(tmp.data)
tmp.data2=mean(all_data.raw_yearly_all_sum_seg_ssp(1:4,:,1,:),4)
tmp.data=sum(tmp.data2,1);
std(tmp.data)

tmp.data2=mean(all_data.raw_yearly_all_sum_seg_ssp(1:6,:,1,:),4)
tmp.data=sum(tmp.data2,1);
std(tmp.data)

plot(tmp.data)

mean(sum(all_data.all_sum_seg_his(1,:),1))


mean(sum(all_data.all_sum_seg_his(1:4,:),1))
mean(all_data.all_sum_seg_his(6,:))

mean(sum(all_data.all_sum_seg_ssp(1:4,:),1))
mean(all_data.all_sum_seg_ssp(6,:))

mean(sum(all_data.all_sum_seg_diff(1:4,:),1))
mean(all_data.all_sum_seg_diff(6,:))


tmp.sum_seg_ssp1= permute(RCM_data.sum_seg_ssp, [1 3 2]);
tmp.sum_seg_ssp2= reshape(tmp.sum_seg_ssp1, [6, 240]);

tmp.sum_seg_ssp3 = sum(tmp.sum_seg_ssp2(1:4,:),1);
plot(tmp.sum_seg_ssp3)
plot(tmp.sum_seg_ssp2(6,:))

save([RCM_info.region{1}, '.mat'], 'all_data')

function [testname_his, error_status] = Func_0023_RCM_CMIP6_testname_his(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=RCM_CMIP6_testname_his(testname);
%
% get the historical testname corresponding to SSP test name
%
%  input:
%  testname             ROMS SSP test name (string)
%
%  output:
%  testname_his         ROMS historical test name (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    18-Apr-2022 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    switch(testname)
        case {'test2117', 'test2127'}
            testname_his='test2117';
        case {'test2118', 'test2128'}
            testname_his='test2118';
        case {'test2119', 'test2129'}
            testname_his='test2119';
        case {'test2120', 'test2130'}
            testname_his='test2120';
        case {'test2121', 'test2131'}
            testname_his='test2121';
        case {'ens2201'}
            testname_his='ens2202';
        case {'ens2203'}
            testname_his='ens2204';
        case {'prob_ens2203'}
            testname_his='prob_ens2204';
           
    end
    error_status=1;
end