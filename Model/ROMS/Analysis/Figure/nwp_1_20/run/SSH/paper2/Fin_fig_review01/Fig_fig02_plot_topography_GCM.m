close all; clear all;  clc;
warning off;

all_testname2  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};

all_region2 ={'NWP'}
% all_region2 ={'AKP4'}

% all_region2 ={'YS'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
all_var2 = {'H'};
scenname='rcp45';
% all_region2 ={'NWP'}
for regionind2=1:length(all_region2)
        close all;
        clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2 scenname

        inputyear1 = [2006:2006]; % % put year which you want to plot [year year ...]
        inputmonth = [12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            dropboxpath='C:\Users\User\Dropbox';
            addpath(genpath([dropboxpath '\source\matlab\Common\export_fig-master']));            
            addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
            addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
            addpath(genpath([dropboxpath '\source\matlab\Common\cptcmap']));
            addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end
        variable =all_var2{1};
        end
        % for snu_desktopd
       
        run('nwp_polygon_point.m');
        regionname=all_region2{regionind2};
        switch(regionname)
            case('NWP')
                refpolygon=nwppolygon;
            case('AKP4') %% Around Korea Peninsula
                refpolygon=akp4polygon;
        end
        lonlat(1)=min(refpolygon(:,1));
        lonlat(2)=max(refpolygon(:,1));
        lonlat(3)=min(refpolygon(:,2));
        lonlat(4)=max(refpolygon(:,2));
 
        load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)

        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\'); % % where figure files will be saved
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        correction_right_fig=[-0.1000,0,0,0];
        correction_upper_fig=[0,-0.1,0,0];
        hold on;
        tifname=strcat(figrawdir, 'fig02','.tif'); %% ~_year_month.jpg
        epsname=strcat(figrawdir, 'fig02','.eps'); %% ~_year_month.jpg

%         f1=figure(1);

%      ---------- RCM-IPSL-L
        testnameind2=1;
        testname=all_testname2{testnameind2};    
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m'];
        run(param_script);
        variable=all_var2{1};
        for yearij=1:length(inputyear1)  % ---------------------- get depth data
            tempyear=inputyear1(yearij);
            yearstr=num2str(tempyear, '%04i');
            tempmonth=inputmonth(1);
            monthstr=num2str(tempmonth, '%02i');

            filedir = strcat(cmip5dir, 'thetao', '\', scenname, '\Omon\', testname, '\'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '\', 'thetao', '*' ]); 
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

            if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(filename, 'lon');
                lat_rho=ncread(filename, 'lat');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho,1);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            end
            data_info = ncinfo(filename, 'thetao'); 

            data = ncread(filename,'thetao',[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
            numx=size(data,1);
            numy=size(data,2);
            lev_bnds= ncread(filename, 'lev_bnds');
            for lonind=1:numx
                for latind=1:numy
                    vert_temp=squeeze(data(lonind,latind,:,1));
                    nanind=find(isnan(vert_temp));                            
                    botind=min(nanind)-1;
                    if (isfinite(vert_temp(end)))
                        botind=length(vert_temp);
                    end
                    if (botind>0)
%                                 bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                        bottomtemp(lonind,latind)=lev_bnds(2,botind);
                    else
                        bottomtemp(lonind,latind)=NaN;
                    end
                    clear nanind botind
                end
            end
            if (exist('mean_data' , 'var') ~= 1)
                mean_data=zeros(size(bottomtemp));
            end
            mean_data=mean_data + (bottomtemp / length(inputyear1));
        end % ---------------------- get depth data
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mask_model_inv = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model_inv(mask_model_inv==0)=NaN;
        mask_model_inv(isfinite(bottomtemp))=NaN;
        
        mean_data= bottomtemp .* mask_model; 
        GCM_IPSL_L_topo=mean_data;
        GCM_IPSL_L_lon=cut_lon_rho;
        GCM_IPSL_L_lat=cut_lat_rho;
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb1=subplot(2,2,1);  % Auto-fitted to the figure.
        pos_sb1=get(sb1, 'pos'); % Get the position.
        pos_sb1 = pos_sb1 + correction_upper_fig;
        delete(sb1); % Delete the subplot axes
        ax1_1=axes;
        set(ax1_1,'pos',pos_sb1);
        pc1_1=m_pcolor(cut_lon_rho',cut_lat_rho', mask_model_inv','parent',ax1_1);
        colormap(ax1_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'parent', ax1_1);  
        hold on

        ax1_2=axes;
        set(ax1_2,'pos',pos_sb1);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc1_2=m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data','parent',ax1_2);
        colormap(ax1_2,gmt_ocean_mod2);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'backcolor', 'none', 'parent', ax1_2);  
        txt_1_2=m_text(120, 20, '(A) IPSL-CM5A-LR', 'FontSize', m_grid_fontsize-5, 'color',[1, 1, 1-0.01]); 

        clear lon_rho mean_data bottomtemp
%      ---------- RCM-IPSL-L

%      ---------- RCM-IPSL-M
        testnameind2=2;
        testname=all_testname2{testnameind2};
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m'];
        run(param_script);
        variable=all_var2{1};
        for yearij=1:length(inputyear1)  % ---------------------- get depth data
            tempyear=inputyear1(yearij);
            yearstr=num2str(tempyear, '%04i');
            tempmonth=inputmonth(1);
            monthstr=num2str(tempmonth, '%02i');

            filedir = strcat(cmip5dir, 'thetao', '\', scenname, '\Omon\', testname, '\'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '\', 'thetao', '*' ]); 
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

            if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(filename, 'lon');
                lat_rho=ncread(filename, 'lat');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho,1);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            end
            data_info = ncinfo(filename, 'thetao'); 

            data = ncread(filename,'thetao',[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
            numx=size(data,1);
            numy=size(data,2);
            lev_bnds= ncread(filename, 'lev_bnds');
            for lonind=1:numx
                for latind=1:numy
                    vert_temp=squeeze(data(lonind,latind,:,1));
                    nanind=find(isnan(vert_temp));                            
                    botind=min(nanind)-1;
                    if (isfinite(vert_temp(end)))
                        botind=length(vert_temp);
                    end
                    if (botind>0)
%                                 bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                        bottomtemp(lonind,latind)=lev_bnds(2,botind);
                    else
                        bottomtemp(lonind,latind)=NaN;
                    end
                    clear nanind botind
                end
            end
            if (exist('mean_data' , 'var') ~= 1)
                mean_data=zeros(size(bottomtemp));
            end
            mean_data=mean_data + (bottomtemp / length(inputyear1));
        end % ---------------------- get depth data
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mask_model_inv = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model_inv(mask_model_inv==0)=NaN;
        mask_model_inv(isfinite(bottomtemp))=NaN;
        
        mean_data= bottomtemp .* mask_model;        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb2=subplot(2,2,2);  % Auto-fitted to the figure.
        pos_sb2=get(sb2, 'pos'); % Get the position.
        pos_sb2 = pos_sb2 + correction_right_fig + correction_upper_fig;
        delete(sb2); % Delete the subplot axes
        ax2_1=axes;
        set(ax2_1,'pos',pos_sb2);
        pc2_1=m_pcolor(cut_lon_rho',cut_lat_rho', mask_model_inv','parent',ax2_1);
        colormap(ax2_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'yticklabels', [], 'parent', ax2_1);  
        hold on

        ax2_2=axes;
        set(ax2_2,'pos',pos_sb2);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc2_2=m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data','parent',ax2_2);
        colormap(ax2_2,gmt_ocean_mod2);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'xticklabels', [], 'xtick',(120:20:160), 'yticklabels', [],'backcolor', 'none', 'parent', ax2_2);  
        txt_2_2=m_text(119, 20, '(B) IPSL-CM5A-MR', 'FontSize', m_grid_fontsize-5, 'color',[1, 1, 1-0.01]); 

        clear lon_rho mean_data bottomtemp
%      ---------- RCM-IPSL-M

%      ---------- RCM-Nor
        testnameind2=3;
        testname=all_testname2{testnameind2};  
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m'];
        run(param_script);
        hold on;
        variable=all_var2{1};
        for yearij=1:length(inputyear1)  % ---------------------- get depth data
            tempyear=inputyear1(yearij);
            yearstr=num2str(tempyear, '%04i');
            tempmonth=inputmonth(1);
            monthstr=num2str(tempmonth, '%02i');

            filedir = strcat(cmip5dir, 'thetao', '\', scenname, '\Omon\', testname, '\'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '\', 'thetao', '*' ]); 
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

            if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(filename, 'lon');
                lat_rho=ncread(filename, 'lat');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho,1);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            end
            data_info = ncinfo(filename, 'thetao'); 

            data = ncread(filename,'thetao',[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
            numx=size(data,1);
            numy=size(data,2);
            lev_bnds= ncread(filename, 'lev_bnds');
            for lonind=1:numx
                for latind=1:numy
                    vert_temp=squeeze(data(lonind,latind,:,1));
                    nanind=find(isnan(vert_temp));                            
                    botind=min(nanind)-1;
                    if (isfinite(vert_temp(end)))
                        botind=length(vert_temp);
                    end
                    if (botind>0)
%                                 bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                        bottomtemp(lonind,latind)=lev_bnds(2,botind);
                    else
                        bottomtemp(lonind,latind)=NaN;
                    end
                    clear nanind botind
                end
            end
            if (exist('mean_data' , 'var') ~= 1)
                mean_data=zeros(size(bottomtemp));
            end
            mean_data=mean_data + (bottomtemp / length(inputyear1));
        end % ---------------------- get depth data
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mask_model_inv = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model_inv(mask_model_inv==0)=NaN;
        mask_model_inv(isfinite(bottomtemp))=NaN;
        
        mean_data= bottomtemp .* mask_model;  
        
        GCM_Nor_topo=mean_data;
        GCM_Nor_lon=cut_lon_rho;
        GCM_Nor_lat=cut_lat_rho;
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb3=subplot(2,2,3);  % Auto-fitted to the figure.
        pos_sb3=get(sb3, 'pos'); % Get the position.
        delete(sb3); % Delete the subplot axes
        ax3_1=axes;
        set(ax3_1,'pos',pos_sb3);
        pc3_1=m_pcolor(cut_lon_rho',cut_lat_rho', mask_model_inv','parent',ax3_1);
        colormap(ax3_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, ...
            'parent', ax3_1, 'xtick',(120:20:160), 'xticklabels',[120 140 160]);  
        hold on

        ax3_2=axes;
        set(ax3_2,'pos',pos_sb3);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc3_2=m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data','parent',ax3_2);
        colormap(ax3_2,gmt_ocean_mod2);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type,  'backcolor', 'none', ...
            'parent', ax3_2, 'xtick',(120:20:160), 'xticklabels',[120 140 160]);  
        txt_3_2=m_text(126, 20, '(C) NorESM1-M', 'FontSize', m_grid_fontsize-5, 'color',[1, 1, 1-0.01]); 

        clear lon_rho mean_data bottomtemp
%      ---------- RCM-Nor

%      ---------- RCM-MPI
        testnameind2=4;
        testname=all_testname2{testnameind2};
        param_script =['C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m'];
        run(param_script);
        variable=all_var2{1};
        for yearij=1:length(inputyear1)  % ---------------------- get depth data
            tempyear=inputyear1(yearij);
            yearstr=num2str(tempyear, '%04i');
            tempmonth=inputmonth(1);
            monthstr=num2str(tempmonth, '%02i');

            filedir = strcat(cmip5dir, 'thetao', '\', scenname, '\Omon\', testname, '\'); % % where data files are
            flag_file_in = false;
            list = dir( [ filedir, '\', 'thetao', '*' ]); 
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

            if (exist('lon_rho' , 'var') ~= 1)
                lon_rho=ncread(filename, 'lon');
                lat_rho=ncread(filename, 'lat');
                [lon_min, lon_max, lat_min, lat_max] = findind_Y(1/20, lonlat(1:4), lon_rho, lat_rho,1);
                cut_lon_rho = lon_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
                cut_lat_rho = lat_rho(lon_min(1):lon_max(1), lat_min(1):lat_max(1));
            end
            data_info = ncinfo(filename, 'thetao'); 

            data = ncread(filename,'thetao',[lon_min(1) lat_min(1) 1 tind], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);  %% cut horizontal area [x,y,z] (wider than target area) 
            numx=size(data,1);
            numy=size(data,2);
            lev_bnds= ncread(filename, 'lev_bnds');
            for lonind=1:numx
                for latind=1:numy
                    vert_temp=squeeze(data(lonind,latind,:,1));
                    nanind=find(isnan(vert_temp));                            
                    botind=min(nanind)-1;
                    if (isfinite(vert_temp(end)))
                        botind=length(vert_temp);
                    end
                    if (botind>0)
%                                 bottomtemp(lonind,latind)=vert_temp(botind)-273.15;
                        bottomtemp(lonind,latind)=lev_bnds(2,botind);
                    else
                        bottomtemp(lonind,latind)=NaN;
                    end
                    clear nanind botind
                end
            end
            if (exist('mean_data' , 'var') ~= 1)
                mean_data=zeros(size(bottomtemp));
            end
            mean_data=mean_data + (bottomtemp / length(inputyear1));
        end % ---------------------- get depth data
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
        mask_model_inv = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model_inv(mask_model_inv==0)=NaN;
        mask_model_inv(isfinite(bottomtemp))=NaN;
        
        mean_data= bottomtemp .* mask_model;        
        
        GCM_MPI_topo=mean_data;
        GCM_MPI_lon=cut_lon_rho;
        GCM_MPI_lat=cut_lat_rho;
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb4=subplot(2,2,4);  % Auto-fitted to the figure.
        pos_sb4=get(sb4, 'pos'); % Get the position.
        pos_sb4 = pos_sb4 + correction_right_fig;

        delete(sb4); % Delete the subplot axes
        ax4_1=axes;
        set(ax4_1,'pos',pos_sb4);
        pc4_1=m_pcolor(cut_lon_rho',cut_lat_rho', mask_model_inv','parent',ax4_1);
        colormap(ax4_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'yticklabels', [], ...
            'parent', ax4_1, 'xtick',(120:20:160), 'xticklabels',[120  140 160]);   
        hold on

        ax4_2=axes;
        set(ax4_2,'pos',pos_sb4);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc4_2=m_pcolor(cut_lon_rho',cut_lat_rho', -mean_data','parent',ax4_2);
        colormap(ax4_2,gmt_ocean_mod2);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        m_grid('fontsize', m_grid_fontsize-5, 'tickdir', m_grid_tickdir_type, 'yticklabels', [], 'backcolor', 'none', ...
            'parent', ax4_2, 'xtick',(120:20:160), 'xticklabels',[120 140 160]);  
        txt_4_2=m_text(125, 20, '(D) MPI-ESM-LR', 'FontSize', m_grid_fontsize-5, 'color',[1, 1, 1-0.01]); 

        
        clear lon_rho mean_data bottomtemp
%      ---------- RCM-MPI
             
        
        h = colorbar;
%         colormap(gmt_ocean_mod2);
        set(h,'fontsize',colorbar_fontsize-5);
        title(h,colorbar_title,'fontsize',colorbar_title_fontsize);
        caxis(colorbar_lev);
        set(h, 'Position', [.9125, pos_sb4(2), 0.0231, pos_sb2(2)+pos_sb2(4)-0.1] + correction_right_fig)  % right, up, width, height
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

%         print -r300
%         saveas(gcf,tifname,'tif');

        print('-dtiff','-r500',tifname)

%         saveas(gcf,epsname,'epsc');
%         
%         set(gcf, 'Units', 'inches');
%         set(gcf,'Position', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 
%         export_fig(gcf,epsname,'-eps');

        hold off;
        close all;
        
        save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\Data02_GCM_topography.mat', ...
            'GCM_IPSL_L_topo', 'GCM_IPSL_L_lon', 'GCM_IPSL_L_lat', ...
            'GCM_Nor_topo', 'GCM_Nor_lon', 'GCM_Nor_lat', ...
            'GCM_MPI_topo', 'GCM_MPI_lon', 'GCM_MPI_lat')
            
end




% figH = figure;
% axLH = gca;
% axRH = axes('color','none');
% mslplot{1}=plot(inputyear,Model.amp_M2(3,:), 'b','parent',axLH);
% mslplot{2}=plot(inputyear,Obs.amp_M2(3,:), 'k','parent',axRH);
% ylabel(axLH,'Model M2 Tidal amplitude (cm)')
% ylabel(axRH,'Obs M2 Tidal amplitude (cm)')
% ax_pos = get(axLH,'position');
% set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
% set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
% %             set(axRH,'color','none','yaxislocation','right');
% 
% %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
% set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
% set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
% set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
% xlabel(axRH, 'Year');
% 
% title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
% % datetick('x','yyyy','keepticks')
% axis tight;
% % ylim(meanplotlev2)
% set(mslplot{1},'LineWidth',2);
% set(mslplot{2},'LineWidth',2);
% grid on
% 
% lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
% 
% %             lgd=legend('Model','TG-UST');
% set(lgd,'FontSize',15);
% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
% set(lgd,'Orientation','horizontal');
% 
% set(gcf,'PaperPosition', [0 0 36 12]) 
% saveas(gcf,pngname,'tif');
% grid off
% close all;