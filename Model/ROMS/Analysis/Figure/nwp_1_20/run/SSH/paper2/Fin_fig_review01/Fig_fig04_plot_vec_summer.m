close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
all_testname2  = {'test53', 'test54', 'test55', 'test56'};

GCM_testnames = {'GCM-IPSL-L', 'GCM-IPSL-M', 'GCM-Nor', 'GCM-MPI'};
RCM_testnames = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
all_region2 ={'AKP4'};

variable = 'vec';
% scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [1993:2005]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        season='summer';
        switch(season)
            case 'all'
                inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]
            case 'summer'
                inputmonth = [6 7 8]; % % put month which you want to plot [month month ...]
            case 'winter'
                inputmonth = [1 2 12]; % % put month which you want to plot [month month ...]
        end
        inputmonth2=inputmonth;
        system_name=computer;
        for folding=1:1
            if (strcmp(system_name,'PCWIN64'))
                % % for windows
                dropboxpath='C:\Users\User\Dropbox';
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
        end
       
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\'); % % where figure files will be saved
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig04_summer','.tif'); %% ~_year_month.jpg

%         f1=figure(1);

%      ---------- CMEMS
%         testnameind2=1;
        testname='CMEMS';    
%         load(['D:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
     
        param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_cmems_', regionname, '.m'];
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-9;
        load([cmemsdir,regionname, '_', testname, '_cmems_land','.mat']);
        cmems_land=cmems_land(1:end-2,:);
        load([cmemsdir, 'AKP4_cmems_vec_1993_2005.mat']);
        load([cmemsdir,regionname,'cmems_gos_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
        
        reshap_cmems_u=reshape(comb_cmems_u, [size(comb_cmems_u,1), size(comb_cmems_u,2), 12, size(comb_cmems_u,3)/12]);
        reshap_cmems_v=reshape(comb_cmems_v, [size(comb_cmems_u,1), size(comb_cmems_u,2), 12, size(comb_cmems_u,3)/12]);
        summer_cmems_u=reshap_cmems_u(:,:,inputmonth2,:);
        summer_cmems_v=reshap_cmems_v(:,:,inputmonth2,:);
        summer_cmems_u2=summer_cmems_u(:,:,:);
        summer_cmems_v2=summer_cmems_v(:,:,:);
        summer_mean_u=mean(summer_cmems_u2,3);
        summer_mean_v=mean(summer_cmems_v2,3);
%         avhrrfilename=['Z:\내 드라이브\Data\Observation\OISST\monthly_kimyy\avhrr_only_monthly_v2_2005.nc'];
%         avhrr_sst = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1) monthij], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
%         avhrr_land(isfinite(avhrr_sst))=NaN;
%         avhrr_land(isnan(avhrr_sst))=1;
%         save([avhrrdir,'avhrr_land.mat'], 'avhrr_land');
        
%         mean_data= mean_data .* mask_model;    

        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb1=subplot(10,4,[2 3 6 7]);  % Auto-fitted to the figure.
        pos_sb1=get(sb1, 'pos'); % Get the position.
        pos_sb1 = pos_sb1 + correction_upper_fig+ correction_large_fig + [-0.0800,0,0,0];
        delete(sb1); % Delete the subplot axes
        ax1_1=axes;
        set(ax1_1,'pos',pos_sb1);
        pc1_1=m_pcolor(cut_lon_rho,cut_lat_rho, cmems_land,'parent',ax1_1);
        colormap(ax1_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                'xticklabels', [], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax1_1);  
        hold on

        ax1_2=axes;
        set(ax1_2,'pos',pos_sb1);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        
        quiv1_2=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                        cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                        summer_mean_u(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                        summer_mean_v(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                        'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth, ...
                        'parent', ax1_2);
        CMEMS_lon = cut_lon_rho;
        CMEMS_lat = cut_lat_rho;
        CMEMS_Svel.u=u_rho;
        CMEMS_Svel.v=v_rho;
%         pc1_2=m_pcolor(cut_lon_rho,cut_lat_rho, mean_data,'parent',ax1_2);
        hold on
        m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_grid_fontsize+2); 
                  
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax1_2);  
        txt_1_2=m_text(118, 49, '(A) CMEMS', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- CMEMS

%      ---------- GCM-IPSL-L
        test_text={'(B) GCM-IPSL-L', '(C) GCM-IPSL-M', '(D) GCM-Nor', '(E) GCM-MPI', '(F) RCM-IPSL-L', '(G) RCM-IPSL-M', '(H) RCM-Nor', '(I) RCM-MPI'};
        for testnameind=1:4
            testname=all_testname1{testnameind};    

            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_CMIP5_', regionname, '.m'];
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([cmip5dir,regionname, '_', testname, '_cmip5_land_vec','.mat']);  % from SSH_paper2_003_plot_decadal_status_raw_historical
            cmip5_land=cmip5_land_vec;
%             if testnameind==3
%                 cmip5_land(end+1,:)=cmip5_land(end,:);
%             elseif testnameind==4
%                 cmip5_land=cmip5_land(1:end-1,:);
%             end

            
            load([cmip5dir,regionname,'_',testname, '_cmip5_', variable,'_summer_',yearstr_min, '_', yearstr_max, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;
            
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,4,[1+8*testnameind 2+8*testnameind 5+8*testnameind 6+8*testnameind]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', cmip5_land','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on
            
            cut_lon_rho=cut_lon_rho+1;
            if testnameind==4
                cut_lat_rho=cut_lat_rho-1;
                if (exist('ref_vec_x_range' , 'var') ~= 1)
                    ref_vec_x_ind = find(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location) == min(abs(cut_lon_rho(:,1)-m_quiver_ref_text_x_location)));
                    ref_vec_y_ind = find(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location) == min(abs(cut_lat_rho(1,:)-m_quiver_ref_text_y_location)))+m_quiver_y_interval;
                    ref_vec_x_range = round(ref_vec_x_ind-(m_quiver_x_interval/2)) : round(ref_vec_x_ind-(m_quiver_x_interval/2));
                    ref_vec_y_range = round(ref_vec_y_ind-(m_quiver_y_interval/2)) : round(ref_vec_y_ind-(m_quiver_y_interval/2));
                end
                u_rho(ref_vec_x_range,ref_vec_y_range+1)=NaN;
                v_rho(ref_vec_x_range,ref_vec_y_range+1)=NaN;
                u_rho(ref_vec_x_range,ref_vec_y_range-1)=m_quiver_ref_u_value;
                v_rho(ref_vec_x_range,ref_vec_y_range-1)=m_quiver_ref_v_value;   
            end
            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            
            quiv{testnameind,2}=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth, ...
                            'parent', ax{testnameind,2});
            GCM_lon{testnameind}=cut_lon_rho;
            GCM_lat{testnameind}=cut_lat_rho;
            GCM_Svel.u{testnameind}=u_rho;
            GCM_Svel.v{testnameind}=v_rho;
                        
            hold on
            m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_grid_fontsize+2); 

            if testnameind==4
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            else
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            end
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
        end
%      ---------- GCM-IPSL-L

drivename='J';
        for testnameind=5:8

%      ---------- RCM-IPSL-L
%         testnameind=5;
            testname=all_testname2{testnameind-4};    

            param_script =['C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_', regionname,'.m'];
            filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([filedir,regionname, '_', testname, '_model_land','.mat']);
            load([filedir,regionname,'_',testname, '_model_summer_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

%             mean_data= mean_data .* mask_model;        
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,4,[3+8*(testnameind-4) 4+8*(testnameind-4) 7+8*(testnameind-4) 8+8*(testnameind-4)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            u_rho=u_rho .*mask_model;
            v_rho=v_rho .*mask_model;
            quiv{testnameind,2}=m_quiver(cut_lon_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            cut_lat_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)', ...
                            u_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            v_rho(1:m_quiver_x_interval:end,1:m_quiver_y_interval:end)' * m_quiver_vector_size, ...
                            'color',m_quiver_vector_color, 'AutoScale','off','LineWidth', m_quiver_LineWidth, ...
                            'parent', ax{testnameind,2});
            RCM_lon{testnameind-4}=cut_lon_rho;
            RCM_lat{testnameind-4}=cut_lat_rho;
            RCM_Svel.u{testnameind-4}=u_rho;
            RCM_Svel.v{testnameind-4}=v_rho;
            hold on
            m_text(m_quiver_ref_text_x_location, m_quiver_ref_text_y_location, m_quiver_ref_text, 'FontSize', m_grid_fontsize+2); 


            if testnameind==8
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            else
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            end
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

        end
      
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)
        
        hold off;
        close all;
        
        save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\Data04_SVel.mat', ...
            'CMEMS_lon', 'CMEMS_lat', 'CMEMS_Svel', ...
            'GCM_lon', 'GCM_lat', 'GCM_Svel', ...
            'RCM_lon', 'RCM_lat', 'RCM_Svel', ...
            'GCM_testnames', 'RCM_testnames')
        
        
end