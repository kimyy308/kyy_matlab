close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
all_testname2  = {'test53', 'test54', 'test55', 'test56'};

all_region2 ={'AKP4'};

variable = 'DET_STD';
% scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [1993:2005]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
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
        tifname=strcat(figrawdir, 'fig06','.tif'); %% ~_year_month.jpg

%         f1=figure(1);

%      ---------- CMEMS
%         testnameind2=1;
        testname='test53';    
%         load(['D:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        cmemsdir='Z:\내 드라이브\Data\Observation\CMEMS\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-9;
        load([cmemsdir,regionname, '_', 'CMEMS', '_cmems_land','.mat']);
        cmems_land(:,89)=cmems_land(:,88);
        cmems_land(:,90)=cmems_land(:,88);

        load([filedir, 'AKP4_cmems_det_std_1993_2005.mat']);
        
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
        pc1_2=m_pcolor(cut_lon_rho,cut_lat_rho, model_sla_std,'parent',ax1_2);
        colormap(ax1_2, parula);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);
        hold on
                       
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax1_2);  
        txt_1_2=m_text(118, 49, '(A) CMEMS', 'FontSize', m_grid_fontsize+4); 
        
        mean_cmems_data = model_sla_std;
        

%      ---------- CMEMS

%      ---------- GCM-IPSL-L
        test_text={'(B) GCM-IPSL-L', '(C) GCM-IPSL-M', '(D) GCM-Nor', '(E) GCM-MPI', '(F) RCM-IPSL-L', '(G) RCM-IPSL-M', '(H) RCM-Nor', '(I) RCM-MPI'};
        for testnameind=1:4
            testname=all_testname1{testnameind};    

            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;
            load([cmip5dir,regionname, '_', testname, '_cmip5_interped_land','.mat']);
            load([cmip5dir,regionname,'_',testname, '_model_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= model_sla_std .* mask_model;  
            ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
            glo_mean_data_ysecs(:,:,testnameind) = mean_data.* ysecsmask;
            mean_data_ysecs = glo_mean_data_ysecs(:,:,testnameind);
            mstd_ysecs(testnameind)=mean(glo_mean_data_ysecs(:), 'omitnan');
            stdstd_ysecs(testnameind)=std(glo_mean_data_ysecs(:), 'omitnan');
            esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);
            glo_mean_data_es(:,:,testnameind) = mean_data.* esmask;
            mean_data_es = glo_mean_data_es(:,:,testnameind);
            mstd_es(testnameind)=mean(glo_mean_data_es(:), 'omitnan');
            stdstd_es(testnameind)=std(glo_mean_data_es(:), 'omitnan');
            
            glo_mask_ysecs(:,:,testnameind)= glo_mean_data_ysecs(:,:,testnameind);
            glo_mask_ysecs(isfinite(glo_mask_ysecs))=1;
            glo_mask_es(:,:,testnameind)= glo_mean_data_es(:,:,testnameind);
            glo_mask_es(isfinite(glo_mask_es))=1;
            
            mean_cmems_data_ysecs = mean_cmems_data.* glo_mask_ysecs(:,:,testnameind);
            mstd_cmems_ysecs(testnameind)=mean(mean_cmems_data_ysecs(:), 'omitnan');
            stdstd_cmems_ysecs(testnameind)=std(mean_cmems_data_ysecs(:), 'omitnan');
            mean_cmems_data_es = mean_cmems_data.* glo_mask_es(:,:,testnameind);
            mstd_cmems_es(testnameind)=mean(mean_cmems_data_es(:), 'omitnan');
            stdstd_cmems_es(testnameind)=std(mean_cmems_data_es(:), 'omitnan');
            
            aaa= mean_data_ysecs(logical(~isnan(mean_data_ysecs).*~isnan(mean_cmems_data_ysecs)));
            bbb= mean_cmems_data_ysecs(logical(~isnan(mean_data_ysecs).*~isnan(mean_cmems_data_ysecs)));
            ccc=corrcoef(aaa,bbb);
            pccs_ysecs(testnameind) = ccc(1,2);
            aaa= mean_data_es(logical(~isnan(mean_data_es).*~isnan(mean_cmems_data_es)));
            bbb= mean_cmems_data_es(logical(~isnan(mean_data_es).*~isnan(mean_cmems_data_es)));
            ccc= corrcoef(aaa,bbb);
            pccs_es(testnameind) = ccc(1,2);
            
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

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},parula);
            caxis(colorbar_lev);
            shading(gca,m_pcolor_shading_method);   
            hold on
            
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

            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            landdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', 'test65', '\run\'); % % where data files are
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([landdir,regionname, '_', 'test65', '_model_interped_land','.mat']);
            load([filedir,regionname,'_',testname, '_model_', 'det_std','_',yearstr_min, '_', yearstr_max, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= model_sla_std .* mask_model;     
            ysecsmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, ysecspolygon);
            mean_data_ysecs = mean_data.* ysecsmask .* glo_mask_ysecs(:,:,testnameind-4);
            mstd_ysecs(testnameind)=mean(mean_data_ysecs(:), 'omitnan');
            stdstd_ysecs(testnameind)=std(mean_data_ysecs(:), 'omitnan');
            esmask=Func_0005_get_mask_from_polygon(cut_lon_rho, cut_lat_rho, espolygon);
            mean_data_es = mean_data.* esmask .* glo_mask_es(:,:,testnameind-4);
            mstd_es(testnameind)=mean(mean_data_es(:), 'omitnan');
            stdstd_es(testnameind)=std(mean_data_es(:), 'omitnan');
            maxstd_es(testnameind)=max(mean_data_es(:));
            mean_cmems_data_ysecs = mean_cmems_data.* glo_mask_ysecs(:,:,testnameind-4);
            mean_cmems_data_es = mean_cmems_data.* glo_mask_es(:,:,testnameind-4);
            
            aaa= mean_data_ysecs(logical(~isnan(mean_data_ysecs).*~isnan(mean_cmems_data_ysecs)));
            bbb= mean_cmems_data_ysecs(logical(~isnan(mean_data_ysecs).*~isnan(mean_cmems_data_ysecs)));
            ccc=corrcoef(aaa,bbb);
            pccs_ysecs(testnameind) = ccc(1,2);
            aaa= mean_data_es(logical(~isnan(mean_data_es).*~isnan(mean_cmems_data_es)));
            bbb= mean_cmems_data_es(logical(~isnan(mean_data_es).*~isnan(mean_cmems_data_es)));
            ccc= corrcoef(aaa,bbb);
            pccs_es(testnameind) = ccc(1,2);
            
            
            
            
            mean(mean_data(:), 'omitnan')
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
            pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},parula);
            caxis(colorbar_lev);
            shading(gca,m_pcolor_shading_method);   
            hold on

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


%         sb4=subplot(10,4,[33 34 37 38]);  % Auto-fitted to the figure.
%         pos_sb4=get(sb4, 'pos');
%         h = colorbar('southoutside','AxisLocation','out');
        
%         for colbarind=1:2
%             if colbarind==1
%                 h{colbarind} = colorbar(ax1_1); 
%                  caxis(ax1_1,colorbar_lev);
%             else
%                 h{colbarind} = colorbar(ax1_2);
%                 caxis(ax1_2,colorbar_lev);
%             end
%         
%         set(h{colbarind},'fontsize',m_grid_fontsize+5);
%         h_title{colbarind}=title(h{colbarind},colorbar_title,'fontsize',m_grid_fontsize+5);
%        
% %         set(h,...
% %         'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC','30^oC'},...
% %                 'YTick',[0 5 10 15 20 25 30]);
%         
% %         set(h_title, 'Position', [pos_sb{4}(1)+330.1+ pos_sb{4}(3)*2, pos_sb{4}(2)+2.5 ])  % right, up
% %         set(h, 'Position', [pos_sb{4}(1)+0.052, pos_sb{4}(2)-0.05, pos_sb{4}(3)*2-0.05, 0.0231])  % right, up, width, height
%         end
        
   
        h = colorbar('eastoutside','AxisLocation','out');
        caxis(colorbar_lev);
        
        set(h,'fontsize',m_grid_fontsize+7);
        h_title=title(h,colorbar_title,'fontsize',m_grid_fontsize+5);
       
%         set(h,...
%         'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC','30^oC'},...
%                 'YTick',[0 5 10 15 20 25 30]);
        
%         set(h_title, 'Position', [pos_sb{4}(1)+330.1+ pos_sb{4}(3)*2, pos_sb{4}(2)+2.5 ])  % right, up
        set(h, 'Position', [pos_sb{8}(1)+0.312, pos_sb{8}(2)+0, 0.0231, pos_sb{8}(3)*2+0.09])  % right, up, width, height
        
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
end