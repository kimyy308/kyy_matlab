close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR'};
all_testname2  = {'test53', 'test54', 'test55', 'test56'};
GCM_testnames = {'GCM-IPSL-L', 'GCM-IPSL-M', 'GCM-Nor', 'GCM-MPI'};
RCM_testnames = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI'};
all_region2 ={'AKP4'};

variable = 'SST';
scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [1982:2005]; % % put year which you want to plot [year year ...]
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
            figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig\'); % % where figure files will be saved
            cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are
        elseif (strcmp(system_name,'GLNXA64'))
        end

        correction_right_fig=[-0.1000,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0.020,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig03','.tif'); %% ~_year_month.jpg

%         f1=figure(1);

%      ---------- AVHRR
        testnameind2=1;
        testname=all_testname2{testnameind2};    
        load(['D:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',regionname,'sst_rms_and_bias_',num2str(inputyear(1),'%04i'),'_',num2str(inputyear(end)),'.mat']);
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
        mean_data=mean(comb_spatial_meanavhrr,3,'omitnan');
        cut_lon_rho =avhrr_lon;
        cut_lat_rho =avhrr_lat;
        AVHRR_lon=avhrr_lon;
        AVHRR_lat=avhrr_lat;
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
%         mask_model_inv = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
%         mask_model_inv(mask_model_inv==0)=NaN;
%         mask_model_inv=ones(size(mean_data));
%         mask_model_inv(isfinite(mean_data))=NaN;
        load([avhrrdir,'avhrr_land.mat']);
        
%         avhrrfilename=['Z:\내 드라이브\Data\Observation\OISST\monthly_kimyy\avhrr_only_monthly_v2_2005.nc'];
%         avhrr_sst = ncread(avhrrfilename,varname,[avhrr_lon_min(1) avhrr_lat_min(1) monthij], [avhrr_lon_max(1)-avhrr_lon_min(1) avhrr_lat_max(1)-avhrr_lat_min(1) 1]);
%         avhrr_land(isfinite(avhrr_sst))=NaN;
%         avhrr_land(isnan(avhrr_sst))=1;
%         save([avhrrdir,'avhrr_land.mat'], 'avhrr_land');
        
        mean_data= mean_data .* mask_model;
        AVHRR_SST = mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb1=subplot(10,4,[2 3 6 7]);  % Auto-fitted to the figure.
        pos_sb1=get(sb1, 'pos'); % Get the position.
        pos_sb1 = pos_sb1 + correction_upper_fig+ correction_large_fig;
        delete(sb1); % Delete the subplot axes
        ax1_1=axes;
        set(ax1_1,'pos',pos_sb1);
        pc1_1=m_pcolor(cut_lon_rho',cut_lat_rho', avhrr_land','parent',ax1_1);
        colormap(ax1_1,[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax1_1);  
        hold on

        ax1_2=axes;
        set(ax1_2,'pos',pos_sb1);
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc1_2=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax1_2);
        colormap(ax1_2,jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C1_2,h1_2]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax1_2);
                clabel(C1_2,h1_2,'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax1_2);  
        txt_1_2=m_text(118, 49, '(A) AVHRR', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- AVHRR

%      ---------- GCM-IPSL-L
        testnameind=1;
        testname=all_testname1{testnameind};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
        
        load([cmip5dir,regionname, '_', testname, '_cmip5_land','.mat']);
        load([cmip5dir,regionname,'_',testname, '_cmip5_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model;     
        
        GCM_lon{testnameind}=cut_lon_rho;
        GCM_lat{testnameind}=cut_lat_rho;
        GCM_SST{testnameind}=mean_data;
        
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[9 10 13 14]);  % Auto-fitted to the figure.
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
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(B) GCM-IPSL-L', 'FontSize', m_grid_fontsize+4); 

%      ---------- GCM-IPSL-L

%      ---------- GCM-IPSL-M
        testnameind=2;
        testname=all_testname1{testnameind};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([cmip5dir,regionname, '_', testname, '_cmip5_land','.mat']);
        load([cmip5dir,regionname,'_',testname, '_cmip5_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model; 
         GCM_lon{testnameind}=cut_lon_rho;
        GCM_lat{testnameind}=cut_lat_rho;
        GCM_SST{testnameind}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[17 18 21 22]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', cmip5_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(c) GCM-IPSL-M', 'FontSize', m_grid_fontsize+4); 

%      ---------- GCM-IPSL-M

%      ---------- GCM-Nor
        testnameind=3;
        testname=all_testname1{testnameind};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
        avhrrdir='Z:\내 드라이브\Data\Observation\OISST\monthly\';

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([cmip5dir,regionname, '_', testname, '_cmip5_land','.mat']);
        load([cmip5dir,regionname,'_',testname, '_cmip5_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model;  
         GCM_lon{testnameind}=cut_lon_rho;
        GCM_lat{testnameind}=cut_lat_rho;
        GCM_SST{testnameind}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[25 26 29 30]);  % Auto-fitted to the figure.
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
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(d) GCM-Nor', 'FontSize', m_grid_fontsize+4); 

%      ---------- GCM-Nor

%      ---------- GCM-MPI
        testnameind=4;
        testname=all_testname1{testnameind};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([cmip5dir,regionname, '_', testname, '_cmip5_land','.mat']);
        load([cmip5dir,regionname,'_',testname, '_cmip5_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model;
        GCM_lon{testnameind}=cut_lon_rho;
        GCM_lat{testnameind}=cut_lat_rho;
        GCM_SST{testnameind}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[33 34 37 38]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', cmip5_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(e) GCM-MPI', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- GCM-MPI

drivename='J';

%      ---------- RCM-IPSL-L
        testnameind=5;
        testname=all_testname2{testnameind-4};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([filedir,regionname, '_', testname, '_model_land','.mat']);
        load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model;  
        RCM_lon{testnameind-4}=cut_lon_rho;
        RCM_lat{testnameind-4}=cut_lat_rho;
        RCM_SST{testnameind-4}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[11 12 15 16]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(f) RCM-IPSL-L', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- RCM-IPSL-L


%      ---------- RCM-IPSL-M
        testnameind=6;
        testname=all_testname2{testnameind-4};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([filedir,regionname, '_', testname, '_model_land','.mat']);
        load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model;    
        RCM_lon{testnameind-4}=cut_lon_rho;
        RCM_lat{testnameind-4}=cut_lat_rho;
        RCM_SST{testnameind-4}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[19 20 23 24]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(g) RCM-IPSL-M', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- RCM-IPSL-M

%      ---------- RCM-Nor
        testnameind=7;
        testname=all_testname2{testnameind-4};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([filedir,regionname, '_', testname, '_model_land','.mat']);
        load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model; 
        RCM_lon{testnameind-4}=cut_lon_rho;
        RCM_lat{testnameind-4}=cut_lat_rho;
        RCM_SST{testnameind-4}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[27 28 31 32]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(h) RCM-Nor', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- RCM-Nor

%      ---------- RCM-MPI
        testnameind=8;
        testname=all_testname2{testnameind-4};    
     
        param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
        filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

        run(param_script);
        m_grid_fontsize=m_grid_fontsize-12;
 
        load([filedir,regionname, '_', testname, '_model_land','.mat']);
        load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min, '_', yearstr_max, '.mat']);
        mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
        mask_model(mask_model==0)=NaN;
  
        mean_data= mean_data .* mask_model; 
        RCM_lon{testnameind-4}=cut_lon_rho;
        RCM_lat{testnameind-4}=cut_lat_rho;
        RCM_SST{testnameind-4}=mean_data;
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        sb{testnameind}=subplot(10,4,[35 36 39 40]);  % Auto-fitted to the figure.
        pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
        pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
        delete(sb{testnameind}); % Delete the subplot axes
        ax{testnameind,1}=axes;
        set(ax{testnameind,1},'pos',pos_sb{testnameind});
        pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
        colormap(ax{testnameind,1},[0.8 0.8 0.8]);
        shading(gca,m_pcolor_shading_method); 
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type, ...
                'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
        hold on

        ax{testnameind,2}=axes;
        set(ax{testnameind,2},'pos',pos_sb{testnameind});
        m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
        pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
        colormap(ax{testnameind,2},jet);
        caxis(colorbar_lev);
        shading(gca,m_pcolor_shading_method);   
        hold on
        [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                    'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                    'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');
                
        m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140],  'box', m_grid_box_type, ...
                'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
        txt{testnameind,2}=m_text(118, 49, '(i) RCM-MPI', 'FontSize', m_grid_fontsize+4); 

        
%      ---------- RCM-MPI





%         sb4=subplot(10,4,[33 34 37 38]);  % Auto-fitted to the figure.
%         pos_sb4=get(sb4, 'pos');
        h = colorbar('southoutside','AxisLocation','out');
%         colormap(gmt_ocean_mod2);
        set(h,'fontsize',m_grid_fontsize+5);
        h_title=title(h,colorbar_title,'fontsize',m_grid_fontsize+5);
        caxis(colorbar_lev);
%         h.ax.xaxis.set_ticks_position('right')
%         set(h,...
%         'YTickLabel',{'0^oC', '5^oC','10^oC','15^oC','20^oC','25^oC','30^oC'},...
%                 'YTick',[0 5 10 15 20 25 30]);
        
        set(h_title, 'Position', [pos_sb{4}(1)+330.1+ pos_sb{4}(3)*2, pos_sb{4}(2)+2.5 ])  % right, up
        set(h, 'Position', [pos_sb{4}(1)+0.052, pos_sb{4}(2)-0.05, pos_sb{4}(3)*2-0.05, 0.0231])  % right, up, width, height

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
        save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\Data03_SST.mat', ...
            'AVHRR_lon', 'AVHRR_lat', 'AVHRR_SST', ...
            'GCM_lon', 'GCM_lat', 'GCM_SST', ...
            'RCM_lon', 'RCM_lat', 'RCM_SST', ...
            'GCM_testnames', 'RCM_testnames')
        
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