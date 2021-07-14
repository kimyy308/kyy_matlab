close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', 'ens08'};
all_testname2  = {'test57', 'test58', 'test59', 'test60', 'ens08'};
% all_testname2  = {'test53', 'test54', 'test55', 'test56', 'ens03'};
GCM_testnames = {'GCM-IPSL-L', 'GCM-IPSL-M', 'GCM-Nor', 'GCM-MPI', 'GCM-Ens'};
RCM_testnames = {'RCM-IPSL-L', 'RCM-IPSL-M', 'RCM-Nor', 'RCM-MPI', 'RCM-Ens'};
all_region2 ={'AKP4'};

variable = 'SSH';
scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;
        sl_caxis= ([-0.2 0.6]);
        inputyear1 = [2081]; % % put year which you want to plot [year year ...]
        inputyear2 = [2090]; % % put year which you want to plot [year year ...]
        inputyear3 = [2100]; % % put year which you want to plot [year year ...]
        yearstr_min1=num2str(inputyear1(1));
        yearstr_max1=num2str(inputyear1(end));
        yearstr_min2=num2str(inputyear2(1));
        yearstr_max2=num2str(inputyear2(end));
        yearstr_min3=num2str(inputyear3(1));
        yearstr_max3=num2str(inputyear3(end));
        drivename='I';
        
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            error_status=Func_0008_set_dropbox_path(system_name);
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            [error_status, refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
            load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)
        end

        figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r3\'); % % where figure files will be saved
        cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        
        tifname=strcat(figrawdir, 'R2Q5C8_rcp45','.tif'); %% ~_year_month.jpg
        ystr1=num2str(inputyear1);
        ystr2=num2str(inputyear2);
        ystr3=num2str(inputyear3);
%         test_text={['(A) RCM-IPSL-L (', ystr1,')'], ['(B) RCM-IPSL-M (', ystr1,')'], ['(C) RCM-Nor (', ystr1,')'], ['(D) RCM-MPI (', ystr1,')'], ['(E) RCM-Ens (', ystr1,')'], ...
%             ['(F) RCM-IPSL-L (', ystr2,')'], ['(G) RCM-IPSL-M (', ystr2,')'], ['(H) RCM-Nor (', ystr2,')'], ['(I) RCM-MPI (', ystr2,')'], ['(J) RCM-Ens (', ystr2,')'], ...
%             ['(K) RCM-IPSL-L (', ystr3,')'], ['(L) RCM-IPSL-M (', ystr3,')'], ['(M) RCM-Nor (', ystr3,')'], ['(N) RCM-MPI (', ystr3,')'], ['(O) RCM-Ens (', ystr3,')']};
        
        test_text={['IPSL-L (', ystr1,')'], ['IPSL-M (', ystr1,')'], ['Nor (', ystr1,')'], ['MPI (', ystr1,')'], ['Ens (', ystr1,')'], ...
            ['IPSL-L (', ystr2,')'], ['IPSL-M (', ystr2,')'], ['Nor (', ystr2,')'], ['MPI (', ystr2,')'], ['Ens (', ystr2,')'], ...
            ['IPSL-L (', ystr3,')'], ['IPSL-M (', ystr3,')'], ['Nor (', ystr3,')'], ['MPI (', ystr3,')'], ['Ens (', ystr3,')']};

%      ---------- GCM
        for testnameind=1:5
%             testname=all_testname1{testnameind};    
            rcmtestname = all_testname2{testnameind};  
            testname=rcmtestname;
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
%             run(param_script);
            
            filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([filedir,regionname, '_', testname, '_model_land','.mat']);
            load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min1, '_', yearstr_max1, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= mean_data .* mask_model;
            
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[1+12*(testnameind-1) 2+12*(testnameind-1) 7+12*(testnameind-1) 8+12*(testnameind-1)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
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
            colormap(ax{testnameind,2},flip(cool));
            caxis(sl_caxis);
            shading(gca,m_pcolor_shading_method);   
            hold on
            
           if testnameind==5
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            else
                m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [], 'xtick',[120,  130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
           end
           txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4);
%            GCM_lon{testnameind}=lon_glo;
%            GCM_lat{testnameind}=lat_glo;
%            GCM_steric_SLR{testnameind}=mean_data;
        end
%      ---------- GCM


%      ---------- RCM
        for testnameind=6:10
%         testnameind=5;
            testname=all_testname2{testnameind-5};    

            rcmtestname = all_testname2{testnameind-5};  
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
%             run(param_script);
%             m_grid_fontsize=m_grid_fontsize-12;
            
            filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([filedir,regionname, '_', testname, '_model_land','.mat']);
            load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min2, '_', yearstr_max2, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= mean_data .* mask_model;
            
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[3+12*(testnameind-6) 4+12*(testnameind-6) 9+12*(testnameind-6) 10+12*(testnameind-6)]);  % Auto-fitted to the figure.
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
            colormap(ax{testnameind,2},flip(cool));
            caxis(sl_caxis);
            shading(gca,m_pcolor_shading_method);   
            hold on
%             [C{testnameind,2},h{testnameind,2}]=m_contour(lon_glo',lat_glo', mean_data', [10, 15, 20], 'color','k', ...
%                         'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
%                     clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
%                         'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');

            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
%             RCM_lon{testnameind-5}=lon_glo;
%             RCM_lat{testnameind-5}=lat_glo;
%             RCM_steric_SLR{testnameind-5}=mean_data;
        end

%         h = colorbar('westoutside','AxisLocation','out');
%         caxis([0 70]);
%         
%         set(h,'fontsize',m_grid_fontsize+7);
%         h_title=title(h,'(cm)','fontsize',m_grid_fontsize+5);
%         
%         set(h, 'Position', [pos_sb{10}(1)-0.092, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.33])  % right, up, width, height
        
        
        
%      ---------- diff
        for testnameind=11:15
%         testnameind=5;
            testname=all_testname2{testnameind-10};    
            
%             run(param_script);
%             m_grid_fontsize=m_grid_fontsize-12;
            
            filedir = strcat(drivename, ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-9;

            load([filedir,regionname, '_', testname, '_model_land','.mat']);
            load([filedir,regionname,'_',testname, '_model_', variable,'_',yearstr_min3, '_', yearstr_max3, '.mat']);
            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= mean_data .* mask_model;

            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[5+12*(testnameind-11) 6+12*(testnameind-11) 11+12*(testnameind-11) 12+12*(testnameind-11)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig*2;
%             pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(cut_lon_rho',cut_lat_rho', model_land','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], 'YaxisLocation', 'right', ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(cut_lon_rho',cut_lat_rho', mean_data','parent',ax{testnameind,2});
            [byrmap, er_status] = Func_0009_get_colormaps('byr');
            colormap(ax{testnameind,2},flip(cool));
            caxis(sl_caxis);
            shading(gca,m_pcolor_shading_method);   
            hold on

            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
%             DIFF_lon{testnameind-10}=lon_glo;
%             DIFF_lat{testnameind-10}=lat_glo;
%             DIFF_steric_SLR{testnameind-10}=mean_data;
        end

        h2 = colorbar('eastoutside','AxisLocation','out');
        caxis(sl_caxis);
        
        set(h2,'fontsize',m_grid_fontsize+7);
        h_title=title(h2,'(m)','fontsize',m_grid_fontsize+5);
        
        set(h2, 'Position', [pos_sb{15}(1)+0.202, pos_sb{15}(2)+0, 0.0231, pos_sb{15}(3)*2+0.33])  % right, up, width, height
%         set(h, 'Position', [pos_sb{10}(1)-0.132, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.33])  % right, up, width, height

      
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2.3 paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
%         save ('Z:\내 드라이브\research\Ph_D_course\2020_SSH_CMIP5_decadal_variation_around_korea_using_downscaling_ocean_model\data_dryad\DataS01_steric_SLR_rcp45.mat', ...
%             'GCM_lon', 'GCM_lat', 'GCM_steric_SLR', ...
%             'RCM_lon', 'RCM_lat', 'RCM_steric_SLR', ...
%             'DIFF_lon', 'DIFF_lat', 'DIFF_steric_SLR', ...
%             'GCM_testnames', 'RCM_testnames')
end

function land_model_interped = get_land_grid(filename,lonlat,lon_glo,lat_glo)
modelinfo=ncinfo(filename);
lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

lon_west = abs(lon - (lonlat(1)-1));
[min_lon_west, lon_min]=min(lon_west);
lon_east = abs(lon - (lonlat(2)+1));
[min_lon_east, lon_max]=min(lon_east);
lat_south = abs(lat - (lonlat(3)-1));
[min_lat_south, lat_min]=min(lat_south);
lat_north = abs(lat - (lonlat(4)+1));
[min_lat_north, lat_max]=min(lat_north);

lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
mask_rho = ncread(filename,'mask_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
interped_mask_rho = griddata(double(lon), double(lat), mask_rho,double(lon_glo),double(lat_glo));   
land_model_interped=interped_mask_rho;
land_model_interped(interped_mask_rho==0)=1;
interped_mask_rho(interped_mask_rho==0)=NaN;
land_model_interped(isfinite(interped_mask_rho))=NaN;
end