close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', 'ens08'};
all_testname2  = {'test57', 'test58', 'test59', 'test60', 'ens08'};

all_region2 ={'AKP4'};

variable = 'SSH';
scenname='rcp45';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
        yearstr_min=num2str(inputyear(1));
        yearstr_max=num2str(inputyear(end));
        inputmonth = [1:12]; % % put month which you want to plot [month month ...]
        system_name=computer;
        for folding=1:1
            error_status=Func_0008_set_dropbox_path(system_name);
            run('nwp_polygon_point.m');
            regionname=all_region2{regionind2};
            [refpolygon, lonlat] = Func_0007_get_polygon_data_from_regionname(regionname);
            load('C:\Users\User\Dropbox\source\matlab\Common\Figure\gmt_ocean_mod2.mat')  % % set colormap (gmt_ocean, nonwhite)
        end

        figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\paper2\fin_fig_r1\'); % % where figure files will be saved
        cmip5dir = strcat('D:\Data\Model\CMIP5\'); % % where data files are

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        
        tifname=strcat(figrawdir, 'fig09','.tif'); %% ~_year_month.jpg
        
        test_text={'(A) GCM-IPSL-L', '(B) GCM-IPSL-M', '(C) GCM-Nor', '(D) GCM-MPI', '(E) GCM-Ens', ...
            '(F) RCM-IPSL-L', '(G) RCM-IPSL-M', '(H) RCM-Nor', '(I) RCM-MPI', '(J) RCM-Ens' ...
            '(K) DIFF-IPSL-L', '(L) DIFF-IPSL-M', '(M) DIFF-Nor', '(N) DIFF-MPI', '(O) DIFF-Ens'};
%      ---------- GCM
        for testnameind=1:5
            testname=all_testname1{testnameind};    
            rcmtestname = all_testname2{testnameind};  
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
%             filedir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', rcmtestname, '\run\'); % % where data files are
%             interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-12;
%             load([cmip5dir,regionname, '_', testname, '_cmip5_interped_land','.mat']);
             

% get 'rcm_slc', 'mean_rcm_slc', 'rcm_steric_slc', 'mean_rcm_steric_slc', 'rcm_mano_slc', 'mean_rcm_mano_slc',
% get 'glo_slc', 'mean_glo_slc', 'glo_steric_slc', 'mean_glo_steric_slc', 'glo_mano_slc', 'mean_glo_mano_slc'
            gcmcount = 4;
            if testnameind<=gcmcount
                load([rcmdir, rcmtestname,'_',regionname, '_rcm_gcm_glo_interped_slc_hist_rcp.mat']);
                all_slc(:,:,testnameind) = glo_slc.*100;
                all_steric_slc(:,:,testnameind) = glo_steric_slc.*100;
                all_mano_slc(:,:,testnameind) = glo_mano_slc.*100;
                all_mean_slc(testnameind,1) = mean_glo_slc.*100;
                all_mean_steric_slc(testnameind,1) = mean_glo_steric_slc.*100;
                all_mean_mano_slc(testnameind,1) = mean_glo_mano_slc.*100;
%                 load([rcmdir,rcmtestname,'_',regionname,'land_model_only_steric_ssh_',num2str(1976,'%04i'),'_',num2str(2005,'%04i'),'.mat']);
                load([rcmdir,rcmtestname,'_',regionname,'land_steric_ssh_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
                land_model_glo2 = land_model_glo;
            else
                all_slc(:,:,testnameind) = mean(all_slc(:,:,gcmcount-3:gcmcount),3);
                all_steric_slc(:,:,testnameind) = mean(all_steric_slc(:,:,gcmcount-3:gcmcount),3);
                all_mano_slc(:,:,testnameind) = mean(all_mano_slc(:,:,gcmcount-3:gcmcount),3);
                all_mean_slc(testnameind,1) = mean(all_mean_slc(gcmcount-3:gcmcount));
                all_mean_steric_slc(testnameind,1) = mean(all_mean_steric_slc(gcmcount-3:gcmcount));
                all_mean_mano_slc(testnameind,1) = mean(all_mean_mano_slc(gcmcount-3:gcmcount));
            end
            
            filename = strcat('D:\Data\Model\ROMS\nwp_1_20\test53\run\2005\ocean_rst2.nc');
            land_model_interped=get_land_grid(filename,lonlat,lon_glo,lat_glo);
            
            
%             sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');

           
% %             land_model = land_model_interped;
% %             land_model(isfinite(land_model_glo2))=1;
            
%             cut_lon_rho=ncread(interpedfilename,'lon_cmems');
%             cut_lat_rho=ncread(interpedfilename,'lat_cmems');
%             mean_data=sl(:,:,end)-sl(:,:,1);
%             mask_model = double(inpolygon(lon_glo,lat_glo,refpolygon(:,1),refpolygon(:,2)));
%             mask_model(mask_model==0)=NaN;
            
            mean_data = all_slc(:,:,testnameind);
            
            ysecsmask=Func_0005_get_mask_from_polygon(lon_glo, lat_glo, ysecs_khoapolygon);
            esmask=Func_0005_get_mask_from_polygon(lon_glo, lat_glo, es_khoapolygon);
            ysecs_slc=all_slc(:,:,testnameind).*ysecsmask;
            es_slc=all_slc(:,:,testnameind).*esmask;
            [mean_all_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(mean_data, lon_glo, lat_glo);            
            [mean_ysecs_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(ysecs_slc, lon_glo, lat_glo);
            [mean_es_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(es_slc, lon_glo, lat_glo);
           
            ysecs_std_all_slc(testnameind,1) = std(ysecs_slc(:),'omitnan');
            es_std_all_slc(testnameind,1) = std(es_slc(:),'omitnan');
            es_max_all_slc(testnameind,1) = max(es_slc(:));
            es_min_all_slc(testnameind,1) = min(es_slc(:));
%             mean_data= mean_data .* mask_model;        
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[1+12*(testnameind-1) 2+12*(testnameind-1) 7+12*(testnameind-1) 8+12*(testnameind-1)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_glo', lat_glo', land_model_glo2','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  'box', m_grid_box_type, ...
                    'xticklabels', [], 'xtick',[120, 130, 140], ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_glo', lat_glo', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},flip(autumn));
            caxis([30 70]);
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
        end
%      ---------- GCM


%      ---------- RCM
        for testnameind=6:10
%         testnameind=5;
            testname=all_testname2{testnameind-5};    

            rcmtestname = all_testname2{testnameind-5};  
            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
%             filedir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are
            rcmdir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', rcmtestname, '\run\'); % % where data files are
%             filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
%             landdir = strcat('G', ':\Data\Model\ROMS\nwp_1_20\', 'test65', '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-12;
            
            rcmcount = 9;
            if testnameind<=rcmcount
                load([rcmdir, testname,'_',regionname, '_rcm_gcm_glo_interped_slc_hist_rcp.mat']);
                all_slc(:,:,testnameind) = rcm_slc.*100;
                all_steric_slc(:,:,testnameind) = rcm_steric_slc.*100;
                all_mano_slc(:,:,testnameind) = rcm_mano_slc.*100;
                all_mean_slc(testnameind,1) = mean_rcm_slc.*100;
                all_mean_steric_slc(testnameind,1) = mean_rcm_steric_slc.*100;
                all_mean_mano_slc(testnameind,1) = mean_rcm_mano_slc.*100;
            else
                all_slc(:,:,testnameind) = mean(all_slc(:,:,rcmcount-3:rcmcount),3);
                all_steric_slc(:,:,testnameind) = mean(all_steric_slc(:,:,rcmcount-3:rcmcount),3);
                all_mano_slc(:,:,testnameind) = mean(all_mano_slc(:,:,rcmcount-3:rcmcount),3);
                all_mean_slc(testnameind,1) = mean(all_mean_slc(rcmcount-3:rcmcount));
                all_mean_steric_slc(testnameind,1) = mean(all_mean_steric_slc(rcmcount-3:rcmcount));
                all_mean_mano_slc(testnameind,1) = mean(all_mean_mano_slc(rcmcount-3:rcmcount));
            end
           
            mean_data = all_slc(:,:,testnameind);
            
            ysecs_slc=all_slc(:,:,testnameind).*ysecsmask;
            es_slc=all_slc(:,:,testnameind).*esmask;
            [mean_all_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(mean_data, lon_glo, lat_glo);
            [mean_ysecs_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(ysecs_slc, lon_glo, lat_glo);
            [mean_es_slc(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(es_slc, lon_glo, lat_glo);
            ysecs_std_all_slc(testnameind,1) = std(ysecs_slc(:),'omitnan');
            es_std_all_slc(testnameind,1) = std(es_slc(:),'omitnan');
            es_max_all_slc(testnameind,1) = max(es_slc(:));
            es_min_all_slc(testnameind,1) = min(es_slc(:));
            
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[3+12*(testnameind-6) 4+12*(testnameind-6) 9+12*(testnameind-6) 10+12*(testnameind-6)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_glo', lat_glo', land_model_interped','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_glo', lat_glo', mean_data','parent',ax{testnameind,2});
            colormap(ax{testnameind,2},flip(autumn));
            caxis([30 70]);
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

        end

        h = colorbar('westoutside','AxisLocation','out');
        caxis([30 70]);
        
        set(h,'fontsize',m_grid_fontsize+7);
        h_title=title(h,'(cm)','fontsize',m_grid_fontsize+5);
        
%         set(h, 'Position', [pos_sb{10}(1)+0.332, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.09])  % right, up, width, height
        set(h, 'Position', [pos_sb{10}(1)-0.092, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.33])  % right, up, width, height

        
        
%      ---------- diff
        for testnameind=11:15
%         testnameind=5;
            testname=all_testname2{testnameind-10};    
            
            rcmtestname = all_testname2{testnameind-10};  
           
            mean_data = all_slc(:,:,testnameind-5) - all_slc(:,:,testnameind-10);
            ysecs_diff=mean_data.*ysecsmask;
            es_diff=mean_data.*esmask;
            [mean_ysecs_diff(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(ysecs_diff, lon_glo, lat_glo);
            [mean_es_diff(testnameind,1), error_status] = Func_0011_get_area_weighted_mean(es_diff, lon_glo, lat_glo);
            
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,6,[5+12*(testnameind-11) 6+12*(testnameind-11) 11+12*(testnameind-11) 12+12*(testnameind-11)]);  % Auto-fitted to the figure.
            pos_sb{testnameind}=get(sb{testnameind}, 'pos'); % Get the position.
            pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig + correction_right_fig*2;
%             pos_sb{testnameind} = pos_sb{testnameind} + correction_large_fig;
            delete(sb{testnameind}); % Delete the subplot axes
            ax{testnameind,1}=axes;
            set(ax{testnameind,1},'pos',pos_sb{testnameind});
            pc{testnameind,1}=m_pcolor(lon_glo', lat_glo', land_model_interped','parent',ax{testnameind,1});
            colormap(ax{testnameind,1},[0.8 0.8 0.8]);
            shading(gca,m_pcolor_shading_method); 
            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type, 'box', m_grid_box_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120, 130, 140], 'YaxisLocation', 'right', ...
                    'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'parent', ax{testnameind,1});  
            hold on

            ax{testnameind,2}=axes;
            set(ax{testnameind,2},'pos',pos_sb{testnameind});
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            pc{testnameind,2}=m_pcolor(lon_glo', lat_glo', mean_data','parent',ax{testnameind,2});
            [byrmap, er_status] = Func_0009_get_colormaps('byr');
            colormap(ax{testnameind,2},byrmap);
            caxis([-20 20]);
            shading(gca,m_pcolor_shading_method);   
            hold on

            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

        end

        h2 = colorbar('eastoutside','AxisLocation','out');
        caxis([-20 20]);
        
        set(h2,'fontsize',m_grid_fontsize+7);
        h_title=title(h2,'(cm)','fontsize',m_grid_fontsize+5);
        
        set(h2, 'Position', [pos_sb{15}(1)+0.202, pos_sb{15}(2)+0, 0.0231, pos_sb{15}(3)*2+0.33])  % right, up, width, height
%         set(h, 'Position', [pos_sb{10}(1)-0.132, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.33])  % right, up, width, height

      
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width*2.3 paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
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