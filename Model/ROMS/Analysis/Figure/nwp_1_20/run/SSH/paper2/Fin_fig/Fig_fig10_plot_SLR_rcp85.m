close all; clear all;  clc;
warning off;

all_testname1  = {'IPSL-CM5A-LR', 'IPSL-CM5A-MR', 'NorESM1-M', 'MPI-ESM-LR', 'ens10'};
all_testname2  = {'test65', 'test66', 'test67', 'test68', 'ens10'};

all_region2 ={'AKP4'};

variable = 'SSH';
scenname='rcp85';
for regionind2=1:length(all_region2)
        close all;

        inputyear = [2006:2100]; % % put year which you want to plot [year year ...]
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

        correction_right_fig=[-0.1600,0,0,0]; % right, up, width, height
        correction_upper_fig=[0,0,0,0];
        correction_large_fig=[0,0,0,0.020];
        hold on;
        tifname=strcat(figrawdir, 'fig10','.tif'); %% ~_year_month.jpg
        test_text={'(A) GCM-IPSL-L', '(B) GCM-IPSL-M', '(C) GCM-Nor', '(D) GCM-MPI', '(E) GCM-Ens','(F) RCM-IPSL-L', '(G) RCM-IPSL-M', '(H) RCM-Nor', '(I) RCM-MPI', '(J) RCM-Ens'};
%      ---------- GCM
        for testnameind=1:5
            testname=all_testname1{testnameind};    

            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat('D:\Data\Model\CMIP5\zos\', scenname, '\interp\', testname, '\'); % % where data files are
            interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');
            run(param_script);
            m_grid_fontsize=m_grid_fontsize-12;
            load([cmip5dir,regionname, '_', testname, '_cmip5_interped_land','.mat']);
            sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');
            cut_lon_rho=ncread(interpedfilename,'lon_cmems');
            cut_lat_rho=ncread(interpedfilename,'lat_cmems');
            mean_data=sl(:,:,end)-sl(:,:,1);

            mask_model = double(inpolygon(cut_lon_rho,cut_lat_rho,refpolygon(:,1),refpolygon(:,2)));
            mask_model(mask_model==0)=NaN;

            mean_data= mean_data .* mask_model;        
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,4,[1+8*(testnameind-1) 2+8*(testnameind-1) 5+8*(testnameind-1) 6+8*(testnameind-1)]);  % Auto-fitted to the figure.
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
            colormap(ax{testnameind,2},flip(autumn));
            caxis([60 100]);
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

        for testnameind=6:10

%      ---------- RCM-IPSL-L
%         testnameind=5;
            testname=all_testname2{testnameind-5};    

            param_script ='C:\Users\user\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m';
            filedir = strcat('D', ':\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            landdir = strcat('G', ':\Data\Model\ROMS\nwp_1_20\', 'test65', '\run\'); % % where data files are

            run(param_script);
            m_grid_fontsize=m_grid_fontsize-12;
            interpedfilename = strcat(filedir, testname,'_',regionname, 'cmems_interped_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

            load([landdir,regionname, '_', 'test65', '_model_interped_land','.mat']);
            sl=ncread(interpedfilename,'interped_sla_yearly_exp_fit');
            cut_lon_rho=ncread(interpedfilename,'lon_cmems');
            cut_lat_rho=ncread(interpedfilename,'lat_cmems');
            mean_data=sl(:,:,end)-sl(:,:,1);
            
            mean_data= mean_data .* mask_model;        
            m_proj(m_proj_name,'lon',[lonlat(1) lonlat(2)],'lat',[lonlat(3) lonlat(4)]);
            sb{testnameind}=subplot(10,4,[3+8*(testnameind-6) 4+8*(testnameind-6) 7+8*(testnameind-6) 8+8*(testnameind-6)]);  % Auto-fitted to the figure.
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
            colormap(ax{testnameind,2},flip(autumn));
            caxis([60 100]);
            shading(gca,m_pcolor_shading_method);   
            hold on
            [C{testnameind,2},h{testnameind,2}]=m_contour(cut_lon_rho',cut_lat_rho', mean_data', [10, 15, 20], 'color','k', ...
                        'linewidth', 1, 'linestyle', '-', 'parent', ax{testnameind,2});
                    clabel(C{testnameind,2},h{testnameind,2},'FontSize',m_grid_fontsize-2,'Color','k', ...
                        'labelspacing', 50000,'Rotation', 0,'fontweight', 'bold');

            m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
                    'xticklabels', [120, 130, 140], 'xtick',[120,  130, 140], ...
                    'yticklabels', [], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
            txt{testnameind,2}=m_text(118, 49, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 

        end

        h = colorbar('eastoutside','AxisLocation','out');
        caxis([60 100]);
        
        set(h,'fontsize',m_grid_fontsize+7);
        h_title=title(h,'(cm)','fontsize',m_grid_fontsize+5);
        
        set(h, 'Position', [pos_sb{10}(1)+0.332, pos_sb{10}(2)+0, 0.0231, pos_sb{10}(3)*2+0.09])  % right, up, width, height

        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 

%         saveas(gcf,tifname,'tif');
        print('-dtiff','-r500',tifname)

        hold off;
        close all;
        
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