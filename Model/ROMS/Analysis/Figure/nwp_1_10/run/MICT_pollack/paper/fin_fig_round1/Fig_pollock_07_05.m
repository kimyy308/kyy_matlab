close all; clear all;  clc;
warning off;

RCM_info.testname = 'test06';
RCM_info.regionname = 'pollock_egg3';
% earlyyear=[1983:1987];
% lateyear=[1988:1992];
RCM_info.months = [1,2]; % % put month which you want to plot [month month ...]
% allyear =[1983:1992];
% refyear =[1983:1987];
% caxisval=[0 150];
% caxisval_diff=[-70 70];


addpath(genpath(['C:\Users\User\Dropbox\source\matlab\function\']))
[tmp.dropboxpath, tmp.erorr_status] = Func_0008_set_dropbox_path(computer);
addpath(genpath([tmp.dropboxpath, '\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\']))
[cmap.byrmap3, tmp.error_status] = Func_0009_get_colormaps('byr3', tmp.dropboxpath);
[cmap.byrmap, tmp.error_status] = Func_0009_get_colormaps('byr2', tmp.dropboxpath);  
[RCM_grid.refpolygon, RCM_grid.domain, tmp.error_status] = Func_0007_get_polygon_data_from_regionname(RCM_info.regionname);

tmp.param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param2_kyy_', RCM_info.regionname, '.m'];
dirs.figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\'); % % where figure files will be saved
dirs.filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
dirs.savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
dirs.inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
LTRANS_RCM_info.testname='Pollock6';

figs.tifname = [dirs.figrawdir, 'Fig07_05.tif'];

figs.text_config.text_period.xpos=127.15;
figs.text_config.text_period.ypos=36.9;
figs.text_config.text_days.xpos=127.15;
figs.text_config.text_days.ypos=40.75;
% figs.text_legend.text_period.xpos=127.15;
% figs.text_legend.text_period.ypos=0.98;

% l1 = annotation('line'); l1.Units = 'inches'; l1.X=[1,1]*0.25;	l1.Y=[0,1.5]; l1.LineWidth = 0.25; l1.Color = 'r';
% l2 = annotation('line'); l2.Units = 'inches'; l2.X=[1,1]*0.48;	l2.Y=[0,1.5]; l2.LineWidth = 0.25; l2.Color = 'r';
% l3 = annotation('line'); l3.Units = 'inches'; l3.X=[1,1]*(0.25+2.1633);	l3.Y=[0,1.5]; l3.LineWidth = 0.25; l3.Color = 'r';
% l4 = annotation('line'); l4.Units = 'inches'; l4.Y=[1,1]*ylab_h;	l4.X=[0,4.5]; l4.LineWidth = 0.25; l4.Color = 'r';
% l5 = annotation('line'); l5.Units = 'inches'; l5.Y=[1,1]*(ylab_h+ytick_h);	l5.X=[0,4.5]; l5.LineWidth = 0.25; l5.Color = 'r';
% l6 = annotation('line'); l6.Units = 'inches'; l6.Y=[1,1]*(ylab_h+ytick_h+ax_h+t_margin);	l6.X=[0,4.5]; l6.LineWidth = 0.25; l6.Color = 'r';

% print(figs.fig,figs.fig_title,'-dtiff','-r300');
% print(figs.fig,figs.fig_title,'-depsc');


for mm=[1,3,4,6]
    RCM_info.years{mm}=1983:1987;
end
for mm=[2,5]
    RCM_info.years{mm}=1988:1992;
end
RCM_info.temp_chekctime=[15, 15, 15, 30, 30, 30];
RCM_info.months=[1,2];
RCM_info.panel_legend={'A', 'B', 'C', 'D', 'E', 'F'};

% % % positioning of the figures
figs.l_scale = 1;
figs.l_unit = 'km';

figs.ylimit = [-42.4e3 inf];
figs.zlimit = [-4000 20];

figs.n_rows = 2;
figs.n_columns = 3;
% figs.ax_width = 1.5;
figs.ax_width = 0.96;

figs.ax_height = 1.0250;
figs.xlabel_width = 0.2;
figs.xtick_width= 0.2;
figs.left_margin = 0.5;
figs.right_margin = 0.2;
figs.ylabel_height = 0.24;
figs.ytick_height= 0.25;
figs.top_margin = 0.05;
figs.fig_pos  = [0, ...
                 0, ...
                 figs.xlabel_width+figs.n_columns*(figs.xtick_width+figs.ax_width+figs.right_margin)+1.5, ...
                 figs.ylabel_height+figs.n_rows*(figs.ytick_height+figs.ax_height+figs.top_margin)+0.5];
             
             
figs.fig_title = 'Figure_7';
figs.exp_title = {{'Open-Ocean'},{'Shelf &','Control','Slope'},{'No Shelf'},{'Shelf &','Gentle','Slope'},{'Shelf &','Steep','Slope'},{'Shelf &','Cliff'}};
figs.fig = figure('name',figs.fig_title,'PaperUnits','inches', ...
    'paperposition',figs.fig_pos,'position',figs.fig_pos*get(groot,'ScreenPixelsPerInch'),'visible','off');
% figs.fig = figure('name',figs.fig_title,'PaperUnits','inches', ...
%     'paperposition',figs.fig_pos,'position',figs.fig_pos*get(groot,'ScreenPixelsPerInch'));

for mm = 1:figs.n_rows*figs.n_columns
    ii=ceil(mm/figs.n_columns); % # of column
    jj=mod(mm-1,figs.n_columns)+1; % # of row
    figs.subplot(mm).posAx = ...
        [figs.left_margin+figs.xlabel_width+figs.xtick_width+(jj-1)*(figs.ax_width+figs.right_margin+figs.xtick_width), ...
        figs.ylabel_height+figs.ytick_height+(figs.n_rows-ii)*(figs.ax_height+figs.top_margin+figs.ytick_height), ...
        figs.ax_width, ...
        figs.ax_height];
    figs.subplot(mm).ax = axes(figs.fig, 'units', 'inches', 'position',figs.subplot(mm).posAx);
    figs.xtick= [128 130 132];
    figs.ytick= [36 : 41];
    
    run(tmp.param_script);
    param.m_grid_fontsize = param.m_grid_fontsize -8;
    
    for yearij = 1:length(RCM_info.years{mm})
        tmp.tempyear = RCM_info.years{mm}(yearij);
        for monthij = 1:length(RCM_info.months)
            tmp.tempmonth = RCM_info.months(monthij);
            tmp.ncname = [dirs.savedir,RCM_info.testname,'_',RCM_info.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            RCM_data.all_checktime{mm}=ncread(tmp.ncname, 'checktime');
            tmp.ind_checktime=find(RCM_data.all_checktime{mm}==RCM_info.temp_chekctime(mm));
            RCM_data.egg_mask{mm}=squeeze(ncread(tmp.ncname, 'mask_par', [1 1 1 tmp.ind_checktime], [inf inf inf 1]));
            lastday_m=size(RCM_data.egg_mask{mm},3);
            if (yearij == 1 && monthij == 1)
                RCM_data.comb_egg_mask{mm}=RCM_data.egg_mask{mm};
            else
                RCM_data.comb_egg_mask{mm}(:,:,end+1:end+lastday_m)=RCM_data.egg_mask{mm};
            end
        end
    end
    RCM_grid.lon_RCM{mm} = ncread(tmp.ncname, 'lon_rho');
    RCM_grid.lat_RCM{mm} = ncread(tmp.ncname, 'lat_rho');
    
    RCM_data.temp_surf{mm}=ncread(tmp.ncname, 'temp_surf', [1 1 1], [inf inf 1]);
    RCM_data.temp_surf{mm}(isnan(RCM_data.temp_surf{mm}))=50000;
    RCM_data.temp_surf{mm}(RCM_data.temp_surf{mm}<50000)=NaN;
    RCM_grid.land{mm}=RCM_data.temp_surf{mm};
    RCM_grid.land{mm}(RCM_data.temp_surf{mm}==50000)=1;

    RCM_data.sum_egg{mm} = sum(RCM_data.comb_egg_mask{mm},3);

    RCM_data.sum_egg{mm}(RCM_data.sum_egg{mm}==0)=NaN;
    RCM_data.lon_sum_data{mm}=sum(RCM_data.sum_egg{mm}, 1, 'omitnan');
    RCM_grid.lat_1d{mm}=mean(RCM_grid.lat_RCM{mm},1,'omitnan');
    RCM_data.lat_sum_data{mm}=sum(RCM_data.sum_egg{mm}, 2, 'omitnan');
    RCM_grid.lon_1d{mm}=mean(RCM_grid.lon_RCM{mm},2,'omitnan');

    if mod(mm,figs.n_columns)==0
        param.colorbar_lev{mm}= [-70 70];
        param.colorbar=cmap.byrmap3;
        RCM_data.sum_egg{mm}= RCM_data.sum_egg{mm-1}-RCM_data.sum_egg{mm-2};
        RCM_data.lon_sum_data{mm}= RCM_data.lon_sum_data{mm-1}-RCM_data.lon_sum_data{mm-2};
        RCM_data.lat_sum_data{mm}= RCM_data.lat_sum_data{mm-1}-RCM_data.lat_sum_data{mm-2};
    else
        param.colorbar_lev{mm}= [0 150];
        param.colorbar=parula;
    end
    
% % %  set axis position
    
% %     xlim(figs.subplot(mm).ax,figs.ylimit/figs.l_scale); 
% %     ylim(figs.subplot(mm).ax,figs.zlimit/figs.l_scale);
% %     figs.subplot(mm).text_period=text(figs.subplot(mm).ax, ...
% %         figs.text_config.text_period.xpos,figs.text_config.text_period.ypos, ...
% %         figs.exp_title{mm},'verticalalignment','top','horizontalalignment','left','units','normalized','interpreter','tex');
% %     figs.subplot(mm).text_legend=text(figs.subplot(mm).ax, ...
% %         figs.text_legend.text_period.xpos,figs.text_legend.text_period.ypos, ...
% %         ['\bf ',char(double('A')-1+mm)], ...
% %         'verticalalignment','bottom','horizontalalignment','left','units','normalized','fontsize',12);
%     switch ii, case figs.n_rows, xlabel(figs.subplot(mm).ax,['$$y$$ (',figs.l_unit,')'],'interpreter','latex'); end
%     switch jj, case 1,      ylabel(figs.subplot(mm).ax,['z (',figs.l_unit,')']); end
    
% % %     pcolor for lands
    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    figs.subplot(mm).pcol = m_pcolor(RCM_grid.lon_RCM{mm}',RCM_grid.lat_RCM{mm}', RCM_grid.land{mm}','parent',figs.subplot(mm).ax);
    colormap(figs.subplot(mm).ax,[0.8 0.8 0.8]);
    shading(gca, param.m_pcolor_shading_method);

    m_grid_apply(mm, figs, param, 1)
    
    % % % pcolor
    hold on
    figs.subplot(mm).ax2 = axes(figs.fig, 'units', 'inches', 'position',figs.subplot(mm).posAx);

    m_proj(param.m_proj_name,'lon',[RCM_grid.domain(1) RCM_grid.domain(2)],'lat',[RCM_grid.domain(3) RCM_grid.domain(4)]);
    figs.subplot(mm).pcol2=m_pcolor(RCM_grid.lon_RCM{mm}',RCM_grid.lat_RCM{mm}', ...
        RCM_data.sum_egg{mm}','parent',figs.subplot(mm).ax2);
    colormap(figs.subplot(mm).ax2,param.colorbar);
    caxis(param.colorbar_lev{mm});
    shading(gca,param.m_pcolor_shading_method);   
    m_grid_apply(mm, figs, param, 2)
    
    if mod(mm,figs.n_columns)~=0
        figs.txt{mm}=m_text(figs.text_config.text_period.xpos, figs.text_config.text_period.ypos, ...
            [num2str(min(RCM_info.years{mm})), '¡ª', num2str(max(RCM_info.years{mm}))], ...
            'FontSize', param.m_grid_fontsize-2); 
    elseif mm==figs.n_columns
        figs.txt{mm}=m_text(figs.text_config.text_period.xpos, figs.text_config.text_period.ypos, ...
            ['(B)', '-', '(A)'], ...
            'FontSize', param.m_grid_fontsize-2);
    elseif mm==figs.n_columns*2
        figs.txt{mm}=m_text(figs.text_config.text_period.xpos, figs.text_config.text_period.ypos, ...
            ['(E)', '-', '(D)'], ...
            'FontSize', param.m_grid_fontsize-2); 
    end
    figs.txt{mm}=m_text(figs.text_config.text_days.xpos, figs.text_config.text_days.ypos, ...
        [num2str(RCM_info.temp_chekctime(mm)) ,' days'], ...
        'FontSize', param.m_grid_fontsize-2); 

end

% % % barh

figs.ax_width_bar_mer = figs.ax_width/4;
figs.ax_height_bar_mer = figs.ax_height;
figs.ax_width_bar_zonal = figs.ax_width;
figs.ax_height_bar_zonal = figs.ax_height/4;
figs.ax_bar_margin = 0.03;

for mm = 1:figs.n_rows*figs.n_columns
    ii=ceil(mm/figs.n_columns); % # of column
    jj=mod(mm-1,figs.n_columns)+1; % # of row
    figs.subplot_bar_mer(mm).posAx = ...
        [figs.left_margin+figs.ax_bar_margin+ figs.ax_width+figs.xlabel_width+figs.xtick_width + (jj-1)*(figs.ax_width+figs.right_margin+figs.xtick_width), ...
        figs.ylabel_height+figs.ytick_height + (figs.n_rows-ii)*(figs.ax_height+figs.top_margin+figs.ytick_height), ...
        figs.ax_width_bar_mer, ...
        figs.ax_height_bar_mer];
    figs.subplot_bar_mer(mm).ax = axes(figs.fig, 'units', 'inches', 'position',figs.subplot_bar_mer(mm).posAx);
    
    latind=find(RCM_grid.lat_1d{mm}>=RCM_grid.domain(3) & RCM_grid.lat_1d{mm}<=RCM_grid.domain(4));
    figs.subplot_bar_mer(mm).bar = barh(RCM_grid.lat_1d{mm}(latind), RCM_data.lon_sum_data{mm}(latind), 'parent', figs.subplot_bar_mer(mm).ax);
    
    figs.subplot_bar_mer(mm).ax.FontSize = param.m_grid_fontsize-4;
    figs.subplot_bar_mer(mm).ax.YTickLabel = [];
    if mod(mm,figs.n_columns)~=0
        figs.subplot_bar_mer(mm).ax.XLim = [0 1500];
        figs.subplot_bar_mer(mm).ax.XTick = [500 1000];
        figs.subplot_bar_mer(mm).ax.XTickLabel = [500, 1000];
    else
        figs.subplot_bar_mer(mm).bar.FaceColor='r';
    end
    figs.subplot_bar_mer(mm).ax.XTickLabelRotation = 45;
    figs.subplot_bar_mer(mm).ax.XMinorTick = 1;
    figs.subplot_bar_mer(mm).ax.YMinorTick = 1;

    grid(figs.subplot_bar_mer(mm).ax,'on')
%     grid(figs.subplot_bar_mer(mm).ax,'minor')
    
%     figs.txt_label{mm}=text(figs.xlabel_width+figs.xtick_width + (jj-1)*(figs.ax_width+figs.right_margin+figs.xtick_width), ...
%         figs.ax_bar_margin + figs.ax_height+figs.ylabel_height+figs.ytick_height + (figs.n_rows-ii)*(figs.ax_height+figs.top_margin+figs.ytick_height), ...
%         RCM_info.panel_legend{mm}, 'fontsize', param.m_grid_fontsize, 'units', 'normalized');    
    figs.txt_label{mm}=text(figs.subplot(mm).ax2, figs.ax_bar_margin, ...
        figs.ax_bar_margin + figs.ax_height, ...
        RCM_info.panel_legend{mm}, 'fontsize', param.m_grid_fontsize, 'units', 'inches', 'verticalalignment', 'bottom');
end

% l1 = annotation('line'); l1.Units = 'inches'; l1.X=[3.4,4.9]; 
% l1.Y=[0.49,0.49]; l1.LineWidth = 1; l1.Color = 'r'; %% 4.9-3.4=1.5 -> xwidth of pcolor

% l1 = annotation('line'); l1.Units = 'inches'; l1.X=[3.67,4.63]-2*((figs.ax_width+figs.right_margin+figs.xtick_width));	
% l1.Y=[0.49,0.49]; l1.LineWidth = 1; l1.Color = 'r'; %% 4.63-3.67=0.96 -> xwidth of pcolor

% l1 = annotation('line'); l1.Units = 'inches'; l1.X=[0,0.96]+figs.xlabel_width+figs.xtick_width;	
% l1.Y=[0.49,0.49]; l1.LineWidth = 1; l1.Color = 'r'; %% 4.63-3.67=0.96 -> xwidth of pcolor
% 
% l2 = annotation('line'); l2.Units = 'inches'; l2.X=[0 0]+figs.xlabel_width+figs.xtick_width+figs.ax_width;	
% l2.Y=[0.49,1.51]; l2.LineWidth = 1; l2.Color = 'r'; %% 4.63-3.67=0.96 -> xwidth of pcolor

% % colorbar (# of particles)
figs.cbar_ax = axes;
figs.cbar_ax.Visible = 'off';
figs.cbar_ax.CLim=param.colorbar_lev{1};
figs.cbar_NTI = colorbar(figs.cbar_ax);
figs.cbar_NTI.Colormap = parula(256);
figs.cbar_NTI.Units = 'inches';
set(figs.cbar_NTI, 'Position', [figs.left_margin-0.05 ...
    figs.ylabel_height+figs.ytick_height, ...
    figs.ax_width/6, ...
    2*figs.ax_height+figs.ytick_height/2])  % right, up, width, height

figs.ctitle_NTI=title(figs.cbar_NTI,'NTI','fontsize',param.m_grid_fontsize+2);

% % % colorbar  (difference)
% h = colorbar('eastoutside','AxisLocation','out');
% caxis(caxisval_diff);
% set(h,'fontsize',m_grid_fontsize+2);
% set(h, 'Position', [pos_sb{3}(1)+0.242, pos_sb{3}(2)+0.002, 0.0231, pos_sb{2}(3)+0.080])  % right, up, width, height
% h_title=title(h,'NTI','fontsize',m_grid_fontsize+2);

figs.cbar_ax2 = axes;
figs.cbar_ax2.Visible = 'off';
figs.cbar_ax2.CLim=param.colorbar_lev{end};
figs.cbar_NTI2 = colorbar(figs.cbar_ax2);
figs.cbar_NTI2.Colormap = cmap.byrmap3;
figs.cbar_NTI2.Units = 'inches';
set(figs.cbar_NTI2, 'Position', [(figs.n_columns)*(figs.ax_width+figs.right_margin+figs.xtick_width+figs.xlabel_width)+figs.ax_width_bar_mer ...
    figs.ylabel_height+figs.ytick_height, ...
    figs.ax_width/6, ...
    2*figs.ax_height+figs.ytick_height/2])  % right, up, width, height

figs.ctitle_NTI2=title(figs.cbar_NTI2,'NTI','fontsize',param.m_grid_fontsize+2);

% (figs.n_columns)*(figs.ax_width+figs.right_margin+figs.xtick_width)
% figs.n_columns
% figs.left_margin+figs.ax_bar_margin+ figs.ax_width+figs.xlabel_width+figs.xtick_width + (jj-1)*(figs.ax_width+figs.right_margin+figs.xtick_width)


print('-dtiff','-r500',figs.tifname); 
RemoveWhiteSpace([], 'file', figs.tifname);
hold off;
close all;

% [mean_early_egg, error_status] = Func_0011_get_area_weighted_mean(mean_early_egg, lon_rho, lat_rho);




% % % barh
% sb1=subplot(6,12,[5 17]);  % Auto-fitted to the figure.
% pos_sb{RCM_info.testnameind,3}=get(sb1, 'pos'); % Get the position.
% pos_sb{RCM_info.testnameind,3} = pos_sb{RCM_info.testnameind,3} + correction_upper_fig+ correction_large_fig + [-0.0800,0.2,0,-0.2];
% delete(sb1); % Delete the subplot axes
% ax{RCM_info.testnameind,3}=axes;
% set(ax{RCM_info.testnameind,3},'pos',pos_sb{RCM_info.testnameind,3});
% 
% barh(lat_1d, lon_sum_early_data_15d)
% %     ylabel('latitude (degree)')
% xlim([0 1500])
% ylim([36 41]); grid minor;
% set(ax{RCM_info.testnameind,1}, 'fontsize', m_grid_fontsize) 
% txt{RCM_info.testnameind,1}=text(127.15, 36.5, '1983¡ª1987', 'FontSize', m_grid_fontsize-2); 
% txt{RCM_info.testnameind,1}=text(127.15, 37.1, '15 days', 'FontSize', m_grid_fontsize-2); 
% set(ax{RCM_info.testnameind,1}, 'XAxisLocation', 'top')
% set(ax{RCM_info.testnameind,1}, 'XTick',[0 500 1000])
% set(ax{RCM_info.testnameind,1}, 'XTicklabel', {})
% set(ax{RCM_info.testnameind,1}, 'YTick',[36 38 40])
% set(ax{RCM_info.testnameind,1}, 'YTicklabel', {'36^oN', '38^oN', '40^oN'})




% % % colorbar (# of particles)
% h = colorbar('westoutside','AxisLocation','out');
% caxis([caxisval]);
% set(h,'fontsize',m_grid_fontsize+2);
% set(h, 'Position', [pos_sb{2}(1)-0.29, pos_sb{2}(2)+0.002, 0.0281, pos_sb{2}(3)+0.080])  % right, up, width, height
% h_title=title(h,'NTI','fontsize',m_grid_fontsize+2);
% 
% % % colorbar  (difference)
% h = colorbar('eastoutside','AxisLocation','out');
% caxis(caxisval_diff);
% set(h,'fontsize',m_grid_fontsize+2);
% set(h, 'Position', [pos_sb{3}(1)+0.242, pos_sb{3}(2)+0.002, 0.0231, pos_sb{2}(3)+0.080])  % right, up, width, height
% h_title=title(h,'NTI','fontsize',m_grid_fontsize+2);
% 
% % % figure save
% set(gcf, 'PaperUnits', 'points');
% set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y/3]);
% set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height*2]) 
% print('-dtiff','-r500',tifname); 
% RemoveWhiteSpace([], 'file', tifname);
% hold off;
% close all;
% 
% save ('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\data_dryad\Data05_numpar.mat', ...
%             'RCM_grid.lon_RCM', 'RCM_grid.lat_RCM', 'RCM_land', ...
%             'sum_early_egg_15d', 'sum_late_egg_15d', 'sum_diff_egg_15d', ...
%             'sum_early_egg_30d', 'sum_late_egg_30d', 'sum_diff_egg_30d')

        
function m_grid_apply(mm, figs, param, axnum)
    if axnum==1 %% for land pcolor
        if (mod(mm,figs.n_columns)==1) %% first column
            if mm==1 %% upper left
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
            elseif mm==(figs.n_rows-1)*figs.n_columns+1  %% lower left
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
            else
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
            end
        elseif (mm>1 && mm<=figs.n_columns) %% first row (except upper left)
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
        elseif (mm>(figs.n_rows-1)*figs.n_columns+1) %% last row (except lower left)
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax);
        end
    elseif axnum==2   %% for NTI pcolor
        if (mod(mm,figs.n_columns)==1) %% first column
            if mm==1 %% upper left
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
            elseif mm==(figs.n_rows-1)*figs.n_columns+1  %% lower left
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
            else
                m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', figs.ytick, 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
            end
        elseif (mm>1 && mm<=figs.n_columns) %% first row (except upper left)
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', [], 'xtick', figs.xtick, 'XaxisLocation', 'top',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
        elseif (mm>(figs.n_rows-1)*figs.n_columns+1) %% last row (except lower left)
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
            m_grid('fontsize', param.m_grid_fontsize, 'tickdir', param.m_grid_tickdir_type, 'box', param.m_grid_box_type,  ...
                    'xticklabels', figs.xtick, 'xtick', figs.xtick, 'XaxisLocation', 'bottom',...
                    'yticklabels', [], 'ytick', figs.ytick, 'backcolor', 'none',  ...
                    'parent', figs.subplot(mm).ax2);
        end
    end
end