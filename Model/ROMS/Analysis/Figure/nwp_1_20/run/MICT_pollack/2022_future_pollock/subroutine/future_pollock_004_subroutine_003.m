% %  Updated 26-Apr-2021 by Yong-Yub Kim, 

% start-------------------- earlier decadal SST, SSS plot
for varind2=1:length(RCM_info.vars)
%% set variable name & figure directory
    tmp.variable=RCM_info.vars{varind2};
    dirs.figdir=[dirs.figrawdir,'ts', tmp.fs, tmp.regionname, tmp.fs, tmp.variable, tmp.fs, ...
        num2str(min(RCM_info.years)), '_', num2str(max(RCM_info.years)), tmp.fs];
    if (exist(strcat(dirs.figdir) , 'dir') ~= 7)
        mkdir(strcat(dirs.figdir));
    end 

%% set figure file name
    tmp.tifname=strcat(dirs.figdir, tmp.testname, '_yearly_ts_', tmp.variable, '_', ...
        num2str(min(RCM_info.years),'%04i'), '_',num2str(max(RCM_info.years),'%04i'), ...
        '_', RCM_info.season,'.tif'); %% ~_year_month.jpg
    if (exist(tmp.tifname , 'file') ~= 2 || fig_flag==2)      
        run(tmp.param_script);
%% set data file name (mat) for reading
        tmp.matname = [dirs.matdir, tmp.testname, '_', tmp.regionname, '_', tmp.variable,...
            '_mean_', num2str(min(RCM_info.years),'%04i'), '-', num2str(max(RCM_info.years),'%04i'),...
        '_', RCM_info.season, '.mat'];
        if (exist(tmp.matname , 'file') ~= 2)
            disp('please get data from subroutine_004_001 first')
        else
            load(tmp.matname);
        end 
        
        axLH = gca;
        figure.mslplot=plot(RCM_info.years, RCM_data.spa_yearly_mean,'k','parent',axLH);
%                 mslplot{2}=plot(sum_ts_egg_mask, 'Color', 'b', 'LineStyle', '-', 'parent',axLH);
%                 mslplot{2}=plot(sum_ts_egg_mask, 'color', [0.8 0.8 0.8], 'parent',axLH);

        ylabel(axLH,'^oC')
        figure.ax_pos = get(axLH,'position');
%         set(figure.axLH,'yaxislocation','left','position',figure.ax_pos+[0 0.02 -0.01 -0.02]);
        
% % %         all months
%         tmp.cal_arr=1; %calendar
%         for cali=1:length(RCM_info.years)-1
%             tmp.cal_arr=[tmp.cal_arr, tmp.cal_arr(end)+sum(eomday(allyear(cali),RCM_info.months))];
%         end
% % %         yearly
        tmp.cal_arr=1; %calendar
        tmp.datenum=datenum(RCM_info.years(1),06,30);
        for cali=1:length(RCM_info.years)-1
            tmp.cal_arr=[tmp.cal_arr, tmp.cal_arr(end)+sum(eomday(RCM_info.years(cali),6))];
            tmp.datenum=[tmp.datenum, datenum(RCM_info.years(cali+1),06,30)];
        end
        
%         set(figure.axLH,'color','none','yaxislocation','left','xtick', tmp.cal_arr,'xticklabel',datestr(comb_ocean_time(tmp.cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
%         set(figure.axLH,'color','none','yaxislocation','left','xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy') ,'position', figure.ax_pos+[0 0.02 -0.01 -0.02]);
%         set(figure.axLH,'yaxislocation','left' ,'position', figure.ax_pos+[0 0.02 -0.01 -0.02]);
%         set(figure.axLH,'yaxislocation','left','xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy'));
%         set(figure.axLH, 'xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy'));
%         set(axLH, 'xaxislocation','bottom', 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy'));
%         set(gca, 'xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy'));
        set(axLH, 'xaxislocation','bottom', 'XTickMode', 'auto', 'XTickLabelMode', 'auto', 'xtick', RCM_info.years,'xticklabel', datestr(tmp.datenum, 'yyyy'));

        axis tight 
        ylim([5 18])
%         ylim([6 11])
%         ylim([4.5 11])

        set(axLH,'xcolor','k', 'box', 'off', 'FontSize',15);
        xlabel(axLH, 'Year');

%         title(['Spawned eggs', ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i'), ...
%             ', ', num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i')]);
        if min(RCM_info.years) == max(RCM_info.years)
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(max(RCM_info.years),'%04i'),') ');                        
        else
            tmp.titlename = strcat(tmp.variable, ', ', RCM_info.season(1:3), ', ', tmp.abb,',(',num2str(min(RCM_info.years),'%04i'),'-',num2str(max(RCM_info.years),'%04i'),') ');
        end
        title(tmp.titlename,'fontsize',param.m_pcolor_title_fontsize);  %%title
        set(figure.mslplot,'LineWidth',2);
        grid on

        figure.lgd=legend([figure.mslplot], 'temp');

%         half_len=length(sum_ts_egg_mask)/2;
%         egg_half_1=mean(sum_ts_egg_mask(1:half_len));
%         egg_half_2=mean(sum_ts_egg_mask(half_len+1:end));

        set(figure.lgd,'FontSize',15);
        set(figure.lgd,'Position',[0.13 0.88, 0.775, 0.03]);
        set(figure.lgd,'Orientation','horizontal');

        set(gcf,'PaperPosition', [0 0 36 12]) 
        
       
        saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
        close all;
        clear RCM_data.mean
        RCM_grid=rmfield(RCM_grid, 'lon_rho');
    end
end