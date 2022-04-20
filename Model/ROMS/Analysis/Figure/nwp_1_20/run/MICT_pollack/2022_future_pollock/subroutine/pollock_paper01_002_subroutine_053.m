 if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
        egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf curl_mask ...
        comb_mov_dist_lon_mean comb_mov_dist_lat_mean
    for yearij = 1:length(allyear)
        tempyear = allyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
            
            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
            mov_dist_lon_mean = ncread(ncname, 'mov_dist_lon_mean', [1 ind_checktime], [inf 1]);
            mov_dist_lat_mean = ncread(ncname, 'mov_dist_lat_mean', [1 ind_checktime], [inf 1]);
            lastday_m=length(mov_dist_lon_mean);
            if (exist('comb_mov_dist_lon_mean')==0)
                comb_ocean_time=ocean_time;
                comb_mov_dist_lon_mean= mov_dist_lon_mean;
                comb_mov_dist_lat_mean= mov_dist_lat_mean;
            else
                comb_ocean_time(end+1:end+lastday_m)=ocean_time;
                comb_mov_dist_lon_mean(end+1:end+lastday_m)=mov_dist_lon_mean;
                comb_mov_dist_lat_mean(end+1:end+lastday_m)=mov_dist_lat_mean;
            end
        end
    end

    half_len=round(length(comb_mov_dist_lat_mean)/2);
    egg_half_1=mean(comb_mov_dist_lat_mean(1:half_len), 'omitnan');
    egg_half_2=mean(comb_mov_dist_lat_mean(half_len+1:end), 'omitnan');
    regime_ts_egg_mask(1:half_len)=egg_half_1;
    regime_ts_egg_mask(half_len+1:length(comb_mov_dist_lat_mean))=egg_half_2;

    axLH = gca;
    mslplot{2}=bar(comb_mov_dist_lat_mean, 'FaceColor', [0.6 0.6 0.6], 'parent',axLH);
    hold on
    mslplot{1}=plot(regime_ts_egg_mask,'k','parent',axLH);
    ylabel(axLH,'degrees (lat)')
    ax_pos = get(axLH,'position');
    set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
    cal_arr=1;
    for cali=1:length(allyear)-1
        cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
    end
    set(axLH,'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);

    axis tight 
    ylim(axLH, [-1,1])
    set(axLH,'ycolor','k', 'box', 'off', 'FontSize',15);
    xlabel(axLH, 'Year');

    set(mslplot{1},'LineWidth',2);
    grid on
    lgd=legend([mslplot{1}], 'moved distance of individuals(lat)');
    set(lgd,'FontSize',15);
    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    set(lgd,'Orientation','horizontal');
    
    title(['latitudinal shift(',num2str(temp_checktime), 'd)',  ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i'), ...
        ', ', num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i')]);
    
    half_len=length(comb_mov_dist_lat_mean)/2;
    egg_half_1=mean(comb_mov_dist_lat_mean(1:half_len), 'omitnan');
    egg_half_2=mean(comb_mov_dist_lat_mean(half_len+1:end), 'omitnan');

    

    set(gcf,'PaperPosition', [0 0 36 12]) 
    saveas(gcf,jpgname,'tif');
    grid off
    hold off
    disp(' ')
    disp([fig_name, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    close all;
end