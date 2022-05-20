 % %  Updated    03-May-2022 by Yong-Yub Kim   % structure

 %% historical plot
    tmp.tifname=strcat(dirs.regime_figdir, tmp.testname_his,'_',tmp.regionname, '_regime_ts_sp_ground_', ...
                num2str(min(RCM_info.years_his),'%04i'),'_',num2str(max(RCM_info.years_his),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 
            
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
    run(tmp.param_script);
    
    %%    initialization
    if(isfield(tmp, 'tlen_his')~=1)
        tmp.tlen_his=0;
        for yearij = 1:length(RCM_info.years_his)
            tmp.tempyear = RCM_info.years_his(yearij);
            for monthij = 1:length(RCM_info.months)
                tmp.tempmonth = RCM_info.months(monthij);
                tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                tmp.tlen_his=tmp.tlen_his + length(ncread(tmp.ncname, 'time'));
            end
        end
        tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
        tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
    end
    tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his);
    tmp.comb_u_rho=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his);
    tmp.comb_v_rho=NaN(tmp.xlen,tmp.ylen,tmp.tlen_his);
    tmp.comb_ocean_time=NaN(tmp.tlen_his);
% %     
% %     
% %     clear wind_curl comb_wind_curl u_rho tmp.comb_u_rho v_rho tmp.comb_v_rho uwind comb_uwind vwind comb_vwind ...
% %         egg_mask tmp.comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 tmp.temp_surf comb_tmp.temp_surf  curl_mask


    for yearij = 1:length(RCM_info.years_his)
        tmp.tempyear = RCM_info.years_his(yearij);
        for monthij = 1:length(RCM_info.months)
            tmp.tempmonth = RCM_info.months(monthij);
            tmp.ncname = [dirs.savedir_his,tmp.testname_his,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
            
            if yearij==1 && monthij==1
                RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
                RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');
                RCM_grid.mean_diff_lon = mean(diff(RCM_grid.lon_rho(:,1)));
                RCM_grid.mean_diff_lat = mean(diff(RCM_grid.lat_rho(1,:)));
%                 RCM_grid.dA = (m_lldist([RCM_grid.mean_diff_lon 0], [0 0])*1e3)^2 * cosd(RCM_grid.lat_rho); %m^2
                RCM_grid.dA = (m_lldist([RCM_grid.mean_diff_lon 0], [0 0]))^2 * cosd(RCM_grid.lat_rho); %km^2
%                 tmp.refpolygon=[0, 0; 0, 39; 360, 39; 360,0];
                tmp.vel_mask = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
            end
            ocean_time=ncread(tmp.ncname, 'time')+datenum(1900,12,31);
            tmp.u_rho=ncread(tmp.ncname, 'u_rho').*tmp.vel_mask;
            tmp.u_rho(tmp.u_rho==0)=NaN;
            tmp.v_rho=ncread(tmp.ncname, 'v_rho').*tmp.vel_mask;
            tmp.v_rho(tmp.v_rho==0)=NaN;
            tmp.egg_mask=ncread(tmp.ncname, 'egg_mask').*tmp.vel_mask;
            tmp.egg_mask(tmp.egg_mask==0)=NaN;
            tmp.egg_mask=tmp.egg_mask.*RCM_grid.dA;
            
            tmp.lastday_m=size(tmp.egg_mask,3);
            if yearij==1 && monthij==1
                tmp.comb_u_rho=tmp.u_rho;
                tmp.comb_v_rho=tmp.v_rho;
                tmp.comb_egg_mask=tmp.egg_mask;
                tmp.comb_ocean_time=ocean_time;
                tmp.endij=tmp.lastday_m;
            else
                tmp.comb_u_rho(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.u_rho;
                tmp.comb_v_rho(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.v_rho;
                tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.egg_mask;
                tmp.comb_ocean_time(tmp.endij+1:tmp.endij+tmp.lastday_m)=ocean_time;
                tmp.endij=tmp.endij+tmp.lastday_m;
            end
        end
    end

    tmp.ts_u_rho=reshape(tmp.comb_u_rho,[size(tmp.comb_u_rho,1)*size(tmp.comb_u_rho,2), size(tmp.comb_u_rho,3)]);
    tmp.mean_ts_u_rho=mean(tmp.ts_u_rho,1,'omitnan');
    tmp.ts_v_rho=reshape(tmp.comb_v_rho,[size(tmp.comb_v_rho,1)*size(tmp.comb_v_rho,2), size(tmp.comb_v_rho,3)]);
    tmp.mean_ts_v_rho=mean(tmp.ts_v_rho,1,'omitnan');
    tmp.ts_egg_mask=reshape(tmp.comb_egg_mask,[size(tmp.comb_egg_mask,1)*size(tmp.comb_egg_mask,2), size(tmp.comb_egg_mask,3)]);
    tmp.sum_ts_egg_mask=sum(tmp.ts_egg_mask,1,'omitnan');

    half_len=round(length(tmp.sum_ts_egg_mask)/2);
    egg_half_1=mean(tmp.sum_ts_egg_mask(1:half_len));
    egg_half_2=mean(tmp.sum_ts_egg_mask(half_len+1:end));
%     regime_ts_egg_mask(1:half_len)=egg_half_1;
%     regime_ts_egg_mask(half_len+1:length(tmp.sum_ts_egg_mask))=egg_half_2;
%     mean_ts_nwv= tmp.mean_ts_v_rho.*cosd(45)-tmp.mean_ts_u_rho.*cosd(45);

    axLH = gca;
    mslplot{2}=bar(tmp.sum_ts_egg_mask, 'FaceColor', [0.6 0.6 0.6], 'parent',axLH);
    hold on
%     mslplot{1}=plot(regime_ts_egg_mask,'k','parent',axLH);
%                 mslplot{2}=plot(tmp.sum_ts_egg_mask, 'Color', 'b', 'LineStyle', '-', 'parent',axLH);
%                 mslplot{2}=plot(tmp.sum_ts_egg_mask, 'color', [0.8 0.8 0.8], 'parent',axLH);

%     ylabel(axLH,'# of eggs')
    ylabel(axLH,'Spawning area (km^2)')    
    ax_pos = get(axLH,'position');
    set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
    cal_arr=1;
    for cali=1:length(RCM_info.years_his)-1
        cal_arr=[cal_arr, cal_arr(end)+sum(eomday(RCM_info.years_his(cali),RCM_info.months))];
    end
    set(axLH,'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(tmp.comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);

    axis tight 

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
    set(axLH,'ycolor','k', 'box', 'off', 'FontSize',15);
%                 set(axRH,'ycolor','b', 'box', 'off', 'FontSize',15);
    xlabel(axLH, 'Year');

    title(['Spawning area', ',', tmp.testname_his, ',', num2str(min(RCM_info.years_his),'%04i'),'-',num2str(max(RCM_info.years_his),'%04i'), ...
        ', ', RCM_info.season]);

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

%     set(mslplot{1},'LineWidth',2);
%                 set(mslplot{2},'LineWidth',2);
    grid on

    lgd=legend([mslplot{1}], 'Spawning area (km^2)');

%     half_len=length(tmp.sum_ts_egg_mask)/2;
%     egg_half_1=mean(tmp.sum_ts_egg_mask(1:half_len));
%     egg_half_2=mean(tmp.sum_ts_egg_mask(half_len+1:end));
%                 txt1=text(5, max(double(mean_ts_nwv))-diff(double([min(mean_ts_nev), max(mean_ts_nwv)]))/32.0 , ...
%                     ['egg# = ', num2str(round(egg_half_1,1)), ' / ', num2str(round(egg_half_2,1))], 'FontSize', 20); 

    set(lgd,'FontSize',15);
    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    set(lgd,'Orientation','horizontal');

    set(gcf,'PaperPosition', [0 0 36 12]) 
    saveas(gcf,tmp.tifname,'tif');RemoveWhiteSpace([], 'file', tmp.tifname);
    grid off
    hold off
    disp(' ')
    disp([num2str(flagi), ' plot is created.'])
    disp(' ')
    disp([' File path is : ',tmp.tifname])
    disp(' ')

    close all;
end


 %% future(ssp) plot
    tmp.tifname=strcat(dirs.regime_figdir, tmp.testname_ssp,'_',tmp.regionname, '_regime_ts_sp_ground_', ...
                num2str(min(RCM_info.years_ssp),'%04i'),'_',num2str(max(RCM_info.years_ssp),'%04i'), 'y_', ...
                RCM_info.season, '.tif'); 
            
if (exist(tmp.tifname , 'file') ~= 2 || flags.fig_switch(flagi)==2)
    run(tmp.param_script);
    
    %%    initialization
    if(isfield(tmp, 'tlen_ssp')~=1)
        tmp.tlen_ssp=0;
        for yearij = 1:length(RCM_info.years_ssp)
            tmp.tempyear = RCM_info.years_ssp(yearij);
            for monthij = 1:length(RCM_info.months)
                tmp.tempmonth = RCM_info.months(monthij);
                tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
                tmp.tlen_ssp=tmp.tlen_ssp + length(ncread(tmp.ncname, 'time'));
            end
        end
        tmp.xlen=size(ncread(tmp.ncname,'lon_rho'),1);
        tmp.ylen=size(ncread(tmp.ncname,'lon_rho'),2);
    end
    tmp.comb_egg_mask=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp);
    tmp.comb_u_rho=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp);
    tmp.comb_v_rho=NaN(tmp.xlen,tmp.ylen,tmp.tlen_ssp);
    tmp.comb_ocean_time=NaN(tmp.tlen_ssp);
% %     
% %     
% %     clear wind_curl comb_wind_curl u_rho tmp.comb_u_rho v_rho tmp.comb_v_rho uwind comb_uwind vwind comb_vwind ...
% %         egg_mask tmp.comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 tmp.temp_surf comb_tmp.temp_surf  curl_mask


    for yearij = 1:length(RCM_info.years_ssp)
        tmp.tempyear = RCM_info.years_ssp(yearij);
        for monthij = 1:length(RCM_info.months)
            tmp.tempmonth = RCM_info.months(monthij);
            tmp.ncname = [dirs.savedir_ssp,tmp.testname_ssp,'_',tmp.regionname,'model_pollock_',num2str(tmp.tempyear,'%04i'),'_',num2str(tmp.tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
            
            if yearij==1 && monthij==1
                RCM_grid.lon_rho = ncread(tmp.ncname, 'lon_rho');
                RCM_grid.lat_rho = ncread(tmp.ncname, 'lat_rho');
                RCM_grid.mean_diff_lon = mean(diff(RCM_grid.lon_rho(:,1)));
                RCM_grid.mean_diff_lat = mean(diff(RCM_grid.lat_rho(1,:)));
%                 RCM_grid.dA = (m_lldist([RCM_grid.mean_diff_lon 0], [0 0])*1e3)^2 * cosd(RCM_grid.lat_rho); %m^2
                RCM_grid.dA = (m_lldist([RCM_grid.mean_diff_lon 0], [0 0]))^2 * cosd(RCM_grid.lat_rho); %km^2
%                 tmp.refpolygon=[0, 0; 0, 39; 360, 39; 360,0];
                tmp.vel_mask = double(inpolygon(RCM_grid.lon_rho,RCM_grid.lat_rho,tmp.refpolygon(:,1),tmp.refpolygon(:,2)));
            end
            ocean_time=ncread(tmp.ncname, 'time')+datenum(1900,12,31);
            tmp.u_rho=ncread(tmp.ncname, 'u_rho').*tmp.vel_mask;
            tmp.u_rho(tmp.u_rho==0)=NaN;
            tmp.v_rho=ncread(tmp.ncname, 'v_rho').*tmp.vel_mask;
            tmp.v_rho(tmp.v_rho==0)=NaN;
            tmp.egg_mask=ncread(tmp.ncname, 'egg_mask').*tmp.vel_mask;
            tmp.egg_mask(tmp.egg_mask==0)=NaN;
            tmp.egg_mask=ncread(tmp.ncname, 'egg_mask').*RCM_grid.dA;
            
            tmp.lastday_m=size(tmp.egg_mask,3);
            if yearij==1 && monthij==1
                tmp.comb_u_rho=tmp.u_rho;
                tmp.comb_v_rho=tmp.v_rho;
                tmp.comb_egg_mask=tmp.egg_mask;
                tmp.comb_ocean_time=ocean_time;
                tmp.endij=tmp.lastday_m;
            else
                tmp.comb_u_rho(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.u_rho;
                tmp.comb_v_rho(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.v_rho;
                tmp.comb_egg_mask(:,:,tmp.endij+1:tmp.endij+tmp.lastday_m)=tmp.egg_mask;
                tmp.comb_ocean_time(tmp.endij+1:tmp.endij+tmp.lastday_m)=ocean_time;
                tmp.endij=tmp.endij+tmp.lastday_m;
            end
        end
    end

    tmp.ts_u_rho=reshape(tmp.comb_u_rho,[size(tmp.comb_u_rho,1)*size(tmp.comb_u_rho,2), size(tmp.comb_u_rho,3)]);
    tmp.mean_ts_u_rho=mean(tmp.ts_u_rho,1,'omitnan');
    tmp.ts_v_rho=reshape(tmp.comb_v_rho,[size(tmp.comb_v_rho,1)*size(tmp.comb_v_rho,2), size(tmp.comb_v_rho,3)]);
    tmp.mean_ts_v_rho=mean(tmp.ts_v_rho,1,'omitnan');
    tmp.ts_egg_mask=reshape(tmp.comb_egg_mask,[size(tmp.comb_egg_mask,1)*size(tmp.comb_egg_mask,2), size(tmp.comb_egg_mask,3)]);
    tmp.sum_ts_egg_mask=sum(tmp.ts_egg_mask,1,'omitnan');

    half_len=round(length(tmp.sum_ts_egg_mask)/2);
    egg_half_1=mean(tmp.sum_ts_egg_mask(1:half_len));
    egg_half_2=mean(tmp.sum_ts_egg_mask(half_len+1:end));
%     regime_ts_egg_mask(1:half_len)=egg_half_1;
%     regime_ts_egg_mask(half_len+1:length(tmp.sum_ts_egg_mask))=egg_half_2;
%     mean_ts_nwv= tmp.mean_ts_v_rho.*cosd(45)-tmp.mean_ts_u_rho.*cosd(45);

    axLH = gca;
    mslplot{2}=bar(tmp.sum_ts_egg_mask, 'FaceColor', [0.6 0.6 0.6], 'parent',axLH);
    hold on
%     mslplot{1}=plot(regime_ts_egg_mask,'k','parent',axLH);
%                 mslplot{2}=plot(tmp.sum_ts_egg_mask, 'Color', 'b', 'LineStyle', '-', 'parent',axLH);
%                 mslplot{2}=plot(tmp.sum_ts_egg_mask, 'color', [0.8 0.8 0.8], 'parent',axLH);

%     ylabel(axLH,'# of eggs')
    ylabel(axLH,'Spawning area (km^2)')    
    ax_pos = get(axLH,'position');
    set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
    cal_arr=1;
    for cali=1:length(RCM_info.years_ssp)-1
        cal_arr=[cal_arr, cal_arr(end)+sum(eomday(RCM_info.years_ssp(cali),RCM_info.months))];
    end
    set(axLH,'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(tmp.comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);

    axis tight 

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
    set(axLH,'ycolor','k', 'box', 'off', 'FontSize',15);
%                 set(axRH,'ycolor','b', 'box', 'off', 'FontSize',15);
    xlabel(axLH, 'Year');

    title(['Spawning area', ',', tmp.testname_his, ',', num2str(min(RCM_info.years_ssp),'%04i'),'-',num2str(max(RCM_info.years_ssp),'%04i'), ...
        ', ', RCM_info.season]);

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

%     set(mslplot{1},'LineWidth',2);
%                 set(mslplot{2},'LineWidth',2);
    grid on

    lgd=legend([mslplot{1}], 'Spawning area (km^2)');

%     half_len=length(tmp.sum_ts_egg_mask)/2;
%     egg_half_1=mean(tmp.sum_ts_egg_mask(1:half_len));
%     egg_half_2=mean(tmp.sum_ts_egg_mask(half_len+1:end));
%                 txt1=text(5, max(double(mean_ts_nwv))-diff(double([min(mean_ts_nev), max(mean_ts_nwv)]))/32.0 , ...
%                     ['egg# = ', num2str(round(egg_half_1,1)), ' / ', num2str(round(egg_half_2,1))], 'FontSize', 20); 

    set(lgd,'FontSize',15);
    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    set(lgd,'Orientation','horizontal');

    set(gcf,'PaperPosition', [0 0 36 12]) 
    saveas(gcf,tmp.tifname,'tif'); RemoveWhiteSpace([], 'file', tmp.tifname);
    grid off
    hold off
    disp(' ')
    disp([num2str(flagi), ' plot is created.'])
    disp(' ')
    disp([' File path is : ',tmp.tifname])
    disp(' ')

    close all;
end