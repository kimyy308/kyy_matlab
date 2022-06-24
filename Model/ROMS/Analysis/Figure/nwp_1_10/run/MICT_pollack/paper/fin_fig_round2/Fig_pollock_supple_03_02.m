clearvars '*' -except regionind testnameind all_region all_testname
        % % % 
        dropboxpath='C:\users\user/Dropbox/';
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\MICT_pollack\paper\subroutine\'))
        addpath(genpath('C:\Users\User\Dropbox\source\matlab\function\'))
        [dropboxpath, erorr_status] = Func_0008_set_dropbox_path(computer);
        shadlev = [-2 2];
        trend_shadlev = [0 4];
        conlev  = 0:5:35;
        dl=1/10;
        
        [byrmap3, error_status] = Func_0009_get_colormaps('byr3', dropboxpath);
        [byrmap, error_status] = Func_0009_get_colormaps('byr2', dropboxpath);        
        yrmap = byrmap(129:256,:);
        
        % for snu_desktop
        testname='test06'
        inputyear = [1988:1992]; % % put year which you want to plot [year year ...]
        
        allyear=[1983:1992];
        
%         refyear = [1983];
        refyear =[1983:1987];

        inputmonth = [1,2]; % % put month which you want to plot [month month ...]
        checktime=[15,30];
        varname ='zeta';
        variable='zeta';
        regionname='pollock_egg3';
        run('nwp_polygon_point.m');
        
% % %         switch region
        [refpolygon, lonlat, error_status] = Func_0007_get_polygon_data_from_regionname(regionname);

        param_script =['C:\users\user/Dropbox/source/matlab/Model/ROMS/Analysis/Figure/nwp_1_10/run/fig_param/fig_param_kyy_', regionname, '.m'];
        figrawdir =strcat('D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\', testname, '\DA\'); % % where figure files will be saved
        filedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\'); % % where data files are
        savedir = strcat('D:\Data\Model\ROMS\nwp_1_10\test06\DA\pollock_6\');
        inputdir = ['/home/auto/MAMS/Data/01_NWP_1_10/Input/'];
        mkdir(savedir);
        LTRANS_testname='Pollock6';

    run(param_script);
    clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
        egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask
    for yearij = 1:length(allyear)
        tempyear = allyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
            if(exist('curl_mask')==0)
                lon_rho = ncread(ncname, 'lon_rho');
                lat_rho = ncread(ncname, 'lat_rho');
                vel_polygon=[134, 40; 130, 42.5; 142 52; 142 46; 139.5 46.5];
                vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
                SK_EEZ_mask = double(inpolygon(lon_rho,lat_rho,SK_EEZ_polygon(:,1),SK_EEZ_polygon(:,2)));
                NK_EEZ_mask = double(inpolygon(lon_rho,lat_rho,NK_EEZ_polygon(:,1),NK_EEZ_polygon(:,2)));
                RU_PR_EEZ_mask = double(inpolygon(lon_rho,lat_rho,RU_PR_EEZ_polygon(:,1),RU_PR_EEZ_polygon(:,2)));
                RU_SH_EEZ_mask = double(inpolygon(lon_rho,lat_rho,RU_SH_EEZ_polygon(:,1),RU_SH_EEZ_polygon(:,2)));
                JP_EEZ_mask = double(inpolygon(lon_rho,lat_rho,JP_EEZ_polygon(:,1),JP_EEZ_polygon(:,2)));
            end
            ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
            u_rho=ncread(ncname, 'u_rho').*vel_mask;
            u_rho(u_rho==0)=NaN;
            v_rho=ncread(ncname, 'v_rho').*vel_mask;
            v_rho(v_rho==0)=NaN;
%                         egg_mask=ncread(ncname, 'egg_mask').*vel_mask;
%                         egg_mask(egg_mask==0)=NaN;
            SK_EEZ_egg_mask=ncread(ncname, 'egg_mask').*SK_EEZ_mask;
            NK_EEZ_egg_mask=ncread(ncname, 'egg_mask').*NK_EEZ_mask;
            RU_PR_EEZ_egg_mask=ncread(ncname, 'egg_mask').*RU_PR_EEZ_mask;
            RU_SH_EEZ_egg_mask=ncread(ncname, 'egg_mask').*RU_SH_EEZ_mask;
            JP_EEZ_egg_mask=ncread(ncname, 'egg_mask').*JP_EEZ_mask;
            egg_mask =ncread(ncname, 'egg_mask');
            lastday_m=size(u_rho,3);
            if (exist('comb_u_rho')==0)
                comb_u_rho=u_rho;
                comb_v_rho=v_rho;
                comb_SK_EEZ_egg_mask=SK_EEZ_egg_mask;
                comb_NK_EEZ_egg_mask=NK_EEZ_egg_mask;
                comb_RU_PR_EEZ_egg_mask=RU_PR_EEZ_egg_mask;
                comb_RU_SH_EEZ_egg_mask=RU_SH_EEZ_egg_mask;
                comb_JP_EEZ_egg_mask=JP_EEZ_egg_mask;
                comb_egg_mask=egg_mask;
                comb_ocean_time=ocean_time;
            else
                comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
                comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
                comb_SK_EEZ_egg_mask(:,:,end+1:end+lastday_m)=SK_EEZ_egg_mask;
                comb_NK_EEZ_egg_mask(:,:,end+1:end+lastday_m)=NK_EEZ_egg_mask;
                comb_RU_PR_EEZ_egg_mask(:,:,end+1:end+lastday_m)=RU_PR_EEZ_egg_mask;
                comb_RU_SH_EEZ_egg_mask(:,:,end+1:end+lastday_m)=RU_SH_EEZ_egg_mask;
                comb_JP_EEZ_egg_mask(:,:,end+1:end+lastday_m)=JP_EEZ_egg_mask;
                comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
                comb_ocean_time(end+1:end+lastday_m)=ocean_time;
            end
        end
    end

    ts_u_rho=reshape(comb_u_rho,[size(comb_u_rho,1)*size(comb_u_rho,2), size(comb_u_rho,3)]);
    mean_ts_u_rho=mean(ts_u_rho,1,'omitnan');
    ts_v_rho=reshape(comb_v_rho,[size(comb_v_rho,1)*size(comb_v_rho,2), size(comb_v_rho,3)]);
    mean_ts_v_rho=mean(ts_v_rho,1,'omitnan');

    ts_SK_EEZ_egg_mask=reshape(comb_SK_EEZ_egg_mask,[size(comb_SK_EEZ_egg_mask,1)*size(comb_SK_EEZ_egg_mask,2), size(comb_SK_EEZ_egg_mask,3)]);
    sum_ts_SK_EEZ_egg_mask=sum(ts_SK_EEZ_egg_mask,1,'omitnan');
    ts_NK_EEZ_egg_mask=reshape(comb_NK_EEZ_egg_mask,[size(comb_NK_EEZ_egg_mask,1)*size(comb_NK_EEZ_egg_mask,2), size(comb_NK_EEZ_egg_mask,3)]);
    sum_ts_NK_EEZ_egg_mask=sum(ts_NK_EEZ_egg_mask,1,'omitnan');
    ts_RU_PR_EEZ_egg_mask=reshape(comb_RU_PR_EEZ_egg_mask,[size(comb_RU_PR_EEZ_egg_mask,1)*size(comb_RU_PR_EEZ_egg_mask,2), size(comb_RU_PR_EEZ_egg_mask,3)]);
    sum_ts_RU_PR_EEZ_egg_mask=sum(ts_RU_PR_EEZ_egg_mask,1,'omitnan');
    ts_RU_SH_EEZ_egg_mask=reshape(comb_RU_SH_EEZ_egg_mask,[size(comb_RU_SH_EEZ_egg_mask,1)*size(comb_RU_SH_EEZ_egg_mask,2), size(comb_RU_SH_EEZ_egg_mask,3)]);
    sum_ts_RU_SH_EEZ_egg_mask=sum(ts_RU_SH_EEZ_egg_mask,1,'omitnan');
    ts_JP_EEZ_egg_mask=reshape(comb_JP_EEZ_egg_mask,[size(comb_JP_EEZ_egg_mask,1)*size(comb_JP_EEZ_egg_mask,2), size(comb_JP_EEZ_egg_mask,3)]);
    sum_ts_JP_EEZ_egg_mask=sum(ts_JP_EEZ_egg_mask,1,'omitnan');
    ts_egg_mask=reshape(comb_egg_mask,[size(comb_egg_mask,1)*size(comb_egg_mask,2), size(comb_egg_mask,3)]);
    sum_ts_egg_mask=sum(ts_egg_mask,1,'omitnan');

%                 half_len=round(length(sum_ts_egg_mask)/2);
%                 egg_half_1=mean(sum_ts_egg_mask(1:half_len));
%                 egg_half_2=mean(sum_ts_egg_mask(half_len+1:end));
%                 regime_ts_egg_mask(1:half_len)=egg_half_1;
%                 regime_ts_egg_mask(half_len+1:length(sum_ts_egg_mask))=egg_half_2;

    mean_ts_nwv= mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45);

    axLH = gca;
%     mslplot{1}=plot(sum_ts_SK_EEZ_egg_mask, 'Color', 'k', 'parent',axLH);
%     hold on
%     mslplot{2}=plot(sum_ts_NK_EEZ_egg_mask, 'Color', 'r', 'parent',axLH);
%     mslplot{3}=plot(sum_ts_RU_PR_EEZ_egg_mask, 'Color', 'b', 'parent',axLH);
%     mslplot{4}=plot(sum_ts_RU_SH_EEZ_egg_mask, 'Color', 'm', 'parent',axLH);
%     mslplot{5}=plot(sum_ts_JP_EEZ_egg_mask, 'Color', 'g', 'parent',axLH);
    mslplot{1}=plot(sum_ts_egg_mask, 'Color', 'k', 'parent',axLH);

%                 mslplot{1}=plot(regime_ts_egg_mask,'k','parent',axLH);
%                 mslplot{2}=plot(sum_ts_egg_mask, 'Color', 'b', 'LineStyle', '-', 'parent',axLH);
%                 mslplot{2}=plot(sum_ts_egg_mask, 'color', [0.8 0.8 0.8], 'parent',axLH);

    ylabel(axLH,'# of released eggs')
    ax_pos = get(axLH,'position');
    set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
    cal_arr=1;
    for cali=1:length(allyear)-1
        cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth))];
    end
    set(axLH,'color','none','yaxislocation','left','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);

    axis tight 

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
    set(axLH,'ycolor','k', 'box', 'off', 'FontSize',15);
%                 set(axRH,'ycolor','b', 'box', 'off', 'FontSize',15);
    xlabel(axLH, 'Year (Jan and Feb only)');
% 
%     title(['Spawned eggs', ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i'), ...
%         ', ', num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i')]);

%                 set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

    set(mslplot{1},'LineWidth',2);
%     set(mslplot{2},'LineWidth',2);
%     set(mslplot{3},'LineWidth',2);
%     set(mslplot{4},'LineWidth',2);
%     set(mslplot{5},'LineWidth',2);
%                 set(mslplot{2},'LineWidth',2);
    grid on

%     lgd=legend([mslplot{1}, mslplot{2}, mslplot{3}, mslplot{4}, mslplot{5}], 'SK', 'NK', 'RU-PR', 'RU-SH', 'JP');
    lgd=legend([mslplot{1}], '# of released eggs');

%                 half_len=length(sum_ts_egg_mask)/2;
%                 egg_half_1=mean(sum_ts_egg_mask(1:half_len));
%                 egg_half_2=mean(sum_ts_egg_mask(half_len+1:end));

%                 txt1=text(5, max(double(mean_ts_nwv))-diff(double([min(mean_ts_nev), max(mean_ts_nwv)]))/32.0 , ...
%                     ['egg# = ', num2str(round(egg_half_1,1)), ' / ', num2str(round(egg_half_2,1))], 'FontSize', 20); 

    set(lgd,'FontSize',15);
    set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
    set(lgd,'Orientation','horizontal');

    set(gcf,'PaperPosition', [0 0 36 12]) 
    jpgname = 'D:\Research\Ph_D_course\2021_Particle tracking of the walleye pollock\02_review\Review_round_02\figure_review_round_02\supple_fig_01_02.tif';
    saveas(gcf,jpgname,'tif');
    grid off
    hold off
    disp(' ')
%     disp([fig_name, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    close all;