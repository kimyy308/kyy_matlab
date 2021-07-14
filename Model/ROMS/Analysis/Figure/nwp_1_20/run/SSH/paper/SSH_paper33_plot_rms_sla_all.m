close all; clear all;  clc;
warning off;
% all_region2 ={'NWP','NWP2'}
% all_region2 ={'ES', 'YS', 'SS', 'NWP2'}
% all_testname2 = {'ens02','test53','test55','test56'};
% all_testname2 = {'test11', 'test12'};
all_testname2 = {'test11', 'test12'};
testname = 'test11';
testname2 = 'test12';
sodatestname = 'SODA_3_4_2';
% all_region2 ={'NWP','ES', 'SS', 'YS', 'AKP'}
% all_region2 ={'NWP','AKP2', 'YS', 'NES', 'SES'}

% all_region2 ={'NWP'};
all_region2 ={'NWP'};
% all_region2 ={'AKP4','YSECS', 'ECS2', 'ES', 'YS', 'NES', 'SES'}



% start-------------------- configuration
for regionind2=1:length(all_region2)
    for folding=1:1

    close all;
    clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 testname testname2 sodatestname
    % % % 
    % % % Read Model SST
    % % % interp
    % % % get RMS
    % % % get BIAS
    system_name=computer;
    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        dropboxpath='C:\Users\KYY\Dropbox';
        addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
        addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
        addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
        addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
        addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
    elseif (strcmp(system_name,'GLNXA64'))
        dropboxpath='/home/kimyy/Dropbox';
        addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
        addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
        addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
        addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
    end

%     regionind2 = 2;

    shadlev = [0 35];
    rms_shadlev = [0 4];
    %     trendlev = [-3 3];  %% trend lev
    trendlev = [-10 10];  %% trend lev
    abstrendlev =[2 7];
    reltrendlev =[-5 5];
    conlev  = 0:5:35;
    meanplotlev =[-0.2 0.2];
    trendplotlev = [0 7];
    trenddifflev = [-10 10];
    sshlev =[-0.3 0.3];
    sshdifflev = [0 20];

    % for snu_desktopd
    % testname=all_testname2{testnameind2}    % % need to change
    inputyear = [1994:2014]; % % put year which you want to plot [year year ...]
    inputmonth = [2 5 8 11]; % % put month which you want to plot [month month ...]

    varname ='zeta';
    variable='SSH';
    run('nwp_polygon_point.m');
    regionname=all_region2{regionind2};
    switch(regionname)
        case('NWP') %% North western Pacific
            lonlat = [115, 164, 15, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('NWP2') %% North western Pacific
            lonlat = [115, 145, 25, 52];  %% whole data area
            refpolygon(1,1)=lonlat(1);
            refpolygon(2,1)=lonlat(2);
            refpolygon(1,2)=lonlat(3);
            refpolygon(2,2)=lonlat(4);
        case('ES') %% East Sea
            refpolygon=espolygon;
        case('NES') %% Northern East Sea
            refpolygon=nespolygon;
        case('SES') %% Southern East Sea
            refpolygon=sespolygon;
        case('SS') %% South Sea
            refpolygon=sspolygon;
        case('YS') %% Yellow Sea
            refpolygon=yspolygon;
        case('ECS') %% East China Sea
            refpolygon=ecspolygon;
        case('ECS2') %% East China Sea
            refpolygon=ecs2polygon;
        case('YSECS') %% East China Sea
            refpolygon=ysecspolygon;
        case('AKP') %% Around Korea Peninsula
            refpolygon=akppolygon;
        case('AKP2') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        case('AKP3') %% Around Korea Peninsula
            refpolygon=akp3polygon;
        case('AKP4') %% Around Korea Peninsula
            refpolygon=akp4polygon;
        case('CA') %% Around Korea Peninsula
            refpolygon=capolygon;
        case('EKB') %% Around Korea Peninsula
            refpolygon=akp2polygon;
        otherwise
            ('?')
    end
    lonlat(1)=min(refpolygon(:,1));
    lonlat(2)=max(refpolygon(:,1));
    lonlat(3)=min(refpolygon(:,2));
    lonlat(4)=max(refpolygon(:,2));

    % % % for EKB
    % regionname='EKB';
    % lonlat = [127, 129.5, 38, 40.5];

    %         load(['G:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
    %             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
    filename = ['E:\Data\Model\ROMS\nwp_1_10\',testname,'\run\',testname,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

    confilename = ['E:\Data\Model\ROMS\nwp_1_10\',testname2,'\run\',testname2,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

    sodafilename =  ['E:\Data\Reanalysis\SODA\', sodatestname, '\',sodatestname,'_',regionname, ...
                '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];



    valnum=0;
    run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
    wrmap = bwrmap(51:100,:);


    if (strcmp(system_name,'PCWIN64'))
        % % for windows
        figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\4th_year\figure\all\',testname,'\',regionname,'\'); % % where figure files will be saved
    %             param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_EKB_RMS.m'
        param_script =['C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_', regionname, '.m']
        filedir = strcat('E:\Data\Model\ROMS\nwp_1_10\', testname, '\run\'); % % where data files are
    elseif (strcmp(system_name,'GLNXA64'))
    end
    %         std(mean(mean(comb_data_filtered,'omitnan'),'omitnan'))

    run(param_script);

    figdir=[figrawdir,'Trend\SSH\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir,regionname);

    cmems_trend=ncread(filename, 'cmems_trend');
    cmems_mask=ones(size(cmems_trend));
    cmems_mask(isnan(cmems_trend))=NaN;

    lon_rho=ncread(filename,'lon_rho');
    lat_rho=ncread(filename,'lat_rho');

    end
    % end-------------------- configuration

    lon_rho=ncread(filename,'lon_rho');
    lat_rho=ncread(filename,'lat_rho');

    cmems_sla = ncread(filename, 'cmems_sla');
    cmems_2y_lowpass_sla = ncread(filename, 'cmems_2y_lowpass');
    xlen=size(cmems_sla,1); ylen = size(cmems_sla,2); tlen = size(cmems_sla,3);
    % % 
    % % for t=1:tlen
    % %     cmems_sla(:,:,t)=cmems_sla(:,:,t).*cmems_mask;
    % %     cmems_2y_lowpass_sla(:,:,t)=cmems_2y_lowpass_sla(:,:,t).*cmems_mask;
    % % end
    % % cmems_sla_2d=reshape(cmems_sla,[xlen*ylen, tlen]);
    % % m_cmems_sla=mean(cmems_sla_2d,1,'omitnan');
    % % cmems_2y_lowpass_sla_2d=reshape(cmems_2y_lowpass_sla,[xlen*ylen, tlen]);
    % % cmems_2y_lowpass_sla_mean=mean(cmems_2y_lowpass_sla_2d,1,'omitnan');
    % % cmems_sla_3d=reshape(cmems_sla,[xlen*ylen*12, tlen/12]);
    % % cmems_sla_yearly_mean=mean(cmems_sla_3d,1,'omitnan');
    % % cmems_sla_4d=reshape(cmems_sla,[xlen*ylen, 12, tlen/12]);
    % % cmems_sla_seasonal_mean=mean(mean(cmems_sla_4d,1,'omitnan'),3,'omitnan');
    % % for t=1:tlen
    % %     if mod(t,12)==0
    % %         tt=12;
    % %     else
    % %         tt=mod(t,12);
    % %     end
    % %     cmems_sla_seasonal_filtered(t)=m_cmems_sla(t)-cmems_sla_seasonal_mean(tt);
    % % end
    % % 
    % % cmems_sla = ncread(filename,'cmems_sla');
    % % for sla_i=1:size(cmems_sla,1)
    % %     for sla_j=1:size(cmems_sla,2)
    % %         cmems_sla_mean(sla_i,sla_j)=mean(cmems_sla(sla_i,sla_j,:),'omitnan');
    % %     end
    % % end
    % % 
    % % interped_ssh=ncread(filename,'interped_ssh');
    % % for sla_i=1:size(cmems_sla,1)
    % %     for sla_j=1:size(cmems_sla,2)
    % %         interped_sla_mean(sla_i,sla_j)=mean(interped_ssh(sla_i,sla_j,:),'omitnan');
    % %         interped_sla(sla_i,sla_j,:)=interped_ssh(sla_i,sla_j,:)-(interped_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
    % %     end
    % % end
    % % interped_sla_divided=reshape(interped_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
    % % clim_interped_sla=mean(interped_sla_divided,4,'omitnan');
    % % for t=1:length(inputyear)
    % %     interped_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped_sla);
    % % end
    % % 
    % % interped2_ssh = ncread(filename2, 'interped_ssh');
    % % 
    % % for sla_i=1:size(cmems_sla,1)
    % %     for sla_j=1:size(cmems_sla,2)
    % %         interped2_sla_mean(sla_i,sla_j)=mean(interped2_ssh(sla_i,sla_j,:),'omitnan');
    % %         interped2_sla(sla_i,sla_j,:)=interped2_ssh(sla_i,sla_j,:)-(interped2_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
    % %     end
    % % end
    % % interped2_sla_divided=reshape(interped2_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
    % % clim_interped2_sla=mean(interped2_sla_divided,4,'omitnan');
    % % for t=1:length(inputyear)
    % %     interped2_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(interped2_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_interped2_sla);
    % % end
    % % 
    % % 
    % % soda_ssh = ncread(sodafilename, 'interped_ssh');
    % % 
    % % for sla_i=1:size(cmems_sla,1)
    % %     for sla_j=1:size(cmems_sla,2)
    % %         soda_sla_mean(sla_i,sla_j)=mean(soda_ssh(sla_i,sla_j,:),'omitnan');
    % %         soda_sla(sla_i,sla_j,:)=soda_ssh(sla_i,sla_j,:)-(soda_sla_mean(sla_i,sla_j)-cmems_sla_mean(sla_i,sla_j));
    % %     end
    % % end
    % % 
    % % soda_sla_divided=reshape(soda_sla,[size(cmems_sla,1), size(cmems_sla,2), 12, length(inputyear)]);
    % % clim_soda_sla=mean(soda_sla_divided,4,'omitnan');
    % % for t=1:length(inputyear)
    % %     soda_sla_filtered(:,:,(t-1)*12+1:(t-1)*12+12)=single(soda_sla(:,:,(t-1)*12+1:(t-1)*12+12)-clim_soda_sla);
    % % end






    % start-------------------- rms_ lowpassed get
    for folding=1:1
            nyears=[1:5];
            nc_varname_prefixes={'interped_sla'};
            nc_titlename_prefixes={'lowpass'};

    %         for nc_varnameij=1:length(nc_varname_prefixes)
    %             nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
    %             nc_varname=[nc_varname_prefix];
    %             eval(['corr=ncread(filename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(confilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_con(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(sodafilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_soda(1,nc_varnameij)=mean(corr(:),'omitnan');
    %         end
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, '_', num2str(nyear),'y_lowpass'];
    %                 jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
    %                     '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    %                 if (exist(jpgname , 'file') ~= 2)



                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(1,nc_varnameij,:,:,:)=sq_diff_all;

                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
    %                     eval(['corr=ncread(confilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_con(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');
    %                     eval(['corr=ncread(sodafilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_soda(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

    %                     eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{1};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mrms,2)));  

    %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    %                     saveas(gcf,jpgname,'tif');

                        disp(' ')
%                         disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
    %                     disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
    %                 end
                end
            end
    end
% end-------------------- rms_ lowpassed get


% % % start-------------------- RMS perfomance by lowpass filter window
    for folding=1:1
         jpgname=strcat(outfile, '_', testname, '_',regionname, 'rms_lp_perfomance', '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(0:length(nyears),mrms(:,1),'r')
            hold on
            msl2plot=plot(0:length(nyears),mrms_con(:,1),'g')
            msl3plot=plot(0:length(nyears),mrms_soda(:,1),'b')

            xlabel('Low pass filter window(year)')
            ylabel('RMS')
            title([regionname, ', RMS performance '])
    %         datetick('x','yyyy','keepticks')
            axis tight;
            ylim([0 7])
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCM-SAT','RCM-GLO','GCM');
            set(lgd,'FontSize',10);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(gcf,'PaperPosition', [0 0 24 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            set(gca,'XTick', (0:1:length(nyears)))
            saveas(gcf,jpgname,'jpg');
            grid off
            end
        close all;
    end
    % % % end-------------------- RMS perfomance by lowpass filter window
    
    
    % % % start-------------------- RMS perfomance by lowpass filter window
% (seasonal)
    for folding=1:1
        for monthij=1:length(inputmonth)
            tempmonth=inputmonth(monthij);
            jpgname=strcat(outfile, '_', testname, '_',regionname, '_', num2str(tempmonth, '%02i'),'_rms_lp_perfomance', '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2)
                if (exist('mrms_seasonal' , 'var') ~= 1)
                    size_window=size(sq_diff_all_rcm_sat,1);
                    size_var=size(sq_diff_all_rcm_sat,2);
                    size_x=size(sq_diff_all_rcm_sat,3);
                    size_y=size(sq_diff_all_rcm_sat,4);
                    size_t=size(sq_diff_all_rcm_sat,5);
                    sq_diff_all_rcm_sat_divided = ...
                        reshape(sq_diff_all_rcm_sat, [size_window, size_var, size_x*size_y, size_t]);
                    sq_diff_all_rcm_glo_divided = ...
                        reshape(sq_diff_all_rcm_glo, [size_window, size_var, size_x*size_y, size_t]);
                    sq_diff_all_gcm_soda_divided = ...
                        reshape(sq_diff_all_gcm_soda, [size_window, size_var, size_x*size_y, size_t]);

                    mrms_rcm_sat_seasonal_spatial=mean(sq_diff_all_rcm_sat_divided,3,'omitnan');
                    mrms_rcm_glo_seasonal_spatial=mean(sq_diff_all_rcm_glo_divided,3,'omitnan');
                    mrms_gcm_soda_seasonal_spatial=mean(sq_diff_all_gcm_soda_divided,3,'omitnan');

                    mrms_rcm_sat_seasonal_divided=reshape(mrms_rcm_sat_seasonal_spatial, [size_window, size_var, size_t/12, 12]);
                    mrms_rcm_glo_seasonal_divided=reshape(mrms_rcm_glo_seasonal_spatial, [size_window, size_var, size_t/12, 12]);
                    mrms_gcm_soda_seasonal_divided=reshape(mrms_gcm_soda_seasonal_spatial, [size_window, size_var, size_t/12, 12]);
                  
                    for windowij=1:size_window
                        for varij=1:size_var
                            mrms_rcm_sat_seasonal(windowij,varij,1,tempmonth)=sqrt(mean(mrms_rcm_sat_seasonal_divided(windowij,varij,:,tempmonth),3,'omitnan'));
                            mrms_rcm_glo_seasonal(windowij,varij,1,tempmonth)=sqrt(mean(mrms_rcm_glo_seasonal_divided(windowij,varij,:,tempmonth),3,'omitnan'));
                            mrms_gcm_soda_seasonal(windowij,varij,1,tempmonth)=sqrt(mean(mrms_gcm_soda_seasonal_divided(windowij,varij,:,tempmonth),3,'omitnan'));
                        end
                    end
                end
%                 mslplot=plot(0:length(nyears),mrms(:,1),'r')
%                 hold on
%                 msl2plot=plot(0:length(nyears),mrms_con(:,1),'g')
%                 msl3plot=plot(0:length(nyears),mrms_soda(:,1),'b')

                mslplot=plot(0:length(nyears),mrms_rcm_sat_seasonal(:,varij,1,tempmonth),'r')
                hold on
                msl2plot=plot(0:length(nyears),mrms_rcm_glo_seasonal(:,varij,1,tempmonth),'g')
                msl3plot=plot(0:length(nyears),mrms_gcm_soda_seasonal(:,varij,1,tempmonth),'b')


                xlabel('Low pass filter window(year)')
                ylabel('RMS')
                title([regionname, ', RMS performance ,', calendarname{tempmonth}(1:3)])
        %         datetick('x','yyyy','keepticks')
                axis tight;
                ylim([0 10])
                set(mslplot,'LineWidth',2);
                set(msl2plot,'LineWidth',2);
                set(msl3plot,'LineWidth',2);
                set(gca,'FontSize',15);
                grid on
                lgd=legend('RCM-SAT','RCM-GLO','GCM');
                set(lgd,'FontSize',10);
                % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                set(gcf,'PaperPosition', [0 0 24 12]) 
                set(lgd,'Orientation','horizontal');
                hold off
                set(gca,'XTick', (0:1:length(nyears)))
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;

            end
        end
    end
    
end
% % % end-------------------- RMS perfomance by lowpass filter window
% (seasonal)




% start-------------------- rms_ lowpassed corrected get
    for folding=1:1
            nyears=[1:5];
            nc_varname_prefixes={'corrected_interped_sla'};
            nc_titlename_prefixes={'lowpass'};

    %         for nc_varnameij=1:length(nc_varname_prefixes)
    %             nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
    %             nc_varname=[nc_varname_prefix];
    %             eval(['corr=ncread(filename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(confilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_con(1,nc_varnameij)=mean(corr(:),'omitnan');
    %             eval(['corr=ncread(sodafilename,','''','corr_', nc_varname,'''',');']);
    %             corr=corr.*cmems_mask;
    %             mcorr_soda(1,nc_varnameij)=mean(corr(:),'omitnan');
    %         end
            for nyeari=1:length(nyears)
                nyear=nyears(nyeari);
                for nc_varnameij=1:length(nc_varname_prefixes)
                    nc_varname_prefix=nc_varname_prefixes{nc_varnameij};
                    nc_varname=[nc_varname_prefix, '_', num2str(nyear),'y_lowpass'];
    %                 jpgname=strcat(outfile, '_', testname,'_',regionname, '_corr_', nc_varname, '_',num2str(min(inputyear),'%04i'), ...
    %                     '_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    %                 if (exist(jpgname , 'file') ~= 2)



                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(1,nc_varnameij,:,:,:)=sq_diff_all;

                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_', '''', ...
                             ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname, '''', ',' '''', '\run\','''', ',', '''', testname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_sat(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Model\ROMS\nwp_1_10\', '''', ',', '''', testname2, '''', ',' '''', '\run\','''', ',', '''', testname2, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_con(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_rcm_glo(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
                        eval(['fname=[','''','E:\Data\Reanalysis\SODA\', '''', ',', '''', sodatestname, '''', ',' '''', '\','''', ',', '''', sodatestname, '''', ',', '''', '_', '''', ',', ...
                            '''', regionname, '''', ',', '''', '_', '''', ',', '''', 'rms_corrected_interped_sla_', '''', ',', '''', num2str(nyear), '''', ...
                            ',', '''', 'y_lowpass_', '''', ',', '''', num2str(min(inputyear)), '''', ',', '''', '_', '''', ',', '''', num2str(max(inputyear)), '''', ',', ...
                            '''', '.mat', '''', '];']); 
                        load(fname)
                        mrms_soda(nyeari+1,nc_varnameij)=mean(mean_rms(:), 'omitnan');
                        for tij=1:size(sq_diff_all,3)
                            sq_diff_all(:,:,tij)=sq_diff_all(:,:,tij).*cmems_mask;
                        end
                        sq_diff_all_gcm_soda(nyeari+1,nc_varnameij,:,:,:)=sq_diff_all;
                        
    %                     eval(['corr=ncread(confilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_con(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');
    %                     eval(['corr=ncread(sodafilename,','''','rms_', nc_varname,'''',');']);
    %                     corr=corr.*cmems_mask;
    %                     mcorr_soda(nyeari+1,nc_varnameij)=mean(corr(:),'omitnan');

                        lon_cmems=ncread(filename,'lon_cmems');
                        lat_cmems=ncread(filename,'lat_cmems');

    %                     eval(['mcorr=mean(corr_',nc_varname,'(:),', '''', 'omitnan', '''',');']);
                        tit_prefix=nc_titlename_prefixes{1};
                        titlename = strcat(tit_prefix,',',testname,',(',num2str(min(inputyear),'%04i'),'-', ...
                            num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mrms,2)));  

    %                     title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

                        set(gcf, 'PaperUnits', 'points');
                        set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                        set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    %                     saveas(gcf,jpgname,'tif');

                        disp(' ')
%                         disp(['corr_', nc_varname, ' plot is created.'])
                        disp(' ')
    %                     disp([' File path is : ',jpgname])
                        disp(' ')

                        hold off
                        close all;
    %                 end
                end
            end
    end
% end-------------------- rms_ lowpassed corrected get


% % % start-------------------- RMS perfomance by lowpass filter window
% (corrected)
    for folding=1:1
         jpgname=strcat(outfile, '_', testname, '_',regionname, 'corrected_rms_lp_perfomance', '.jpg'); %% ~_year_month.jpg
%             if (exist(jpgname , 'file') ~= 2)
            mslplot=plot(0:length(nyears),mrms_soda(:,1),'b')
            hold on
            msl2plot=plot(0:length(nyears),mrms_con(:,1),'g')
            msl3plot=plot(0:length(nyears),mrms(:,1),'r')

            xlabel('Low pass filter window(year)')
            ylabel('RMS')
            title([regionname, ', RMS performance '])
    %         datetick('x','yyyy','keepticks')
            axis tight;
            ylim([0 7])
            set(mslplot,'LineWidth',2);
            set(msl2plot,'LineWidth',2);
            set(msl3plot,'LineWidth',2);
            set(gca,'FontSize',15);
            grid on
            lgd=legend('RCM-SAT','RCM-GLO','GCM');
            set(lgd,'FontSize',10);
            % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(gcf,'PaperPosition', [0 0 24 12]) 
            set(lgd,'Orientation','horizontal');
            hold off
            set(gca,'XTick', (0:1:length(nyears)))
            saveas(gcf,jpgname,'jpg');
            grid off
%             end
        close all;
    end
    % % % end-------------------- RMS perfomance by lowpass filter window
    % (corrected)




