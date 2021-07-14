close all; clear all; clc;

dropboxpath='C:\Users\User\Dropbox';
% addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
% addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
addpath(genpath([dropboxpath '\source\matlab\Common\t_tide_v1.3beta']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
% addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
% addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path

inputyear=1989:2005;

% % set info of station
for folding=1:1
station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

num_sta=length(station.name);

station.name{1}='Anheung';
station.name{2}='Gunsan';
station.name{3}='Mokpo';
station.name{4}='Heuksando';
station.name{5}='Chujado';
station.name{6}='Wando';
station.name{7}='Yeosu';
station.name{8}='Tongyeong';
station.name{9}='Gadeokdo';
station.name{10}='Busan';
station.name{11}='Ieodo';
station.name{12}='Jeju';
station.name{13}='Seogwipo';
station.name{14}='Geomundo';
station.name{15}='Ulsan';
station.name{16}='Pohang';
station.name{17}='Mukho';
station.name{18}='Sokcho';
station.name{19}='Ulleungdo';

load('Z:\내 드라이브\MEPL\project\SSH\4th_year\09_backup\UST\99_재생산자료\01_매시별\tide_all.mat')

Obs.lat(1)=36.67;
Obs.lat(2)=35.97;
Obs.lat(3)=34.77;
Obs.lat(4)=34.68;
Obs.lat(5)=33.95;
Obs.lat(6)=34.30;
Obs.lat(7)=34.73;
Obs.lat(8)=34.82;
Obs.lat(9)=35.02;
Obs.lat(10)=35.08;
Obs.lat(11)=NaN;
Obs.lat(12)=33.52;
Obs.lat(13)=33.23;
Obs.lat(14)=34.02;
Obs.lat(15)=35.5;
Obs.lat(16)=36.03;
Obs.lat(17)=37.55;
Obs.lat(18)=38.2;
Obs.lat(19)=37.48;

end
% % get sum of 4 major tidal constitute (Observation)

% % get py data
for yeari=1:length(inputyear)
    for stai=1:num_sta
        if (stai~=11)
            data_1y_ind=find(Obs.data{stai}(:,1)==inputyear(yeari));
            Obs.data_1y{stai,yeari}=Obs.data{stai}(data_1y_ind,6);

            [tname,tfreq,tcon,tout]=t_tide(Obs.data_1y{stai,yeari},...
                   'interval',1, ...                     % hourly data
                   'start',datenum(inputyear(yeari),1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
                   'latitude',Obs.lat(stai),...               % Latitude of Obs
                   'rayleigh',1, 'error','wboot', 'output', 'none');                       % Use SNR=1 for synthesis. 

            Obs.tcon{stai,yeari}=tcon;

            tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio

            tide_info.name{1}='M2  ';
            tide_info.name{2}='S2  ';
            tide_info.name{3}='K1  ';
            tide_info.name{4}='O1  ';

            num_tide_all=size(tname,1);
            num_tide_tgt=length(tide_info.name);
            for coni=1:num_tide_all
                for tide_namei=1:num_tide_tgt
                    if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                        tide_info.index(tide_namei)=coni;
                    end
                end
            end

            Obs.amp_4con(stai,yeari)= sum(tcon(tide_info.index(:),1));
            Obs.amp_M2(stai,yeari)=tcon(tide_info.index(1),1);
            Obs.amp_S2(stai,yeari)=tcon(tide_info.index(2),1);
            Obs.amp_K1(stai,yeari)=tcon(tide_info.index(3),1);
            Obs.amp_O1(stai,yeari)=tcon(tide_info.index(4),1);

        else
            Obs.data_1y{stai,yeari}=NaN;
            Obs.amp_4con(stai,yeari)= NaN;
            Obs.amp_M2(stai,yeari)=NaN;
            Obs.amp_S2(stai,yeari)=NaN;
            Obs.amp_K1(stai,yeari)=NaN;
            Obs.amp_O1(stai,yeari)=NaN;
        end  
    end
    Obs.year(yeari)=inputyear(yeari);
end

all_testname2 = {'test53', 'test54', 'test55', 'test56'};

for testnameind2=1:length(all_testname2)
    
    % % %         flag configuration
    for folding=1:1
        fig_flags{1,1}='Mokpo M2 time series';

    end
    for flagi=1:8
        fig_flags{flagi,2}=0;
    end

    fig_flags{1,2}=2;
    fig_flags{2,2}=2;
    fig_flags{3,2}=2;
    fig_flags{4,2}=2;

    testname=all_testname2{testnameind2};
% % get py 
    for yeari=1:length(inputyear)

        filename=['J:\Data\Model\ROMS\nwp_1_20\', testname, '\run\',num2str(inputyear(yeari)),'\sta.nc'];
        Model.zeta=ncread(filename, 'zeta')*100;
        Model.lon=ncread(filename, 'lon_rho');
        Model.lat=ncread(filename, 'lat_rho');
        Model.ocean_time=ncread(filename, 'ocean_time');
        
        if length(Model.zeta)>=24*365
            for stai=1:num_sta    
                [tname,tfreq,tcon,tout]=t_tide(Model.zeta(stai,:),...
                       'interval',1, ...                     % hourly data
                       'start',datenum(inputyear(yeari),1,1,0,0,0), ...               % start time is datestr(tuk_time(1))
                       'latitude',Model.lat(stai),...               % Latitude of Model
                       'rayleigh',1, 'error','wboot', 'output', 'none');                       % Use SNR=1 for synthesis. 

               Model.tcon{stai,yeari}=tcon;

                tsnr=(tcon(:,1)./tcon(:,2)).^2;  % signal to noise ratio

                num_tide_all=size(tname,1);
                num_tide_tgt=length(tide_info.name);
                for coni=1:num_tide_all
                    for tide_namei=1:num_tide_tgt
                        if (strcmp(tide_info.name{tide_namei}, tname(coni,:))==1)
                            tide_info.index(tide_namei)=coni;
                        end
                    end
                end

                Model.amp_4con(stai,yeari)= sum(tcon(tide_info.index(:),1));
                Model.amp_M2(stai,yeari)=tcon(tide_info.index(1),1);
                Model.amp_S2(stai,yeari)=tcon(tide_info.index(2),1);
                Model.amp_K1(stai,yeari)=tcon(tide_info.index(3),1);
                Model.amp_O1(stai,yeari)=tcon(tide_info.index(4),1);
            end
        else
            for stai=1:num_sta    
                Model.tcon{stai,yeari}=NaN;
                Model.amp_4con(stai,yeari)=NaN;
                Model.amp_M2(stai,yeari)=NaN;
                Model.amp_S2(stai,yeari)=NaN;
                Model.amp_K1(stai,yeari)=NaN;
                Model.amp_O1(stai,yeari)=NaN;
            end
        end
    end
    
%     Obs.M2_diff=Obs.amp_M2-Obs.amp_M2;
%     Obs.S2_diff=Obs.amp_S2-Obs.amp_S2;
%     Obs.K1_diff=Obs.amp_K1-Obs.amp_K1;
%     Obs.O1_diff=Obs.amp_O1-Obs.amp_O1;
%     Model.M2_diff=Model_fy.amp_M2-Model_py.amp_M2;
%     Model.S2_diff=Model_fy.amp_S2-Model_py.amp_S2;
%     Model.K1_diff=Model_fy.amp_K1-Model_py.amp_K1;
%     Model.O1_diff=Model_fy.amp_O1-Model_py.amp_O1;

    figrawdir =strcat('Z:\내 드라이브\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\station\'); % % where figure files will be saved

    figdir=[figrawdir,'CLIM\'];
    if (exist(strcat(figdir) , 'dir') ~= 7)
        mkdir(strcat(figdir));
    end 
    outfile = strcat(figdir);
        
    fig_flag=fig_flags{1,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_Mokpo_M2_time_series',num2str(min(inputyear)), '-', num2str(max(inputyear)), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)
            figH = figure;
            axLH = gca;
            axRH = axes('color','none');
            mslplot{1}=plot(inputyear,Model.amp_M2(3,:), 'b','parent',axLH);
            mslplot{2}=plot(inputyear,Obs.amp_M2(3,:), 'k','parent',axRH);
            ylabel(axLH,'Model M2 Tidal amplitude (cm)')
            ylabel(axRH,'Obs M2 Tidal amplitude (cm)')
            ax_pos = get(axLH,'position');
            set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
            set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%             set(axRH,'color','none','yaxislocation','right');

%             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
            set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
            set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
            xlabel(axRH, 'Year');

            title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);
            grid on

            lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
            
%             lgd=legend('Model','TG-UST');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');

            set(gcf,'PaperPosition', [0 0 36 12]) 
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
    
    fig_flag=fig_flags{2,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_Mokpo_S2_time_series',num2str(min(inputyear)), '-', num2str(max(inputyear)), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)
            figH = figure;
            axLH = gca;
            axRH = axes('color','none');
            mslplot{1}=plot(inputyear,Model.amp_S2(3,:), 'b','parent',axLH);
            mslplot{2}=plot(inputyear,Obs.amp_S2(3,:), 'k','parent',axRH);
            ylabel(axLH,'Model S2 Tidal amplitude (cm)')
            ylabel(axRH,'Obs S2 Tidal amplitude (cm)')
            ax_pos = get(axLH,'position');
            set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
            set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%             set(axRH,'color','none','yaxislocation','right');

%             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
            set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
            set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
            xlabel(axRH, 'Year');

            title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);
            grid on

            lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
            
%             lgd=legend('Model','TG-UST');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');

            set(gcf,'PaperPosition', [0 0 36 12]) 
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
    
    fig_flag=fig_flags{3,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_Mokpo_K1_time_series',num2str(min(inputyear)), '-', num2str(max(inputyear)), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)
            figH = figure;
            axLH = gca;
            axRH = axes('color','none');
            mslplot{1}=plot(inputyear,Model.amp_K1(3,:), 'b','parent',axLH);
            mslplot{2}=plot(inputyear,Obs.amp_K1(3,:), 'k','parent',axRH);
            ylabel(axLH,'Model K1 Tidal amplitude (cm)')
            ylabel(axRH,'Obs K1 Tidal amplitude (cm)')
            ax_pos = get(axLH,'position');
            set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
            set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%             set(axRH,'color','none','yaxislocation','right');

%             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
            set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
            set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
            xlabel(axRH, 'Year');

            title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);
            grid on

            lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
            
%             lgd=legend('Model','TG-UST');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');

            set(gcf,'PaperPosition', [0 0 36 12]) 
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
    
    fig_flag=fig_flags{4,2};
    while (fig_flag)
        pngname=strcat(outfile, '_', testname,'_Mokpo_O1_time_series',num2str(min(inputyear)), '-', num2str(max(inputyear)), '.tif'); %% ~_year_month.jpg
        if (exist(pngname , 'file') ~= 2 || fig_flag==2)
            figH = figure;
            axLH = gca;
            axRH = axes('color','none');
            mslplot{1}=plot(inputyear,Model.amp_O1(3,:), 'b','parent',axLH);
            mslplot{2}=plot(inputyear,Obs.amp_O1(3,:), 'k','parent',axRH);
            ylabel(axLH,'Model O1 Tidal amplitude (cm)')
            ylabel(axRH,'Obs O1 Tidal amplitude (cm)')
            ax_pos = get(axLH,'position');
            set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
            set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%             set(axRH,'color','none','yaxislocation','right');

%             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
            set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
            set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
            set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
            xlabel(axRH, 'Year');

            title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
            % datetick('x','yyyy','keepticks')
            axis tight;
            % ylim(meanplotlev2)
            set(mslplot{1},'LineWidth',2);
            set(mslplot{2},'LineWidth',2);
            grid on

            lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
            
%             lgd=legend('Model','TG-UST');
            set(lgd,'FontSize',15);
            set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
            set(lgd,'Orientation','horizontal');

            set(gcf,'PaperPosition', [0 0 36 12]) 
            saveas(gcf,pngname,'tif');
            grid off
            close all;
        end
        fig_flag=0;
    end
    
 

end