close all; clear all;  clc;
warning off;

all_testname2 = {'test53', 'test54', 'test55', 'test56', 'ens03'};
% all_testname2 = {'test53'};
  
all_region2 ={'AKP4'}

% all_region2 ={'YS'};

% all_var2 = {'SST', 'SSS', 'SSH', 'BT'};
all_var2 = {'SSH'};

% all_region2 ={'NWP'}
for testnameind2=1:length(all_testname2)
        close all;
%         clearvars '*' -except regionind2 testnameind2 all_region2 all_testname2 all_var2
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
            addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Roms_tools\Preprocessing_tools']));
            addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
        elseif (strcmp(system_name,'GLNXA64'))
            dropboxpath='/home/kimyy/Dropbox';
            addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
            addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
            addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
            addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
        end

%         shadlev = [0 35];
%         rms_shadlev = [0 4];
%     %     trendlev = [-3 3];  %% trend lev
%         trendlev = [-10 10];  %% trend lev
%         abstrendlev =[4 7];
%         reltrendlev =[-5 5];
%         conlev  = 0:5:35;
        meanplotlev =[-0.2 0.2];
        meanplotlev2 =[-0.4 0.4];
%         trendplotlev = [3 7];
%         sshlev =[-0.7 1.3];
%         sshdifflev = [40 70];

        % for snu_desktopd
        testname=all_testname2{testnameind2}    % % need to change
        inputyear = [1977:2005]; % % put year which you want to plot [year year ...]
        inputmonth = [1 2 3 4 5 6 7 8 9 10 11 12]; % % put month which you want to plot [month month ...]

        varname ='zeta'
%         varname ='zos';

        run('nwp_polygon_point.m');
%         regionname=all_region2{regionind2};
        regionname='AKP4';

        for folding=1:1
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
                case('ECS2') %% East China Sea2
                    refpolygon=ecs2polygon;
                case('YSECS') %% YS & East China Sea
                    refpolygon=ysecspolygon;
                case('AKP') %% Around Korea Peninsula
                    refpolygon=akppolygon;
                case('AKP2') %% Around Korea Peninsula
                    refpolygon=akp2polygon;
                case('AKP3') %% Around Korea Peninsula
                    refpolygon=akp3polygon;
                case('AKP4') %% Around Korea Peninsula
                    refpolygon=akp4polygon;
                case('CA') %% Coastal Area around korea peninsula
                    refpolygon=capolygon;
                case('EKB') %% Coastal Area around korea peninsula
                    refpolygon=ekbpolygon;
                case('BOH') %% Coastal Area around korea peninsula
                    refpolygon=bohpolygon;
                case('TEST') %% for debugging
                    refpolygon=testpolygon;
                otherwise
                    ('?')
            end
                lonlat(1)=min(refpolygon(:,1));
                lonlat(2)=max(refpolygon(:,1));
                lonlat(3)=min(refpolygon(:,2));
                lonlat(4)=max(refpolygon(:,2));
        end
        
        % % %         flag configuration
        for folding=1:1
            fig_flags{1,1}='msl time series (seasonal filtered), (trend, corr)';
            fig_flags{2,1}='msl time series (trend, corr)';
            fig_flags{3,1}='tidal station map';
            fig_flags{4,1}='moving averaged msl time series (trend, corr)';


            for flagi=1:5
                fig_flags{flagi,2}=0;
            end
            fig_flags{1,2}=0;
            fig_flags{2,2}=0;
            fig_flags{3,2}=2;
            fig_flags{4,2}=2;
            fig_flags{5,2}=2;
        end
        
        
        % % % for EKB
        % regionname='EKB';
        % lonlat = [127, 129.5, 38, 40.5];

%         load(['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%             'ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);
%         filename = ['G:\Data\Model\ROMS\nwp_1_20\',testname,'\run\',testname,'_',regionname, ...
%                     '_ssh_analysis_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.nc'];

        valnum=0;
        run('C:\Users\kyy\Dropbox\source\matlab\Common\Figure\bwr_map')  % % set colormap (blue and white)
        wrmap = bwrmap(51:100,:);
        
        scenname='historical';
        system_name=computer;
        if (strcmp(system_name,'PCWIN64'))
            % % for windows
            figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\SSH\5th_year\figure\nwp_1_20\',testname,'\',regionname,'\'); % % where figure files will be saved
            param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_20\run\fig_param\fig_param_kyy_EKB_RMS.m'
            filedir = strcat('H:\Data\Model\ROMS\nwp_1_20\', testname, '\run\'); % % where data files are
            cmemsdir='E:\Data\Observation\CMEMS\';
        elseif (strcmp(system_name,'GLNXA64'))
        end
        
        run(param_script);
        
        filename = strcat(filedir, testname,'_',regionname, '_ssh_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.nc');

        tidal_st= [ 117.717, 39.000; ...  % TANGGU (Tianjin)
                    119.450, 34.750; ...  % LIANYUNGANG
                    119.550, 35.383; ...  % SHIJIUSUO (RIzhao)
                    119.600, 39.900; ...  % QINHUANGDAO
                    121.383, 37.533; ...  % YANTAI
                    121.617, 32.133; ...  % LUSI (Qidong, Jiangsu)
                    121.683, 38.867; ...  % DALIAN
                    126.592, 37.452; ...  % INCHEON
                    129.036, 35.096; ...  % BUSAN
                    130.408, 33.619; ...  % HAKATA (Fukuoka)
                    130.914, 37.491; ...  % ULLEUNG
                    140.858, 43.209];      % OSHORO II (Sapporo)
        
                
        tidal_name = { 'TANGGU (Tianjin)', ...
                       'LIANYUNGANG', ...
                       'SHIJIUSUO_(RIzhao)', ...
                       'QINHUANGDAO', ...
                       'YANTAI', ...
                       'LUSI_(Qidong)', ...
                       'DALIAN', ...
                       'INCHEON', ...
                       'BUSAN', ...
                       'HAKATA (Fukuoka)', ...
                       'ULLEUNG', ...
                       'OSHORO_II_(Sapporo)'};
                   
        PSMSL_dir = 'E:\Data\Observation\PSMSL';
        PSMSL_filename = { [PSMSL_dir, '\', '1403_TANGGU.txt'], ...
                           [PSMSL_dir, '\', '1405_LIANYUNGANG.txt'], ...
                           [PSMSL_dir, '\', '1404_SHIJIUSUO.txt'], ...
                           [PSMSL_dir, '\', '0614_QINHUANGDAO.txt'], ...
                           [PSMSL_dir, '\', '0731_YANTAI.txt'], ...
                           [PSMSL_dir, '\', '0979_LUSI.txt'], ...
                           [PSMSL_dir, '\', '0723_DALIAN.txt'], ...
                           [PSMSL_dir, '\', '0956_INCHEON.txt'], ...
                           [PSMSL_dir, '\', '0955_BUSAN.txt'], ...
                           [PSMSL_dir, '\', '1094_HAKATA.txt'], ...
                           [PSMSL_dir, '\', '1490_ULLEUNG.txt'], ...
                           [PSMSL_dir, '\', '1027_OSHORO_II.txt']};
        
        
        rlrfilename = 'E:\Data\Observation\PSMSL\rlr_monthly\filelist.txt';
        formatSpec = '%6s%12s%13s%42s%5s%5s%s%[^\n\r]';
        fileID = fopen(rlrfilename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
        dataArray{4} = strtrim(dataArray{4});
        dataArray{7} = strtrim(dataArray{7});
        fclose(fileID);
        raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
        for col=1:length(dataArray)-1
            raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
        end
        numericData = NaN(size(dataArray{1},1),size(dataArray,2));
        for col=[1,2,3,5,6]
            % 입력 셀형 배열의 텍스트를 숫자로 변환합니다. 숫자형이 아닌 텍스트를 NaN으로 바꿨습니다.
            rawData = dataArray{col};
            for row=1:size(rawData, 1)
                % 숫자형이 아닌 접두사 및 접미사를 검색하고 제거하는 정규 표현식을 만듭니다.
                regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                try
                    result = regexp(rawData(row), regexstr, 'names');
                    numbers = result.numbers;

                    % 천 단위가 아닌 위치에서 쉼표를 검색했습니다.
                    invalidThousandsSeparator = false;
                    if numbers.contains(',')
                        thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                        if isempty(regexp(numbers, thousandsRegExp, 'once'))
                            numbers = NaN;
                            invalidThousandsSeparator = true;
                        end
                    end
                    % 숫자형 텍스트를 숫자로 변환합니다.
                    if ~invalidThousandsSeparator
                        numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                        numericData(row, col) = numbers{1};
                        raw{row, col} = numbers{1};
                    end
                catch
                    raw{row, col} = rawData{row};
                end
            end
        end
        rawNumericColumns = raw(:, [1,2,3,5,6]);
        rawStringColumns = string(raw(:, [4,7]));
        filelist = raw;
        clearvars rlrfilename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns;               
                       
        listlength=size(filelist,1);
        for listi=1:listlength
            sta_file(listi)=filelist{listi,1};
            lon_file(listi)=filelist{listi,3};
            lat_file(listi)=filelist{listi,2};
        end
        
        valid_sta = double(inpolygon(lon_file,lat_file,refpolygon(:,1),refpolygon(:,2)));
        valid_sta_num=sta_file(find(valid_sta==1));
        valid_sta_ind=find(valid_sta==1);
        valid_listlength=length(valid_sta_ind);
        
        for val_listi=1:valid_listlength
            valid_tempname{val_listi}=filelist{valid_sta_ind(val_listi),4};
            temp_staname=char(valid_tempname{val_listi});
            valid_name2{val_listi}=strtrim(temp_staname(1:end-2));
            PSMSL_lon2(val_listi)=lon_file(valid_sta_ind(val_listi));
            PSMSL_lat2(val_listi)=lat_file(valid_sta_ind(val_listi));
            
            PSMSL_data_dir=['E:\Data\Observation\PSMSL\rlr_monthly\data'];
            rlrfilename = [PSMSL_data_dir, '\', num2str(valid_sta_num(val_listi)), '.rlrdata'];
            formatSpec = '%12s%7s%s%[^\n\r]';
            fileID = fopen(rlrfilename,'r');
            dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string',  'ReturnOnError', false);
            dataArray{3} = strtrim(dataArray{3});
            fclose(fileID);
            raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
            for col=1:length(dataArray)-1
                raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
            end
            numericData = NaN(size(dataArray{1},1),size(dataArray,2));
            for col=[1,2]
                % 입력 셀형 배열의 텍스트를 숫자로 변환합니다. 숫자형이 아닌 텍스트를 NaN으로 바꿨습니다.
                rawData = dataArray{col};
                for row=1:size(rawData, 1)
                    % 숫자형이 아닌 접두사 및 접미사를 검색하고 제거하는 정규 표현식을 만듭니다.
                    regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                    try
                        result = regexp(rawData(row), regexstr, 'names');
                        numbers = result.numbers;

                        % 천 단위가 아닌 위치에서 쉼표를 검색했습니다.
                        invalidThousandsSeparator = false;
                        if numbers.contains(',')
                            thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                            if isempty(regexp(numbers, thousandsRegExp, 'once'))
                                numbers = NaN;
                                invalidThousandsSeparator = true;
                            end
                        end
                        % 숫자형 텍스트를 숫자로 변환합니다.
                        if ~invalidThousandsSeparator
                            numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                            numericData(row, col) = numbers{1};
                            raw{row, col} = numbers{1};
                        end
                    catch
                        raw{row, col} = rawData{row};
                    end
                end
            end
            rawNumericColumns = raw(:, [1,2]);
            rawStringColumns = string(raw(:, 3));
            R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); 
            rawNumericColumns(R) = {NaN}; 
            PSMSL_raw{val_listi} = raw;
            clearvars rlrfilename formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp rawNumericColumns rawStringColumns R;
        end
        
        load('E:\Data\Observation\PSMSL\rlr_monthly\PSMSL_time.mat')
        if min(inputyear)>1977
            PSMSL_time=PSMSL_time((inputyear(1)-1977)*12+1:348);
        end

        for ind_sta=1:length(PSMSL_raw)
            for PSMSL_tind=1:length(PSMSL_time)
%             find(PSMSL_raw{5}(:,1)==2.005875000000000e+03)
                tind=find(cell2mat(PSMSL_raw{ind_sta}(:,1))==PSMSL_time(PSMSL_tind));
                PSMSL_obs2{ind_sta}(PSMSL_tind,1)=floor(PSMSL_time(PSMSL_tind));
                PSMSL_month=mod(PSMSL_tind-1,12)+1;
                PSMSL_obs2{ind_sta}(PSMSL_tind,2)=PSMSL_month;
                if isempty(tind)
                    PSMSL_obs2{ind_sta}(PSMSL_tind,3)=NaN;
                else
                    if cell2mat(PSMSL_raw{ind_sta}(tind,2))==-99999
                        PSMSL_obs2{ind_sta}(PSMSL_tind,3)=NaN;
                    else
                        PSMSL_obs2{ind_sta}(PSMSL_tind,3)=cell2mat(PSMSL_raw{ind_sta}(tind,2)) ./ 1000; % mm -> m
                    end
                end
            end
            if(sum(isfinite(PSMSL_obs2{ind_sta}(:,3)))>=12*12)
                time_flag(ind_sta)=1;
            else
                time_flag(ind_sta)=0;
            end
        end
        
        ind_time_flag=find(time_flag==1);
        
        for ind_sta=1:sum(time_flag)
            PSMSL_obs{ind_sta}=PSMSL_obs2{ind_time_flag(ind_sta)};
            PSMSL_lon(ind_sta)=PSMSL_lon2(ind_time_flag(ind_sta));
            PSMSL_lat(ind_sta)=PSMSL_lat2(ind_time_flag(ind_sta));
            valid_name{ind_sta}=valid_name2{(ind_time_flag(ind_sta))};
        end
        
        m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
        m_gshhs_i('color',m_gshhs_line_color);
        m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
        m_line(PSMSL_lon(:),PSMSL_lat(:),'marker','o','color','r','linewi',2,...
          'linest','none','markersize',8,'markerfacecolor','r');
        m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

        
        info=ncinfo(filename);
        lon_rho=ncread(filename, 'lon_rho');
        lat_rho=ncread(filename, 'lat_rho');
        raw_ssh=ncread(filename, 'raw_ssh');
        
        model_mask=ones(size(lon_rho));
        model_mask(isnan(raw_ssh(:,:,1)))=NaN;
        
        for ind_sta=1:length(PSMSL_obs)
            for i=1:size(lon_rho,1)
                for j=1:size(lat_rho,2)
                    dist(i,j)=m_lldist([lon_rho(i,j), PSMSL_lon(ind_sta)], [lat_rho(i,j), PSMSL_lat(ind_sta)]) * model_mask(i,j);
                end
            end
            ind_sta_model(ind_sta)=find(dist(:)==min(dist(:)));
            ind_sta_model_lon(ind_sta)=mod(ind_sta_model(ind_sta),size(lon_rho,1));
            ind_sta_model_lat(ind_sta)=floor(ind_sta_model(ind_sta)/size(lon_rho,1))+1;
            model_msl(:,ind_sta)=raw_ssh(ind_sta_model_lon(ind_sta),ind_sta_model_lat(ind_sta),:);
        end
%         lon_rho(ind_sta_model(1))
%         lat_rho(ind_sta_model(1))
%         lon_rho(ind_sta_model_lon(1), ind_sta_model_lat(1))
%         lat_rho(ind_sta_model_lon(1), ind_sta_model_lat(1))
        clim_model_msl=reshape(model_msl,[12, size(model_msl,1)/12, length(PSMSL_obs)]);
        clim_mean_msl=squeeze(mean(clim_model_msl,2));
        for ind_sta=1:length(PSMSL_obs)
            clim_obs_msl(:,:,ind_sta)=reshape(PSMSL_obs{ind_sta}(:,3),[12, size(model_msl,1)/12]);
        end
        clim_mean_obs_msl=squeeze(mean(clim_obs_msl,2,'omitnan'));
        
        trendtime=inputyear(1):1/12:inputyear(end)+1-1/12;
        trend(1:ind_sta)=NaN;
        for ind_sta=1:length(PSMSL_obs)
            for ind_time=1:size(model_msl,1)
%                 mod(ind_time-1,12)+1
                clim_filtered_msl(ind_time,ind_sta)=model_msl(ind_time,ind_sta)-clim_mean_msl(mod(ind_time-1,12)+1,ind_sta);
                clim_filtered_obs(ind_time,ind_sta)=PSMSL_obs{ind_sta}(ind_time,3)- clim_mean_obs_msl(mod(ind_time-1,12)+1,ind_sta);
                PSMSL_msl(ind_time,ind_sta)=PSMSL_obs{ind_sta}(ind_time,3);
            end
            p=polyfit(trendtime,squeeze(clim_filtered_msl(:,ind_sta))',1);
            trend(ind_sta)=p(1) * 1000.0;
        end
%         plot(squeeze(clim_filtered_msl(:,1))*100)
        
        figdir=[figrawdir,'Tidal\'];
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end 
        
        for i =1:length(inputyear) 
            tempyear=inputyear(i);
            for month=1:12
                xData((12*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month,'%02i'),'-01',]);
            end
        end
        
       
% % %         msl time series (seasonal filtered)  
        fig_flag=fig_flags{1,2};
        while (fig_flag)
            for ind_sta=1:length(PSMSL_obs)
                outfile = strcat(figdir,valid_name{ind_sta});
                jpgname=strcat(outfile, '_', testname, '_msl_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    %                 for varind=1:length(inputyear)*12
    %                     msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
    %                 end
                    msl_filt=squeeze(clim_filtered_msl(:,ind_sta))';
                    msl_filt_obs=squeeze(clim_filtered_obs(:,ind_sta))';

            %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                    msl_filt=msl_filt-mean(msl_filt); 
                    msl_filt_obs=msl_filt_obs-mean(msl_filt_obs,'omitnan');
                    
                    p=polyfit(xData,msl_filt,1);
                    msl2=xData*p(1)+p(2);
                    p2=polyfit(xData(find(isfinite(msl_filt_obs)==1)),msl_filt_obs(find(isfinite(msl_filt_obs)==1)),1);
                    msl3=xData*p2(1)+p2(2);
                    
                    mslplot=plot(xData,msl_filt,'k');
                    hold on
                    mslplot3=plot(xData,msl_filt_obs,'b');
                    mslplot2=plot(xData,msl2,'Color','r');
                    xlabel('year')
                    ylabel('Mean SSH (m)')
                    title([tidal_name{ind_sta}, ', Mean SSH(',num2str(min(inputyear),'%04i'), ...
                        '-',num2str(max(inputyear),'%04i'),'), ',num2str(round(trend(ind_sta),2)), ' mm/y'])
                    datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev)
                    msl_corr=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
                    txt1=text(xData(end-40), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(msl_corr(1,2),2))], 'FontSize', 20); 

                    set(gcf,'PaperPosition', [0 0 36 12]) 

                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',20);
                    
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                    pause(2)
                end
            end
            fig_flag=0;
        end
        
% % %         msl time series 
        fig_flag=fig_flags{2,2};
        while (fig_flag)
            for ind_sta=1:length(PSMSL_obs)
                outfile = strcat(figdir,tidal_name{ind_sta});
                jpgname=strcat(outfile, '_', testname, '_msl_nonfilt_', ...
                    num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
                if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    %                 for varind=1:length(inputyear)*12
    %                     msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
    %                 end
                    msl_filt=squeeze(model_msl(:,ind_sta))';
                    msl_filt_obs=squeeze(PSMSL_msl(:,ind_sta))';

            %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                    msl_filt=msl_filt-mean(msl_filt); 
                    msl_filt_obs=msl_filt_obs-mean(msl_filt_obs, 'omitnan');
                    
                    p=polyfit(xData,msl_filt,1);
                    msl2=xData*p(1)+p(2);
                    mslplot=plot(xData,msl_filt,'k');
                    hold on
                    mslplot3=plot(xData,msl_filt_obs,'b');
                    mslplot2=plot(xData,msl2,'Color','r');
                    xlabel('year')
                    ylabel('Mean SSH (m)')
                    title([tidal_name{ind_sta}, ', Mean SSH(',num2str(min(inputyear),'%04i'), ...
                        '-',num2str(max(inputyear),'%04i'),'), ',num2str(round(trend(ind_sta),2)), ' mm/y'])
                    datetick('x','yyyy','keepticks')
                    axis tight;
                    ylim(meanplotlev2)
%                     msl_corr=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
%                     txt1=text(xData(end-40), min(meanplotlev2)+diff(meanplotlev2)/16.0 ,['R = ', num2str(round(msl_corr(1,2),2))], 'FontSize', 20); 

                    set(gcf,'PaperPosition', [0 0 36 12]) 

                    set(mslplot,'LineWidth',2);
                    set(gca,'FontSize',20);
                    
                    grid on
                    hold off
                    saveas(gcf,jpgname,'jpg');
                    grid off
                    close all;
                    pause(1)
                end
            end
            fig_flag=0;
        end

% % %         msl trends (seasonal filtered)  
        fig_flag=fig_flags{3,2};
        while (fig_flag)
            outfile = strcat(figdir);
            jpgname=strcat(outfile, testname, '_trends_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
%                 for varind=1:length(inputyear)*12
%                     msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%                 end
                
                for ind_sta=1:length(PSMSL_obs)
                    msl_filt=squeeze(clim_filtered_msl(:,ind_sta))'*1000.0;
                    msl_filt_obs=squeeze(clim_filtered_obs(:,ind_sta))'*1000.0;

            %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                    msl_filt=msl_filt-mean(msl_filt); 
                    msl_filt_obs=msl_filt_obs-mean(msl_filt_obs,'omitnan');

                    p=polyfit(trendtime(find(isfinite(msl_filt_obs)==1)),msl_filt(find(isfinite(msl_filt_obs)==1)),1);
                    msl2=xData*p(1)+p(2);
                    p2=polyfit(trendtime(find(isfinite(msl_filt_obs)==1)),msl_filt_obs(find(isfinite(msl_filt_obs)==1)),1);
                    msl3=xData*p2(1)+p2(2);

                    trend(ind_sta)=p(1);
                    trend2(ind_sta)=p2(1);
                end
                mslplot=plot(trend,'k');
                hold on
%                 mslplot3=plot(xData,msl_filt_obs,'b');
                mslplot2=plot(trend2,'Color','r');
                xlabel('Station#')
                ylabel('Trends (mm/yr)')
                title([testname, ', SSH trends(',num2str(min(inputyear),'%04i'), ...
                    '-',num2str(max(inputyear),'%04i'),'), '])
%                 datetick('x','yyyy','keepticks')
                axis tight;
                ylim([-1 8])
                
                set(gcf,'PaperPosition', [0 0 36 12]) 

                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);

                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
                pause(1)
            end
            
            
            fig_flag=0;
        end
        
        % % %         msl corr (seasonal filtered)  
        fig_flag=fig_flags{4,2};
        while (fig_flag)
            outfile = strcat(figdir);
            jpgname=strcat(outfile, testname, '_corr_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
%                 for varind=1:length(inputyear)*12
%                     msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%                 end
                
                for ind_sta=1:length(PSMSL_obs)
                    msl_filt=squeeze(clim_filtered_msl(:,ind_sta))';
                    msl_filt_obs=squeeze(clim_filtered_obs(:,ind_sta))';

            %         msl=squeeze(mean(mean(comb_data_filtered,1,'omitnan'),2,'omitnan'));
                    msl_filt=msl_filt-mean(msl_filt); 
                    msl_filt_obs=msl_filt_obs-mean(msl_filt_obs,'omitnan');
                    
                    msl_corr_temp=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
                    msl_corr(ind_sta)=msl_corr_temp(1,2);

                    p=polyfit(xData,msl_filt,1);
                    msl2=xData*p(1)+p(2);
                    p2=polyfit(xData(find(isfinite(msl_filt_obs)==1)),msl_filt_obs(find(isfinite(msl_filt_obs)==1)),1);
                    msl3=xData*p2(1)+p2(2);
                end
                mslplot=plot(msl_corr,'k');
                hold on
%                 mslplot3=plot(xData,msl_filt_obs,'b');
                xlabel('Station#')
                ylabel('Corr coef')
                title([testname, ', SSH corr(',num2str(min(inputyear),'%04i'), ...
                    '-',num2str(max(inputyear),'%04i'),'), '])
%                 datetick('x','yyyy','keepticks')
                axis tight;
                ylim([-0.1 0.7])
%                 msl_corr=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
%                 txt1=text(xData(end-40), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(msl_corr(1,2),2))], 'FontSize', 20); 

                set(gcf,'PaperPosition', [0 0 36 12]) 

                set(mslplot,'LineWidth',2);
                set(gca,'FontSize',20);

                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
                pause(1)
            end
            
            fig_flag=0;
        end
        
        
         % % %         all msl time series (seasonal filtered)  
        fig_flag=fig_flags{5,2};
        while (fig_flag)
            outfile = strcat(figdir);
            jpgname=strcat(outfile, testname, '_all_msl_', ...
                num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
            if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
%                 for varind=1:length(inputyear)*12
%                     msl_filt(varind)=squeeze(mean(mean(ncread(filename,'ssh_filtered',[1 1 varind], [inf inf 1]),1,'omitnan'),2,'omitnan'));
%                 end 
                hold on
                clear mslplot mslplot2
                for ind_sta=1:length(PSMSL_obs)
                    msl_filt=squeeze(clim_filtered_msl(:,ind_sta))';
                    msl_filt_obs=squeeze(clim_filtered_obs(:,ind_sta))';
                    

                    msl_filt=msl_filt-mean(msl_filt); 
                    msl_filt_obs=msl_filt_obs-mean(msl_filt_obs,'omitnan');
                    
                    cmap = get(groot,'defaultaxescolororder');
                    cmap_b = rgb2hsv(cmap);
                    cmap_b(:,2) = cmap_b(:,2)*.3;
                    cmap_b(:,3) = cmap_b(:,3)*.3+.7;
                    cmap_b = hsv2rgb(cmap_b);
                    
                    mslplot{ind_sta}=plot(xData,msl_filt,'Color', cmap(1,:));
                    mslplot2{ind_sta}=plot(xData,msl_filt_obs,'Color', cmap_b(1,:));
                    
                    
                

%                     msl_corr_temp=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
%                     msl_corr(ind_sta)=msl_corr_temp(1,2);
% 
%                     p=polyfit(xData,msl_filt,1);
%                     msl2=xData*p(1)+p(2);
%                     p2=polyfit(xData(find(isfinite(msl_filt_obs)==1)),msl_filt_obs(find(isfinite(msl_filt_obs)==1)),1);
%                     msl3=xData*p2(1)+p2(2);
                end
                mslplot{ind_sta+1}=plot(xData,mean(clim_filtered_msl,2)','Color', cmap(1,:));
                mslplot2{ind_sta+1}=plot(xData,mean(clim_filtered_obs,2)','Color', cmap_b(1,:));
                set(mslplot{ind_sta+1},'LineWidth',2);
                set(mslplot2{ind_sta+1},'LineWidth',2);
%                 mslplot3=plot(xData,msl_filt_obs,'b');
                xlabel('year')
                ylabel('Mean SSH (m)')
                title([testname, ', Mean SSH(',num2str(min(inputyear),'%04i'), ...
                        '-',num2str(max(inputyear),'%04i'),'), '])
                datetick('x','yyyy','keepticks')
                axis tight;

                ylim(meanplotlev)
%                 msl_corr=corrcoef(msl_filt(isfinite(msl_filt_obs)), msl_filt_obs(isfinite(msl_filt_obs)));
%                 txt1=text(xData(end-40), min(meanplotlev)+diff(meanplotlev)/16.0 ,['R = ', num2str(round(msl_corr(1,2),2))], 'FontSize', 20); 

                set(gcf,'PaperPosition', [0 0 36 12]) 

                set(gca,'FontSize',20);

                grid on
                hold off
                saveas(gcf,jpgname,'jpg');
                grid off
                close all;
                pause(1)
                
                

            end
            
            fig_flag=0;
        end


end
