if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                run(param_script);
                clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
                    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask mean_ts_aoi
                for yearij = 1:length(allyear)
                    tempyear = allyear(yearij);
                    for monthij = 1:length(inputmonth)
                        tempmonth = inputmonth(monthij);
                        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
                        if(exist('curl_mask')==0)
                            lon_rho = ncread(ncname, 'lon_rho');
                            lat_rho = ncread(ncname, 'lat_rho');
                            vel_polygon=[128, 37.5; 128, 39; 130, 39; 130, 37.5];
                            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
                        end
                        ocean_time=ncread(ncname, 'time')+datenum(1900,12,31);
                        u_rho=ncread(ncname, 'u_rho').*vel_mask;
                        u_rho(u_rho==0)=NaN;
                        v_rho=ncread(ncname, 'v_rho').*vel_mask;
                        v_rho(v_rho==0)=NaN;
                        uwind=ncread(ncname, 'uwind').*vel_mask;
                        uwind(uwind==0)=NaN;
                        vwind=ncread(ncname, 'vwind').*vel_mask;
                        vwind(vwind==0)=NaN;
                        egg_mask=ncread(ncname, 'mask_15day').*vel_mask;
                        egg_mask(egg_mask==0)=NaN;
                        
                        lastday_m=size(u_rho,3);
                        if (exist('comb_u_rho')==0)
                            comb_u_rho=u_rho;
                            comb_v_rho=v_rho;
                            comb_uwind=uwind;
                            comb_vwind=vwind;
                            comb_egg_mask=egg_mask;
                            comb_ocean_time=ocean_time;
                        else
                            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
                            comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
                            comb_uwind(:,:,end+1:end+lastday_m)=uwind;
                            comb_vwind(:,:,end+1:end+lastday_m)=vwind;
                            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
                            comb_ocean_time(end+1:end+lastday_m)=ocean_time;
                        end
                    end
                end
                
                ts_u_rho=reshape(comb_u_rho,[size(comb_u_rho,1)*size(comb_u_rho,2), size(comb_u_rho,3)]);
                mean_ts_u_rho=mean(ts_u_rho,1,'omitnan');
                ts_v_rho=reshape(comb_v_rho,[size(comb_v_rho,1)*size(comb_v_rho,2), size(comb_v_rho,3)]);
                mean_ts_v_rho=mean(ts_v_rho,1,'omitnan');
                ts_uwind=reshape(comb_uwind,[size(comb_uwind,1)*size(comb_uwind,2), size(comb_uwind,3)]);
                mean_ts_uwind=mean(ts_uwind,1,'omitnan');
                ts_vwind=reshape(comb_vwind,[size(comb_vwind,1)*size(comb_vwind,2), size(comb_vwind,3)]);
                mean_ts_vwind=mean(ts_vwind,1,'omitnan');
                ts_egg_mask=reshape(comb_egg_mask,[size(comb_egg_mask,1)*size(comb_egg_mask,2), size(comb_egg_mask,3)]);
                sum_ts_egg_mask=sum(ts_egg_mask,1,'omitnan');
                
                mean_ts_nwv= mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45);

                filename = 'Z:\내 드라이브\Data\Observation\NOAA\AOI\AOI_Index.txt';
                startRow = 1;
                formatSpec = '%5s%5s%s%[^\n\r]';
                fileID = fopen(filename,'r');
                dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
                fclose(fileID);
                raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
                for col=1:length(dataArray)-1
                    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
                end
                numericData = NaN(size(dataArray{1},1),size(dataArray,2));
                for col=[1,2,3]
                    rawData = dataArray{col};
                    for row=1:size(rawData, 1)
                        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
                        try
                            result = regexp(rawData(row), regexstr, 'names');
                            numbers = result.numbers;
                            invalidThousandsSeparator = false;
                            if numbers.contains(',')
                                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                                    numbers = NaN;
                                    invalidThousandsSeparator = true;
                                end
                            end
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
                AOI_year = cell2mat(raw(:, 1));
                AOI_month = cell2mat(raw(:, 2));
                AOI_value = cell2mat(raw(:, 3));
                clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;
                
% %                 10 year lowpass filter
                AOI_value=lowpass(AOI_value, 1/120, 1, 'Steepness', 0.99);  % phase correction is applied automatically, but butterworth isn't.

%                 bar(AOI_year+AOI_month/12, AOI_value_filt)
%                 plot(AOI_year+AOI_month/12, AOI_value)
%                 plot(AOI_year(1:end-27)+AOI_month(1:end-27)/12, AOI_value_filt(28:end))
                

                AOI_year=reshape(AOI_year, [12, length(AOI_year)/12]);
                AOI_month=reshape(AOI_month, [12, length(AOI_month)/12]);
                AOI_value=reshape(AOI_value, [12, length(AOI_value)/12]);
                AOI_year_start=find(AOI_year(1,:)==min(allyear));
                AOI_year_end=find(AOI_year(1,:)==max(allyear));
                
%                 plot(AOI_value(:,AOI_year_start:AOI_year_end))
                cal_arr=1;
                for cali=1:length(allyear)
                    if cali==length(allyear)
                        for monthij=2:length(inputmonth)
                            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
                        end
                    else
                        for monthij=1:length(inputmonth)
                            cal_arr=[cal_arr, cal_arr(end)+sum(eomday(allyear(cali),inputmonth(monthij)))];
                        end
                    end
                end
                
                mean_ts_aoi=NaN(1,length(mean_ts_nwv));
                cal_ind=1;
                for cali=1:length(allyear)-1
                    for monthij=1:length(inputmonth)
                        mean_ts_aoi(cal_arr(cal_ind):cal_arr(cal_ind+1)-1)=AOI_value(inputmonth(monthij),AOI_year_start+cali-1);
                        cal_ind=cal_ind+1;
                    end
                end
                cali=length(allyear);
                for monthij=1:length(inputmonth)-1
                    mean_ts_aoi(cal_arr(cal_ind):cal_arr(cal_ind+1)-1)=AOI_value(inputmonth(monthij),AOI_year_start+cali-1);
                    cal_ind=cal_ind+1;
                end
                monthij=length(inputmonth);
                mean_ts_aoi(cal_arr(cal_ind):end)=AOI_value(inputmonth(monthij),AOI_year_start+cali-1);

                
                
%                 mean_ts_aoi=AOI_value(inputmonth,AOI_year_start:AOI_year_end);
                
                axLH = gca;
                axRH = axes('color','none');
                mslplot{1}=plot(-mean_ts_uwind,'m','parent',axLH);
%                 mslplot{2}=bar(mean_ts_aoi, 'k*','parent',axRH, 'Markersize', 10);
                mslplot{2}=bar(mean_ts_aoi, 'FaceColor', [0.9 0.9 0.9],'parent',axRH);
  
                ylabel(axLH,'Easterly speed, m/s')
                ylabel(axRH,'AOI')
                ax_pos = get(axLH,'position');
                set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
%                 set(axRH,'color','none','yaxislocation','right','xtick', comb_ocean_time(1:30*length(inputmonth):end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
%                 set(axRH,'color','none','yaxislocation','right','xtick', 1:60*length(inputmonth):length(comb_ocean_time),'xticklabel',datestr(comb_ocean_time(1:60*length(inputmonth):end), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
                
                set(axRH,'color','none','yaxislocation','right','xtick', cal_arr,'xticklabel',datestr(comb_ocean_time(cal_arr), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
                axis tight 

                set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
                set(axLH,'ycolor','m', 'box', 'off', 'FontSize',15);
                set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
                xlabel(axRH, 'Year');

                title(['Easterly', ',', 'AOI', ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i')]);
%                 datetick(axLH, 'x','yymm','keepticks');
%                 datetick(axRH, 'x','yymm','keepticks');
%                 datetick(axLH, 'x','yymm','keeplimits');
%                 datetick(axRH, 'x','yymm','keeplimits');
%                 datetick(axLH, 'x','yymm');
%                 datetick(axRH, 'x','yymm');
                set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');

                set(mslplot{1},'LineWidth',2);
                set(mslplot{2},'LineWidth',2);
                grid on
                
                
                lgd=legend([mslplot{1} mslplot{2}], 'Easterly', 'AOI');

%                 txt1=text(5, max(double(mean_ts_nev))-diff(double([min(mean_ts_nev), max(mean_ts_nev)]))/32.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 
                    
                half_len=length(-mean_ts_uwind)/2;
          
%                 wind_half_1=mean(-mean_ts_uwind(1:half_len));
%                 wind_half_2=mean(-mean_ts_uwind(half_len+1:end));
%                 txt1=text(5, max(double(mean_ts_aoi))-diff(double([min(mean_ts_aoi), max(mean_ts_aoi)]))/32.0 , ...
%                     ['easterly = ', num2str(round(wind_half_1,2)), ' / ', num2str(round(wind_half_2,2))], 'FontSize', 20); 

%                 aoi_half_1=mean(mean_ts_aoi(1:half_len), 'omitnan');
%                 aoi_half_2=mean(mean_ts_aoi(half_len+1:end), 'omitnan');
%                 txt2=text(5, max(double(mean_ts_aoi))-diff(double([min(mean_ts_aoi), max(mean_ts_aoi)]))/8.0 , ...
%                     ['aoi = ', num2str(round(aoi_half_1,2)), ' / ', num2str(round(aoi_half_2,2))], 'FontSize', 20); 


                set(lgd,'FontSize',15);
                set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
                set(lgd,'Orientation','horizontal');

                set(gcf,'PaperPosition', [0 0 36 12]) 
                saveas(gcf,jpgname,'tif');
                grid off
                
                disp(' ')
                disp([fig_name, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')
                
                close all;
            end