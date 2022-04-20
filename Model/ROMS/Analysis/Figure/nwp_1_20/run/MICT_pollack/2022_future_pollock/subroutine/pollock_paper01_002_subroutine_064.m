if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
        egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask mean_ts_EAWMI EAWMI_*
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
            lastday_m=size(u_rho,3);
            if (exist('comb_u_rho')==0)
                comb_ocean_time=ocean_time;
            else
                comb_ocean_time(end+1:end+lastday_m)=ocean_time;
            end
        end
    end
% 

    opts = spreadsheetImportOptions("NumVariables", 13);
    opts.Sheet = "EAWMI";
    opts.DataRange = "A2:M51";
    opts.VariableNames = ["yy", "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    climate2021S1 = readtable("Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\climate_2021.xlsx", opts, "UseExcel", false);
    climate2021S1 = table2array(climate2021S1);
    clear opts

    EAWMI_year = climate2021S1(:,1);
%     EAWMI_month = 1:12;
    EAWMI_value = climate2021S1(:,2:13)';
    EAWMI_value= EAWMI_value(:);
    EAWMI_year = repmat(EAWMI_year, [1 12])';
    EAWMI_year = EAWMI_year(:);
    EAWMI_month = 1:12;
    EAWMI_month = repmat(EAWMI_month, [length(EAWMI_year)/12 1])';
    EAWMI_month = EAWMI_month(:);
%     clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp;



% % % % % %                 10 year lowpass filter
% % % %     EAWMI_value=lowpass(EAWMI_value, 1/120, 1, 'Steepness', 0.99);  % phase correction is applied automatically, but butterworth isn't.

    
%                 bar(EAWMI_year+EAWMI_month/12, EAWMI_value_filt)
%                 plot(EAWMI_year+EAWMI_month/12, EAWMI_value)
%                 plot(EAWMI_year(1:end-27)+EAWMI_month(1:end-27)/12, EAWMI_value_filt(28:end))

    EAWMI_year=reshape(EAWMI_year, [12, length(EAWMI_year)/12]);
    EAWMI_month=reshape(EAWMI_month, [12, length(EAWMI_month)/12]);
    EAWMI_value=reshape(EAWMI_value, [12, length(EAWMI_value)/12]);
    EAWMI_year_start=find(EAWMI_year(1,:)==min(allyear));
    EAWMI_year_end=find(EAWMI_year(1,:)==max(allyear));

%                 plot(EAWMI_value(:,EAWMI_year_start:EAWMI_year_end))
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
    
    EAWMI_value_yearly = mean(EAWMI_value(inputmonth, :), 1);
    EAWMI_value_yearly=lowpass(EAWMI_value_yearly, 1/10, 1, 'Steepness', 0.7);  % phase correction is applied automatically, but butterworth isn't.
    EAWMI_value_yearly = EAWMI_value_yearly(allyear-1969);
    EAWMI_year_yearly = mean(EAWMI_year(inputmonth, allyear-1969), 1);
    EAWMI_month_yearly = mean(EAWMI_month(inputmonth, allyear-1969), 1);
%     EAWMI_value_yearly = mean(EAWMI_value(allyear-1969), 1);

    axLH = gca;
    mslplot{1}=bar(EAWMI_year_yearly, EAWMI_value_yearly, 'FaceColor', [0.8 0.8 0.8], 'parent',axLH);
    hold on
    half_len=round(length(EAWMI_value_yearly)/2);
    egg_half_1=mean(EAWMI_value_yearly(1:half_len), 'omitnan');
    egg_half_2=mean(EAWMI_value_yearly(half_len+1:end), 'omitnan');
    regime_ts_egg_mask(1:half_len)=egg_half_1;
    regime_ts_egg_mask(half_len+1:length(EAWMI_value_yearly))=egg_half_2;
    mslplot{2}=plot(EAWMI_year_yearly, regime_ts_egg_mask,'k','parent',axLH);

    ylabel(axLH,'EAWMI')

    axis tight 
    ylim([32 55])
    xlabel(axLH, 'Year');
    title([ 'EAWMI', ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i')]);
    set(mslplot{2},'LineWidth',4);
    grid on

    lgd=legend([ mslplot{1}], 'EAWMI');

    set(gca, 'FontSize',20);
    set(lgd,'FontSize',15);
    set(lgd,'Position',[0.13 0.85, 0.775, 0.03]);
    set(lgd,'Orientation','horizontal');

    set(gcf,'PaperPosition', [0 0 36 12]) 
    saveas(gcf,jpgname,'tif');
    grid off

    disp(' ')
    disp([fig_name, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    hold off 
    close all;
end