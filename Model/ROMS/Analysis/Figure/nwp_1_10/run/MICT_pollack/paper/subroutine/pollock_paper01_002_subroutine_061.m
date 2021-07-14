% %  Updated    14-Jun-2021 by Yong-Yub Kim   % # of eggs -> # of indivisuals
if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    ind=1;
    clear ref_data comb_ref_data
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            ref_data=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
            lastday_m=size(ref_data,3);
            if (exist('comb_ref_data')==0)
                comb_ref_data=ref_data;
            else
                comb_ref_data(:,:,end+1:end+lastday_m)=ref_data;
            end
        end
    end
    sum_ref_data=sum(comb_ref_data,3, 'omitnan');
    lon_sum_ref_data=sum(sum_ref_data, 1, 'omitnan');
%                 mean_ref_data(mean_ref_data==0)=NaN;

    clear data comb_data
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            data=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
            lastday_m=size(data,3);
            if (exist('comb_data')==0)
                comb_data=data;
            else
                comb_data(:,:,end+1:end+lastday_m)=data;
            end
        end
    end

    sum_later_data=sum(comb_data,3, 'omitnan');
    lon_sum_later_data=sum(sum_later_data, 1, 'omitnan');
%                 mean_later_data(mean_later_data==0)=NaN;
    lon_sum_data=lon_sum_later_data-lon_sum_ref_data;
    lon_sum_data(lon_sum_data==0)=NaN;
    lon_sum_percent=lon_sum_data./lon_sum_ref_data .*100.0;
    
    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
%                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
%                 mean_data = sum(comb_egg_mask,3)/size(comb_egg_mask,3);
%     
%     mean_data = mean(comb_egg_mask,3);

    testnameind=1;
%     lon_sum_data=sum(sum_data, 1, 'omitnan');
    lat_1d=mean(lat_rho,1,'omitnan');
    barh(lat_1d, lon_sum_percent)
    ylabel('latitude (degree)')
    xlabel('percentage of change')
    xlim([-150 150])
%     ylim([37 41])
    set(gca, 'fontsize', m_pcolor_title_fontsize)
    titlename = strcat('particle change ratio, ',testname, ',(', ...
        num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
        num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 
    
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    grid on
    grid minor
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_ver paper_position_hor paper_position_width paper_position_height]) 

    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([fig_flags{1,1}, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    close all;
end