% %  'latitudinal particle density (elapsed ??days, horizontal bar)';
% %  Updated    14-Jun-2021 by Yong-Yub Kim   % # of eggs -> sum of probabilities

if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    ind=1;
    clear egg_mask comb_egg_mask
    for yearij = 1:length(inputyear)
        tempyear = inputyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])

            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
            
            lastday_m=size(egg_mask,3);
            if (exist('comb_egg_mask')==0)
                comb_egg_mask=egg_mask;
            else
                comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
            end
        end
    end
    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
%                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
%                 mean_data = sum(comb_egg_mask,3)/size(comb_egg_mask,3);
    sum_data = sum(comb_egg_mask,3, 'omitnan');
    sum_data(sum_data==0)=NaN;
    testnameind=1;
    lon_sum_data=sum(sum_data, 1, 'omitnan');
    lon_sum_data_nan=lon_sum_data;
    lon_sum_data_nan(lon_sum_data==0)=NaN;
    lat_1d=mean(lat_rho,1,'omitnan');
    [a,b]=find(lon_sum_data==max(lon_sum_data))
    max(lon_sum_data)
    max_lat=lat_1d(a,b)
    cumsum_lon_sum_data=cumsum(lon_sum_data);
%     quantile(cumsum_lon_sum_data,0.5)
    lat_1d_nan=lat_1d;
    lat_1d_nan(isnan(lon_sum_data_nan))=NaN;
    min(lat_1d_nan)
    max(lat_1d_nan)
%     barh(lat_1d, cumsum_lon_sum_data)
    barh(lat_1d, lon_sum_data)
    ylabel('latitude (degree)')
    xlabel('# of individuals')
    xlim([0 1500])
%     ylim([36 42.5])
    set(gca, 'fontsize', m_pcolor_title_fontsize)
    titlename = strcat(num2str(temp_checktime, '%02i'), 'd,meridional particle existence ratio, ',testname, ',(', ...
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