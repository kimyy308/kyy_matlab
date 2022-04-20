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

            egg_mask=ncread(ncname, 'egg_mask');
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
    
    mean_data = mean(comb_egg_mask,3);
    mask_model2 = double(inpolygon(lon_rho,lat_rho,refpolygon(:,1),refpolygon(:,2)));
    mask_model2(mask_model2==0)=NaN;
    mean_data = mean_data .* mask_model2;
    
    sum_data = sum(comb_egg_mask,3, 'omitnan');
    sum_data(sum_data==0)=NaN;
    sum_data = sum_data .* mask_model2;

    testnameind=1;
    lon_sum_data=sum(sum_data, 1, 'omitnan');
    lat_1d=mean(lat_rho,1,'omitnan');
    hAxes=axes;
    hbarplot=barh(lat_1d, lon_sum_data)
    ylabel('latitude (degree)')
    xlabel('# of particles')
    xlim([0 25000])
    ylim([37 51])
    set(gca, 'fontsize', m_pcolor_title_fontsize)
    titlename = strcat('meridional egg #, ',testname, ',(', ...
        num2str(min(inputyear),'%04i'),'-', num2str(max(inputyear),'%04i'), ',',  ...
        num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 
    hAxes.XAxis.Exponent = 0;
    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
    grid on
    grid minor
    set(gcf, 'PaperUnits', 'points');
    set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
    set(gcf,'PaperPosition', [paper_position_ver paper_position_hor paper_position_width paper_position_height]) 

    saveas(gcf,jpgname,'tif'); RemoveWhiteSpace([], 'file', jpgname);

    disp(' ')
    disp([fig_flags{1,1}, ' plot is created.'])
    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')
    close all;
end