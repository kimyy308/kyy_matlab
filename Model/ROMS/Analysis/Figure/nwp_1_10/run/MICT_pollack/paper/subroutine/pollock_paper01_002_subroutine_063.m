% %  'particle latitude boxplot (elapsed ?? days)';
% %  Updated    23-Jun-2021 by Yong-Yub Kim   

if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
    run(param_script);
    
    temp_checktime = 15;    
    clear ref_egg_mask ref_comb_egg_mask
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])

            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            ref_egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
            
            lastday_m=size(ref_egg_mask,3);
            if (exist('ref_comb_egg_mask')==0)
                ref_comb_egg_mask=ref_egg_mask;
            else
                ref_comb_egg_mask(:,:,end+1:end+lastday_m)=ref_egg_mask;
            end
        end
    end
    
    ref_sum_data = sum(ref_comb_egg_mask,3, 'omitnan');
    ref_sum_data(ref_sum_data==0)=NaN;
    ref_lon_sum_data=sum(ref_sum_data, 1, 'omitnan');
    
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
    
    sum_data = sum(comb_egg_mask,3, 'omitnan');
    sum_data(sum_data==0)=NaN;
    lon_sum_data=sum(sum_data, 1, 'omitnan');
    
    lon_rho = ncread(ncname, 'lon_rho');
    lat_rho = ncread(ncname, 'lat_rho');
    lat_1d=mean(lat_rho,1,'omitnan');    
    
    clear lat_pars ref_lat_pars
    for lati=1:length(lat_1d)
        if exist('lat_pars','var')~=1 && lon_sum_data(lati)>0
            lat_pars=repelem(lat_1d(lati), lon_sum_data(lati));
        elseif exist('lat_pars','var')==1 && lon_sum_data(lati)>0
            lat_pars(end+1:end+lon_sum_data(lati))=repelem(lat_1d(lati), lon_sum_data(lati));
        end
        
        if exist('ref_lat_pars','var')~=1 && ref_lon_sum_data(lati)>0
            ref_lat_pars=repelem(lat_1d(lati), ref_lon_sum_data(lati));
        elseif exist('ref_lat_pars','var')==1 && ref_lon_sum_data(lati)>0
            ref_lat_pars(end+1:end+ref_lon_sum_data(lati))=repelem(lat_1d(lati), ref_lon_sum_data(lati));
        end
    end
    
% % % % % % 30d

    temp_checktime = 30;    
    clear ref_egg_mask ref_comb_egg_mask
    for yearij = 1:length(refyear)
        tempyear = refyear(yearij);
        for monthij = 1:length(inputmonth)
            tempmonth = inputmonth(monthij);
            ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])

            all_checktime=ncread(ncname, 'checktime');
            ind_checktime=find(all_checktime==temp_checktime);
            ref_egg_mask=squeeze(ncread(ncname, 'mask_par', [1 1 1 ind_checktime], [inf inf inf 1]));
            
            lastday_m=size(ref_egg_mask,3);
            if (exist('ref_comb_egg_mask')==0)
                ref_comb_egg_mask=ref_egg_mask;
            else
                ref_comb_egg_mask(:,:,end+1:end+lastday_m)=ref_egg_mask;
            end
        end
    end
    
    ref_sum_data = sum(ref_comb_egg_mask,3, 'omitnan');
    ref_sum_data(ref_sum_data==0)=NaN;
    ref_lon_sum_data_30d=sum(ref_sum_data, 1, 'omitnan');
    
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
    
    sum_data = sum(comb_egg_mask,3, 'omitnan');
    sum_data(sum_data==0)=NaN;
    lon_sum_data_30d=sum(sum_data, 1, 'omitnan');
    
    clear lat_pars_30d ref_lat_pars_30d
    for lati=1:length(lat_1d)
        if exist('lat_pars_30d','var')~=1 && lon_sum_data_30d(lati)>0
            lat_pars_30d=repelem(lat_1d(lati), lon_sum_data_30d(lati));
        elseif exist('lat_pars_30d','var')==1 && lon_sum_data_30d(lati)>0
            lat_pars_30d(end+1:end+lon_sum_data_30d(lati))=repelem(lat_1d(lati), lon_sum_data_30d(lati));
        end
        
        if exist('ref_lat_pars_30d','var')~=1 && ref_lon_sum_data_30d(lati)>0
            ref_lat_pars_30d=repelem(lat_1d(lati), ref_lon_sum_data_30d(lati));
        elseif exist('ref_lat_pars_30d','var')==1 && ref_lon_sum_data_30d(lati)>0
            ref_lat_pars_30d(end+1:end+ref_lon_sum_data_30d(lati))=repelem(lat_1d(lati), ref_lon_sum_data_30d(lati));
        end
    end

    
clear tot_lat_pars
tot_lat_pars=[ref_lat_pars'; lat_pars'; ref_lat_pars_30d'; lat_pars_30d'];
g_lat_pars1= repmat({'Early (15d)'},length(ref_lat_pars),1);
g_lat_pars2 = repmat({'Late (15d)'},length(lat_pars),1);
g_lat_pars3= repmat({'Early (30d)'},length(ref_lat_pars_30d),1);
g_lat_pars4 = repmat({'Late (30d)'},length(lat_pars_30d),1);
g_lat_pars = [g_lat_pars1; g_lat_pars2; g_lat_pars3; g_lat_pars4];
boxplot(tot_lat_pars, g_lat_pars, 'BoxStyle', 'outline', 'Colors', 'brbr', 'ColorGroup', g_lat_pars, ...
    'Whisker', 1.5)
ylabel('Latitude');
[latstat.mean, latstat.median, latstat.min, latstat.max, latstat.numel, latstat.meanci, latstat.std, latstat.gname, latstat.predci] = ...
    grpstats(tot_lat_pars, g_lat_pars, {'mean', 'median', 'min', 'max', 'numel', 'meanci',  'std', 'gname', 'predci'}, 'alpha', 0.05)
    
 [h,p]=ttest2(ref_lat_pars, lat_pars)
 [h,p]=ttest2(ref_lat_pars_30d, lat_pars_30d)
    set(gca, 'fontsize', m_pcolor_title_fontsize)
    titlename = strcat( 'particle latitude, ',testname, ',(', ...
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