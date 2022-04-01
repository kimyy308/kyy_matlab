disp('subroutine SSH_2p_2nd_010_sub_003_get_RMS_bndy_GLORYS') 

tmp.matsavefile=[dirs.matdir, filesep, tmp.testname, '_', tmp.variable, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(tmp.matsavefile)

tmp.gcmmatsavefile=[dirs.glorysdir, filesep, 'GLORYS', '_', tmp.variable_GCM, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
load(tmp.gcmmatsavefile)

% abc=1;
RCM_grid.n_north=length(RCM_grid.mask_rho_north);
RCM_grid.n_south=length(RCM_grid.mask_rho_south);
RCM_grid.n_west=length(RCM_grid.mask_rho_west);
RCM_grid.n_east=length(RCM_grid.mask_rho_east);

RCM_all=RCM_data.([tmp.variable, '_stddepth_north']);
RCM_all(end+1:end+RCM_grid.n_west,:,:,:)=RCM_data.([tmp.variable, '_stddepth_east']);
RCM_all(end+1:end+RCM_grid.n_south,:,:,:)=RCM_data.([tmp.variable, '_stddepth_south']);
RCM_all(end+1:end+RCM_grid.n_west,:,:,:)=RCM_data.([tmp.variable, '_stddepth_west']);

RCM_all_bndy.comb_mean=squeeze(mean(RCM_all,3, 'omitnan'));
RCM_all_bndy.comb_mean=mean(RCM_all_bndy.comb_mean,3, 'omitnan');

% pcolor(RCM_all_mean'); shading flat; colorbar;
% pcolor(RCM_data.salt_stddepth_south(:,:,1,1)'); shading flat; colorbar;
% pcolor(RCM_data.salt_stddepth_north(:,:,1,1)'); shading flat; colorbar;

GCM_all=GCM_data.([tmp.variable_GCM, '_stddepth_north']);
GCM_all(end+1:end+size(GCM_data.so_stddepth_east,1),:,:,:)=GCM_data.([tmp.variable_GCM, '_stddepth_east']); 
GCM_all(end+1:end+size(GCM_data.so_stddepth_south,1),:,:,:)=GCM_data.([tmp.variable_GCM, '_stddepth_south']);
GCM_all(end+1:end+size(GCM_data.so_stddepth_west,1),:,:,:)=GCM_data.([tmp.variable_GCM, '_stddepth_west']);

GCM_all_mean=squeeze(mean(GCM_all,3, 'omitnan'));
GCM_all_mean=mean(GCM_all_mean,3, 'omitnan');
% pcolor(GCM_all_mean'); shading flat; colorbar;


for bndyij=1:length(RCM_info.bndy_directions)
    tmp.direction=RCM_info.bndy_directions{bndyij};
    %% get horizontal distance(km)
    RCM_grid.hor_len=length(RCM_grid.(['mask_rho_', tmp.direction]));
    RCM_grid.(['hor_dist_', tmp.direction])= ...
        m_lldist([RCM_grid.(['lon_rho_', tmp.direction])(1), RCM_grid.(['lon_rho_', tmp.direction])(2)], ...
       [RCM_grid.(['lat_rho_', tmp.direction])(1), RCM_grid.(['lat_rho_', tmp.direction])(2)]);
   RCM_grid.stddepth_thick = diff([0, ...
        (-RCM_grid.stddepth(2:end)+-RCM_grid.stddepth(1:end-1))/2.,...
        abs(RCM_grid.stddepth(end))+abs(RCM_grid.stddepth(end)-RCM_grid.stddepth(end-1))/2]); 
   RCM_grid.(['area_', tmp.direction]) = repmat(RCM_grid.stddepth_thick, [RCM_grid.hor_len 1]) ...
       .*RCM_grid.(['hor_dist_', tmp.direction]);
   RCM_grid.(['stddepth_',  tmp.direction]) = repmat(RCM_grid.stddepth, [RCM_grid.hor_len 1]) ...
       .*RCM_grid.(['hor_dist_', tmp.direction]);
   RCM_grid.(['stddepth_',  tmp.direction]) = repmat(RCM_grid.stddepth, [RCM_grid.hor_len 1]) ...
       .*RCM_grid.(['hor_dist_', tmp.direction]);
   
%    RCM_grid.stddepth_thick=diff([0; (RCM_vert_grid.depth(2:end)+RCM_vert_grid.depth(1:end-1))/2; 3750]);
   
    tmp.data=reshape(GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction]), ... 
        [size(GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction]),1), ...
        size(GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction]),2), ...
        size(GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction]),3) ...
        *size(GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction]),4)]);
    GCM_bndy_mean.([tmp.direction])=mean(tmp.data, 3);
    
    switch tmp.direction
        case {'north'}
            [GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, GCM_grid.(['lon_rho_', tmp.direction]));
            [RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, RCM_grid.(['lon_rho_', tmp.direction]));
        case {'south'}
            [GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, GCM_grid.(['lon_rho_', tmp.direction]));
            [RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, RCM_grid.(['lon_rho_', tmp.direction]));
        case {'west'}
            [GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, GCM_grid.(['lat_rho_', tmp.direction]));
            [RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, RCM_grid.(['lat_rho_', tmp.direction]));            
        case {'east'}
            [GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, GCM_grid.(['lat_rho_', tmp.direction]));
            [RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction])] = meshgrid(RCM_grid.stddepth, RCM_grid.(['lat_rho_', tmp.direction]));            
    end
    GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_', tmp.direction, '_interped']) = ...
        griddata(GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction]), GCM_bndy_mean.([tmp.direction]), ...
            RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction]));
    
    yearlen=size(RCM_data.([tmp.variable, '_stddepth_north']),3);
    monthlen=size(RCM_data.([tmp.variable, '_stddepth_north']),4);
    
    for ti=1:size(tmp.data,3)
        yearij=floor((ti-1)/monthlen)+1;
        monthij=ti-((yearij-1)*monthlen);
        GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction, '_interped'])(:,:,ti)= ...
            griddata(GCM_grid.(['stddepth_',  tmp.direction]), GCM_grid.(['lonlat_',  tmp.direction]), tmp.data(:,:,ti), ...
             RCM_grid.(['stddepth_',  tmp.direction]), RCM_grid.(['lonlat_',  tmp.direction]));
        RCM_eval.(['DIFF_', tmp.direction])(:,:,ti)= ...
            (RCM_data.([tmp.variable, '_stddepth_', tmp.direction])(:,:,yearij,monthij) ...
            - GCM_data.([tmp.variable_GCM, '_stddepth_', tmp.direction, '_interped'])(:,:,ti));
    end
        RCM_eval.(['MS_', tmp.direction])= ...
            (RCM_eval.(['DIFF_', tmp.direction])) .^2./size(tmp.data,3);
        RCM_eval.(['RMS_', tmp.direction])=sum(RCM_eval.(['MS_', tmp.direction]), 3);
        RCM_eval.(['RMS_', tmp.direction]) = sqrt(RCM_eval.(['RMS_', tmp.direction]));
        RCM_eval.(['BIAS_', tmp.direction])= mean(RCM_eval.(['DIFF_', tmp.direction]), 3);
        
         pcolor(RCM_eval.BIAS_north'); shading flat; colorbar;
        
    tmp.data=reshape(RCM_data.([tmp.variable, '_stddepth_', tmp.direction]), ... 
        [size(RCM_data.([tmp.variable, '_stddepth_', tmp.direction]),1), ...
        size(RCM_data.([tmp.variable, '_stddepth_', tmp.direction]),2), ...
        size(RCM_data.([tmp.variable, '_stddepth_', tmp.direction]),3) ...
        *size(RCM_data.([tmp.variable, '_stddepth_', tmp.direction]),4)]);
    RCM_bndy_mean.([tmp.variable, '_stddepth_', tmp.direction])=mean(tmp.data, 3);
end

GCM_all_bndy.comb_mean=GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_north_interped']);
GCM_all_bndy.comb_mean(end+1:end+size(GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_east_interped']),1),:,:,:)= ...
    GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_east_interped']); 
GCM_all_bndy.comb_mean(end+1:end+size(GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_south_interped']),1),:,:,:)= ...
    GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_south_interped']);
GCM_all_bndy.comb_mean(end+1:end+size(GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_west_interped']),1),:,:,:)= ...
    GCM_bndy_mean.([tmp.variable_GCM, '_stddepth_west_interped']);


RCM_all_bndy.comb_area=  RCM_grid.area_north;
RCM_all_bndy.comb_area(end+1:end+RCM_grid.n_west,:)=RCM_grid.area_west;
RCM_all_bndy.comb_area(end+1:end+RCM_grid.n_south,:)=RCM_grid.area_south;
RCM_all_bndy.comb_area(end+1:end+RCM_grid.n_east,:)=RCM_grid.area_east;


RCM_all=RCM_data.([tmp.variable, '_stddepth_north']);
RCM_all(end+1:end+RCM_grid.n_west,:,:,:)=RCM_data.([tmp.variable, '_stddepth_east']);
RCM_all(end+1:end+RCM_grid.n_south,:,:,:)=RCM_data.([tmp.variable, '_stddepth_south']);
RCM_all(end+1:end+RCM_grid.n_west,:,:,:)=RCM_data.([tmp.variable, '_stddepth_west']);

% switch RCM_info.bndy_directions{bndyij}
%     case {'north'}
% 
%     case {'south'}
% 
%     case {'west'}
% 
%     case {'east'}
% 
% end
RCM_all_bndy.comb_area2= RCM_all_bndy.comb_area(logical(~isnan(GCM_all_bndy.comb_mean).*~isnan(RCM_all_bndy.comb_mean)));
GCM_all_bndy.comb_mean2=GCM_all_bndy.comb_mean(logical(~isnan(GCM_all_bndy.comb_mean).*~isnan(RCM_all_bndy.comb_mean)));
RCM_all_bndy.comb_mean2=RCM_all_bndy.comb_mean(logical(~isnan(GCM_all_bndy.comb_mean).*~isnan(RCM_all_bndy.comb_mean)));

GCM_all_bndy.areamean=sum(GCM_all_bndy.comb_mean2.*RCM_all_bndy.comb_area2)/sum(RCM_all_bndy.comb_area2)
RCM_all_bndy.areamean=sum(RCM_all_bndy.comb_mean2.*RCM_all_bndy.comb_area2)/sum(RCM_all_bndy.comb_area2)

corrcoef(GCM_all_bndy.comb_mean2, RCM_all_bndy.comb_mean2)

% pcolor(RCM_grid.lonlat_north, RCM_grid.stddepth_north, RCM_bndy_mean.salt_stddepth_north); shading flat; colorbar;
% pcolor(RCM_grid.lonlat_north, RCM_grid.stddepth_north, GCM_bndy_mean.so_stddepth_north_interped); shading flat; colorbar;


tmp.matsavefile=[dirs.matdir, filesep, tmp.testname, '_', tmp.variable, '_bndy_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
save(tmp.matsavefile, 'RCM_all', 'RCM_all_bndy', 'RCM_bndy_mean', 'RCM_grid', 'RCM_eval')

tmp.gcmmatsavefile=[dirs.glorysdir, filesep, 'GLORYS', '_', tmp.variable_GCM, '_bndy_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
save(tmp.gcmmatsavefile, 'GCM_all', 'GCM_all_bndy', 'GCM_bndy_mean', 'GCM_grid')

