disp('subroutine SSH_2p_2nd_010_sub_001_get_bndy') 

tmp.lap_time_j=tic;

tmp.elapsed=toc(tmp.lap_time_j);
tmp.totlap= (length(RCM_info.years)*length(RCM_info.months));
tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-1+1)/1);


for bndyij=1:length(RCM_info.bndy_directions)
    for yearij=1:length(RCM_info.years)
        tmp.tempyear=RCM_info.years(yearij);
        tmp.yearstr=num2str(tmp.tempyear, '%04i');
        for monthij=1:length(RCM_info.months)
            tmp.tempmonth=RCM_info.months(monthij);
            tmp.monthstr=num2str(tmp.tempmonth, '%02i');
            tmp.bndyfilename=[dirs.glorysdir, 'mercatorglorys12v1_gl12_mean_', tmp.yearstr, tmp.monthstr, '.nc'];

            if (isfield(RCM_grid, 'lon_rho') ~= 1)
%                 for gridi=1:length(RCM_grid.gridname)
% %                     RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.(['filename_', RCM_grid.gridname{gridi}]), RCM_grid.gridname{gridi});
%                     RCM_grid.(RCM_grid.gridname{gridi})=ncread(RCM_grid.filename, RCM_grid.gridname{gridi});
%                 end
%                 ncinfo(tmp.bndyfilename)
                GCM_grid.lon_glo=ncread(tmp.bndyfilename, 'longitude');
                GCM_grid.lat_glo=ncread(tmp.bndyfilename, 'latitude');
                GCM_grid.depth=ncread(tmp.bndyfilename, 'depth');
                GCM_grid.N = length(GCM_grid.depth);
                GCM_grid.depth= [10; -GCM_grid.depth];
                
                [GCM_grid.lon_min, GCM_grid.lon_max, GCM_grid.lat_min, GCM_grid.lat_max] = ...
                    findind_Y(1/12, RCM_grid.domain(1:4), GCM_grid.lon_glo, GCM_grid.lat_glo);
                GCM_grid.lon_count= GCM_grid.lon_max - GCM_grid.lon_min +1;
                GCM_grid.lat_count= GCM_grid.lat_max - GCM_grid.lat_min +1;
                
                
                GCM_grid.lat_rho_north=GCM_grid.lat_glo(GCM_grid.lat_max);
                GCM_grid.lat_rho_east=GCM_grid.lat_glo(GCM_grid.lat_min:GCM_grid.lat_max);
                GCM_grid.lat_rho_south=GCM_grid.lat_glo(GCM_grid.lat_min);
                GCM_grid.lat_rho_west=GCM_grid.lat_glo(GCM_grid.lat_min:GCM_grid.lat_max);
                
                GCM_grid.lon_rho_north=GCM_grid.lon_glo(GCM_grid.lon_min:GCM_grid.lon_max);
                GCM_grid.lon_rho_east=GCM_grid.lon_glo(GCM_grid.lon_max);
                GCM_grid.lon_rho_south=GCM_grid.lon_glo(GCM_grid.lon_min:GCM_grid.lon_max);
                GCM_grid.lon_rho_west=GCM_grid.lon_glo(GCM_grid.lon_min);
                
                GCM_grid.stddepth = RCM_grid.stddepth;
%                 GCM_grid.mask_rho_north=GCM_grid.mask_rho(:,end);
%                 GCM_grid.mask_rho_east=GCM_grid.mask_rho(end,:);
%                 GCM_grid.mask_rho_south=GCM_grid.mask_rho(:,1);
%                 GCM_grid.mask_rho_west=GCM_grid.mask_rho(1,:);
                
                
%                 RCM_grid.theta_s=ncread(tmp.bndyfilename, 'theta_s');
%                 RCM_grid.theta_b=ncread(tmp.bndyfilename, 'theta_b');
%                 RCM_grid.hc=ncread(tmp.bndyfilename, 'hc');
%                 RCM_grid.N=length(ncread(tmp.bndyfilename, 's_rho'));
%                 RCM_grid.Vstretching=ncread(tmp.bndyfilename, 'Vstretching');
%                 RCM_grid.Vtransform=ncread(tmp.bndyfilename, 'Vtransform');
                
%                 cut_lon_rho = ...
%                     RCM_grid.lon_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
%                 cut_lat_rho = ...
%                     RCM_grid.lat_rho(RCM_grid.lon_min(1):RCM_grid.lon_max(1), RCM_grid.lat_min(1):RCM_grid.lat_max(1));
            end
            
%             RCM_data.(['zeta','_', RCM_info.bndy_directions{bndyij}])(:,yearij,monthij)= ncread(tmp.bndyfilename,['zeta','_', RCM_info.bndy_directions{bndyij}], ...
%                     [1 monthij], ...
%                     [inf 1]);
%             ncinfo(tmp.bndyfilename)

% % %             RCM_data.(['z','_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij) = ...
% % %                 zlevs(RCM_grid.Vtransform, RCM_grid.Vstretching, ...
% % %                 RCM_grid.(['h','_', RCM_info.bndy_directions{bndyij}]), squeeze(RCM_data.(['zeta','_', RCM_info.bndy_directions{bndyij}])(:,yearij,monthij)), ...
% % %                 RCM_grid.theta_s, RCM_grid.theta_b, RCM_grid.hc, RCM_grid.N,  'r')';
%             zw = zlevs(h,zeta,theta_s,theta_b,hc,N,Vstretching,Vtransform,'w');

            switch RCM_info.bndy_directions{bndyij}
                case {'north'}
                   GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij)= ...
                       ncread(tmp.bndyfilename, tmp.variable_GCM, ...
                    [GCM_grid.lon_min GCM_grid.lat_max 1 1], [GCM_grid.lon_count 1 inf 1]);  
                case {'south'}
                    GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij)= ...
                       ncread(tmp.bndyfilename, tmp.variable_GCM, ...
                    [GCM_grid.lon_min GCM_grid.lat_min 1 1], [GCM_grid.lon_count 1 inf 1]);  
                case {'west'}
                    GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij)= ...
                       ncread(tmp.bndyfilename, tmp.variable_GCM, ...
                    [GCM_grid.lon_min GCM_grid.lat_min 1 1], [1 GCM_grid.lat_count inf 1]);
                case {'east'}
                     GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij)= ...
                       ncread(tmp.bndyfilename, tmp.variable_GCM, ...
                    [GCM_grid.lon_max GCM_grid.lat_min 1 1], [1 GCM_grid.lat_count inf 1]);           
            end
            
            for indi=1:size(GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}]), 1)
                tmp.data(2:GCM_grid.N+1)=GCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(indi,:,yearij,monthij);
                tmp.data(1)=tmp.data(2);
%                 tmp.data(end+1)=tmp.data(end);
%                 tmp.z(2:RCM_grid.N+1) = RCM_data.(['z_', RCM_info.bndy_directions{bndyij}])(indi,:,yearij,monthij);
%                 tmp.z(1) = -5500;
%                 tmp.z(end+1)=100;
                GCM_data.([tmp.variable_GCM,'_stddepth_', RCM_info.bndy_directions{bndyij}])(indi,:,yearij, monthij) = ...
                    interp1(GCM_grid.depth, tmp.data, GCM_grid.stddepth);
                tmp=rmfield(tmp, 'data');
%                 tmp=rmfield(tmp, 'z');
                %print time
                tmp.elapsed=toc(tmp.lap_time_j);
                tmp.templap= (yearij-1)*length(RCM_info.months)+monthij-1;
                fprintf(repmat('\b',1,sum(tmp.nchar)))  % remove printed time
                tmp.nchar(1)= fprintf(' %.0f sec. (%.0f sec. left)', tmp.elapsed, tmp.elapsed*(tmp.totlap-tmp.templap+1)/tmp.templap);
            end
            
%             RCM_data.([tmp.variable_GCM,'_stddepth_', RCM_info.bndy_directions{bndyij}])(find(RCM_grid.mask_rho_north==0),:,yearij,monthij) = NaN;
            
%             RCM_data.([tmp.variable_GCM,'_stddepth_', RCM_info.bndy_directions{bndyij}])= griddata( ...
%                 tmp.lonlat_2d, RCM_data.(['z','_', RCM_info.bndy_directions{bndyij}]), ...
%                 RCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(:,:,yearij,monthij), ...
%                 tmp.lonlat_2d_std, tmp.stddepth_2d);
%             interp1(RCM_data.(['z','_', RCM_info.bndy_directions{bndyij}])(950,:), ...
%                 RCM_data.([tmp.variable_GCM,'_', RCM_info.bndy_directions{bndyij}])(950,:,yearij,monthij), ...
%                 RCM_grid.stddepth)
            
            
            
            
%             for i = 1 : L
%                 if mod(i,100) == 0
%                     disp(['i = ',num2str(i),'/',num2str(L)])
%                 end
%                 for j = 1 : M
%                     org_depth = z(:,i,j);
%                     org_wdepth = zw(:,i,j);
%                     
%                     org_temp = temp(:,i,j); org_salt = salt(:,i,j);
%                     org_u = u(:,i,j); org_v = v(:,i,j); org_w = w(:,i,j);
%                     
%                     re_temp(:,i,j) = interp1(org_depth, org_temp, stddepth);
%                     re_salt(:,i,j) = interp1(org_depth, org_salt, stddepth);
%                     re_u(:,i,j) = interp1(org_depth, org_u, stddepth);
%                     re_v(:,i,j) = interp1(org_depth, org_v, stddepth);
%                     re_w(:,i,j) = interp1(org_wdepth, org_w, stddepth);
%                 end
%             end


        end
    end
end
abc=1

tmp.matsavefile=['/data2/kimyy/Reanalysis/MYOCEAN', filesep, 'GLORYS', '_', tmp.variable_GCM, '_', ...
    num2str(min(RCM_info.years),'%04i'),'_',num2str(max(RCM_info.years),'%04i'), '_', ...
    RCM_info.season, '.mat'];
save(tmp.matsavefile, 'GCM_data', 'GCM_grid');
% north, east, south, west read (variables)
% vertical grid definition
% interpolation
% 
% reanalysis data read
% interpolation
% get RMSE, get PCCs
