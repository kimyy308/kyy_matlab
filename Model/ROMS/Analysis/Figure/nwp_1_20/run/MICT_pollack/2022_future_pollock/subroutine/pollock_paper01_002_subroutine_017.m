 if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                run(param_script);
                clear wind_curl comb_wind_curl u_rho comb_u_rho v_rho comb_v_rho uwind comb_uwind vwind comb_vwind ...
                    egg_mask comb_egg_mask sp_ground comb_sp_ground wind_curl2 comb_wind_curl2 temp_surf comb_temp_surf  curl_mask
                for yearij = 1:length(allyear)
                    tempyear = allyear(yearij);
                    for monthij = 1:length(inputmonth)
                        tempmonth = inputmonth(monthij);
                        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])  
                        if(exist('curl_mask')==0)
                            lon_rho = ncread(ncname, 'lon_rho');
                            lat_rho = ncread(ncname, 'lat_rho');
                            wstr_curl_polygon=[128, 39.2; 128, 40.2; 129, 40.2; 129, 39.2];
                            curl_mask = double(inpolygon(lon_rho,lat_rho,wstr_curl_polygon(:,1),wstr_curl_polygon(:,2)));
                            
%                             vel_polygon=[128, 38; 128, 39; 129, 39; 129, 38];
%                             vel_polygon=[128, 38; 128, 39.5; 129, 39.5; 129, 38];
%                             vel_polygon=[128, 38; 128, 39.3; 129, 39.3; 129, 38];
%                             vel_polygon=[128, 38; 128, 40; 129, 40; 129, 38];
%                             vel_polygon=[128, 38; 128, 39.5; 130, 39.5; 130, 38];
%                             vel_polygon=[128, 38; 128, 39.5; 130, 39.5; 130, 39; 129, 38];
%                             vel_polygon=[129, 38; 129, 39; 130, 39; 130, 38];
%                             vel_polygon=[129.5, 42; 129.5, 42.5; 132, 42.5; 132, 42];
                            vel_polygon=[130, 40.5; 130, 41.5; 131, 41.5; 131, 40.5];
%                             vel_polygon=pollock_egg2polygon;
                            vel_mask = double(inpolygon(lon_rho,lat_rho,vel_polygon(:,1),vel_polygon(:,2)));
                        end
                        wind_curl=ncread(ncname, 'wstr_curl').*curl_mask;
                        wind_curl(wind_curl==0)=NaN;
                        wind_curl2=ncread(ncname, 'wstr_curl').*vel_mask;
                        wind_curl2(wind_curl2==0)=NaN;
                        u_rho=ncread(ncname, 'u_rho').*vel_mask;
                        u_rho(u_rho==0)=NaN;
                        v_rho=ncread(ncname, 'v_rho').*vel_mask;
                        v_rho(v_rho==0)=NaN;
                        uwind=ncread(ncname, 'uwstr').*vel_mask;
                        uwind(uwind==0)=NaN;
                        vwind=ncread(ncname, 'vwstr').*vel_mask;
                        vwind(vwind==0)=NaN;
                        egg_mask=ncread(ncname, 'mask_15day').*vel_mask;
                        egg_mask(egg_mask==0)=NaN;
                        sp_ground=ncread(ncname, 'egg_mask').*vel_mask;
                        sp_ground(sp_ground==0)=NaN;
                        temp_surf=ncread(ncname, 'temp_surf').*vel_mask;
                        temp_surf(temp_surf==0)=NaN;
                        lastday_m=size(wind_curl,3);
                        if (exist('comb_wind_curl')==0)
                            comb_wind_curl=wind_curl;
                            comb_wind_curl2=wind_curl2;
                            comb_u_rho=u_rho;
                            comb_v_rho=v_rho;
                            comb_uwind=uwind;
                            comb_vwind=vwind;
                            comb_egg_mask=egg_mask;
                            comb_sp_ground=sp_ground;
                            comb_temp_surf=temp_surf;
                        else
                            comb_wind_curl(:,:,end+1:end+lastday_m)=wind_curl;
                            comb_wind_curl2(:,:,end+1:end+lastday_m)=wind_curl2;
                            comb_u_rho(:,:,end+1:end+lastday_m)=u_rho;
                            comb_v_rho(:,:,end+1:end+lastday_m)=v_rho;
                            comb_uwind(:,:,end+1:end+lastday_m)=uwind;
                            comb_vwind(:,:,end+1:end+lastday_m)=vwind;
                            comb_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
                            comb_sp_ground(:,:,end+1:end+lastday_m)=sp_ground;
                            comb_temp_surf(:,:,end+1:end+lastday_m)=temp_surf;
                        end
                    end
                end
                
                ts_wind_curl=reshape(comb_wind_curl,[size(comb_wind_curl,1)*size(comb_wind_curl,2), size(comb_wind_curl,3)]);
                mean_ts_wind_curl=mean(ts_wind_curl,1,'omitnan');
                ts_wind_curl2=reshape(comb_wind_curl2,[size(comb_wind_curl2,1)*size(comb_wind_curl2,2), size(comb_wind_curl2,3)]);
                mean_ts_wind_curl2=mean(ts_wind_curl2,1,'omitnan');
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
                ts_sp_ground=reshape(comb_sp_ground,[size(comb_sp_ground,1)*size(comb_sp_ground,2), size(comb_sp_ground,3)]);
                sum_ts_sp_ground=sum(ts_sp_ground,1,'omitnan');
                ts_temp_surf=reshape(comb_temp_surf,[size(comb_temp_surf,1)*size(comb_temp_surf,2), size(comb_temp_surf,3)]);
                mean_ts_temp_surf=mean(ts_temp_surf,1,'omitnan');
                
                sum_ratio=sum_ts_egg_mask./sum_ts_sp_ground;
                sum_ratio(isnan(sum_ratio))=0;
                sum_ratio(isinf(sum_ratio))=1;
%                 pcolor(lon_rho',lat_rho', wind_curl(:,:,1)'); shading flat      
%                 plot(mean_ts_v_rho)
                corrcoef(mean_ts_wind_curl, mean_ts_v_rho)
                corrcoef(mean_ts_wind_curl, mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45))
                corrcoef(mean_ts_wind_curl2, mean_ts_v_rho)
                corrcoef(mean_ts_wind_curl2, mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45)) % curl & nwv
                corrcoef(-mean_ts_uwind.*cosd(45)+mean_ts_vwind.*cosd(45), mean_ts_v_rho) %southeastery & nv
                corrcoef(-mean_ts_uwind, mean_ts_v_rho)  % eastery & nv
                corrcoef(-mean_ts_uwind, mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45))   % eastery & nwv !!!
                corrcoef(mean_ts_vwind, mean_ts_v_rho.*cosd(45)+mean_ts_u_rho.*cosd(45))   % southery & nev !!!
                
                corrcoef(sum_ts_egg_mask, mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45))  %egg & nwv !!!
                corrcoef(sum_ts_egg_mask, mean_ts_v_rho)  % egg & nv
                corrcoef(sum_ts_egg_mask, mean_ts_u_rho)  % egg & ev
                corrcoef(sum_ts_egg_mask, mean_ts_v_rho.*cosd(45)+mean_ts_u_rho.*cosd(45)) %egg & nev
                
                corrcoef(mean_ts_uwind, mean_ts_wind_curl) % rcurl & westery
                corrcoef(mean_ts_uwind, mean_ts_wind_curl2) % curl & westery
%                 corrcoef(sum_ratio, mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45))
                
                mean_ts_nwvel= mean_ts_v_rho.*cosd(45)-mean_ts_u_rho.*cosd(45);
%                 plot(mean_ts_wind_curl./mean(mean_ts_wind_curl))
%                 hold on
%                 plot(-mean_ts_uwind./mean(-mean_ts_uwind))
% %                 plot(-mean_ts_nwvel./mean(-mean_ts_nwvel))
% %                 plot(-sum_ts_egg_mask./mean(-sum_ts_egg_mask))
%                 hold off
                
                for lon_i=1:size(comb_u_rho,1)
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(-comb_uwind(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:));
                            sp_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            sp_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                comb_egg_mask2=comb_egg_mask;
                comb_egg_mask2(isnan(comb_egg_mask2))=0;
                for lon_i=1:size(comb_u_rho,1)   % egg_nv_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_egg_mask2(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:));
                            egg_nv_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            egg_nv_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % egg_nwv_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_egg_mask2(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:).*cosd(45)-comb_u_rho(lon_i,lat_i,:).*cosd(45));
                            egg_nwv_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            egg_nwv_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % egg_nev_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_egg_mask2(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:).*cosd(45)+comb_u_rho(lon_i,lat_i,:).*cosd(45));
                            egg_nev_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            egg_nev_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % eastery_v_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(-comb_uwind(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:));
                            eastery_v_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            eastery_v_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % southery_u_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_vwind(lon_i,lat_i,:),comb_u_rho(lon_i,lat_i,:));
                            southery_u_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            southery_u_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % curl_v_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_wind_curl2(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:));
                            curl_v_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            curl_v_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % temp_v_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_temp_surf(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:));
                            temp_v_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            temp_v_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % std_curl
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=std(comb_wind_curl2(lon_i,lat_i,:));
                            std_curl(lon_i,lat_i)=tempcorr;
                        else
                            std_curl(lon_i,lat_i)=NaN;
                        end
                    end
                end
                  
                for lon_i=1:size(comb_u_rho,1)  % std_uwind
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=std(comb_uwind(lon_i,lat_i,:));
                            std_uwind(lon_i,lat_i)=tempcorr;
                        else
                            std_uwind(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % std_v_rho
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=std(comb_v_rho(lon_i,lat_i,:));
                            std_v_rho(lon_i,lat_i)=tempcorr;
                        else
                            std_v_rho(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % eastery_curl_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(-comb_uwind(lon_i,lat_i,:),comb_wind_curl2(lon_i,lat_i,:));
                            eastery_curl_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            eastery_curl_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % eastery_nwv_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(-comb_uwind(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:).*cosd(45)-comb_u_rho(lon_i,lat_i,:).*cosd(45));
                            eastery_nwv_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            eastery_nwv_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % temp_nwv_corr
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=corrcoef(comb_temp_surf(lon_i,lat_i,:),comb_v_rho(lon_i,lat_i,:).*cosd(45)-comb_u_rho(lon_i,lat_i,:).*cosd(45));
                            temp_nwv_corr(lon_i,lat_i)=tempcorr(1,2);
                        else
                            temp_nwv_corr(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
                for lon_i=1:size(comb_u_rho,1)  % std_nwv
                    for lat_i=1:size(comb_u_rho,2)
                        if (isfinite(comb_v_rho(lon_i,lat_i,1)))
                            tempcorr=std(comb_v_rho(lon_i,lat_i,:).*cosd(45)-comb_u_rho(lon_i,lat_i,:).*cosd(45));
                            std_nwv(lon_i,lat_i)=tempcorr;
                        else
                            std_nwv(lon_i,lat_i)=NaN;
                        end
                    end
                end
                
%                 sp_corr=reshape(sp_corr,[61,105]);
%                 pcolor(lon_rho', lat_rho', sp_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', egg_nv_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', egg_nwv_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', egg_nev_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', eastery_v_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', southery_u_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', curl_v_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', temp_v_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', std_curl'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', std_uwind'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', std_v_rho'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', std_nwv'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', eastery_curl_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', eastery_nwv_corr'); shading flat; colorbar;
%                 pcolor(lon_rho', lat_rho', temp_nwv_corr'); shading flat; colorbar;

%                 mean_data(mean_data==0)=NaN;
                testnameind=1;
                close all;
            end