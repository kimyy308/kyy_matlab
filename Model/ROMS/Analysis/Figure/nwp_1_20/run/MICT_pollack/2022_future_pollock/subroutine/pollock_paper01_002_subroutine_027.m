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
                            vel_polygon=[130.5, 41.5; 130.5, 42.5; 132, 42.5; 132, 41.5];
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
                
                mean_ts_nev= mean_ts_v_rho.*cosd(45)+mean_ts_u_rho.*cosd(45);

                
                axLH = gca;
                axRH = axes('color','none');
%                 mslplot{1}=plot(comb_ocean_time,sum_ts_egg_mask,'b','parent',axLH);
%                 mslplot{2}=plot(comb_ocean_time,mean_ts_nwv, 'k','parent',axRH);
                mslplot{1}=plot(mean_ts_vwind,'k','parent',axLH);
                mslplot{2}=plot(mean_ts_nev, 'r','parent',axRH);
                ylabel(axLH,'Southerly speed, m/s')
                ylabel(axRH,'Northeastward Velocity, m/s')
                ax_pos = get(axLH,'position');
                set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
%                 set(axRH,'color','none','yaxislocation','right','xtick', comb_ocean_time(1:30*length(inputmonth):end), 'position', ax_pos+[0 0.02 -0.01 -0.02]);
                set(axRH,'color','none','yaxislocation','right','xtick', 1:60*length(inputmonth):length(comb_ocean_time),'xticklabel',datestr(comb_ocean_time(1:60*length(inputmonth):end), 'yyyy'),'position', ax_pos+[0 0.02 -0.01 -0.02]);
                axis tight 

                set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
                set(axLH,'ycolor','k', 'box', 'off', 'FontSize',15);
                set(axRH,'ycolor','r', 'box', 'off', 'FontSize',15);
                xlabel(axRH, 'Year');

                title(['Southerly', ',', 'NEV', ',', num2str(min(allyear),'%04i'),'-',num2str(max(allyear),'%04i')]);
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
                
                
                lgd=legend([mslplot{1} mslplot{2}], 'Southerly', 'Northeast V');
                tempcorr=corrcoef(mean_ts_vwind, mean_ts_nev);  %egg & nwv !!!

                txt1=text(5, max(double(mean_ts_nev))-diff(double([min(mean_ts_nev), max(mean_ts_nev)]))/32.0 ,['R = ', num2str(round(tempcorr(1,2),2))], 'FontSize', 20); 

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