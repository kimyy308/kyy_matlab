if (exist(jpgname , 'file') ~= 2 || fig_flag==2)
                clear mean_data comb_egg_mask
                run(param_script);
                ind=1;
                for yearij = 1:length(allyear)
                    tempyear = allyear(yearij);
                    clear yearly_egg_mask
                    for monthij = 1:length(inputmonth)
                        tempmonth = inputmonth(monthij);
                        ncname = [savedir,testname,'_',regionname,'model_pollock_',num2str(tempyear,'%04i'),'_',num2str(tempmonth,'%02i'),'.nc'];
                        disp([num2str(yearij), 'y_',num2str(monthij),'m'])

                        egg_mask=ncread(ncname, 'egg_mask');
                        lastday_m=size(egg_mask,3);
                        if (exist('yearly_egg_mask')==0)
                            yearly_egg_mask=egg_mask;
                        else
                            yearly_egg_mask(:,:,end+1:end+lastday_m)=egg_mask;
                        end
                    end
                    yearly_egg_mask=mean(yearly_egg_mask,3);
                    yearly_egg_mask(yearly_egg_mask==0)=NaN;
                    if (exist('comb_egg_mask')==0)
                        comb_egg_mask=yearly_egg_mask;
                    else
                        comb_egg_mask(:,:,end+1)=yearly_egg_mask;
                    end
                end
                lon_rho = ncread(ncname, 'lon_rho');
                lat_rho = ncread(ncname, 'lat_rho');
%                 pcolor(sum(comb_egg_mask,3)'/size(comb_egg_mask,3)); shading flat; colorbar
                mean_data = squeeze(mean(comb_egg_mask,1,'omitnan'));
                mean_data(mean_data==0)=NaN;
                testnameind=1;
                
                latmin=min(find(isfinite(mean(mean_data,2))));
                latmax=max(find(isfinite(mean(mean_data,2))));
                
                [mesh_allyear, mesh_lat]=meshgrid(allyear, squeeze(mean(lat_rho,1)));
                pcolor(mesh_allyear(latmin:latmax,:)', mesh_lat(latmin:latmax,:)', mean_data(latmin:latmax,:)');
                colormap(yrmap);
                caxis([0, 1.0]);
                shading(gca,m_pcolor_shading_method);   
                
                col_bar{testnameind,2}=colorbar;

%                if testnameind==4
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [120,  130, 140], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
%                 else
%                     m_grid('fontsize', m_grid_fontsize, 'tickdir', m_grid_tickdir_type,  ...
%                         'xticklabels', [], 'xtick',[120,  130, 140], ...
%                         'yticklabels', [32, 36, 40, 44, 48], 'ytick',[32, 36, 40, 44, 48], 'backcolor', 'none', 'parent', ax{testnameind,2});  
%                 end
%                 txt{testnameind,2}=m_text(119, 50, test_text{testnameind}, 'FontSize', m_grid_fontsize+4); 
%                 txt{testnameind,3}=m_text(119, 48, test_text2{testnameind}, 'FontSize', m_grid_fontsize+4); 

                titlename = strcat('sp ground, ',testname, ',(', ...
                    num2str(min(allyear),'%04i'),'-', num2str(max(allyear),'%04i'), ',',  ...
                    num2str(min(inputmonth),'%02i'),'-', num2str(max(inputmonth),'%02i'), ')'); 

                title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title
                
                set(gcf, 'PaperUnits', 'points');
                set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
                set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

                saveas(gcf,jpgname,'tif');

                disp(' ')
                disp([fig_name, ' plot is created.'])
                disp(' ')
                disp([' File path is : ',jpgname])
                disp(' ')
                close all;
            end