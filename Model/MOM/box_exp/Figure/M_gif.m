clear all

topog_field=true;
vec_surf_field=true;
vec_100m_field=true;
vec_200m_field=true;
ssh_field=true;
temp_surf_field=true;
temp_100m_field=true;
temp_200m_field=true;
temp_129E_field=true;
temp_130E_field=true;
temp_37N_field=true;
temp_38N_field=true;
temp_39N_field=true;
temp_40N_field=true;
temp_41N_field=true;
vort_surf_field=true;
h_flux_vert_int_field=true;

for j=0:9
  readname='D:\need_to_presentation/oman_restore_04_18/197';
  readname2='12/ocean_snap_197';
  readname3='12.nc';
  rname=strcat(readname,num2str(j),readname2,num2str(j),readname3);
  %nc2 = netcdf.open(rname)
  u2_temp=ncread(rname,'u');
  v2_temp=ncread(rname,'v');
  temp2_temp=ncread(rname,'temp');
  eta_t2_temp=ncread(rname,'eta_t');
  lon2=ncread(rname,'xt_ocean');
  lat2=ncread(rname,'yt_ocean');
  level2=ncread(rname,'st_ocean');
  ilevel2 = -level2;
  for i=1:12
    u2(1:length(lon2),1:length(lat2),1:length(level2),j*12+i) = u2_temp(1:length(lon2),1:length(lat2),1:length(level2),i);
    v2(1:length(lon2),1:length(lat2),1:length(level2),j*12+i) = v2_temp(1:length(lon2),1:length(lat2),1:length(level2),i);
    temp2(1:length(lon2),1:length(lat2),1:length(level2),j*12+i) = temp2_temp(1:length(lon2),1:length(lat2),1:length(level2),i);
    eta_t2(1:length(lon2),1:length(lat2),j*12+i) =eta_t2_temp(1:length(lon2),1:length(lat2),i);
  end
end

for j=0:9
  readname='D:\need_to_presentation/oman_restore_04_18/197';
  readname2='12/ocean_diag_197';
  readname3='12.nc';
  rname=strcat(readname,num2str(j),readname2,num2str(j),readname3);
  w2_temp=ncread(rname,'wt');
  for i=1:12
    w2(1:length(lon2),1:length(lat2),1:length(level2),j*12+i) = w2_temp(1:length(lon2),1:length(lat2),1:length(level2),i);
  end
end

ilevel2=-level2(:);
[lat,lon] = meshgrid(lat2,lon2);
[lat3,ilevel3] = meshgrid(ilevel2(1:20),lat2);
[lon4,ilevel4] = meshgrid(ilevel2(1:15),lon2(1:25));

for i=1:120
    for j=1:110
        for k=1:30
            for l=1:120
                if (u2(i,j,k,l)< -1.0e+5)
                    u2(i,j,k,l)=NaN;
                end
                if (v2(i,j,k,l)< -1.0e+5)
                    v2(i,j,k,l)=NaN;
                end
                if (temp2(i,j,k,l)< -1.0e+5)
                    temp2(i,j,k,l)=NaN;
                end
                if (eta_t2(i,j,l)< -1.0e+5)
                    eta_t2(i,j,l)=NaN;
                end
            end
        end
    end
end

for i=1:120
    u2(80:84,14:17,1:30,i)=0.2;
    v2(80:84,14:17,1:30,i)=0.00001;
    v2(11,5:9,16,i)=0.2;
    w2(11,5:9,16,i)=0.00001;
end

%%
%% temperature(100m) & surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 100m Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 100m Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,10,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


%%
%% 100m temperature & vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_100m_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('100m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('100m Vector & Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_100m_field/197',num2str(year),num2str(month),'.png');
%         gifname = strcat('100m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('100m Vector & Temperature (%s year %2s month)',year,month);  
    end
    pcolor(lon2,lat2,(temp2(:,:,10,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
%     set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,10,l)*2,v2(1:3:120,1:3:110,10,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end

%%
%% 150m temperature & vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_150m_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('150m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('150m Vector & Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_150m_field/197',num2str(year),num2str(month),'.png');
%         gifname = strcat('150m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('150m Vector & Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,15,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
%     set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,15,l)*2,v2(1:3:120,1:3:110,15,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end





%%
%% 200m temperature & vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_200m_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('200m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('200m Vector & Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_200m_field/197',num2str(year),num2str(month),'.png');
%         gifname = strcat('200m Vector & Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('200m Vector & Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,20,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
%     set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,20,l)*2,v2(1:3:120,1:3:110,20,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


%%
%% 129E Temperature field
%%
ilevel2=-level2(:)
for l=1:120         
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/temp_129E_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Meridional Vector & Vertical Temperature along 129°E(', num2str(year), ' year ,', num2str(month), ' month)');
%         gifname = sprintf('Meridional Vector & Vertical Temperature along 129°E (%s year %2s month)',year,month);
        gifname = sprintf('Meridional Vector & Vertical Temperature along 129°E');
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/temp_129E_field/197',num2str(year),num2str(month),'.png');
%         gifname = strcat('Meridional Vector & Vertical Temperature along 129°E (', num2str(year), ' year ,', num2str(month), ' month)');
%         gifname = sprintf('Meridional Vector & Vertical Temperature along 129°E (%s year %2s month)',year,month);
        gifname = sprintf('Meridional Vector & Vertical Temperature along 129°E');
    end
    pcolor(lat2,ilevel2(1:20)/20,(reshape(temp2(11,1:110,1:20,l),110,20))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees North (°N)','fontsize',15);
    ylabel('Depth (M)','fontsize',15);
    %set(gcf,'Position',[300 300 1500 600]);   %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
    set(gcf,'PaperPosition',[0 0 30 12]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
    %axis equal;
    set(gca,'fontsize',25);
    quiver(ilevel3(1:5:110,1:20),lat3(1:5:110,1:20)/20,(reshape(v2(11,1:5:110,1:20,l),22,20))*3,(reshape(w2(11,1:5:110,1:20,l),22,20))*2,'-k','filled','AutoScale','off');
    hold off;
    set(gca,'yTickLabel',{-160:40:0})
    text(34.55,-8.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


%%
%% 37N Temperature field
%%

for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/temp_37N_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Meridional Velocity & Vertical Temperature along 37°N(', num2str(year), ' year ,', num2str(month), ' month)');
%         gifname = sprintf('Meridional Vector & Vertical Temperature along 37°N (%s year %2s month)',year,month);
        gifname = sprintf('Meridional Vector & Vertical Temperature along 37°N');
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/temp_37N_field/197',num2str(year),num2str(month),'.png');
%         gifname = strcat('Meridional Velocity & Vertical Temperature along 37°N (', num2str(year), ' year ,', num2str(month), ' month)');
%         gifname = sprintf('Meridional Vector & Vertical Temperature along 37°N (%s year %2s month)',year,month);
        gifname = sprintf('Meridional Vector & Vertical Temperature along 37°N');
    end
    pcolor(ilevel4,lon4,reshape(temp2(1:25,30,1:15,l),25,15));
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees North (°N)','fontsize',15);
    ylabel('Depth (M)','fontsize',15);
    %set(gcf,'Position',[300 300 1500 600]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
    set(gcf,'PaperPosition',[0 0 30 12]);  %왼쪽 여백, 위쪽 여백, 가로길이, 세로 길이
    %axis equal;
    set(gca,'fontsize',25);
   [C,h]=contour(ilevel4,lon4,reshape(v2(1:25,30,1:15,l),25,15),'ShowText','on');
    h.LineColor='black';
    hold off;
    saveas(gcf,filename,'png');
end

% [C,h]=contour(ilevel4,lon4,reshape(v2(1:25,30,1:15,l),25,15),'ShowText','on');
% h.LineColor='black';
% quiver(lat,lon,u2(:,:,1,1),v2(1:3:120,1:3:110,1,1),2);
% quiver(lat,lon,u2(:,:,1,2),v2(1:3:120,1:3:110,1,2),2);



%%
%% Sea Surface Height & vector field
%%
for l=120:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/ssh_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Sea Surface Height(%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/ssh_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Sea Surface Height (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(eta_t2(:,:,l))');
    hold on;
    shading interp;
%      caxis([-0.15,0.15]);
    c=colorbar;
    c.Label.String= 'SSH (m)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end




%%
%% temperature(surface) & surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,1,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end



%%
%% temperature(150m) & surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_150_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 150m Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_150_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 150m Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,15,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


%%
%% temperature(200m) & surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_200_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 200m Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_200_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 200m Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,20,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end






%%
%% temperature(70m) & surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_70_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 70m Temperature (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_70_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & 70m Temperature (%s year %2s month)',year,month);
    end
    pcolor(lon2,lat2,(temp2(:,:,7,l))');
    hold on;
    shading interp;
    caxis([0.2,20]);
    c=colorbar;
    c.Label.String= 'Temp (°C)';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap jet;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


%%
%% surface vector field
%%
for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_only_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vec_surf_only_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector (%s year %2s month)',year,month);
    end
    quiver(lon(1:3:120,1:3:110),lat(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end



%%
%% surface vorticity field
%%
rname='D:\need_to_presentation/oman_restore_04_18/vort.nc';
vort5=ncread(rname,'SURF_VORT');
lon5=ncread(rname,'XU_OCEAN');
lat5=ncread(rname,'YU_OCEAN');
level5=ncread(rname,'ST_OCEAN');
ilevel5 = -level5;

ilevel5=-level5(:);
[lat6,lon6] = meshgrid(lat5,lon5);

for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vort_surf_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Relative Vorticity (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vort_surf_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Relative Vorticity (%s year %2s month)',year,month);
    end
    pcolor(lon5,lat5,(vort5(:,:,1,l))');
    hold on;
    shading interp;
    caxis([-3e-5,3e-5]);
    c=colorbar;
    c.Label.String= 'Vorticity';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap bluewhitered;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon6(1:3:120,1:3:110),lat6(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end





%%
%% surface vorticity field + red dot
%%
rname='D:\need_to_presentation/oman_restore_04_18/vort.nc';
vort5=ncread(rname,'SURF_VORT');
lon5=ncread(rname,'XU_OCEAN');
lat5=ncread(rname,'YU_OCEAN');
level5=ncread(rname,'ST_OCEAN');
ilevel5 = -level5;

ilevel5=-level5(:);
[lat6,lon6] = meshgrid(lat5,lon5);

for l=1:120
    year=num2str(floor((l-1)/12));
    month=num2str(l-floor((l-1)/12)*12);
    if (str2num(month)<=9)
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vort_surf_field/197',num2str(year),'0',num2str(month));
%         gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Relative Vorticity (%s year %2s month)',year,month);
    else
        filename= strcat('D:\need_to_presentation/oman_restore_04_18/MATLAB_gif/vort_surf_field/197',num2str(year),num2str(month),'.png');
 %       gifname = strcat('Surface Vector & 100m Temperature (', num2str(year), ' year ,', num2str(month), ' month)');
        gifname = sprintf('Surface Vector & Relative Vorticity (%s year %2s month)',year,month);
    end
    pcolor(lon5,lat5,(vort5(:,:,1,l))');
    hold on;
    shading interp;
    caxis([-3e-5,3e-5]);
    c=colorbar;
    c.Label.String= 'Vorticity';
    c.Label.FontSize= 12;
    % c.Limits(1)=0.2;
    % c.Limits(2)=20.;
    colormap bluewhitered;
    title(gifname,'fontsize',15);
    xlabel('Degrees East (°E)','fontsize',15);
    ylabel('Degrees North (°N)','fontsize',15);
    %set(gcf,'Position',[800 300 800 600]);  %left, down, right, up
    set(gcf,'PaperPosition',[0 0 20 15]);  %left, down, right, up
    axis equal;
    set(gca,'fontsize',15);
    quiver(lon6(1:3:120,1:3:110),lat6(1:3:120,1:3:110),u2(1:3:120,1:3:110,1,l)*2,v2(1:3:120,1:3:110,1,l)*2,'k','AutoScale','off');
    plot(5,20,'r.','MarkerSize',20)
    hold off;
    text(135.9,35.2,'0.2m/s')
    saveas(gcf,filename,'png');
end


