%%% plot observation data
clearvars; close all;

line_gap = 2;
% coast_fname = 'd:\bongan\4_data\coastline\Korea\coast_sin.dat';
% dum2 = load( coast_fname ); % coastline
% coa_lon = dum2(:,1);
% coa_lat = dum2(:,2);

for fig_title = {'1804_SST','1804_SSS','1807_SST','1807_SSS'}
    switch fig_title{1}(1:4)
        case '1804'
            load('1804_CTD_info.mat'); data_ctd = data18;
            title_str = {' 2018.04.21',' ~2018.04.29'};
            lat_lim = [31.0, 37.0];
            lon_lim = [123.0, 131.00];
        case '1807'
            load('1807_CTD_info.mat');
            title_str = {' 2018.07.17',' ~2018.07.25'};
            lat_lim = [34.5, 38.5];
            lon_lim = [127.0, 133.00];
    end

    x=data_ctd(:,2); y=data_ctd(:,1);
    switch fig_title{1}(end-2:end)
        case 'SST'
            v=data_ctd(:,9);
            c_label = 'SST (\circC)'; lev_step = 2; c_lim = [9.5,19.0];
            if( strcmp(fig_title{1}(1:4),'1807') ), c_lim = c_lim+8.2; end
        case 'SSS'
            v=data_ctd(:,10);
            c_label = 'SSS (P_s)'; lev_step = 0.4; c_lim = [32.45,34.65];
    end
    x_=[x; 126;   127.5; 126.8; 128.2; 130.3; 129.5; 131.89];
    y_=[y; 34.25; 34.75;  32.2;  33.7; 37.67; 37.01;  37.00];
    v_=[v; NaN;     NaN;   NaN;   NaN;   NaN;   NaN;    NaN];

    lon=lon_lim(1):0.01:lon_lim(2);
    lat=lat_lim(1):0.01:lat_lim(2);
    [X,Y]=meshgrid(lon,lat);

    V=griddata(x_,y_,v_,X,Y,'linear');

    figure('name',fig_title{1});
    pcolor(X,Y,V); shading flat; cax=colorbar; hold on; daspect([1 cosd(mean(lat_lim)) 1]);
    contour(X,Y,V,'showtext','on','textstep',lev_step,'linecolor','k');
    scatter(x,y,30,v,'markerfacecolor','flat','markeredgecolor','k');

    geoshow('landareas.shp','facecolor','none');
    % geoshow( coa_lat, coa_lon, 'DisplayType', 'polygon', 'facecolor', [0.7 0.7 0.7], 'EdgeColor', 'k' ); axis tight;
    xlim(lon_lim); ylim(lat_lim);

    caxis(c_lim);
    set(cax,'fontsize',14); ylabel(cax,c_label);
    set(gca,'position',get(gca,'position')+[0 0 0 0]);
    set(gca,'fontsize',15,'linewidth',2,'tickdir','both','ticklength',[.01 .01]);
    lon_tick = -180:line_gap:180;
    lon_tick = lon_tick( lon_tick >= lon_lim(1) & lon_tick <= lon_lim(2) );
    lat_tick = -90:line_gap:90;
    lat_tick = lat_tick( lat_tick >= lat_lim(1) & lat_tick <= lat_lim(2) );
    for ii = 1 : length(lon_tick)
        lon_label{ii} = [num2str(lon_tick(ii)),'\circE'];
    end
    for ii = 1 : length(lat_tick)
        lat_label{ii} = [num2str(lat_tick(ii)),'\circN'];
    end
    set(gca,'xtick',lon_tick,'ytick',lat_tick,'xticklabel',lon_label,'yticklabel',lat_label);
    text(lon_lim(1),lat_lim(2),title_str,'fontsize',15,'verticalalignment','top','horizontalalignment','left');

end