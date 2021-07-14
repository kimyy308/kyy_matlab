% % xlabel('Longitude (^oE)','fontsize',27);
% % ylabel('Latitude (^oN)','fontsize',27);

xlim([min(lon2(lon_min_ind:lon_max_ind)) max(lon2(lon_min_ind:lon_max_ind))]);
ylim([min(lat2(lat_min_ind:lat_max_ind)) max(lat2(lat_min_ind:lat_max_ind))]);
xdivvar=4.0;
ydivvar=6.0;
lontick_min_ind=10; %%128.5;
lontick_max_ind=50; %%132.5;
lattick_min_ind=30;  %%37.0;
lattick_max_ind=90;  %%43.0;
% xtickvar=lon2(lontick_min_ind):(lon2(lontick_max_ind)-lon2(lontick_min_ind))/xdivvar:lon2(lontick_max_ind);
% ytickvar=lat2(lattick_min_ind):(lat2(lattick_max_ind)-lat2(lattick_min_ind))/ydivvar:lat2(lattick_max_ind);
xtickvar=lon2(lontick_min_ind):(lon2(lontick_max_ind)-lon2(lontick_min_ind))/xdivvar:lon2(lontick_max_ind);
ytickvar=lat2(lattick_min_ind):(lat2(lattick_max_ind)-lat2(lattick_min_ind))/ydivvar:lat2(lattick_max_ind);
set(gca,...
    'XTickLabel',{'129^oE','','131^oE','','133^oE'},...
    'XTick',xtickvar);
set(gca,...
    'YTickLabel',{'37^oN','38^oN','39^oN','40^oN','41^oN', '42^oN', '43^oN'},...
        'YTick',ytickvar);
set(gca, 'TickDir','out');
% set(gca,...
%     'YTickLabel',{'36^oN','40^oN','44^oN'},...   ---> 36 40 44 36 40 (error)
%         'YTick',[36 38 40 42 44]);
% set(gcf,'PaperPosition',[0 0 25 15]);  %left, down, right, up
% axis equal;
% set(gca,'fontsize',15);
set(gca,'fontsize',20);
set(gca, 'color', [0.8,0.8,0.8]);


% if (pcolindex==1)
% % %     for horizontal vector
% %     tex=text(135.7,35.1,'0.2m/s');
% %     set(tex,'fontsize',7);
% %     fig=gcf;
%     tex=text(131.9,35.7,'0.2m/s');
%     set(tex,'fontsize',7);
%     fig=gcf;
% else
%     tex=text(135.9,35.2,'0.2m/s');
%     set(tex,'fontsize',10);
%     fig=gcf;
% end
% fig.InvertHardcopy='off';
% saveas(gcf,figname,'tiff');