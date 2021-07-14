% % xlabel('Longitude (^oE)','fontsize',27);
% % ylabel('Latitude (^oN)','fontsize',27);


set(gca,...
    'XTickLabel',{'130^oE','135^oE','140^oE'},...
    'XTick',[130 135 140]);
set(gca,...
    'YTickLabel',{'36^oN','38^oN','40^oN','42^oN','44^oN'},...
        'YTick',[36 38 40 42 44]);
% set(gca,...
%     'YTickLabel',{'36^oN','40^oN','44^oN'},...   ---> 36 40 44 36 40 (error)
%         'YTick',[36 38 40 42 44]);
set(gcf,'PaperPosition',[0 0 25 15]);  %left, down, right, up
axis equal;
% set(gca,'fontsize',15);
set(gca,'fontsize',12);
set(gca, 'color', [0.8,0.8,0.8]);


if (pcolindex==1)
% %     for horizontal vector
%     tex=text(135.7,35.1,'0.2m/s');
%     set(tex,'fontsize',7);
%     fig=gcf;
    tex=text(131.9,35.7,'0.2m/s');
    set(tex,'fontsize',7);
    fig=gcf;
else
    tex=text(135.9,35.2,'0.2m/s');
    set(tex,'fontsize',10);
    fig=gcf;
end
fig.InvertHardcopy='off';
saveas(gcf,figname,'tiff');