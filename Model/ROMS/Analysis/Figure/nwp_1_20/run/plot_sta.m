clear all; clc; close all;

zeta_1980=ncread('E:\Data\Model\ROMS\nwp_1_20\test46\run\1980\sta.nc','zeta');
zeta_1984=ncread('E:\Data\Model\ROMS\nwp_1_20\test46\run\1984\sta.nc','zeta');
zeta_1986=ncread('E:\Data\Model\ROMS\nwp_1_20\test46\run\1986\sta.nc','zeta');
zeta_1987=ncread('E:\Data\Model\ROMS\nwp_1_20\test46\run\1987\sta.nc','zeta');

% plot(zeta_1980(1,:));
% plot(zeta_1987(1,:));
plot(zeta_1987(1,:)-zeta_1980(1,1:8761));

mean(zeta_1987(:,:)-zeta_1980(:,1:8761),2);
mean(zeta_1987(:,:),2)-mean(zeta_1984(:,:),2);

mean(zeta_1987(:,:),2)-mean(zeta_1986(:,:),2);

% mean(zeta_1984(:,:),2)-mean(zeta_1980(:,:),2);
% mean(zeta_1987(:,:),2)-mean(zeta_1980(:,:),2);

% % % % % taiwan strait, onshore transport
% % es1plot=plot(xData, aiwankur(1:12*year),'k');
% % hold on
% % es2plot=plot(xData, (-o_intruy(1:12*year)),'r');
% % datetick('x','yyyy')
% % % es3plot=plot(soyat(1:60),'b');
% % set(gca,'YLim',[-4 4]);
% % % set(gca,'XLim',[1 60]);
% % set(es1plot,'LineWidth',4);
% % set(es2plot,'LineWidth',4);
% % % set(es3plot,'LineWidth',4);
% % title('Monthly mean transports of the model results','fontsize',17);
% % xlabel('Time (month)','color','k','FontSize',17,'fontweight','bold');
% % ylabel('Volume Transport (sv)','color','k','FontSize',17,'fontweight','bold');
% % lgd=legend('Taiwan','Onshore');
% % set(lgd,'FontSize',15);
% % set(lgd,'Position',[0.13 0.87, 0.775, 0.03]);
% % set(lgd,'Orientation','horizontal');
% % set(gca,'FontSize',13);
% % grid on;
% % saveas(gcf,[figdir, 'trans_taiwan_model.png'],'png');
% % hold off;
