clc;clear all;close all

location=7;

high=['G:\grd2_EP\2013_500\viscosity500.txt'];
hdata=load(high);
hdata=hdata(:,location);
sdata1=load('smagorinsky100_spinup2.txt');
sdata1=sdata1(:,location);
sdata2=load('smagorinsky100_spinup3.txt');
sdata2=sdata2(:,location);
sdata3=load('smagorinsky100_spinup4.txt');
sdata3=sdata3(:,location);
data_leng=length(sdata3);
figure(1)
figure('position',[500 300 1100 550],'PaperUnits','inches','PaperPosition',[0 0 11 5]);
set(gca,'Position',[0.08 0.12 .89 0.8]);
hold on
plot([1:1:data_leng],sdata1(:,1),'r','linewidth',1.5);
plot([1:1:data_leng],sdata2(:,1),'m','linewidth',1.5);
plot([1:1:data_leng],sdata3(:,1),'g','linewidth',1.5);
plot([1:1:data_leng],hdata(:,1),'k','linewidth',1.5);

 set(gca,'fontsize',17,'box','on','xgrid','on','gridlinestyle','--');
 legend_str{1}='Smagorinsky-spinup2'; 
 legend_str{2}='Smagorinsky-spinup3'; 
 legend_str{3}='Smagorinsky-spinup4'; 
 legend_str{4}='visc.500';

  
le1=legend(legend_str,'Location', 'NorthWest');
xlim([1 12])
  set(le1,'fontsize',13)
  location_name={'Korea Strait' 'ruku' 'Taiwan' 'kuroshio' 'tsugaru' 'Soya' 'Yellow Sea'};
  title([location_name{location},' volume transport'],'fontsize',17);
   xlabel('Month','fontsize',18);ylabel('Transport(Sv)','fontsize',18);
 saveas(gcf,[location_name{location},'transport'],'tif')
