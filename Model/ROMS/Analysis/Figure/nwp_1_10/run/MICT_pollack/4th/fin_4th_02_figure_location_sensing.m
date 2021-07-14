close all; clear all; clc;


testname='nwp_1_10_EnOI';
% testname='avg_ens_10km_mean';

% if (strcmp(testname,'avg_ens_10km_mean')==1)
    filename=['E:\Data\Model\LTRANS\output\test06_DA_3rd\nwp_1_10_EnOI_comb_res_location_1983_2019.mat'];
% else
%     filename=[testname,'_comb_res_location.mat'];
% end

load(filename);

for nyear=1:length(totyear)
    for nmonth=1:length(totmonth)
        year=totyear(nyear);
        month=totmonth(nmonth);
        hold on;
        tt=1:t-1;
        ttt=1:24:t-25;
%         egg1plot = plot(tt/24,squeeze(wb_egg_in(nyear,nmonth,:)), 'color', [0/255, 114/255, 189/255]);
% % % %         % egg2plot = plot(tt/24,squeeze(skc_egg_in(1,1,:)), 'color', [217/255, 83/255, 25/255]);
%         egg2plot = plot(tt/24,squeeze(skc_egg_in(nyear,nmonth,:)), 'color', 'k');
%         egg3plot = plot(tt/24,squeeze(nkc_egg_in(nyear,nmonth,:)), 'color', [237/255 117/255 32/255]);
% % % %         egg4plot = plot(tt/24,squeeze(off_egg_in(nyear,nmonth,:)), 'color', [0.8 0.8 0.8]);
%         egg4plot = plot(tt/24,squeeze(off_egg_in(nyear,nmonth,:)), 'color', 'g');
%         set(egg1plot,'LineWidth',4);
%         set(egg2plot,'LineWidth',3);
%         set(egg3plot,'LineWidth',3);
%         set(egg4plot,'LineWidth',3);
        
        eggcolor(1,:)=[0/255, 114/255, 189/255];
        eggcolor(2,:)=[255/255, 255/255, 255/255];
        eggcolor(3,:)=[237/255 117/255 32/255];
        eggcolor(4,:)=[0/255 255/255 0/255];
        
        egg_in_all(nyear,nmonth,:,1)=wb_egg_in(nyear,nmonth,:);
        egg_in_all(nyear,nmonth,:,2)=skc_egg_in(nyear,nmonth,:);
        egg_in_all(nyear,nmonth,:,3)=nkc_egg_in(nyear,nmonth,:);
        egg_in_all(nyear,nmonth,:,4)=off_egg_in(nyear,nmonth,:);
%         plotbar=bar(tt/24,squeeze(egg_in_all(nyear,nmonth,:,:)),'stacked')
        plotbar=bar(ttt/24,squeeze(egg_in_all(nyear,nmonth,1:24:size(egg_in_all,3),:)),'stacked');

        line_15d = line('XData',[15 15], 'YData',[0 63], 'Color','k','LineStyle','--');
        line_30d = line('XData',[30 30], 'YData',[0 63], 'Color','k','LineStyle','--');

        set(gca,'YLim',[0 64]);
        set(line_15d,'LineWidth',6);
        set(line_30d,'LineWidth',6);


        title(['Egg location, ',num2str(year,'%04i'),', ',num2str(month,'%02i')],'fontsize',17);
        xlabel('Days since spawning','color','k','FontSize',17,'fontweight','bold');
        ylabel('The number of Egg','color','k','FontSize',17,'fontweight','bold');
        lgd=legend('Wonsan Bay','SK Coast', 'NK Coast', 'offshore');
        set(lgd,'FontSize',13);
        set(gcf,'PaperPosition', [0 0 36 18]) 
        set(lgd,'Orientation','horizontal');
        set(gca,'FontSize',15);
        grid on;
        grid minor;
        figdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\',testname,'\LTRANS\location\'); % % where figure files will be saved
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end
        hold off;
        drawnow;
        
        saveas(gcf,[figdir, 'ltrans_location_',num2str(year,'%04i'),'_',num2str(month,'%02i'),'.png'],'png');
        close all;
        
        textfile=[figdir, 'ltrans_location_',num2str(year,'%04i'),'_',num2str(month,'%02i'), '.dat'];
        fid=fopen(textfile,'w')
        fprintf(fid, 'Time(hour)  Wonsan bay  South       North       Offshore');
        formatSpec = '%12d %12.2f %12.2f %12.2f %12.2f\n';
        time=tt(25:end)-24;
        fprintf(fid,'\n');
        for iii=1:length(time)
            tempmat(1)=time(iii);
            tempmat(2:5)=egg_in_all(nyear,nmonth,iii,:);
            fprintf(fid,formatSpec,tempmat);
        end
        fclose(fid)
    end
    
        year=totyear(nyear);
%         month=totmonth(nmonth);
        hold on;
        mean_wb_egg_in=mean(wb_egg_in(nyear,:,:),2);
        mean_skc_egg_in=mean(skc_egg_in(nyear,:,:),2);
        mean_nkc_egg_in=mean(nkc_egg_in(nyear,:,:),2);
        mean_off_egg_in=mean(off_egg_in(nyear,:,:),2);
        tt=1:t-1;
%         egg1plot = plot(tt/24,squeeze(mean_wb_egg_in(1,1,:)), 'color', [0/255, 114/255, 189/255]);
%         % egg2plot = plot(tt/24,squeeze(skc_egg_in(1,1,:)), 'color', [217/255, 83/255, 25/255]);
%         egg2plot = plot(tt/24,squeeze(mean_skc_egg_in(1,1,:)), 'color', 'k');
%         egg3plot = plot(tt/24,squeeze(mean_nkc_egg_in(1,1,:)), 'color', [237/255 117/255 32/255]);
% % %         egg4plot = plot(tt/24,squeeze(off_egg_in(nyear,nmonth,:)), 'color', [0.8 0.8 0.8]);
%         egg4plot = plot(tt/24,squeeze(mean_off_egg_in(1,1,:)), 'color', 'g');
%         set(egg1plot,'LineWidth',4);
%         set(egg2plot,'LineWidth',3);
%         set(egg3plot,'LineWidth',3);
%         set(egg4plot,'LineWidth',3);
        
        mean_egg_in_all(1,1,:,1)=mean_wb_egg_in(1,1,:);
        mean_egg_in_all(1,1,:,2)=mean_skc_egg_in(1,1,:);
        mean_egg_in_all(1,1,:,3)=mean_nkc_egg_in(1,1,:);
        mean_egg_in_all(1,1,:,4)=mean_off_egg_in(1,1,:);
%         bar(tt/24,squeeze(mean_egg_in_all(1,1,:,:)),'stacked');
        plotbar=bar(ttt/24,squeeze(mean_egg_in_all(1,1,1:24:size(mean_egg_in_all,3),:)),'stacked');

        line_15d = line('XData',[15 15], 'YData',[0 63], 'Color','k','LineStyle','--');
        line_30d = line('XData',[30 30], 'YData',[0 63], 'Color','k','LineStyle','--');

        set(gca,'YLim',[0 64]);
        set(line_15d,'LineWidth',6);
        set(line_30d,'LineWidth',6);
        
        title(['Egg location, ',num2str(year,'%04i'),', mean'],'fontsize',17);
        xlabel('Days since spawning','color','k','FontSize',17,'fontweight','bold');
        ylabel('The number of Egg','color','k','FontSize',17,'fontweight','bold');
        lgd=legend('Wonsan Bay','SK Coast', 'NK Coast', 'offshore');
        set(lgd,'FontSize',13);
        set(gcf,'PaperPosition', [0 0 36 18]) 
        set(lgd,'Orientation','horizontal');
        set(gca,'FontSize',15);
        grid on;
        grid minor;
        figdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\4th_year\Figure\',testname,'\LTRANS\location\'); % % where figure files will be saved
        if (exist(strcat(figdir) , 'dir') ~= 7)
            mkdir(strcat(figdir));
        end
        hold off;
        drawnow;
        textfile=[figdir, 'mean_ltrans_location_',num2str(year,'%04i'),'.dat'];
        fid=fopen(textfile,'w')
        fprintf(fid, 'Time(hour)  Wonsan bay  South       North       Offshore');
        formatSpec = '%12d %12.2f %12.2f %12.2f %12.2f\n';
        time=tt(25:end)-24;
        fprintf(fid,'\n');
        for iii=1:length(time)
            tempmat(1)=time(iii);
            tempmat(2:5)=mean_egg_in_all(1,1,iii,:);
            fprintf(fid,formatSpec,tempmat);
        end
        fclose(fid)
        saveas(gcf,[figdir, 'mean_ltrans_location_',num2str(year,'%04i'),'.png'],'png');
        close all;
end

% xData((12*(i-1))+month) = datenum([num2str(month,'%02i'),'-01-',tempyear]);
% datetick('x','yy','keeplimits')
% axis tight;

% set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);

mean(mean(wb_egg_in(1:8,:,2159)))
mean(mean(wb_egg_in(9:18,:,2159)))
% plot(mean(wb_egg_in(:,:,2159),2))