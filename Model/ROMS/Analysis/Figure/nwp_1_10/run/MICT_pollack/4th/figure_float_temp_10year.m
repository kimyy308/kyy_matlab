close all;clear all;clc
%% read data
% ltransdir = 'E:\Data\Model\ROMS\nwp_1_10\test06\DA\ltrans\';
% figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10_EnOI\LTRANS\temp\'
ltransdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\avg_ens_10km_mean\LTRANS\';
figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\avg_ens_10km_mean\LTRANS\temp\'

% totyear = [1981 1986];
totyear = [1981 1991 2001];

totmonth = 1:3;
numpar = 96 - 32;


for nmonth = 1:length(totmonth)
    hold on
    for nyear = 1:length(totyear)
        tempyear = totyear(nyear);
%         year = num2str(tempyear,'%04i');
    %     for j = [1 2 3];
        
    
        tempmonth = totmonth(nmonth);
        month = num2str(tempmonth,'%2.2d');
        
        for yind_10=1:10
            year = num2str(tempyear+yind_10-1,'%04i');
            filepath = [ltransdir, year, '/'];
            filename = [filepath,'output_',year,'_',month,'_01_0090.nc'];
            time = ncread(filename,'model_time');
            tempaa = ncread(filename,'temperature');
            rawaa=tempaa(1:numpar,121:end);
            aa=tempaa(1:numpar,121:end);  %% to remove initial shock -> 121 ~ end
            aa(aa<-2) = -2;
            for par = 1:numpar
                hour = 1;
                while (hour < size(aa,2))
                    diffaa=diff(aa(par,:));
                    if(abs(aa(par,hour)>20) || aa(par,hour) <-10)
                        aa(par,:)=NaN;
                    end
                    if abs(diffaa(hour))>10
                        aa(par,hour+1)=aa(par,hour);
    %                     diffbb=diff(aa(par,:));
    %                     diffaa(1:hour)=diffbb(1:hour);
    %                     hour = 1;
                        disp(['noise corrected', num2str(par)])
                    end
                    if (hour > 300)
                        if (mean(diff(aa(par,hour-300:hour)))==0)
                            aa(par,:)=NaN;
                        end
                    end
                    hour=hour+1;
                end
            end
            mean_aa2(yind_10,:,:) = mean(aa,'omitnan');
        end
        mean_aa=squeeze(mean(mean_aa2));
         
        leng = length(aa(1,:));
        t(1) = datenum(totyear(1),tempmonth,6,0,0,0);
        for i = 2:leng
            t(i) = t(i-1) + 1/24;
        end
        ftime = datestr(t);
        % figure
        % hold on
        % for i = 1:leng
        %     plot(aa(i,:))
        % end

        c1 = 2*ones(size(aa,2),1);
        c4 = 5*ones(size(aa,2),1);
        c2 = 10*ones(size(aa,2),1);
        c3 = 13*ones(size(aa,2),1);

        grid on
        %         for i = 1:length(aa(:,1))
        %             plot(t,aa(i,:),'Color',[0 0 0]+0.05*5,'linewidth',1,'markersize',30);
        %         end
        if nyear==1
            plot(t,mean_aa,'b','linewidth',6,'markersize',30);
%             plot(t,mean_aa,'b-','linewidth',6,'markersize',30);
        elseif nyear==2
            plot(t,mean_aa,'k','linewidth',6,'markersize',30);
%             plot(t,mean_aa,'k-','linewidth',6,'markersize',30);
        else
            plot(t,mean_aa,'g','linewidth',6,'markersize',30);
%             plot(t,mean_aa,'g-','linewidth',6,'markersize',30);
        end
        
        
    end
    
    savepath = strcat(figdir,'/');
    if (exist(savepath , 'dir') ~= 7)
        mkdir(savepath);
    end 
    plot(t,c1,'r-','linewidth',2,'markersize',30);
    plot(t,c2,'r-','linewidth',2,'markersize',30,'linestyle','--');
    plot(t,c3,'r-','linewidth',2,'markersize',30,'linestyle','--');
    plot(t,c4,'r-','linewidth',2,'markersize',30);
    line('XData',[t(24*24) t(24*24)],'YData',[-8 25],'color',[0.8 0.8 0.8],'linewidth',2,'linestyle','--');

    ylim([-8 25])
    set(gcf,'Position',[200 100 1200 600])
    % set(gca,'xtick',t(1):5:t(end))
    datetick('x','mm/dd','keepticks')
    xlim([t(1) t(end)])
    legend([num2str(totyear(1)),'-',num2str(totyear(1)+9)],[num2str(totyear(2)),'-',num2str(totyear(2)+9)],[num2str(totyear(3)),'-',num2str(totyear(end)+9)],'location','north','orientation','horizontal')
    ylabel('^oC','fontsize',20);
    set(gca,'fontsize',20)
    titlename = [month,' Particle Temperature Time Series'];
    title(titlename,'fontsize',30)
    set(gcf,'PaperPositionMode','auto')
    savename = [savepath,'ltrans_temp_10year_',num2str(totyear(1),'%04i'),'_',num2str(totyear(end),'%04i'),'_',month,'.png'];
    saveas(gcf,savename,'png');
    hold off
    close all;
end