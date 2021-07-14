close all;clear all;clc
%% read data



ltransdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10_EnOI\LTRANS\';
figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\nwp_1_10_EnOI\LTRANS\temp\'
% ltransdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\avg_ens_10km_mean\LTRANS\';
% figdir = 'D:\OneDrive - 서울대학교\MEPL\project\MICT_pollack\3rd_year\figure\avg_ens_10km_mean\LTRANS\temp\'

totyear = [1982 1987 2007];
% totyear = [1983 1987 1991 1995];

% mean_range = mean(diff(totyear)); 
mean_range = 5; 
totday = 30;

totmonth = 1:3;
numpar = 96 - 32;
 

for include_std=0:0
    for nmonth = 1:length(totmonth)
        hold on 
        for nyear = 1:length(totyear)

            tempyear = totyear(nyear);
%             heatmap1=1;
%             heatmap2=1-(nyear-1)*1/(length(totyear)-1);
%             heatmap3=1-(nyear-1)*1/(length(totyear)-1);
            jet_cbar=jet(round((max(totyear)-min(totyear))/mean_range+1));

            heatmap1=jet_cbar((tempyear-totyear(1))/mean_range+1,1);
            heatmap2=jet_cbar((tempyear-totyear(1))/mean_range+1,2);
            heatmap3=jet_cbar((tempyear-totyear(1))/mean_range+1,3);
    %         year = num2str(tempyear,'%04i');
        %     for j = [1 2 3];
            tempmonth = totmonth(nmonth);
            month = num2str(tempmonth,'%2.2d');

            for yind=1:mean_range
                year = num2str(tempyear+yind-1,'%04i');
                filepath = [ltransdir, year, '/'];
                filename = [filepath,'output_',year,'_',month,'_01_0090.nc'];
                time = ncread(filename,'model_time');
                tempaa = ncread(filename,'temperature');
                rawaa=tempaa(1:numpar,121:end);  
                aa=tempaa(1:numpar,121:end);%% remove noise during first 5 days
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
                mean_aa2(yind,:,:) = mean(aa,1,'omitnan');
                comb_aa2(yind,:,:) = aa;
            end
            std_aa=squeeze(std(reshape(comb_aa2,[mean_range*numpar size(aa,2)]),'omitnan'));
            mean_aa=squeeze(mean(mean_aa2));

%             leng = length(aa(1,:));
            leng = 24 * totday - 120;
            mean_aa=mean_aa(1:leng);
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
% 
%             c1 = 2*ones(size(aa,2),1);
%             c4 = 5*ones(size(aa,2),1);
%             c2 = 10*ones(size(aa,2),1);
%             c3 = 13*ones(size(aa,2),1);

            c1 = 2*ones(leng,1);
            c4 = 5*ones(leng,1);
            c2 = 10*ones(leng,1);
            c3 = 13*ones(leng,1);
            
            grid on
            %         for i = 1:length(aa(:,1))
            %             plot(t,aa(i,:),'Color',[0 0 0]+0.05*5,'linewidth',1,'markersize',30);
            %         end
            tinterval=24;
            
%             mean(mean_aa(1:680))

            if nyear==1
    %             p1=plot(t,mean_aa,'b','linewidth',6,'markersize',30);
                disp(num2str(nyear))
                disp(num2str(nmonth))
                disp(mean(mean_aa(1:tinterval:end)))
                p(nyear)=plot(t,mean_aa,'color','b','linewidth',6,'markersize',30);
                if (include_std==1)
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)-std_aa(1:tinterval:end)','b','linewidth',1,'linestyle',':');
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)+std_aa(1:tinterval:end)','b','linewidth',1,'linestyle',':');
                end
    %             errorbar(t(1:tinterval:end),mean_aa(1:tinterval:end),std_aa(1:tinterval:end),'b','linewidth',2,'markersize',30);
            elseif nyear==2
                disp(num2str(nyear))
                disp(num2str(nmonth))
                disp(mean(mean_aa(1:tinterval:end)))
    %             p2=plot(t,mean_aa,'k','linewidth',6,'markersize',30);
                p(nyear)=plot(t,mean_aa,'color','g','linewidth',6,'markersize',30);
                if (include_std==1)
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)-std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)+std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                end
    %             errorbar(t(1:tinterval:end),mean_aa(1:tinterval:end),std_aa(1:tinterval:end),'color',[heatmap1 heatmap2 heatmap3],'linewidth',2,'markersize',30);
             elseif nyear==3
                 disp(num2str(nyear))
                disp(num2str(nmonth))
                disp(mean(mean_aa(1:tinterval:end)))
    %             p2=plot(t,mean_aa,'k','linewidth',6,'markersize',30);
                p(nyear)=plot(t,mean_aa,'color','r','linewidth',6,'markersize',30);
                if (include_std==1)
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)-std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)+std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                end
    %             errorbar(t(1:tinterval:end),mean_aa(1:tinterval:end),std_aa(1:tinterval:end),'color',[heatmap1
            else
                disp(num2str(nyear))
                disp(num2str(nmonth))
                disp(mean(mean_aa(1:tinterval:end)))
    %             p3=plot(t,mean_aa,'g','linewidth',6,'markersize',30);
                p(nyear)=plot(t,mean_aa,'color',[heatmap1 heatmap2 heatmap3],'linewidth',6,'markersize',30);
                if (include_std==1)
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)-std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                    plot(t(1:tinterval:end),mean_aa(1:tinterval:end)+std_aa(1:tinterval:end)','color',[heatmap1 heatmap2 heatmap3],'linewidth',1,'linestyle',':');
                end
    %             errorbar(t(1:tinterval:end),mean_aa(1:tinterval:end),std_aa(1:tinterval:end),'color',[heatmap1 heatmap2 heatmap3],'linewidth',2,'markersize',30);
            end
            tempplot_legend{nyear}=[num2str(totyear(nyear)),'-',num2str(totyear(nyear)+mean_range-1)];
        end

        savepath = strcat(figdir,'/');
        if (exist(savepath , 'dir') ~= 7)
            mkdir(savepath);
        end 
        plot(t,c1,'k-','linewidth',2,'markersize',30);
        plot(t,c2,'k-','linewidth',2,'markersize',30,'linestyle','--');
        plot(t,c3,'k-','linewidth',2,'markersize',30,'linestyle','--');
        plot(t,c4,'k-','linewidth',2,'markersize',30);
        ymax = 13.0;  ymin = 0.0;
        line('XData',[t(24*24) t(24*24)],'YData',[ymin ymax],'color',[0.8 0.8 0.8],'linewidth',2,'linestyle','--');
        ylim([ymin ymax])
        set(gcf,'Position',[200 100 1200 600])
        % set(gca,'xtick',t(1):5:t(end))
        datetick('x','mm/dd','keepticks')
        xlim([t(1) t(end)])
    %     legend([p1 p2 p3],{[num2str(totyear(1)),'-',num2str(totyear(1)+4)],[num2str(totyear(2)),'-',num2str(totyear(2)+4)],[num2str(totyear(3)),'-',num2str(totyear(end)+4)]},'location','north','orientation','horizontal')
        lgd=legend(tempplot_legend,'location','north','orientation','horizontal')
        if (include_std==1)
            legend([p(1:end)],tempplot_legend,'location','north','orientation','horizontal')
        end
        ylabel('temperature(^oC)','fontsize',30);
        xlabel('date (mm/dd)','fontsize',30);
        set(gca,'fontsize',30)
        set(lgd,'fontsize',10)
        titlename = [month,' Particle Temperature Time Series'];
        title(titlename,'fontsize',30)
        set(gcf,'PaperPositionMode','auto')
        if (include_std==1)
            savename = [savepath,'jet_ltrans_temp_std_',num2str(mean_range,'%02i'),'year_',num2str(totyear(1),'%04i'),'_',num2str(totyear(end),'%04i'),'_',month,'_', num2str(totday,'%02i'),'.png'];
        else
            savename = [savepath,'jet_ltrans_temp_',num2str(mean_range,'%02i'),'year_',num2str(totyear(1),'%04i'),'_',num2str(totyear(end),'%04i'),'_',month,'_', num2str(totday,'%02i'),'.png'];
        end
        saveas(gcf,savename,'png');

        hold off
        close all;
    end
end