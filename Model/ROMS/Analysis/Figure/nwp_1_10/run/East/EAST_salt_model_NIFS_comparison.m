close all; clear all;  clc;
warning off;

close all;
clearvars '*' -except regionind2 all_region2
% % % 

system_name=computer;
if (strcmp(system_name,'PCWIN64'))
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox';
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
    addpath(genpath([matlabroot,'\toolbox\matlab\imagesci\'])); %% add new netcdf path
elseif (strcmp(system_name,'GLNXA64'))
    dropboxpath='/home/kimyy/Dropbox';
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
end

shadlev = [0 35];
rms_shadlev = [0 4];
%     trendlev = [-3 3];  %% trend lev
trendlev = [-0.5 0.5];  %% trend lev
conlev  = 0:5:35;
meanplotlev =[25 35];

% for snu_desktopd
testname='test08'   % % need to change
inputyear = [2009:2017]; % % put year which you want to plot [year year ...]
inputmonth = [2 4 6 8 10 12]; % % put month which you want to plot [month month ...]

varname ='salt'
run('nwp_polygon_point.m');


% % % for EKB
% regionname='EKB';
% lonlat = [127, 129.5, 38, 40.5];

load(['E:\Data\Observation\NIFS\xls\',testname,'_nwp_1_10_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'),'.mat']);

% valnum=0;
% % %     valid cell number
%  for vi=1:size(comb_spatial_meanressh,1)
%      for vj=1:size(comb_spatial_meanressh,2)
%          if (isnan(comb_spatial_meanressh(vi,vj,1))~=1)
%             valnum=valnum+1;
%          end
%      end
%  end

if (strcmp(system_name,'PCWIN64'))
    % % for windows
    figrawdir =strcat('D:\OneDrive - 서울대학교\MEPL\project\EAST\3rd_year\figure\nwp_1_10\',testname,'\NIFS\'); % % where figure files will be saved
    param_script ='C:\Users\kyy\Dropbox\source\matlab\Model\ROMS\Analysis\Figure\nwp_1_10\run\fig_param\fig_param_kyy_ECS_CWD.m'
    filedir = strcat('G:\nwp_1_10\output\', testname, '\run\'); % % where data files are
    NIFSdir='E:\Data\Observation\NIFS\xls\';
    run(param_script);
elseif (strcmp(system_name,'GLNXA64'))
end

figdir=[figrawdir,'Trend\'];
if (exist(strcat(figdir) , 'dir') ~= 7)
    mkdir(strcat(figdir));
end 
outfile = strcat(figdir);

NIFS6.lat(NIFS6.lat==0)=NaN;
NIFS6.lon(NIFS6.lon==0)=NaN;
maxlat=max(max(max(max(max(NIFS6.lat)))));
minlat=min(min(min(min(min(NIFS6.lat)))));
maxlon=max(max(max(max(max(NIFS6.lon)))));
minlon=min(min(min(min(min(NIFS6.lon)))));

refpolygon=[maxlon,maxlat; ...
            maxlon, minlat; ...
            minlon, minlat; ...
            minlon,maxlat];
lonlat(1)=min(refpolygon(:,1));
lonlat(2)=max(refpolygon(:,1));
lonlat(3)=min(refpolygon(:,2));
lonlat(4)=max(refpolygon(:,2));

% bwr10000(1:5000,1)=0.0002:0.0002:1;
% bwr10000(1:5000,2)=0.0002:0.0002:1;
% bwr10000(1:5000,3)=1;
% bwr10000(5001:10000,1)=1;
% bwr10000(5001:10000,2)=1:-0.0002:0.0002;
% bwr10000(5001:10000,3)=1:-0.0002:0.0002;

% NIFS surface salt trend plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    for i=1:size(NIFS6.salt,2)  %% line
        for j=1:size(NIFS6.salt,3)  %% sta
            if (isnan(NIFS6.salt_surf_comb_trend(i,j))==1)
            else
                if (NIFS6.salt_surf_comb_trend(i,j)<=min(trendlev))
                    colind=1;
                elseif (NIFS6.salt_surf_comb_trend(i,j)>=max(trendlev))
                    colind=100;
                else
                    colind=round((NIFS6.salt_surf_comb_trend(i,j)-min(trendlev))/diff(trendlev) *100.0);
                end
                round(NIFS6.salt_surf_comb_trend(i,j) ,3 );
                if (sum(sum(~isnan(NIFS6.salt(:,i,j,:,1)))) < sum(sum(isnan(NIFS6.salt(:,i,j,:,1)))))
                    NIFS6.salt_surf_comb_trend(i,j)=NaN;
                end
                h1=m_plot(NIFS6.lon(1,i,j,1,1),NIFS6.lat(1,i,j,1,1),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:)); 
            end
        end
    end
    mean_trend=mean(mean(NIFS6.salt_surf_comb_trend,'omitnan'),'omitnan');
%     std(std(NIFS6.salt_surf_comb_trend,'omitnan'),'omitnan')
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
     hold off;

    titlename = strcat('Salt trend, ','NIFS',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend,4)));  %% + glacier contribution

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'/y','fontsize',colorbar_title_fontsize);
    caxis(trendlev);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, 'NIFS_', 'surface_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;

    
% Model surface salt trend plot

    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    for i=1:size(MODEL6.salt,2)  %% line
        for j=1:size(MODEL6.salt,3)  %% sta
            if (isnan(MODEL6.salt_surf_comb_trend(i,j))==1)
            else
                if (MODEL6.salt_surf_comb_trend(i,j)<=min(trendlev))
                    colind=1;
                elseif (MODEL6.salt_surf_comb_trend(i,j)>=max(trendlev))
                    colind=100;
                else
                    colind=round((MODEL6.salt_surf_comb_trend(i,j)-min(trendlev))/diff(trendlev) *100.0);
                end
                round(MODEL6.salt_surf_comb_trend(i,j) ,3 );
                if (sum(sum(~isnan(MODEL6.salt(:,i,j,:,1)))) < sum(sum(isnan(MODEL6.salt(:,i,j,:,1)))))
                    MODEL6.salt_surf_comb_trend(i,j)=NaN;
                end
                h1=m_plot(NIFS6.lon(1,i,j,1,1),NIFS6.lat(1,i,j,1,1),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:)); 
            end
        end
    end
    mean_trend=mean(mean(MODEL6.salt_surf_comb_trend,'omitnan'),'omitnan');
%     std(std(MODEL6.salt_surf_comb_trend,'omitnan'),'omitnan')
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
     hold off;

    titlename = strcat('Salt trend, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend,4)));  %% + glacier contribution

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'/y','fontsize',colorbar_title_fontsize);
    caxis(trendlev);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, testname,'_', 'surface_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;
    
    
% NIFS 75m salt trend plot
    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    for i=1:size(NIFS6.salt,2)  %% line
        for j=1:size(NIFS6.salt,3)  %% sta
            if (isnan(NIFS6.salt75_comb_trend(i,j))==1)
            else
                if (NIFS6.salt75_comb_trend(i,j)<=min(trendlev))
                    colind=1;
                elseif (NIFS6.salt75_comb_trend(i,j)>=max(trendlev))
                    colind=100;
                else
                    colind=round((NIFS6.salt75_comb_trend(i,j)-min(trendlev))/diff(trendlev) *100.0);
                end
                round(NIFS6.salt75_comb_trend(i,j) ,3 );
                if (sum(sum(~isnan(NIFS6.salt(:,i,j,:,1)))) < sum(sum(isnan(NIFS6.salt(:,i,j,:,1)))))
                    NIFS6.salt75_comb_trend(i,j)=NaN;
                end
                h1=m_plot(NIFS6.lon(1,i,j,1,1),NIFS6.lat(1,i,j,1,1),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:)); 
            end
        end
    end
    mean_trend=mean(mean(NIFS6.salt75_comb_trend,'omitnan'),'omitnan');
%     std(std(NIFS6.salt75_comb_trend,'omitnan'),'omitnan')
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
     hold off;

    titlename = strcat('Salt trend, ','NIFS',',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend,4)));  %% + glacier contribution

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'/y','fontsize',colorbar_title_fontsize);
    caxis(trendlev);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, 'NIFS_', '75m_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;

    
% Model 75m salt trend plot

    m_proj(m_proj_name,'lon',[lonlat(1)-0.5 lonlat(2)+0.5],'lat',[lonlat(3)-0.5 lonlat(4)+0.5]);
    hold on;
    for i=1:size(MODEL6.salt,2)  %% line
        for j=1:size(MODEL6.salt,3)  %% sta
            if (isnan(MODEL6.salt75_comb_trend(i,j))==1)
            else
                if (MODEL6.salt75_comb_trend(i,j)<=min(trendlev))
                    colind=1;
                elseif (MODEL6.salt75_comb_trend(i,j)>=max(trendlev))
                    colind=100;
                else
                    colind=round((MODEL6.salt75_comb_trend(i,j)-min(trendlev))/diff(trendlev) *100.0);
                end
                round(MODEL6.salt75_comb_trend(i,j) ,3 );
                if (sum(sum(~isnan(MODEL6.salt(:,i,j,:,1)))) < sum(sum(isnan(MODEL6.salt(:,i,j,:,1)))))
                    MODEL6.salt75_comb_trend(i,j)=NaN;
                end
                h1=m_plot(NIFS6.lon(1,i,j,1,1),NIFS6.lat(1,i,j,1,1),'marker','o','color','k', ...
                  'markerfacecolor',bwrmap(colind,:)); 
            end
        end
    end
    mean_trend=mean(mean(MODEL6.salt75_comb_trend,'omitnan'),'omitnan');
%     std(std(MODEL6.salt75_comb_trend,'omitnan'),'omitnan')
    m_gshhs_i('color',m_gshhs_line_color);
    m_gshhs_i('patch',m_gshhs_land_color);   % gray colored land
     hold off;

    titlename = strcat('Salt trend, ',testname,',(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ','M=',num2str(round(mean_trend,4)));  %% + glacier contribution

    title(titlename,'fontsize',m_pcolor_title_fontsize);  %%title

    % set colorbar 
    h = colorbar;
    colormap(bwrmap);
    set(h,'fontsize',colorbar_fontsize);
    title(h,'/y','fontsize',colorbar_title_fontsize);
    caxis(trendlev);

    % set grid
    m_grid('fontsize', m_grid_fontsize, 'box', m_grid_box_type, 'tickdir', m_grid_tickdir_type);

%     set(gcf, 'PaperUnits', 'points');
%     set(gcf, 'PaperSize', [hor_paper_size_x, hor_paper_size_y]);
%     set(gcf,'PaperPosition', [paper_position_hor paper_position_ver paper_position_width paper_position_height]) 

    jpgname=strcat(outfile, testname,'_', '75m_salt_trend_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.tif'); %% ~_year_month.jpg
    saveas(gcf,jpgname,'tif');

    disp(' ')
    disp([' File path is : ',jpgname])
    disp(' ')

    hold off
    close all;    
    
    

for i =1:length(inputyear) 
    tempyear=inputyear(i);
    for month=1:6
        xData((6*(i-1))+month) = datenum([num2str(tempyear),'-',num2str(month*2,'%02i'),'-15',]);
    end
end

for numline=1:size(NIFS6.line,2)
    if (isnan(sum(NIFS6.salt_surf_line_comb(:,numline),'omitnan'))==1)

    else
        lineline=unique(NIFS6.line(:,numline,:,:,:),'stable');
        lineline2=lineline(~isnan(lineline));
        lineline3=lineline2(lineline2>0);

        idx=isnan(NIFS6.salt_surf_line_comb(:,numline)');
        p=polyfit(xData(~idx),NIFS6.salt_surf_line_comb(~idx,numline)',1);
        line_trend2=xData*p(1)+p(2);
        msaltplot=plot(xData,MODEL6.salt_surf_line_comb(:,numline)','color','k');
        hold on
        msaltplot2=plot(xData,NIFS6.salt_surf_line_comb(:,numline)','o','color','r','markerfacecolor','k','markersize',2);

%         msaltplot2=plot(xData,line_trend2,'Color','r')
        idx2=isnan(MODEL6.salt_surf_line_comb(:,numline)');
        idx3=idx+idx2;
        idx3(idx3==2)=1;
        constant_cor=corrcoef(NIFS6.salt_surf_line_comb(~idx3,numline)',MODEL6.salt_surf_line_comb(~idx3,numline)');
        
        salt_bias=MODEL6.salt_surf_line_comb(~idx3,numline)' - NIFS6.salt_surf_line_comb(~idx3,numline)';
        mean_salt_bias=mean(salt_bias);
        mean_salt_RMS=mean(abs(salt_bias));
        hold off
        jpgname=strcat(outfile,testname,'_',num2str(lineline3), '_surf_line_mean_salt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('year')
        ylabel('Mean Salinity along NIFS line')
        title([num2str(lineline3), ', surf mean salt(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(NIFS6.salt_surf_line_comb_trend(numline),4)), ' /y'])
        datetick('x','yymmm','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(msaltplot,'LineWidth',2);
        set(msaltplot2,'LineWidth',2);
        set(gca,'FontSize',20);
        set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
        txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
        txt2=text(xData(22), min(meanplotlev)+1 ,['RMSD = ', num2str(round(mean_salt_RMS,2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
        txt3=text(xData(53), min(meanplotlev)+1 ,['bias = ', num2str(round(mean_salt_bias,2))], 'FontSize', m_quiver_ref_text_fontsize); 
        lgd=legend(['Model(',testname,')'],'NIFS');
        
        set(txt1,'FontSize',20);
        set(txt2,'FontSize',20);
        set(txt3,'FontSize',20);
        set(lgd,'FontSize',20);
        set(lgd,'Position',[0.7, 0.2, 0.2, 0.03]);
        set(lgd,'Orientation','vertical');
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [1000, 400]);
        set(gcf,'PaperPosition', [0 0 1000 400]) 

        
        grid on
        grid minor
        hold off
        saveas(gcf,jpgname,'jpg');
        grid off
    end
end

close all;


for numline=1:size(NIFS6.line,2)
    if (isnan(sum(NIFS6.salt75_line_comb(:,numline),'omitnan'))==1)

    else
        lineline=unique(NIFS6.line(:,numline,:,:,:),'stable');
        lineline2=lineline(~isnan(lineline));
        lineline3=lineline2(lineline2>0);

        idx=isnan(NIFS6.salt75_line_comb(:,numline)');
        p=polyfit(xData(~idx),NIFS6.salt75_line_comb(~idx,numline)',1);
        line_trend2=xData*p(1)+p(2);
        msaltplot=plot(xData,MODEL6.salt75_line_comb(:,numline)','color','k');
        hold on
        msaltplot2=plot(xData,NIFS6.salt75_line_comb(:,numline)','o','color','r','markerfacecolor','k','markersize',2);

%         msaltplot2=plot(xData,line_trend2,'Color','r')
        idx2=isnan(MODEL6.salt75_line_comb(:,numline)');
        idx3=idx+idx2;
        idx3(idx3==2)=1;
        constant_cor=corrcoef(NIFS6.salt75_line_comb(~idx3,numline)',MODEL6.salt75_line_comb(~idx3,numline)');
        
        salt_bias=MODEL6.salt75_line_comb(~idx3,numline)' - NIFS6.salt75_line_comb(~idx3,numline)';
        mean_salt_bias=mean(salt_bias);
        mean_salt_RMS=mean(abs(salt_bias));
        hold off
        jpgname=strcat(outfile,testname,'_',num2str(lineline3), '0_75m_line_mean_salt_',num2str(min(inputyear),'%04i'),'_',num2str(max(inputyear),'%04i'), '.jpg'); %% ~_year_month.jpg
        xlabel('year')
        ylabel('Mean Salinity along NIFS line')
        title([num2str(lineline3), ', 0-75m mean salt(',num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i'),'), ',num2str(round(NIFS6.salt75_line_comb_trend(numline),4)), ' /y'])
        datetick('x','yymmm','keepticks')
        axis tight;
        ylim(meanplotlev)
        set(msaltplot,'LineWidth',2);
        set(msaltplot2,'LineWidth',2);
        set(gca,'FontSize',20);
        set(gca,'YTick',min(meanplotlev):1:max(meanplotlev))
        txt1=text(xData(3), min(meanplotlev)+1 ,['R = ', num2str(round(constant_cor(1,2),2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
        txt2=text(xData(22), min(meanplotlev)+1 ,['RMSD = ', num2str(round(mean_salt_RMS,2)), ', '], 'FontSize', m_quiver_ref_text_fontsize); 
        txt3=text(xData(53), min(meanplotlev)+1 ,['bias = ', num2str(round(mean_salt_bias,2))], 'FontSize', m_quiver_ref_text_fontsize); 
        lgd=legend(['Model(',testname,')'],'NIFS');
        
        set(txt1,'FontSize',20);
        set(txt2,'FontSize',20);
        set(txt3,'FontSize',20);
        set(lgd,'FontSize',20);
        set(lgd,'Position',[0.7, 0.2, 0.2, 0.03]);
        set(lgd,'Orientation','vertical');
        
        set(gcf, 'PaperUnits', 'points');
        set(gcf, 'PaperSize', [1000, 400]);
        set(gcf,'PaperPosition', [0 0 1000 400])
        
        grid on
        grid minor
        hold off
        saveas(gcf,jpgname,'jpg');
        grid off
    end
end

close all;