function status=scatter_isodepth_v_ROMS(testname, outfile, filedir, section, section2, year, inputmonth, refvalue, param_script, option)
% clear all;close all;
%==========================================================================
vert_param;
var='vert_temp';
run(param_script);
%==========================================================================


sectionname = [num2str(section(1)),'_', num2str(section(2)), '_',num2str(section(3)), '_', num2str(section(4))];
fname = [filedir, 'isodepth\', 'ROMS_', testname, '_isodepth_', num2str(refvalue), '_section_', sectionname, '_', num2str(year(1)), '_', num2str(year(end)), '.mat'];
load(fname);
fname = [filedir, 'isodepth\', 'ROMS_', testname, '_', 'vert_v', '_section_', sectionname, '_', num2str(year(1)), '_', num2str(year(end)), '.mat'];
load(fname);

abs_section_left= abs(x(1,:)-section2(1));
abs_section_right= abs(x(1,:)-section2(2));
min_section_left = min(abs_section_left);
min_section_right = min(abs_section_right);

section_left = find(abs_section_left == min_section_left);
section_right = find( abs_section_right== min_section_right);

if (section_left>1)
    if (section2(1)-x(1,section_left) < 0)
        section_left=section_left-1;
    end
end
if (section_right<size(x,2))
    if (section2(2)-x(1,section_right) > 0)
        section_right = section_right+1;
    end
end

% % 
% % abc=(1:144)';
% % abcd=repmat(abc,1,25);
% % abcde=reshape(abcd, [12,size(abcd,1)/12, size(abcd,2)]);
% % abcdef=abcde(:,2,:);

maxmonth = 120;

isodepth = isodepth(1:maxmonth,:);
sectiondata = sectiondata(1:maxmonth,:,:);

seasonal_isodepth = reshape(isodepth,[12,size(isodepth,1)/12, size(isodepth,2)]);  %% [clim_month year x] 
mean_seasonal_isodepth =  mean(seasonal_isodepth(:,:,section_left:section_right),3);
% mean_isodepth_feb = mean(seasonal_isodepth(2,:,section_left:section_right),3);
% mean_isodepth_aug = mean(seasonal_isodepth(8,:,section_left:section_right),3);

seasonal_surf_vel = reshape(sectiondata(:,40,:),[12,size(sectiondata,1)/12, size(sectiondata,3)]);
mean_seasonal_surf_vel = mean(seasonal_surf_vel(:,:,section_left:section_right),3);
% mean_surf_vel_feb = mean(seasonal_surf_vel(2,:,section_left:section_right),3);
% mean_surf_vel_aug = mean(seasonal_surf_vel(8,:,section_left:section_right),3);


mean_isodepth = mean(isodepth(:,section_left:section_right),2);
mean_surf_vel = mean(sectiondata(:,size(sectiondata,2),section_left:section_right),3);
max_surf_vel = max(sectiondata(:,size(sectiondata,2),section_left:section_right), [], 3);
% mean_slope = (isodepth(:,section_right)-isodepth(:,section_left)) /  (m_lldist([section2(1) section2(2)], [section(3) section(4)]) *1000.0);
for i= 1:maxmonth
    mean_slope(i,:)=  polyfit(x(1,section_left:section_right),isodepth(i,section_left:section_right),1);
end
% mean_slope(:,:)=  polyfit(repmat(x(1,section_left:section_right),120,1),isodepth(:,section_left:section_right),1);
mean_seasonal_slope = reshape(mean_slope(:,1),[12,size(mean_slope,1)/12]);

marker_size=20;

msi_fa(1:10) = mean_seasonal_isodepth(2,:);
msi_fa(11:20) = mean_seasonal_isodepth(8,:);
mssv_fa(1:10) = mean_seasonal_surf_vel(2,:);
mssv_fa(11:20) = mean_seasonal_surf_vel(8,:);
mss_fa(1:10) = mean_seasonal_slope(2,:);
mss_fa(11:20) = mean_seasonal_slope(8,:);
% tempmsj = jet(2);
msj(1:10,:) = repmat(bwrmap(1,:),10,1);
msj(11:20,:) = repmat(bwrmap(100,:),10,1);


bwr_all(1:12:109,:) = repmat(bwrmap(1,:),10,1);
bwr_all(2:12:110,:) = repmat(bwrmap(1,:),10,1);
bwr_all(3:12:111,:) = repmat(bwrmap(15,:),10,1);
bwr_all(4:12:112,:) = repmat(bwrmap(30,:),10,1);
bwr_all(5:12:113,:) = repmat(bwrmap(70,:),10,1);
bwr_all(6:12:114,:) = repmat(bwrmap(85,:),10,1);
bwr_all(7:12:115,:) = repmat(bwrmap(100,:),10,1);
bwr_all(8:12:116,:) = repmat(bwrmap(100,:),10,1);
bwr_all(9:12:117,:) = repmat(bwrmap(85,:),10,1);
bwr_all(10:12:118,:) = repmat(bwrmap(70,:),10,1);
bwr_all(11:12:119,:) = repmat(bwrmap(15,:),10,1);
bwr_all(12:12:120,:) = repmat(bwrmap(30,:),10,1);



switch(option)
    case(1)
%         sc=scatter(mean_isodepth, mean_surf_vel,marker_size,jet(size(mean_isodepth,1)), 'filled');
        sc=scatter(mean_isodepth, mean_surf_vel,marker_size,bwr_all, 'filled');
        lsl=lsline;
        r = corrcoef(squeeze(mean_isodepth), squeeze(mean_surf_vel), 'Alpha',0.05)
    case(2)
        scatter(mean_slope(:,1), mean_surf_vel,marker_size,jet(size(mean_isodepth,1)), 'filled');
        r = corrcoef(squeeze(mean_slope(:,1)), squeeze(mean_surf_vel))
        lsl=lsline;
    case(3)
        sc=scatter(msi_fa, mssv_fa,marker_size,msj, 'filled');
        r = corrcoef(msi_fa, mssv_fa)
        lsl=lsline;
    case(4)
        sc=scatter(mss_fa, mssv_fa,marker_size,msj, 'filled');
        r = corrcoef(mss_fa, mssv_fa)
        lsl=lsline;
    case(5)
        monthind=str2num(outfile(end));
        sc=scatter(mean_seasonal_isodepth(monthind,:), mean_seasonal_surf_vel(monthind,:), ...
            marker_size,jet(length(mean_seasonal_isodepth(monthind,:))), 'filled');
        r = corrcoef(mean_seasonal_isodepth(monthind,:), mean_seasonal_surf_vel(monthind,:));
        lsl=lsline;
    case(6)
        monthind=str2num(outfile(end));
        sc=scatter(mean_seasonal_slope(monthind,:), mean_seasonal_surf_vel(monthind,:), ...
            marker_size,jet(length(mean_seasonal_isodepth(monthind,:))), 'filled');
        r = corrcoef(mean_seasonal_slope(monthind,:), mean_seasonal_surf_vel(monthind,:));
        lsl=lsline;
end



% sc=scatter(mean_seasonal_isodepth(2,:), mean_seasonal_surf_vel(2,:),marker_size,jet(size(mean_seasonal_isodepth,2)), 'filled');
% r3 = corrcoef(mean_seasonal_isodepth(2,:), mean_seasonal_surf_vel(2,:))
% lsl=lsline;

% scatter(mean_seasonal_isodepth(8,:), mean_seasonal_surf_vel(8,:),marker_size,jet(size(mean_seasonal_isodepth,2)), 'filled');
% r4 = corrcoef(mean_seasonal_isodepth(8,:), mean_seasonal_surf_vel(8,:))





set(gca,'box', vert_grid_box, 'linewidth',vert_grid_box_linewidth, 'tickdir',vert_grid_tickdir_type,'fontsize',vert_grid_fontsize)

if (option==2 || option ==4 || option ==6)
    xlabel('slope','FontSize',vert_grid_fontsize,'fontweight','bold')
else
    xlabel('10^oC depth (m)','FontSize',vert_grid_fontsize,'fontweight','bold')
end
ylabel('Surface velocity (m/s)','FontSize',vert_grid_fontsize,'fontweight','bold')
h = colorbar;
if (option == 1 || option==3 || option ==4)
    colormap(bwrmap);
%     caxis(['feb', 'aug']);
    h.Ticks=[0, 0.15 0.3 0.7 0.85 1];
    h.TickLabels=['Jan, Feb'; 'Dec, Mar'; 'Nov, Apr'; 'May, Oct'; 'Jun, Sep'; 'Jul, Aug'];
else
    colormap jet;
    caxis([year(1) 2000+maxmonth/12]);
end
set(h,'fontsize',colorbar_fontsize-3);
title(h,'Month','fontsize',colorbar_title_fontsize);

gcaaa=gca;
textypos =  gcaaa.YLim(2)-(gcaaa.YLim(2)-gcaaa.YLim(1))/20;
textxpos =  gcaaa.XLim(2)-(gcaaa.XLim(2)-gcaaa.XLim(1))/4.5;
text1 = text(textxpos, textypos, ['r = ', num2str(r(1,2),'%3.2f')], 'FontSize', colorbar_title_fontsize); 

jpgname=strcat(outfile, '_', testname, '_', num2str(year(1),'%04i'), '_', num2str(2000+maxmonth/12,'%04i'), '.jpg'); %% ~_year_month.jpg
saveas(gcf,jpgname,'jpg');
status=1;
return