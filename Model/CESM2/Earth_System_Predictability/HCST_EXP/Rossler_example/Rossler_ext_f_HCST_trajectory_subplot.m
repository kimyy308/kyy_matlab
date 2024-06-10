close all; clear all; clc;

dir.lroot = '/Volumes/kyy_raid/kimyy/Model/Rossler_predictability_example';

dir.hcst = [dir.lroot, '/', 'ensPRED'];
dir.ctl = [dir.lroot, '/', 'CTL'];

cfg.ensnum=50;
cfg.t_ini_min=5000;
% cfg.t_ini_max=95;
cfg.t_ini_max=10500;
cfg.pred_len=1000;
cfg.t_ini_min_CTL=5000;
cfg.t_ini_max_CTL=19900;
tstep=10;
ext_flag=0;

fname=['/Volumes/kyy_raid/kimyy/Model/Rossler_predictability_example', '/', ...
    'rossler_temp.mat'];

if exist(fname)==0

    %% read control run
    for t_ini=cfg.t_ini_min_CTL:tstep:cfg.t_ini_max_CTL
        str.t_ini=num2str(t_ini, '%08i');
        dir.ctl_tmp=[dir.ctl, '/', 'out'];
        fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
    %     CTL_X.(['p',str.t_ini])=single(ncread(fname_ctl, 'X'));
    %     CTL_Y.(['p',str.t_ini])=single(ncread(fname_ctl, 'Y'));
    %     CTL_Z.(['p',str.t_ini])=single(ncread(fname_ctl, 'Z'));
        ncid=netcdf.open(fname_ctl, 'NOWRITE');
        xvarid=netcdf.inqVarID(ncid,'X');
        yvarid=netcdf.inqVarID(ncid,'Y');
        zvarid=netcdf.inqVarID(ncid,'Z');
        CTL_X.(['p',str.t_ini])=netcdf.getVar(ncid,xvarid);
        CTL_Y.(['p',str.t_ini])=netcdf.getVar(ncid,yvarid);
        CTL_Z.(['p',str.t_ini])=netcdf.getVar(ncid,zvarid);
        disp(num2str(t_ini))
        netcdf.close(ncid)
    end
    
    %% data combine for CTL
    for t_ini=cfg.t_ini_min_CTL:tstep:cfg.t_ini_max_CTL
        str.t_ini=num2str(t_ini, '%08i');
    %     t_plot=t_ini:0.01:t_ini+0.99;
    %     plot(t_plot, CTL_Z.(['p',str.t_ini]), 'b');
    
    %     data_comb.t((t_ini-cfg.t_ini_min_CTL)+1:(t_ini-cfg.t_ini_min_CTL)+100)=t_plot;
        data_comb.X((t_ini-cfg.t_ini_min_CTL)+1:(t_ini-cfg.t_ini_min_CTL)+10)=CTL_X.(['p',str.t_ini]);
        data_comb.Y((t_ini-cfg.t_ini_min_CTL)+1:(t_ini-cfg.t_ini_min_CTL)+10)=CTL_Y.(['p',str.t_ini]);
        data_comb.Z((t_ini-cfg.t_ini_min_CTL)+1:(t_ini-cfg.t_ini_min_CTL)+10)=CTL_Z.(['p',str.t_ini]);
    end


    %% read ensemble predictions
    for t_ini=cfg.t_ini_min:tstep:cfg.t_ini_max
        str.t_ini=num2str(t_ini, '%08i');
        for ensnum=1:cfg.ensnum
            str_ens=num2str(ensnum, '%02i');
    %         disp('pred:', str.t_ini, 'ens:',str_ens);
            dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
    
            for t_pred=t_ini:tstep:t_ini+cfg.pred_len
                str.t_hcst=num2str(t_pred, '%08i');
                fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_hcst, '.nc'];
                %% raw
                ncid=netcdf.open(fname_hcst, 'NOWRITE');
                xvarid=netcdf.inqVarID(ncid,'X');
    %             yvarid=netcdf.inqVarID(ncid,'Y');
    %             zvarid=netcdf.inqVarID(ncid,'Z');
                HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(netcdf.getVar(ncid,xvarid));
                netcdf.close(ncid)
            end
            disp([str.t_ini, ', ens : ', str_ens]);
        end
    end

    %% additional hindcast read for trajectory plot
    t_ini=cfg.t_ini_min;
    str.t_ini=num2str(cfg.t_ini_min, '%08i');
    for ensnum=1:cfg.ensnum
        str_ens=num2str(ensnum, '%02i');
    %         disp('pred:', str.t_ini, 'ens:',str_ens);
        dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
    
        for t_pred=cfg.t_ini_min:tstep:cfg.t_ini_max_CTL
            str.t_hcst=num2str(t_pred, '%08i');
            fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_hcst, '.nc'];
    
            ncid=netcdf.open(fname_hcst, 'NOWRITE');
            xvarid=netcdf.inqVarID(ncid,'X');
            yvarid=netcdf.inqVarID(ncid,'Y');
            zvarid=netcdf.inqVarID(ncid,'Z');
            HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(netcdf.getVar(ncid,xvarid));
            HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(netcdf.getVar(ncid,yvarid));
            HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(netcdf.getVar(ncid,zvarid));
    %         disp(num2str(t_ini))
            netcdf.close(ncid)
         
        end
        disp([str.t_ini, ', ens : ', str_ens]);
    end

    save(fname, '-v7.3');
else
    load(fname);
end


% p to data
for t_ini=cfg.t_ini_min
    str.t_ini=num2str(t_ini, '%08i');
    for ensi=1:cfg.ensnum
        str.ens=num2str(ensi, '%02i');
        for ly=tstep:tstep:14500
            str.t_pred=num2str(t_ini+ly,'%08i');
            hcst_comb.X.(['i',str.t_ini])(ensi,(ly-tstep)+1:ly) = ...
                HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);

            hcst_comb.Z.(['i',str.t_ini])(ensi,(ly-tstep)+1:ly) = ...
                HCST_Z.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);
        end
    end
end

title('x')
hold on
for ensi=1:cfg.ensnum
    plot(0.01:0.01:145,hcst_comb.X.i00005000(ensi,:), 'b')
end
plot(0.01:0.01:145,mean(hcst_comb.X.i00005000,1),'k', 'linewidth',3)
hold off
xlabel('tau')
ylabel('x')
set(gca,'fontsize',20)

title('z')
hold on
for ensi=1:cfg.ensnum
    plot(0.01:0.01:145,hcst_comb.Z.i00005000(ensi,:), 'b')
end
plot(0.01:0.01:145,mean(hcst_comb.Z.i00005000,1),'k', 'linewidth',3)
hold off
xlabel('tau')
ylabel('z')
set(gca,'fontsize',20)

corrcoef(data_comb.X(1:14500), mean(hcst_comb.X.i00005000(:,1:14500),1))

%% assign data by lead time
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    for lt_i=1:cfg.pred_len/tstep
        lt=(lt_i-1)*tstep;
%         if lt==0
%             lt=1;
%         end
        for t_ini=cfg.t_ini_min:tstep:cfg.t_ini_max
            str.t_ini=num2str(t_ini, '%08i');
            t_p=t_ini+lt;
            str.t_p=num2str(t_p, '%08i');
            HCST_X.data_lt(ensi,lt_i,(t_ini-cfg.t_ini_min+tstep)/tstep) = ...
                HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])(1);
        end
    end
end
HCST_X.data_lt_em=squeeze(mean(HCST_X.data_lt,1));




%% corr coeficient calculation -------------------------------------------------------

% corrcoef(data_comb.X(901:100:3401),  HCST_X.data_lt(1,1,6:31))  % (lt=1, 60~85)
% corrcoef(data_comb.X(900:100:3400),  HCST_X.data_lt_em(500,1:26))  % (lt=500, 60~85)

hold off
plot(data_comb.Z(1000:10:6500))
hold on
plot(HCST_X.data_lt_em(100,:))

hold on
for i=1:50
    plot(HCST_X.(['ens',num2str(i,'%02i')]).i00005000.p00010500)
end


% for HCST
tis=[1, 100:100:cfg.pred_len*100];

for tti=1:length(tis)
    ti=tis(tti);
    tmp.tmin=(cfg.t_ini_min-cfg.t_ini_min_CTL+cfg.pred_len)*100+mod(ti,100);
    tmp.tmax=(cfg.t_ini_max-cfg.t_ini_min_CTL)*100+mod(ti,100);
    tmp.ts=1+cfg.pred_len-floor(ti/100); % 1 + 10 ~ 1 + 11 - 10
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1-floor(ti/100); % 95 - 55 + 1 ~ 95 -55 + 1 - 10

    [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        HCST_X.data_lt_em(tti,tmp.ts:tmp.te));  % (lt=1000, 65~95)
    HCST_X.corr_em(tti)=tmp.corr(1,2);
    HCST_X.corr_em_p(tti)=tmp.corr_p(1,2);
    HCST_X.sig_t_em(tti)=std(HCST_X.data_lt_em(tti,tmp.ts:tmp.te));
    tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        squeeze(HCST_X.data_lt_em(tti,tmp.ts:tmp.te)));  
    HCST_X.cov_em(tti)=tmp.cov(1,2);
    for ensi=1:cfg.ensnum
        [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            HCST_X.data_lt(ensi,tti,tmp.ts:tmp.te));  
        HCST_X.corr(ensi,tti)=tmp.corr(1,2);
        HCST_X.corr_p(ensi,tti)=tmp.corr_p(1,2);
        HCST_X.sig_t_indv(ensi,tti)=std(HCST_X.data_lt(ensi,tti,tmp.ts:tmp.te));
        tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            squeeze(HCST_X.data_lt(ensi,tti,tmp.ts:tmp.te)));  
        HCST_X.cov_indv(ensi,tti)=tmp.cov(1,2);
    end
    HCST_X.corr_med(tti)=median(HCST_X.corr(:,tti));
    HCST_X.corr_mean(tti)=mean(HCST_X.corr(:,tti));
    HCST_X.sig_t_med(tti)=squeeze(median(HCST_X.sig_t_indv(:,tti),1));
    HCST_X.sig_t_mean(tti)=squeeze(mean(HCST_X.sig_t_indv(:,tti),1));
    HCST_X.cov_med(tti)=median(HCST_X.cov_indv(:,tti));
    HCST_X.cov_mean(tti)=mean(HCST_X.cov_indv(:,tti),1);
end





%% data combine for HCST (ensmean) (tau=1:50)
data_comb_hcst.X(1:5000)=0;
data_comb_hcst.Y(1:5000)=0;
data_comb_hcst.Z(1:5000)=0;
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    t_ini=cfg.t_ini_min;
    str.t_ini=num2str(t_ini, '%04i');
    for lt=1:50
        t_p=t_ini+floor((lt-1));
        str.t_p=num2str(t_p, '%04i');
        data_comb_hcst.X((lt-1)*100+1:lt*100)=data_comb_hcst.X((lt-1)*100+1:lt*100) + ...
            HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
        data_comb_hcst.Y((lt-1)*100+1:lt*100)=data_comb_hcst.Y((lt-1)*100+1:lt*100) + ...
            HCST_Y.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
        data_comb_hcst.Z((lt-1)*100+1:lt*100)=data_comb_hcst.Z((lt-1)*100+1:lt*100) + ...
            HCST_Z.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
    end
end
data_comb_hcst.X=data_comb_hcst.X/cfg.ensnum;
data_comb_hcst.Y=data_comb_hcst.Y/cfg.ensnum;
data_comb_hcst.Z=data_comb_hcst.Z/cfg.ensnum;


%% data combine for HCST(ind) (tau=1:50)

for ensi=1:1
    str.ens=num2str(ensi, '%02i');
    t_ini=cfg.t_ini_min;
    str.t_ini=num2str(t_ini, '%04i');
    for lt=1:50
        t_p=t_ini+floor((lt-1));
        str.t_p=num2str(t_p, '%04i');
        data_comb_hcst.X_ind((lt-1)*100+1:lt*100)= HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
        data_comb_hcst.Y_ind((lt-1)*100+1:lt*100)= HCST_Y.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
        data_comb_hcst.Z_ind((lt-1)*100+1:lt*100)= HCST_Z.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])';
    end
end


fig_cfg.fig_size = [0,0,13,14.5]; %% paper size (original)
fig_cfg.cb_size = [1, 13, 11, 0.3];
loc_column_first=1;
loc_row_first=7;
fig_h = figure('name', 'fig_lorenz','PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

%% plot using surf (CTL, OBS)
ax_m=subplot(2,2,1);
fig_cfg.ax_size = [loc_column_first, loc_row_first, 5, 5];
period=1:14000;
fig_weig=1;
% ax_m=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m, 'Parent', fig_h)
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m,'fontsize',12);
grid on
hold on
plots.X=data_comb.X(period+400);
plots.Y=data_comb.Y(period+400);
plots.Z=data_comb.Z(period+400);
tmp.parula=parula((100*2)+(max(period))); % from light blue, until darker yellow
if isfield(plots, 'C0')==1
plots=rmfield(plots,'C0');
end
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
title1=title('(a) Observation', 'fontsize', 20);

sc1=scatter3(ax_m,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);

view(20,30)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-25 25])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-25 25]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 50]);
cb = colorbar(ax_m,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [810, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


%% plot using surf (ensemble mean)
% ax_m=subplot(2,2,2);
fig_cfg.ax_size = [loc_column_first + 6, loc_row_first, 5, 5];
period=1:50000;
fig_weig=1;
% ax_m=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m, 'Parent', fig_h)
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m,'fontsize',12);
grid on
hold on
plots.X=data_comb_hcst.X(period);
plots.Y=data_comb_hcst.Y(period);
plots.Z=data_comb_hcst.Z(period);
tmp.parula=parula((cfg.pred_len+2)*(max(period)/10)); % from light blue, until darker yellow
plots=rmfield(plots,'C0');
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line

title1=title('(b) Ensemble Mean', 'fontsize', 20);
sc1=scatter3(ax_m,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);
view(45,30)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-25 25])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-25 25]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 50]);
cb = colorbar(ax_m,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [810, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


%% plot using surf (Individual)
ax_m=subplot(2,2,3);
fig_cfg.ax_size = [loc_column_first, loc_row_first-6, 5, 5];
period=1:5000;
fig_weig=1;
% ax_m=axes('fontsize',14, 'fontname','freeserif'); 
set(ax_m, 'Parent', fig_h)
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);
set(ax_m,'fontsize',12);
grid on
hold on
plots.X=data_comb_hcst.X_ind(period);
plots.Y=data_comb_hcst.Y_ind(period);
plots.Z=data_comb_hcst.Z_ind(period);
tmp.parula=parula((cfg.pred_len+2)*(max(period)/10)); % from light blue, until darker yellow
plots=rmfield(plots,'C0');
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf(ax_m, [plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
title1=title('(c) Individual', 'fontsize', 20);
sc1=scatter3(ax_m,plots.X(1), plots.Y(1), plots.Z(1), 'MarkerFaceColor', plots.C0(1,1,:), ...
    'MarkerEdgeColor', 'none', 'LineWidth', 3);
view(45,30)
caxis([0 max(period)/100]);
xlabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlim([-20 20])
ylabel('$$ y(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
ylim([-25 25]);
zlabel('$$ z(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
zlim([0 50]);
cb = colorbar(ax_m,'units','inches','Location', 'northoutside', 'position',fig_cfg.cb_size);
set(cb,'fontsize',15,'fontname','freeserif','TickDir','both');
cb_title=title(cb,'$$ \tau $$','fontsize', 22, 'Position', [810, 0, 0]); % hor, ver, ?
set(cb_title, 'interpreter', 'latex');


% q3=norminv(0.75);
% q95=norminv(0.975);
% w95=(q95-q3)/(2*q3);

ax_m=subplot(2,2,4);
fig_cfg.ax_size = [loc_column_first+6, loc_row_first-6, 5, 5];
set(ax_m, 'Parent', fig_h)
set(ax_m,'Units','inches','Position',fig_cfg.ax_size);

% boxplot(HCST_X.corr(:,1), 1, 'BoxStyle', 'outline', 'Colors', 'krrrbbbbbbbbbbbbbbbbbbbbb', 'ColorGroup', RCM_fig.tot_trend_rand_name, ...
%     'Whisker', w95, 'symbol','');

sc_1=scatter(ax_m, 1:11, HCST_X.corr_em);
set(sc_1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k')
set(sc_1, 'SizeData', 50.*fig_weig)

hold on
bp1=boxplot(ax_m, HCST_X.corr(:,1:11), 0:10, 'BoxStyle', 'outline', 'symbol', '', 'Whisker', 2);
% set(bp1, 'Widths', 0.4)
set(bp1, {'linew'}, {1.5})

% bp1=boxchart(double(HCST_X.corr(:,1:11)),  'symbol', '');

ylim([-0.2 1])
title1=title('(d) ACC(\tau)', 'fontsize', 20);
grid on

ylabel('ACC for $$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlabel('lead $$ \tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)

set(ax_m,'fontsize',15);
set(ax_m,'YAxisLocation', 'right');

% box_vars = findall(fig_h,'Tag','Box');
% hLegend = legend(box_vars([3,2,4]), {'Group 1','Group 2','Group 3'});



%% 95% significance low & high range
conf = 0.95;
n=926;
alpha = 1 - conf;
pLo = alpha/2;
pUp = 1 - alpha/2;
crit = tinv([pLo pUp], n-1);
% xbar = mean(r); % = 0
xbar = 0; % = 0
r_crit=sqrt((crit.^2)./(n-2+(crit).^2));

line(ax_m, 0:12, repmat(r_crit(1), 1,13), 'color', 'r', 'LineStyle', '--');

print(gcf, ['lorenz_trajectory_',num2str(ext_flag), '.png'], '-dpng');















