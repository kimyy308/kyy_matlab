close all; clear all; clc;

dir.lroot = ['~/Desktop/backup/Research/Postdoc/03_IBS/', ...
    '2022_predictability_assimilation_run/Lorenz_predictability_example/Model_src-v1'];

dir.hcst = [dir.lroot, '/', 'ensPRED'];
dir.hcst_lqini = [dir.lroot, '/', 'ensPRED_lqini']; % low quality initialized
dir.hcst_unini = [dir.lroot, '/', 'ensPRED_unini']; % low quality initialized
dir.ctl = [dir.lroot, '/', 'ext_f_CTL'];
dir.ctl = [dir.lroot, '/', 'CTL'];

cfg.ensnum=100;
cfg.t_ini_min=55;
% cfg.t_ini_max=95;
cfg.t_ini_max=990;
cfg.pred_len=10;
cfg.t_ini_min_CTL=51;
% cfg.t_ini_max_CTL=100;
cfg.t_ini_max_CTL=990;

ext_factor = 2;
ext_flag=0;

%% read control run
for t_ini=51:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'out'];
    fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
    if ext_flag==0
%% raw
    CTL_X.(['p',str.t_ini])=ncread(fname_ctl, 'X');
    CTL_Y.(['p',str.t_ini])=ncread(fname_ctl, 'Y');
    CTL_Z.(['p',str.t_ini])=ncread(fname_ctl, 'Z');
    else
%% with external forcing
    CTL_X.(['p',str.t_ini])=ncread(fname_ctl, 'X') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    CTL_Y.(['p',str.t_ini])=ncread(fname_ctl, 'Y') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    CTL_Z.(['p',str.t_ini])=ncread(fname_ctl, 'Z') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    end
end

for t_ini=51:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'out'];
    fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
    CTL_X_raw.(['p',str.t_ini])=ncread(fname_ctl, 'X');
    CTL_X_ext.(['p',str.t_ini])=0 + (t_ini-51)/ext_factor+[0.01:0.01:1]';
end
for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
    data_comb_raw.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X_raw.(['p',str.t_ini]);
    data_comb_ext.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X_ext.(['p',str.t_ini]);
end
% std(data_comb_ext.X)
% std(data_comb_raw.X)
% std(data_comb_raw.X)/std(data_comb_ext.X)

%% read ensemble predictions
% for t_ini=55:85
for t_ini=cfg.t_ini_min:cfg.t_ini_max
    str.t_ini=num2str(t_ini, '%04i');
    for ensnum=1:cfg.ensnum
        str_ens=num2str(ensnum, '%02i');
%         disp('pred:', str.t_ini, 'ens:',str_ens);
        dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
        dir.hcst_lqini_tmp=[dir.hcst_lqini, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
        dir.hcst_unini_tmp=[dir.hcst_unini, '/', '0055', '/', 'ens',str_ens, '/', 'out'];

        for t_pred=1:cfg.pred_len
            str.t_hcst=num2str(t_ini+t_pred-1, '%04i');
            fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_hcst, '.nc'];
            fname_hcst_lqini=[dir.hcst_lqini_tmp, '/', 'out.', str.t_hcst, '.nc'];
            fname_hcst_unini=[dir.hcst_unini_tmp, '/', 'out.', str.t_hcst, '.nc'];
            if ext_flag==0
            %% raw
                HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'X');
                HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Y');
                HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Z');

            %% raw_unini
            HCST_unini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_unini, 'X');
% % %             HCST_unini_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_unini, 'Y');
% % %             HCST_unini_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_unini, 'Z'); 
            else
            %% with ext (ini)
            HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'X') + (t_ini-51+t_pred-1)/ext_factor+[0.01:0.01:1]';
            HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Y') + (t_ini-51+t_pred-1)/ext_factor+[0.01:0.01:1]'; 
            HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst, 'Z') + (t_ini-51+t_pred-1)/ext_factor+[0.01:0.01:1]';
% % %             HCST_lqini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'X');
% % %             HCST_lqini_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'Y');
% % %             HCST_lqini_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_lqini, 'Z');

            %% with ext (unini)
            HCST_unini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=ncread(fname_hcst_unini, 'X') + (t_ini-51+t_pred-1)/ext_factor+[0.01:0.01:1]';

            end
            
                    end
        disp([str.t_ini, ', ens : ', str_ens]);
    end
end

hold on


for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    t_plot=t_ini:0.01:t_ini+0.99;
    plot(t_plot, CTL_Z.(['p',str.t_ini]), 'b');

    data_comb.t((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=t_plot;
    data_comb.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X.(['p',str.t_ini]);
    data_comb.Y((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_Y.(['p',str.t_ini]);
    data_comb.Z((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_Z.(['p',str.t_ini]);
end
hold off

for t_ini=cfg.t_ini_min:cfg.pred_len:cfg.t_ini_max
    str.t_ini=num2str(t_ini, '%04i');
    for ensi=1:cfg.ensnum
        str.ens=num2str(ensi, '%02i');
        for ly=1:cfg.pred_len
            str.t_pred=num2str(t_ini+ly-1,'%04i');
            hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini])((ly-1)*100+1:(ly-1)*100+100) = ...
                HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);
            hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini])((ly-1)*100+1:(ly-1)*100+100) = ...
                HCST_unini_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_pred]);
        end
    end
end

%% assign data by lead time
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    for lt=1:cfg.pred_len*100
        for t_ini=cfg.t_ini_min:cfg.t_ini_max
            str.t_ini=num2str(t_ini, '%04i');
            t_p=t_ini+floor((lt-1)/100);
            str.t_p=num2str(t_p, '%04i');
            HCST_X.data_lt(ensi,lt,t_ini-54) = ...
                HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])(lt-(t_p-t_ini)*100);
            HCST_unini_X.data_lt(ensi,lt,t_ini-54) = ...
                HCST_unini_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])(lt-(t_p-t_ini)*100);
        end
    end
end
HCST_X.data_lt_em=squeeze(mean(HCST_X.data_lt,1));
HCST_unini_X.data_lt_em=squeeze(mean(HCST_unini_X.data_lt,1));


% % % % % %% corr, not same period
% % % % % for lt=1:500
% % % % %     tmp.corr=corrcoef(HCST_X.data_lt_em(1,:), HCST_X.data_lt_em(lt,:));
% % % % %     HCST_X.autocorr(lt)=tmp.corr(1,2);
% % % % %     HCST_X.sig_t(lt)=std(HCST_X.data_lt_em(lt,:));
% % % % %     
% % % % %     tmp.corr=corrcoef(HCST_unini_X.data_lt_em(1,:), HCST_unini_X.data_lt_em(lt,:));
% % % % %     HCST_unini_X.autocorr(lt)=tmp.corr(1,2);
% % % % %     HCST_unini_X.sig_t(lt)=std(HCST_unini_X.data_lt_em(lt,:));
% % % % % 
% % % % %     tmp.corr=corrcoef(data_comb.X(1:4500), data_comb.X(1+lt-1:4500+lt-1));
% % % % %     CTL_X.autocorr(lt)=tmp.corr(1,2);
% % % % %     CTL_X.sig_t(lt)=std(data_comb.X(1+lt-1:4500+lt-1));
% % % % %     for ensi=1:20
% % % % %         HCST_X.sig_t_indv(ensi,lt)=std(HCST_X.data_lt(ensi,lt,:));
% % % % %         HCST_unini_X.sig_t_indv(ensi,lt)=std(HCST_unini_X.data_lt(ensi,lt,:));
% % % % %     end
% % % % %     HCST_X.sig_t_med(lt)=median(HCST_X.sig_t_indv(:,lt),1);
% % % % %     HCST_unini_X.sig_t_med(lt)=median(HCST_unini_X.sig_t_indv(:,lt),1);
% % % % %     HCST_X.sig_t_mean(lt)=mean(HCST_X.sig_t_indv(:,lt),1);
% % % % %     HCST_unini_X.sig_t_mean(lt)=mean(HCST_unini_X.sig_t_indv(:,lt),1);
% % % % % end
% % % % % plot(HCST_X.autocorr)
% % % % % plot(HCST_X.sig_t)
% % % % % HCST_X.sig_t(1).*0.37;
% % % % % 
% % % % % plot(HCST_unini_X.autocorr);
% % % % % plot(HCST_unini_X.sig_t);
% % % % % % 37% = e-folding time scale
% % % % % plot(CTL_X.autocorr)
% % % % % plot(CTL_X.sig_t)


% % % % % plot(CTL_X.sig_t)
% % % % % hold on
% % % % % plot(HCST_X.sig_t)
% % % % % plot(HCST_unini_X.sig_t);
% % % % % plot(HCST_X.sig_t_med)
% % % % % plot(HCST_unini_X.sig_t_med);
% % % % % hold off
% % % % % legend
% % % % % set(gca, 'fontsize', 15.*fig_weig)


% LT = trenddecomp(data_comb.Z);
% plot(LT)










%% plot plot plot plot plot


close all;
clear plot_legend
subplot(2,1,1)
set(gcf,'PaperPosition',[0 0 30 45]);
fig_weig=2;

% (cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1;
% (cfg.t_ini_max-cfg.t_ini_min_CTL)*100;

%% normal plot
plot_legend(1)=plot(data_comb.t((cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1:(cfg.t_ini_max-cfg.t_ini_min_CTL)*100), ...
    data_comb.X((cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1:(cfg.t_ini_max-cfg.t_ini_min_CTL)*100), 'k', 'linewidth', 4.*fig_weig);
hold on

% %% using plot
% hcst_ensm.X.i0065=zeros(1,500);
% for ensi=1:20
%     str.ens=num2str(ensi, '%02i');
%     plot(65:0.01:69.99, hcst_comb.X.(['ens',str.ens]).i0065,'color', ...
%         [19 0 255 30]/255, 'linewidth', 2); %blue, transparency: 30/255;
%     hcst_ensm.X.i0065=hcst_ensm.X.i0065+hcst_comb.X.(['ens',str.ens]).i0065;
% end
% hcst_ensm.X.i0065=hcst_ensm.X.i0065./20;
% plot(65:0.01:69.99, hcst_ensm.X.i0065,'color', [19 0 255]/255, 'linewidth', 4);  % blue




%% using surf for gradient color plot (tstart +10 ~ tstart +10 + pred_len)
str.t_ini=num2str(cfg.t_ini_min+10, '%04i');
hcst_ensm.X.(['i',str.t_ini])=zeros(1,100*cfg.pred_len);
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    tmp.parula=parula((cfg.pred_len+2)*100); % from light blue, until darker yellow
    plots.X=cfg.t_ini_min+10:0.01:cfg.t_ini_min+10+cfg.pred_len-0.01;
    plots.Y=hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini]);
    plots.Z=1:length(plots.X);
    for ti=1:length(plots.X)
        plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
        plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
        plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
    end
%     plots.C0(1:length(plots.X),1:2,4)=30/255; RGBa cannot be applied
%     colormap([19 0 255]/255);
    objsurf=surf([plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 2.*fig_weig, 'EdgeAlpha', 0.2);           
     view(2);   % Default 2-D view)
%     objsurf.CData=[19 0 255]/255;
%     objsurf.AlphaData=1;
   
    hcst_ensm.X.(['i',str.t_ini])=hcst_ensm.X.(['i',str.t_ini])+hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini]);
end
hcst_ensm.X.(['i',str.t_ini])=hcst_ensm.X.(['i',str.t_ini])./cfg.ensnum;

% plot(65:0.01:69.99, hcst_ensm.X.i0065,'color', [19 0 255]/255, 'linewidth', 4);  % blue

%% ensemble mean plot using surf (tstart +10 ~ tstart +10 + pred_len)
plots.X=cfg.t_ini_min+10:0.01:cfg.t_ini_min+10+cfg.pred_len-0.01;
plots.Y=hcst_ensm.X.(['i',str.t_ini]);
plots.Z=1:length(plots.X);
for ti=1:length(plots.X)
    plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
    plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
    plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
end
 objsurf=surf([plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
        plots.C0, ...  % Reshape and replicate data
     'FaceColor', 'none', ...    % Don't bother filling faces with color
     'EdgeColor', 'interp', ...  % Use interpolated color for edges
     'LineWidth', 4.*fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
 view(2);   % Default 2-D view)

clear plots



%% unini plot (tstart + 10 + pred_len + 10 ~ tstart + 10 + pred_len + 10 + pred_len)
% plot(data_comb.t, data_comb.X, 'k', 'linewidth', 3)
% hold on
str.t_ini=num2str(cfg.t_ini_min+10+cfg.pred_len+10, '%04i');

hcst_ensm_unini.X.(['i',str.t_ini])=zeros(1,cfg.pred_len*100);
for ensi=1:cfg.ensnum
    str.ens=num2str(ensi, '%02i');
    plot(cfg.t_ini_min+10+cfg.pred_len+10:0.01:cfg.t_ini_min+10+cfg.pred_len+10+cfg.pred_len-0.01, ...
        hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini]),'color', [255 0 0 30]/255, 'linewidth', 2.*fig_weig); % red, 30/255 transparency
    hcst_ensm_unini.X.(['i',str.t_ini])=hcst_ensm_unini.X.(['i',str.t_ini])+hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini]);
end
hcst_ensm_unini.X.(['i',str.t_ini])=hcst_ensm_unini.X.(['i',str.t_ini])./cfg.ensnum;
plot_legend(4)=plot(cfg.t_ini_min+10+cfg.pred_len+10:0.01:cfg.t_ini_min+10+cfg.pred_len+10+cfg.pred_len-0.01, hcst_ensm_unini.X.(['i',str.t_ini]),'color',[255 0 0]/255, 'linewidth', 4.*fig_weig); % dark red


%% normal scatter plot, lt=1 (tstart + 1 : tstart + 10)
for ensi=1:cfg.ensnum
    s=scatter(cfg.t_ini_min+1:cfg.t_ini_min+10,squeeze(HCST_X.data_lt(ensi,1,2:11)),'filled');
    s.Marker='p'; %pentagram
%     s.MarkerEdgeColor=[143 189 136]/255; % green
    s.MarkerEdgeColor=[tmp.parula(1+100,:)]; 
    s.MarkerFaceColor = s.MarkerEdgeColor;
    s.SizeData=100.*fig_weig;
    s.MarkerFaceAlpha=1; %0.8
    s.MarkerEdgeAlpha=s.MarkerFaceAlpha;
    s.AlphaDataMode='manual';
end
plot_legend(2)=scatter(cfg.t_ini_min+1:cfg.t_ini_min+10,squeeze(HCST_X.data_lt_em(1,2:11)),'filled');
plot_legend(2).Marker='s'; %square
% s.MarkerEdgeColor=[0 128 0]/255; % dark green
plot_legend(2).MarkerEdgeColor=[tmp.parula(1+100,:)]; 
plot_legend(2).MarkerFaceColor = plot_legend(2).MarkerEdgeColor;
plot_legend(2).MarkerFaceAlpha=1;
plot_legend(2).MarkerEdgeAlpha=plot_legend(2).MarkerFaceAlpha;
plot_legend(2).SizeData=200.*fig_weig;
plot_legend(2).AlphaDataMode='manual';



%% normal scatter plot, lt=500 (tstart + 10 + pred_len ~ tstart + 10 + pred_len + 10 -1)
for ensi=1:cfg.ensnum
    s=scatter(cfg.t_ini_min+10+cfg.pred_len:cfg.t_ini_min+10+cfg.pred_len+10-1, ...
        squeeze(HCST_X.data_lt(ensi,end,1+10:1+10+cfg.pred_len-1)),'filled');
    s.Marker='p';
%     s.MarkerEdgeColor=[255 208 0]/255; % yellow
    s.MarkerEdgeColor=tmp.parula(cfg.pred_len*100+100,:); 
    s.MarkerFaceColor = s.MarkerEdgeColor;
    s.SizeData=100.*fig_weig;
    s.MarkerFaceAlpha=1;  %0.8
    s.MarkerEdgeAlpha=s.MarkerFaceAlpha;
    s.AlphaDataMode='manual';
end
plot_legend(3)=scatter(cfg.t_ini_min+10+cfg.pred_len:cfg.t_ini_min+10+cfg.pred_len+10-1, ...
    squeeze(HCST_X.data_lt_em(end,1+10:1+10+cfg.pred_len-1)),'filled');
plot_legend(3).Marker='s';
% s.MarkerEdgeColor=[255 140 0]/255; % dark orange
plot_legend(3).MarkerEdgeColor=[tmp.parula(cfg.pred_len*100+100,:)]; 
plot_legend(3).MarkerFaceColor = plot_legend(3).MarkerEdgeColor;
plot_legend(3).MarkerFaceAlpha=1;
plot_legend(3).MarkerEdgeAlpha=plot_legend(3).MarkerFaceAlpha;
plot_legend(3).SizeData=200.*fig_weig;
plot_legend(3).AlphaDataMode='manual';
plot_legend(3).AlphaDataMapping='none';

hold off
set(gca,'fontsize',15.*fig_weig);
grid on;
s=text(56,45,['Lorenz Model: ', 'X', ' versus ', '\tau'], 'fontsize', 22.*fig_weig, ...
    'FontName', 'Arial', 'FontWeight', 'Bold', 'FontAngle', 'italic');
% s=text(56,35,'Lorenz Model:', 'fontsize', 22, ...
%     'FontName', 'Sans Serif', 'FontWeight', 'Bold');
% s=text(59.9,34.9, '$$ x $$', 'Interpreter', 'latex', 'fontsize', 27, 'FontWeight', 'Bold' );
% s=text(60.4,35,'versus', 'fontsize', 22, ...
%     'FontName', 'Sans Serif');
% s=text(62.0,34.9, '$$ t $$', 'Interpreter', 'latex', 'fontsize', 27, 'FontWeight', 'Bold' );

if ext_flag==1
    ylabel('$$ x(\tau) + f(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
else
    ylabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
end
xlabel('$$ \tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)

pos=get(gca,'Position');
outerpos=get(gca,'OuterPosition');

% legend([plot_legend], {'CTL  ', 'INI[LT=0]  ', ['INI[LT=', num2str(cfg.pred_len), ']  '], 'unINI  '}, ...
%     'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 22.*fig_weig);
% legend([plot_legend], ...
%     {'$$ C  $$', '$$ I[lt=0] $$', ['$$ I[lt=', num2str(cfg.pred_len), ']  $$'], '$$ UI $$'}, ...
%     'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 22.*fig_weig, ...
%     'interpreter', 'latex');%% ver. 1
legend([plot_legend], ...
    {'$$ O  \hspace{5mm} $$ ', '$$ X_{I}[l\tau=0] \hspace{5mm} $$', ['$$ X_{I}[l\tau=', num2str(cfg.pred_len), '] \hspace{5mm} $$'], '$$ X_{UI} $$'}, ...
    'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 35.*fig_weig, ...
    'interpreter', 'latex'); %% fontsize 22*w = default

% print(gcf, 'abc.png', '-dpng');






%% corr coeficient calculation -------------------------------------------------------

% corrcoef(data_comb.X(901:100:3401),  HCST_X.data_lt(1,1,6:31))  % (lt=1, 60~85)
% corrcoef(data_comb.X(900:100:3400),  HCST_X.data_lt_em(500,1:26))  % (lt=500, 60~85)

% for HCST
for ti=1:cfg.pred_len*100
    tmp.tmin=(cfg.t_ini_min-cfg.t_ini_min_CTL+cfg.pred_len)*100+mod(ti,100);
    tmp.tmax=(cfg.t_ini_max-cfg.t_ini_min_CTL)*100+mod(ti,100);
    tmp.ts=1+cfg.pred_len-floor(ti/100); % 1 + 10 ~ 1 + 11 - 10
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1-floor(ti/100); % 95 - 55 + 1 ~ 95 -55 + 1 - 10

    [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        HCST_X.data_lt_em(ti,tmp.ts:tmp.te));  % (lt=1000, 65~95)
    HCST_X.corr_em(ti)=tmp.corr(1,2);
    HCST_X.corr_em_p(ti)=tmp.corr_p(1,2);
    HCST_X.sig_t_em(ti)=std(HCST_X.data_lt_em(ti,tmp.ts:tmp.te));
    tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        squeeze(HCST_X.data_lt_em(ti,tmp.ts:tmp.te)));  
    HCST_X.cov_em(ti)=tmp.cov(1,2);
    for ensi=1:cfg.ensnum
        [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            HCST_X.data_lt(ensi,ti,tmp.ts:tmp.te));  
        HCST_X.corr(ensi,ti)=tmp.corr(1,2);
        HCST_X.corr_p(ensi,ti)=tmp.corr_p(1,2);
        HCST_X.sig_t_indv(ensi,ti)=std(HCST_X.data_lt(ensi,ti,tmp.ts:tmp.te));
        tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            squeeze(HCST_X.data_lt(ensi,ti,tmp.ts:tmp.te)));  
        HCST_X.cov_indv(ensi,ti)=tmp.cov(1,2);
    end
    HCST_X.corr_med(ti)=median(HCST_X.corr(:,ti));
    HCST_X.corr_mean(ti)=mean(HCST_X.corr(:,ti));
    HCST_X.sig_t_med(ti)=squeeze(median(HCST_X.sig_t_indv(:,ti),1));
    HCST_X.sig_t_mean(ti)=squeeze(mean(HCST_X.sig_t_indv(:,ti),1));
    HCST_X.cov_med(ti)=median(HCST_X.cov_indv(:,ti));
    HCST_X.cov_mean(ti)=mean(HCST_X.cov_indv(:,ti),1);
end
% for HCST_unini
for ti=1:cfg.pred_len*100
    tmp.tmin=(cfg.t_ini_min-cfg.t_ini_min_CTL+cfg.pred_len)*100+mod(ti,100);
    tmp.tmax=(cfg.t_ini_max-cfg.t_ini_min_CTL)*100+mod(ti,100);
    tmp.ts=1+cfg.pred_len-floor(ti/100);
    tmp.te=cfg.t_ini_max-cfg.t_ini_min+1-floor(ti/100);

    [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        HCST_unini_X.data_lt_em(ti,tmp.ts:tmp.te));  % (lt=500, 60~85)
    HCST_unini_X.corr_em(ti)=tmp.corr(1,2);
    HCST_unini_X.corr_em_p(ti)=tmp.corr_p(1,2);
    HCST_unini_X.sig_t_em(ti)=std(HCST_unini_X.data_lt_em(ti,tmp.ts:tmp.te));
    tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
        squeeze(HCST_unini_X.data_lt_em(ti,tmp.ts:tmp.te)));  % (lt=500, 60~85)
    HCST_unini_X.cov_em(ti)=tmp.cov(1,2);
    for ensi=1:cfg.ensnum
        [tmp.corr, tmp.corr_p]=corrcoef(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            HCST_unini_X.data_lt(ensi,ti,tmp.ts:tmp.te));  % (lt=500, 60~85)
        HCST_unini_X.corr(ensi,ti)=tmp.corr(1,2);
        HCST_unini_X.corr_p(ensi,ti)=tmp.corr_p(1,2);
        HCST_unini_X.sig_t_indv(ensi,ti)=std(HCST_unini_X.data_lt(ensi,ti,tmp.ts:tmp.te));
        tmp.cov=cov(data_comb.X(tmp.tmin:100:tmp.tmax),  ...
            squeeze(HCST_unini_X.data_lt(ensi,ti,tmp.ts:tmp.te)));  % (lt=500, 60~85)
        HCST_unini_X.cov_indv(ensi,ti)=tmp.cov(1,2);
    end
    HCST_unini_X.corr_med(ti)=median(HCST_unini_X.corr(:,ti));
    HCST_unini_X.corr_mean(ti)=mean(HCST_unini_X.corr(:,ti));
    HCST_unini_X.sig_t_med(ti)=squeeze(median(HCST_unini_X.sig_t_indv(:,ti),1));
    HCST_unini_X.sig_t_mean(ti)=squeeze(mean(HCST_unini_X.sig_t_indv(:,ti),1));
    HCST_unini_X.cov_med(ti)=median(HCST_unini_X.cov_indv(:,ti));
    HCST_unini_X.cov_mean(ti)=mean(HCST_unini_X.cov_indv(:,ti),1);
end

% % % % % %% all plot
% % % % % plot(HCST_X.corr_em)
% % % % % hold on
% % % % % plot(HCST_unini_X.corr_em)
% % % % % plot(HCST_X.corr_med)
% % % % % plot(HCST_unini_X.corr_med)
% % % % % hold off



% assign corr values
clear plots_hc_em plots_hc_med plots_hc_mean plots_hc_unini_em plots_hc_unini_em plots_hc_unini_med plots_hc_unini_mean
clear plots_hc plots_hc_unini plots_hc_cov_indv  plots_hc_unini_cov_indv plots_hc_sig_t_indv plots_hc_unini_sig_t_indv
clear plots_hc_sig_t_em plots_hc_sig_t_mean plots_hc_cov_em plots_hc_cov_mean plots_hc_unini_sig_t_em plots_hc_unini_sig_t_mean plots_hc_unini_cov_em plots_hc_unini_cov_mean
tt=[1, 1*100:100:cfg.pred_len*100];
for ti=1:length(tt)
    ttt=tt(ti);
    plots_hc_em(ti)=HCST_X.corr_em(ttt);
    plots_hc_med(ti)=HCST_X.corr_med(ttt);
    plots_hc_mean(ti)=HCST_X.corr_mean(ttt);
    
    plots_hc_unini_em(ti)=HCST_unini_X.corr_em(ttt);
    plots_hc_unini_med(ti)=HCST_unini_X.corr_med(ttt);
    plots_hc_unini_mean(ti)=HCST_unini_X.corr_mean(ttt);
    
    
    for ensi=1:cfg.ensnum
        plots_hc(ensi,ti)=HCST_X.corr(ensi,ttt);
        plots_hc_unini(ensi,ti)=HCST_unini_X.corr(ensi,ttt);
        plots_hc_cov_indv(ensi,ti)=HCST_X.cov_indv(ensi,ttt);
        plots_hc_unini_cov_indv(ensi,ti)=HCST_unini_X.cov_indv(ensi,ttt);
        plots_hc_sig_t_indv(ensi,ti)=HCST_X.sig_t_indv(ensi,ttt);
        plots_hc_unini_sig_t_indv(ensi,ti)=HCST_unini_X.sig_t_indv(ensi,ttt);
    end
    
    plots_hc_sig_t_em(ti)=HCST_X.sig_t_em(ttt);
    plots_hc_sig_t_mean(ti)=HCST_X.sig_t_mean(ttt);
    plots_hc_cov_em(ti)=HCST_X.cov_em(ttt);
    plots_hc_cov_mean(ti)=HCST_X.cov_mean(ttt);
    
    plots_hc_unini_sig_t_em(ti)=HCST_unini_X.sig_t_em(ttt);
    plots_hc_unini_sig_t_mean(ti)=HCST_unini_X.sig_t_mean(ttt);
    plots_hc_unini_cov_em(ti)=HCST_unini_X.cov_em(ttt);
    plots_hc_unini_cov_mean(ti)=HCST_unini_X.cov_mean(ttt);

end



% % % % overestimation factor by number of ensemble members
% % % 
% % % % overestimation rate, mean(sig(idv))/sig(em)
% % % % plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
% % % %     'k', 'parent', axRH, 'linewidth', 2)
% % % 
% % % for ti=1:length(tt)
% % %     [a,b]=sort(plots_hc(:,ti));
% % %     ind_med_1=b(cfg.ensnum/2);
% % %     ind_med_2=b(cfg.ensnum/2+1);
% % %     plots_hc_med_corresponding_cov(ti)=(plots_hc_cov_indv(ind_med_1,ti) + plots_hc_cov_indv(ind_med_2,ti))/2;
% % %     plots_hc_med_corresponding_sig_t(ti)=(plots_hc_sig_t_indv(ind_med_1,ti) + plots_hc_sig_t_indv(ind_med_2,ti))/2;
% % % 
% % %     [a,b]=sort(plots_hc_unini(:,ti));
% % %     ind_med_1=b(10);
% % %     ind_med_2=b(11);
% % %     plots_hc_unini_med_corresponding_cov(ti)=(plots_hc_unini_cov_indv(ind_med_1,ti) + plots_hc_unini_cov_indv(ind_med_2,ti))/2;
% % %     plots_hc_unini_med_corresponding_sig_t(ti)=(plots_hc_unini_sig_t_indv(ind_med_1,ti) + plots_hc_unini_sig_t_indv(ind_med_2,ti))/2;
% % % end
% % % 
% % % [plots_hc_cov_em./plots_hc_med_corresponding_cov]
% % % [plots_hc_med_corresponding_sig_t./plots_hc_sig_t_em]
% % % 
% % % plot_legend(7)=plot([tt, cfg.pred_len*100+100], [plots_hc_em./plots_hc_mean, plots_hc_unini_em(end)./plots_hc_unini_mean(end)], ...
% % %     'k', 'parent', axRH, 'linewidth', 2.*fig_weig);
% % % hold on
% % % % plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
% % % %     'k--', 'parent', axRH, 'linewidth', 2)
% % % 
% % % plot_legend(8)=plot([tt, cfg.pred_len*100+100], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
% % %     'k--', 'linewidth', 2.*fig_weig);
% % % set(plot_legend(8),'parent', axRH)  %% it works!!!!!!!!!!!!!!!!!!!!!!!


% % plots_hc_sig_t_em(end) <-> plots_hc_sig_t_idv(:,end)
% % plots_hc_cov_em(end) <-> plots_hc_cov_idv(:,end)  mean(plots_hc_cov_idv(:,end)) median(plots_hc_cov_idv(:,end))

% plots_hc_cov_em(end)./plots_hc_sig_t_em(end)
% plot(plots_hc_cov_idv(:,end)./plots_hc_sig_t_idv(:,end))



%% axis define
tmp.t_dummy=cfg.pred_len*100+100;
% close all;
clear plot_legend

% figH = figure;
axRH = axes('color','none');
axLH = axes('color','none');
ax_pos = get(axLH,'position');

subplot(2,1,2,axRH)
subplot(2,1,2,axLH)


% HCST_em
hold on
for ti=length(tt):-1:1  % for color legend, last should be blue
    ttt=tt(ti);
    plot_legend(1)=scatter(ttt, plots_hc_em(ti), 'filled', 'parent', axLH);
    plot_legend(1).Marker='s';
    % s.MarkerEdgeColor=[255 140 0]/255; % dark orange
    plot_legend(1).MarkerEdgeColor=[tmp.parula(ttt+100,:)]; 
    plot_legend(1).MarkerFaceColor = plot_legend(1).MarkerEdgeColor;
    plot_legend(1).MarkerFaceAlpha=1;
    plot_legend(1).MarkerEdgeAlpha=plot_legend(1).MarkerFaceAlpha;
    plot_legend(1).SizeData=400.*fig_weig;
    plot_legend(1).AlphaDataMode='manual';
    plot_legend(1).AlphaDataMapping='none';
end

% HCST_members
for ti=length(tt):-1:1
    ttt=tt(ti);
    for ensi=1:cfg.ensnum
        plot_legend(2)=scatter(ttt, plots_hc(ensi,ti), 'filled', 'parent', axLH);
        plot_legend(2).Marker='p';
        % s.MarkerEdgeColor=[255 140 0]/255; % dark orange
        plot_legend(2).MarkerEdgeColor=[tmp.parula(ttt+100,:)]; 
        plot_legend(2).MarkerFaceColor = plot_legend(2).MarkerEdgeColor;
        plot_legend(2).MarkerFaceAlpha=0.3;
        plot_legend(2).MarkerEdgeAlpha=plot_legend(2).MarkerFaceAlpha;
        plot_legend(2).SizeData=200.*fig_weig;
        plot_legend(2).AlphaDataMode='manual';
        plot_legend(2).AlphaDataMapping='none';
    end
end

% HCST_members_median
for ti=length(tt):-1:1
    ttt=tt(ti);
    plot_legend(5)=scatter(ttt, plots_hc_med(ti), 'filled', 'parent', axLH);
    plot_legend(5).Marker='o';
    % s.MarkerEdgeColor=[255 140 0]/255; % dark orange
    plot_legend(5).MarkerEdgeColor=[tmp.parula(ttt+100,:)]; 
    plot_legend(5).MarkerFaceColor = plot_legend(5).MarkerEdgeColor;
    plot_legend(5).MarkerFaceAlpha=1;
    plot_legend(5).MarkerEdgeAlpha=plot_legend(5).MarkerFaceAlpha;
    plot_legend(5).SizeData=400.*fig_weig;
    plot_legend(5).AlphaDataMode='manual';
    plot_legend(5).AlphaDataMapping='none';
end


% HCST_unini_em
plot_legend(3)=scatter(tmp.t_dummy, plots_hc_unini_em(6), 'filled', 'parent', axLH);
plot_legend(3).Marker='s';
plot_legend(3).MarkerEdgeColor=[255 0 0]/255;  % red
plot_legend(3).MarkerFaceColor = plot_legend(3).MarkerEdgeColor;
plot_legend(3).MarkerFaceAlpha=1;
plot_legend(3).MarkerEdgeAlpha=plot_legend(3).MarkerFaceAlpha;
plot_legend(3).SizeData=400.*fig_weig;
plot_legend(3).AlphaDataMode='manual';
plot_legend(3).AlphaDataMapping='none';


% HCST_unini_members
for ensi=1:cfg.ensnum
    plot_legend(4)=scatter(tmp.t_dummy, plots_hc_unini(ensi,6), 'filled', 'parent', axLH);
    plot_legend(4).Marker='p';
    plot_legend(4).MarkerEdgeColor=[255 0 0]/255; 
    plot_legend(4).MarkerFaceColor = plot_legend(4).MarkerEdgeColor;
    plot_legend(4).MarkerFaceAlpha=0.3;
    plot_legend(4).MarkerEdgeAlpha=plot_legend(4).MarkerFaceAlpha;
    plot_legend(4).SizeData=200.*fig_weig;
    plot_legend(4).AlphaDataMode='manual';
    plot_legend(4).AlphaDataMapping='none';
end


% HCST_unini_members_median
ttt=tt(ti);
plot_legend(6)=scatter(tmp.t_dummy, plots_hc_unini_med(ti), 'filled', 'parent', axLH);
plot_legend(6).Marker='o';
% s.MarkerEdgeColor=[255 140 0]/255; % dark orange
plot_legend(6).MarkerEdgeColor=[255 0 0]/255;  % red
plot_legend(6).MarkerFaceColor = plot_legend(6).MarkerEdgeColor;
plot_legend(6).MarkerFaceAlpha=1;
plot_legend(6).MarkerEdgeAlpha=plot_legend(6).MarkerFaceAlpha;
plot_legend(6).SizeData=300.*fig_weig;
plot_legend(6).AlphaDataMode='manual';
plot_legend(6).AlphaDataMapping='none';

for tttti=1:cfg.pred_len
    tmp.xticks{tttti}=num2str(tttti-1);
end
tmp.xticks{cfg.pred_len+1}=num2str(tttti);
tmp.xticks{cfg.pred_len+2}='$$ X_{UI} $$';

% set(axLH, 'yaxislocation', 'left', ...
%     'position',ax_pos+[0 0.02 -0.01 -0.02], ...
%     'xlim', [-20 tmp.t_dummy+20], ...
%     'xtick', [tt tmp.t_dummy], ...
%     'xticklabel', tmp.xticks, ...
%     'Ycolor', 'b', 'box', 'off', ...
%     'FontSize', 22.*fig_weig);
% axLH.TickLabelInterpreter = 'latex';
% ylabel(axLH, '$$ r $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% xlabel(axLH,'$$ lt $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)

set(axLH, 'yaxislocation', 'left', ...
    'position',ax_pos+[0 0.02 -0.01 -0.02], ...
    'xlim', [-20 tmp.t_dummy+20], ...
    'xtick', [tt tmp.t_dummy], ...
    'xticklabel', tmp.xticks, ...
    'Ycolor', 'b', 'box', 'off', ...
    'FontSize', 22.*fig_weig);
axLH.TickLabelInterpreter = 'latex';
ylabel(axLH, '$$ r $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
xlabel(axLH,'$$ l\tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)


% overestimation rate, mean(sig(idv))/sig(em)
% plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
%     'k', 'parent', axRH, 'linewidth', 2)

for ti=1:length(tt)
    [a,b]=sort(plots_hc(:,ti));
    ind_med_1=b(cfg.ensnum/2);
    ind_med_2=b(cfg.ensnum/2+1);
    plots_hc_med_corresponding_cov(ti)=(plots_hc_cov_indv(ind_med_1,ti) + plots_hc_cov_indv(ind_med_2,ti))/2;
    plots_hc_med_corresponding_sig_t(ti)=(plots_hc_sig_t_indv(ind_med_1,ti) + plots_hc_sig_t_indv(ind_med_2,ti))/2;

    [a,b]=sort(plots_hc_unini(:,ti));
    ind_med_1=b(cfg.ensnum/2);
    ind_med_2=b(cfg.ensnum/2+1);
    plots_hc_unini_med_corresponding_cov(ti)=(plots_hc_unini_cov_indv(ind_med_1,ti) + plots_hc_unini_cov_indv(ind_med_2,ti))/2;
    plots_hc_unini_med_corresponding_sig_t(ti)=(plots_hc_unini_sig_t_indv(ind_med_1,ti) + plots_hc_unini_sig_t_indv(ind_med_2,ti))/2;
end

[plots_hc_cov_em./plots_hc_med_corresponding_cov];
[plots_hc_med_corresponding_sig_t./plots_hc_sig_t_em];

[plots_hc_unini_cov_em./plots_hc_unini_med_corresponding_cov];
[plots_hc_unini_med_corresponding_sig_t./plots_hc_unini_sig_t_em];



% % mean version
% plot_legend(7)=plot([tt, tmp.t_dummy], [plots_hc_em./plots_hc_mean, plots_hc_unini_em(end)./plots_hc_unini_mean(end)], ...
%     'k', 'parent', axRH, 'linewidth', 2.*fig_weig);
% % hold on
% % plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
% %     'k--', 'parent', axRH, 'linewidth', 2)
% plot_legend(8)=plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
%     'k--', 'linewidth', 2.*fig_weig);
% set(plot_legend(8),'parent', axRH)  %% it works!!!!!!!!!!!!!!!!!!!!!!!

% median version
plot_legend(7)=plot([tt, tmp.t_dummy], [plots_hc_em./plots_hc_med, plots_hc_unini_em(end)./plots_hc_unini_med(end)], ...
    'k', 'parent', axRH, 'linewidth', 2.*fig_weig);
% hold on
% plot([tt, tmp.t_dummy], [plots_hc_sig_t_mean./plots_hc_sig_t_em, plots_hc_unini_sig_t_mean(end)./plots_hc_unini_sig_t_em(end)], ...
%     'k--', 'parent', axRH, 'linewidth', 2)
plot_legend(8)=plot([tt, tmp.t_dummy], [plots_hc_med_corresponding_sig_t./plots_hc_sig_t_em, plots_hc_unini_med_corresponding_sig_t(end)./plots_hc_unini_sig_t_em(end)], ...
    'k-.', 'linewidth', 2.*fig_weig);
set(plot_legend(8),'parent', axRH)  %% it works!!!!!!!!!!!!!!!!!!!!!!!

plot_legend(9)=plot([tt, tmp.t_dummy], [plots_hc_cov_em./plots_hc_med_corresponding_cov, plots_hc_unini_cov_em(end)./plots_hc_unini_med_corresponding_cov(end)], ...
    'k:', 'linewidth', 2.*fig_weig);
set(plot_legend(9),'parent', axRH)  %% it works!!!!!!!!!!!!!!!!!!!!!!!

% set(axRH, 'color', 'none', ...
%     'yaxislocation', 'right', ...
%     'position', ax_pos+[0 0.02 -0.01 -0.02], ...
%     'xlim', [-20 tmp.t_dummy+20], ...
%     'xtick', [tt tmp.t_dummy], ...
%     'xticklabel', tmp.xticks, ...
%     'ycolor', 'k', 'box', 'off', ...
%     'FontSize', 22.*fig_weig);
% axRH.TickLabelInterpreter = 'latex';

set(axRH, 'color', 'none', ...
    'yaxislocation', 'right', ...
    'xlim', [-20 tmp.t_dummy+20], ...
    'xtick', [tt tmp.t_dummy], ...
    'xticklabel', tmp.xticks, ...
    'ycolor', 'k', 'box', 'off', ...
    'FontSize', 22.*fig_weig);
axRH.TickLabelInterpreter = 'latex';

% ylabel(axRH, '\alpha', 'fontsize', 25.*fig_weig)
ylabel(axRH, 'Overestimaion Factor', 'fontsize', 25.*fig_weig)

xlabel(axRH,'$$ l\tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
grid(axRH)
% legend([plot_legend(7:9)], {'$$ r_{C, \overline{Y}}/\overline{r_{C,Y}} $$', ...
%     '$$ \overline{\sigma_Y}/\sigma_{\overline{Y}} $$', ...
%     '$$ \overline{S_{C,Y}}/S_{C,\overline{Y}} $$'}, ...
%     'location', 'northeast', 'interpreter', 'latex', ...
%     'Orientation', 'vertical', 'Fontsize', 25.*fig_weig);  %% for mean

% legend([plot_legend(7:9)], {'$$ r_{C, \overline{Y}}/Md(r_{C,Y}) $$', ...
%     '$$ \sigma_{Y_{Md(r)}}/\sigma_{\overline{Y}} $$', ...
%     '$$ S_{C,\overline{Y}}/S_{C,Y_{Md(r)}} $$'}, ...
%     'location', 'northeast', 'interpreter', 'latex', ...
%     'Orientation', 'vertical', 'Fontsize', 25.*fig_weig);  %% for median

% legend([plot_legend(2), plot_legend(1), plot_legend(5), plot_legend(7:9)], ...
%     { '$$ r_{C,Y} $$', ...
%     '$$ r_{C, \overline{Y}} $$', ...
%     '$$ Md(r_{C,Y}) $$', ...
%     '$$ r_{C, \overline{Y}}/Md(r_{C,Y}) $$', ...
%     '$$ \sigma_{Y_{Md(r)}}/\sigma_{\overline{Y}} $$', ...
%     '$$ S_{C,\overline{Y}}/S_{C,Y_{Md(r)}} $$'}, ...
%     'location', 'northeast', 'interpreter', 'latex', ...
%     'Orientation', 'vertical', 'Fontsize', 25.*fig_weig, ...
%     'NumColumns', 2);  %% for median

% [hLg]=legend([plot_legend(2), plot_legend(7), plot_legend(1), plot_legend(8), plot_legend(5), plot_legend(9)], ...
%     { '$$ r_{C,X} $$', ...
%     '$$ r_{C, \overline{X}}/Md(r_{C,X}) \hspace{5mm} $$', ...
%     '$$ r_{C, \overline{X}} $$', ...
%     '$$ \sigma_{Y_{Md(r)}}/\sigma_{\overline{X}} \hspace{5mm} $$', ...
%     '$$ Md(r_{C,X}) $$', ...
%     '$$ S_{C,\overline{X}}/S_{C,Y_{Md(r)}} $$'}, ...
%     'location', 'northoutside', 'interpreter', 'latex', ...
%     'Orientation', 'vertical', 'Fontsize', 35.*fig_weig, ...
%     'NumColumns', 3);  %% for median

[hLg]=legend([plot_legend(2), plot_legend(7), plot_legend(1), plot_legend(8), plot_legend(5), plot_legend(9)], ...
    { '$$ r_{O,X} $$', ...
    '$$ r_{O, \overline{X}}/Md(r_{O,X}) \hspace{5mm} $$', ...
    '$$ r_{O, \overline{X}} $$', ...
    '$$ \sigma_{X_{Md(r)}}/\sigma_{\overline{X}} \hspace{5mm} $$', ...
    '$$ Md(r_{O,X}) $$', ...
    '$$ S_{O,\overline{X}}/S_{O,X_{Md(r)}} $$'}, ...
    'location', 'northoutside', 'interpreter', 'latex', ...
    'Orientation', 'vertical', 'Fontsize', 35.*fig_weig, ...
    'NumColumns', 3);  %% for median

% icons = findobj(icons,'Type','patch');
% icons = findobj(icons,'Marker','none','-xor');
% set(icons(1:3),'MarkerSize',15);
% set(hLg, 'NumColumns', 3);

ax_pos_LH = get(axLH,'position');
set(axLH, 'position',ax_pos_LH+[0 0.02 -0.01 0.05]);

% set(gcf,'PaperPosition',[0 0 30 20]);
print(gcf, ['ghi_ext_flag',num2str(ext_flag), '.png'], '-dpng');



%% -----

% plots_hc_sig_t_mean(:)./plots_hc_sig_t_em(:);
% plots_hc_unini_sig_t_mean(:)./plots_hc_unini_sig_t_em(:);



% % % % % %% CTL doublearrow
% % % % % % arrow_x = mean([outerpos(1)+outerpos(3) pos(1)+pos(3)]);
% % % % % arrow_x = [outerpos(1)+outerpos(3) pos(1)+pos(3)];
% % % % % arrow_x = arrow_x(1).*2/3 + arrow_x(2).*1/3;
% % % % % % arrow_y = [pos(2) pos(2)+pos(4)]; % full y range
% % % % % yl=ylim;
% % % % % dy=(pos(4))/(yl(2)-yl(1));
% % % % % % CTL_mean=mean(data_comb.X(501:3400));
% % % % % CTL_mean=data_comb.X(3400);
% % % % % CTL_std=std(data_comb.X(501:3400));
% % % % % CTL_min=CTL_mean-CTL_std/2;
% % % % % CTL_max=CTL_mean+CTL_std/2;
% % % % % pos_CTL_min=pos(2)+(CTL_min-yl(1)).*dy;
% % % % % pos_CTL_max=pos(2)+(CTL_max-yl(1)).*dy;
% % % % % arrow_y = [pos_CTL_min pos_CTL_max]; % full y range
% % % % % 
% % % % % % annotation('textarrow',[arrow_x arrow_x], arrow_y,'String','CTL', 'linewidth', 4, 'color', 'k')
% % % % % annotation('doublearrow',[arrow_x arrow_x], arrow_y, 'linewidth', 4.*fig_weig, 'color', 'k')
% % % % % 
% % % % % %% unini doublearrow
% % % % % arrow_x = [outerpos(1)+outerpos(3) pos(1)+pos(3)];
% % % % % arrow_x = arrow_x(1).*1/3 + arrow_x(2).*2/3;
% % % % % % unini_mean=mean(HCST_unini_X.data_lt_em(500,1:25)); % lt=500, min (t=55~80)
% % % % % unini_mean=HCST_unini_X.data_lt_em(500,25); % lt=500, min (t=55~80)
% % % % % unini_std=std(HCST_unini_X.data_lt_em(500,1:25));
% % % % % unini_min=unini_mean-unini_std/2;
% % % % % unini_max=unini_mean+unini_std/2;
% % % % % pos_unini_min=pos(2)+(unini_min-yl(1)).*dy;
% % % % % pos_unini_max=pos(2)+(unini_max-yl(1)).*dy;
% % % % % arrow_y = [pos_unini_min pos_unini_max]; % full y range
% % % % % 
% % % % % % annotation('textarrow',[arrow_x arrow_x], arrow_y,'String','CTL', 'linewidth', 4, 'color', 'k')
% % % % % annotation('doublearrow',[arrow_x arrow_x], arrow_y, 'linewidth', 4.*fig_weig, 'color', [255 0 0]/255)



% %% normal scatter plot, uninitialized run
% for ensi=1:20
%     s=scatter(71:80,squeeze(HCST_unini_X.data_lt(ensi,1,17:26)),'filled');
%     s.Marker='p';
%     s.MarkerEdgeColor=[255 208 0]/255; % yellow
%     s.MarkerFaceColor=[255 208 0]/255; % yellow
%     s.SizeData=100;
% end
% s=scatter(71:80,squeeze(HCST_unini_X.data_lt_em(1,17:26)),'filled');
% s.Marker='s';
% s.MarkerEdgeColor=[255 140 0]/255; % dark orange
% s.MarkerFaceColor=[255 140 0]/255; % dark orange
% s.SizeData=150;
% hold off
% set(gca,'fontsize',15)


% corrcoef(data_comb.X(3401:3401+500-1), hcst_ensm_unini.X.i0085)

% color sky blue : [135 206 235]/255;
% color blue : [19 0 255]/255;
% color green : [143 189 136]/255; (star)
% color light green : [144 238 143]/255; (circle)
% color dark green : [0 128 0]/255; (square)
% color light red : [255 127 127]/255;
% color dark red : [139 0 0]/255;
% color red : [255 0 0]/255;
% color dark yellow : [255 208 0]/255;
% color dark orange : [255 140 0]/255;














% % % %% external forcing
% % % ext_forc(1:5000)=(1:5000)/250;
% % % data_comb.X_ext=data_comb.X+ext_forc;
% % % %% normal plot
% % % plot(data_comb.t, data_comb.X_ext, 'k', 'linewidth', 3)
% % % hold on
% % % hcst_ensm.X.i0055=zeros(1,500);
% % % for ensi=1:20
% % %     str.ens=num2str(ensi, '%02i');
% % %     plot(55:0.01:59.99, hcst_comb.X.(['ens',str.ens]).i0055 + ext_forc(401:401+500-1),'b');
% % %     hcst_ensm.X.i0055=hcst_ensm.X.i0055+hcst_comb.X.(['ens',str.ens]).i0055;
% % % end
% % % hcst_ensm.X.i0055=hcst_ensm.X.i0055./20;
% % % plot(55:0.01:59.99, hcst_ensm.X.i0055 + ext_forc(401:401+500-1) ,'r','linewidth', 2);
% % % 
% % % %% unini plot
% % % % plot(data_comb.t, data_comb.X, 'k', 'linewidth', 3)
% % % % hold on
% % % hcst_ensm_unini.X.i0085=zeros(1,500);
% % % for ensi=1:20
% % %     str.ens=num2str(ensi, '%02i');
% % %     plot(85:0.01:89.99, hcst_comb_unini.X.(['ens',str.ens]).i0085 + ext_forc(3401:3401+500-1),'b');
% % %     hcst_ensm_unini.X.i0085=hcst_ensm_unini.X.i0085+hcst_comb_unini.X.(['ens',str.ens]).i0085;
% % % end
% % % hcst_ensm_unini.X.i0085=hcst_ensm_unini.X.i0085./20;
% % % plot(85:0.01:89.99, hcst_ensm_unini.X.i0085 + ext_forc(3401:3401+500-1),'r','linewidth', 2);
% % % hold off
% % % 
% % % corrcoef(data_comb.X(3401:3401+500-1), hcst_ensm_unini.X.i0085)
% % % corrcoef(data_comb.X(3401:3401+500-1)+ext_forc(3401:3401+500-1)', hcst_ensm_unini.X.i0085+ext_forc(3401:3401+500-1)')
% % % 
% % % 
% % % 
% % % %% corr ensemble mean
% % % corrcoef(data_comb.X(401:401+500-1), hcst_ensm.X.i0065)
% % % stat.a=cov(data_comb.X(401:401+500-1), hcst_ensm.X.i0065);
% % % stat.b=std(data_comb.X(401:401+500-1));
% % % stat.c=std(hcst_ensm.X.i0065);
% % % for ensi=1:20
% % %     str.ens=num2str(ensi, '%02i');
% % %     tmp.corr=corrcoef(data_comb.X(401:401+500-1), hcst_comb.X.(['ens',str.ens]).i0065);
% % %     corrs_hcst(ensi)=tmp.corr(1,2);
% % %     tmp.cov=cov(data_comb.X(401:401+500-1), hcst_comb.X.(['ens',str.ens]).i0065,0);
% % %     covs_hcst(ensi)=tmp.cov(1,2);
% % %     stds_hcst(ensi)=std(hcst_comb.X.(['ens',str.ens]).i0065);
% % % end
% % % mean(corrs_hcst)
% % % mean(covs_hcst)
% % % mean(stds_hcst)
% % % 
% % % 
% % % corr_ensm=stat.a(1,2)/(stat.b*stat.c) % ensm corr
% % % corr_ensm.*(stat.c./mean(stds_hcst))
% % % 
% % % a=[1,3,5,6,7];
% % % b=[2,5,2,7,8];
% % % mean(a.*b)
% % % mean(a).*mean(b)
% % % 
% % % 
% % % 
% % % 
% % % %% correlation case study for ensemble mean & individual base using random function  -------------------------------
% % % for i=1:50
% % %     for j=1:100
% % %         randvar_LE(i,j)=rand(1)*5+j/100;
% % %         randvar_HCST(i,j)=rand(1)*5+j/80;
% % %     end
% % % end
% % % 
% % % 
% % % for i=1:50
% % %     for j=1:50
% % %         tmp.corr=corrcoef(randvar_LE(i,:), randvar_HCST(j,:));
% % %         randvar_corr(i,j)=tmp.corr(1,2);
% % %     end
% % % end
% % % mean_randvar_LE=mean(randvar_LE,1);
% % % mean_randvar_HCST=mean(randvar_HCST,1);
% % % 
% % % figure(1)
% % % plot(mean_randvar_LE,'b');
% % % hold on
% % % plot(mean_randvar_HCST,'r');
% % % hold off
% % % 
% % % figure(2)
% % % plot(randvar_LE(1,:),'b');
% % % hold on
% % % plot(randvar_HCST(1,:),'r');
% % % hold off
% % % 
% % % 
% % % corrcoef(mean_randvar_LE, mean_randvar_HCST)
% % % corrcoef(randvar_LE(1,:), mean_randvar_HCST)
% % % 
% % % mean(randvar_corr(:))
% % % 
% % % 
% % % 
% % % 
% % % 
% % % %% two axis example
% % % % % % figH = figure;
% % % % % % axLH = gca;
% % % % % % axRH = axes('color','none');
% % % % % % mslplot{1}=plot(inputyear,Model.amp_S2(3,:), 'b','parent',axLH);
% % % % % % mslplot{2}=plot(inputyear,Obs.amp_S2(3,:), 'k','parent',axRH);
% % % % % % ylabel(axLH,'Model S2 Tidal amplitude (cm)')
% % % % % % ylabel(axRH,'Obs S2 Tidal amplitude (cm)')
% % % % % % ax_pos = get(axLH,'position');
% % % % % % set(axLH,'yaxislocation','left','position',ax_pos+[0 0.02 -0.01 -0.02]);
% % % % % % set(axRH,'color','none','yaxislocation','right','xtick', inputyear, 'position', ax_pos+[0 0.02 -0.01 -0.02]);
% % % % % % %             set(axRH,'color','none','yaxislocation','right');
% % % % % % 
% % % % % % %             set(axRH,'xticklabel',[], 'xtick', get(axLH,'xtick'));
% % % % % % set(axLH,'xlim', get(axRH, 'xlim'), 'xtick', get(axRH,'xtick'),'xticklabel',[], 'xaxislocation','top');
% % % % % % set(axLH,'ycolor','b', 'box', 'off', 'FontSize',15);
% % % % % % set(axRH,'ycolor','k', 'box', 'off', 'FontSize',15);
% % % % % % xlabel(axRH, 'Year');
% % % % % % 
% % % % % % title([num2str(min(inputyear),'%04i'),'-',num2str(max(inputyear),'%04i')])
% % % % % % % datetick('x','yyyy','keepticks')
% % % % % % axis tight;
% % % % % % % ylim(meanplotlev2)
% % % % % % set(mslplot{1},'LineWidth',2);
% % % % % % set(mslplot{2},'LineWidth',2);
% % % % % % grid on
% % % % % % 
% % % % % % lgd=legend([mslplot{1} mslplot{2}], 'Model','TG-UST');
% % % % % % 
% % % % % % %             lgd=legend('Model','TG-UST');
% % % % % % set(lgd,'FontSize',15);
% % % % % % set(lgd,'Position',[0.13 0.88, 0.775, 0.03]);
% % % % % % set(lgd,'Orientation','horizontal');
% % % % % % 
% % % % % % set(gcf,'PaperPosition', [0 0 36 12]) 
% % % % % % saveas(gcf,pngname,'tif');
% % % % % % grid off
% % % % % % close all;
% % % 
% % % 
% % % nums=100;
% % % wf=1;
% % % a=(rand(nums,1))*wf;
% % % b=(zeros(nums,1)+0.5)*wf;
% % % c=(rand(nums,1))*wf;
% % % 
% % % sum((a-b).^2)/nums
% % % sum((a-c).^2)/nums
% % % 
% % % sum(abs((a-b)))/nums
% % % sum(abs((a-c)))/nums
% % % 
% % % 
% % % 25*6
% % % 100*3
% % % 
% % % e


% % %         %% Fire %%
% % %         F_1(ky:ly,2:nx-1) = (diffc_F(ky:ly,2:nx-1).*F_t_per(ky:ly,2:nx-1) + diffc_F(ky+1:ly+1,2:nx-1).*F_t_per(ky+1:ly+1,2:nx-1)) ...
% % %             ./ (diffc_F(ky:ly,2:nx-1)+diffc_F(ky+1:ly+1,2:nx-1));
% % % 
% % %         F_2(ky:ly,2:nx-1) = (diffc_F(ky:ly,2:nx-1).*F_t_per(ky:ly,2:nx-1)+diffc_F(ky-1:ly-1,2:nx-1).*F_t_per(ky-1:ly-1,2:nx-1))./ (diffc_F(ky:ly,2:nx-1)+diffc_F(ky-1:ly-1,2:nx-1));
% % %         F_3(ky:ly,2:nx-1) = (diffc_F(ky:ly,2:nx-1).*F_t_per(ky:ly,2:nx-1)+diffc_F(ky:ly,3:nx-1+1).*F_t_per(ky:ly,3:nx-1+1))./ (diffc_F(ky:ly,2:nx-1)+diffc_F(ky:ly,3:nx-1+1));
% % %         F_4(ky:ly,2:nx-1) = (diffc_F(ky:ly,2:nx-1).*F_t_per(ky:ly,2:nx-1)+diffc_F(ky:ly,1:nx-1-1).*F_t_per(ky:ly,1:nx-1-1))./ (diffc_F(ky:ly,2:nx-1)+diffc_F(ky:ly,1:nx-1-1)); 
% % %         
% % %         % carrying capacity of Fire
% % %         CC_F_2d(:,:)  =  sum_V_t_new_h(:,:).*PDFGENUS_F(:,:)*firescaling;  % modified by Axel
% % %         
% % %         % time integration: Fire
% % %         dF(ky:ly,2:nx-1) = gr_F.*F_t_per(ky:ly,2:nx-1).*(1-((F_t_per(ky:ly,2:nx-1))./CC_F_2d(ky:ly,2:nx-1)))+wblrnd(1,0.1,154,286)*sigma.*PDFGENUS_F(ky:ly,2:nx-1);
% % %         dF(isnan(dF)) = 0;
% % % 
% % % 
% % %         F_t_new(ky:ly,2:nx-1,i) = F_t_per(ky:ly,2:nx-1) + ...
% % %             dt*(diffc_F(ky:ly,2:nx-1)./(dy.^2)*2.* (F_1(ky:ly,2:nx-1) - 2*F_t_per(ky:ly,2:nx-1) + F_2(ky:ly,2:nx-1)) ) + dt*(diffc_F(ky:ly,2:nx-1)./(dxmo.^2)*2.*(F_3(ky:ly,2:nx-1) -2*F_t_per(ky:ly,2:nx-1) +F_4(ky:ly,2:nx-1))) + dt*dF(ky:ly,2:nx-1);
% % %         
% % %         
% % %         
% % %         F_t_new(ky:ly,2:nx-1,i) = F_t_per(ky:ly,2:nx-1) + dt*dF(ky:ly,2:nx-1);
% % %         F_t_per(:,:) = F_t_new(:,:,i); % Burned area fraction (unit: %)
% % % 
% % %         F_t_per(F_t_per>1e5)=0;
% % %         F_t_per(F_t_per<0)=0;
% % %         F_t_per(isnan(F_t_per))=0;
