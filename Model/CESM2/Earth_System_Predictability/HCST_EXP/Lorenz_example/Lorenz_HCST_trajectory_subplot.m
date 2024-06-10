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

fname=['/Users/kimyy/Desktop/backup/Research/Postdoc/03_IBS', ...
    '/2022_predictability_assimilation_run/Lorenz_predictability_example/', ...
    'lorenz_temp.mat'];

if exist(fname)==0

%% read control run
for t_ini=51:cfg.t_ini_max_CTL
    str.t_ini=num2str(t_ini, '%04i');
    dir.ctl_tmp=[dir.ctl, '/', 'out'];
    fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
    if ext_flag==0
%% raw
    CTL_X.(['p',str.t_ini])=single(ncread(fname_ctl, 'X'));
    CTL_Y.(['p',str.t_ini])=single(ncread(fname_ctl, 'Y'));
    CTL_Z.(['p',str.t_ini])=single(ncread(fname_ctl, 'Z'));
    else
%% with external forcing
    CTL_X.(['p',str.t_ini])=ncread(fname_ctl, 'X') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    CTL_Y.(['p',str.t_ini])=ncread(fname_ctl, 'Y') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    CTL_Z.(['p',str.t_ini])=ncread(fname_ctl, 'Z') + (t_ini-51)/ext_factor+[0.01:0.01:1]';
    end
end

% for t_ini=51:cfg.t_ini_max_CTL
%     str.t_ini=num2str(t_ini, '%04i');
%     dir.ctl_tmp=[dir.ctl, '/', 'out'];
%     fname_ctl=[dir.ctl_tmp, '/', 'out.', str.t_ini, '.nc'];
%     CTL_X_raw.(['p',str.t_ini])=ncread(fname_ctl, 'X');
%     CTL_X_ext.(['p',str.t_ini])=0 + (t_ini-51)/ext_factor+[0.01:0.01:1]';
% end
% for t_ini=cfg.t_ini_min_CTL:cfg.t_ini_max_CTL
%     data_comb_raw.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X_raw.(['p',str.t_ini]);
%     data_comb_ext.X((t_ini-cfg.t_ini_min_CTL)*100+1:(t_ini-cfg.t_ini_min_CTL+1)*100)=CTL_X_ext.(['p',str.t_ini]);
% end
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
                HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'X'));
                HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'Y'));
                HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'Z'));

            %% raw_unini
            HCST_unini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'X'));
            HCST_unini_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'Y'));
            HCST_unini_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'Z')); 
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
    save(fname, '-v7.3');
else
    load(fname);
end

hold on


%% data combine for CTL
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
    for lt_i=1:11
        lt=(lt_i-1)*100;
        
        if lt==0
            lt=1;
        end
        
        for t_ini=cfg.t_ini_min:cfg.t_ini_max
            str.t_ini=num2str(t_ini, '%04i');
            t_p=t_ini+floor((lt-1)/100);
            str.t_p=num2str(t_p, '%04i');
            HCST_X.data_lt(ensi,lt_i,t_ini-54) = ...
                HCST_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])(lt-(t_p-t_ini)*100);
            HCST_unini_X.data_lt(ensi,lt_i,t_ini-54) = ...
                HCST_unini_X.(['ens',str.ens]).(['i',str.t_ini]).(['p',str.t_p])(lt-(t_p-t_ini)*100);
        end
    end
end
HCST_X.data_lt_em=squeeze(mean(HCST_X.data_lt,1));
HCST_unini_X.data_lt_em=squeeze(mean(HCST_unini_X.data_lt,1));






% % % % % %% plot plot plot plot plot
% % % % % 
% % % % % 
% % % % % close all;
% % % % % clear plot_legend
% % % % % subplot(2,1,1)
% % % % % set(gcf,'PaperPosition',[0 0 30 45]);
% % % % % fig_weig=2;
% % % % % 
% % % % % % (cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1;
% % % % % % (cfg.t_ini_max-cfg.t_ini_min_CTL)*100;
% % % % % 
% % % % % %% normal plot
% % % % % plot_legend(1)=plot(data_comb.t((cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1:(cfg.t_ini_max-cfg.t_ini_min_CTL)*100), ...
% % % % %     data_comb.X((cfg.t_ini_min-cfg.t_ini_min_CTL+1)*100+1:(cfg.t_ini_max-cfg.t_ini_min_CTL)*100), 'k', 'linewidth', 4.*fig_weig);
% % % % % hold on
% % % % % 
% % % % % % %% using plot
% % % % % % hcst_ensm.X.i0065=zeros(1,500);
% % % % % % for ensi=1:20
% % % % % %     str.ens=num2str(ensi, '%02i');
% % % % % %     plot(65:0.01:69.99, hcst_comb.X.(['ens',str.ens]).i0065,'color', ...
% % % % % %         [19 0 255 30]/255, 'linewidth', 2); %blue, transparency: 30/255;
% % % % % %     hcst_ensm.X.i0065=hcst_ensm.X.i0065+hcst_comb.X.(['ens',str.ens]).i0065;
% % % % % % end
% % % % % % hcst_ensm.X.i0065=hcst_ensm.X.i0065./20;
% % % % % % plot(65:0.01:69.99, hcst_ensm.X.i0065,'color', [19 0 255]/255, 'linewidth', 4);  % blue
% % % % % 
% % % % % 
% % % % % 
% % % % % 
% % % % % %% using surf for gradient color plot (tstart +10 ~ tstart +10 + pred_len)
% % % % % str.t_ini=num2str(cfg.t_ini_min+10, '%04i');
% % % % % hcst_ensm.X.(['i',str.t_ini])=zeros(1,100*cfg.pred_len);
% % % % % for ensi=1:cfg.ensnum
% % % % %     str.ens=num2str(ensi, '%02i');
% % % % %     tmp.parula=parula((cfg.pred_len+2)*100); % from light blue, until darker yellow
% % % % %     plots.X=cfg.t_ini_min+10:0.01:cfg.t_ini_min+10+cfg.pred_len-0.01;
% % % % %     plots.Y=hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini]);
% % % % %     plots.Z=1:length(plots.X);
% % % % %     for ti=1:length(plots.X)
% % % % %         plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
% % % % %         plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
% % % % %         plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
% % % % %     end
% % % % % %     plots.C0(1:length(plots.X),1:2,4)=30/255; RGBa cannot be applied
% % % % % %     colormap([19 0 255]/255);
% % % % %     objsurf=surf([plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
% % % % %         plots.C0, ...  % Reshape and replicate data
% % % % %      'FaceColor', 'none', ...    % Don't bother filling faces with color
% % % % %      'EdgeColor', 'interp', ...  % Use interpolated color for edges
% % % % %      'LineWidth', 2.*fig_weig, 'EdgeAlpha', 0.2);           
% % % % %      view(2);   % Default 2-D view)
% % % % % %     objsurf.CData=[19 0 255]/255;
% % % % % %     objsurf.AlphaData=1;
% % % % %    
% % % % %     hcst_ensm.X.(['i',str.t_ini])=hcst_ensm.X.(['i',str.t_ini])+hcst_comb.X.(['ens',str.ens]).(['i',str.t_ini]);
% % % % % end
% % % % % hcst_ensm.X.(['i',str.t_ini])=hcst_ensm.X.(['i',str.t_ini])./cfg.ensnum;
% % % % % 
% % % % % % plot(65:0.01:69.99, hcst_ensm.X.i0065,'color', [19 0 255]/255, 'linewidth', 4);  % blue
% % % % % 
% % % % % %% ensemble mean plot using surf (tstart +10 ~ tstart +10 + pred_len)
% % % % % plots.X=cfg.t_ini_min+10:0.01:cfg.t_ini_min+10+cfg.pred_len-0.01;
% % % % % plots.Y=hcst_ensm.X.(['i',str.t_ini]);
% % % % % plots.Z=1:length(plots.X);
% % % % % for ti=1:length(plots.X)
% % % % %     plots.C0(ti,1:2,1)=tmp.parula(ti+100,1);
% % % % %     plots.C0(ti,1:2,2)=tmp.parula(ti+100,2);
% % % % %     plots.C0(ti,1:2,3)=tmp.parula(ti+100,3);
% % % % % end
% % % % %  objsurf=surf([plots.X(:) plots.X(:)], [plots.Y(:) plots.Y(:)], [plots.Z(:) plots.Z(:)], ...
% % % % %         plots.C0, ...  % Reshape and replicate data
% % % % %      'FaceColor', 'none', ...    % Don't bother filling faces with color
% % % % %      'EdgeColor', 'interp', ...  % Use interpolated color for edges
% % % % %      'LineWidth', 4.*fig_weig, 'EdgeAlpha', 1);            % Make a thicker line
% % % % %  view(2);   % Default 2-D view)
% % % % % 
% % % % % clear plots
% % % % % 
% % % % % 
% % % % % 
% % % % % %% unini plot (tstart + 10 + pred_len + 10 ~ tstart + 10 + pred_len + 10 + pred_len)
% % % % % % plot(data_comb.t, data_comb.X, 'k', 'linewidth', 3)
% % % % % % hold on
% % % % % str.t_ini=num2str(cfg.t_ini_min+10+cfg.pred_len+10, '%04i');
% % % % % 
% % % % % hcst_ensm_unini.X.(['i',str.t_ini])=zeros(1,cfg.pred_len*100);
% % % % % for ensi=1:cfg.ensnum
% % % % %     str.ens=num2str(ensi, '%02i');
% % % % %     plot(cfg.t_ini_min+10+cfg.pred_len+10:0.01:cfg.t_ini_min+10+cfg.pred_len+10+cfg.pred_len-0.01, ...
% % % % %         hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini]),'color', [255 0 0 30]/255, 'linewidth', 2.*fig_weig); % red, 30/255 transparency
% % % % %     hcst_ensm_unini.X.(['i',str.t_ini])=hcst_ensm_unini.X.(['i',str.t_ini])+hcst_comb_unini.X.(['ens',str.ens]).(['i',str.t_ini]);
% % % % % end
% % % % % hcst_ensm_unini.X.(['i',str.t_ini])=hcst_ensm_unini.X.(['i',str.t_ini])./cfg.ensnum;
% % % % % plot_legend(4)=plot(cfg.t_ini_min+10+cfg.pred_len+10:0.01:cfg.t_ini_min+10+cfg.pred_len+10+cfg.pred_len-0.01, hcst_ensm_unini.X.(['i',str.t_ini]),'color',[255 0 0]/255, 'linewidth', 4.*fig_weig); % dark red
% % % % % 
% % % % % 
% % % % % %% normal scatter plot, lt=1 (tstart + 1 : tstart + 10)
% % % % % for ensi=1:cfg.ensnum
% % % % %     s=scatter(cfg.t_ini_min+1:cfg.t_ini_min+10,squeeze(HCST_X.data_lt(ensi,1,2:11)),'filled');
% % % % %     s.Marker='p'; %pentagram
% % % % % %     s.MarkerEdgeColor=[143 189 136]/255; % green
% % % % %     s.MarkerEdgeColor=[tmp.parula(1+100,:)]; 
% % % % %     s.MarkerFaceColor = s.MarkerEdgeColor;
% % % % %     s.SizeData=100.*fig_weig;
% % % % %     s.MarkerFaceAlpha=1; %0.8
% % % % %     s.MarkerEdgeAlpha=s.MarkerFaceAlpha;
% % % % %     s.AlphaDataMode='manual';
% % % % % end
% % % % % plot_legend(2)=scatter(cfg.t_ini_min+1:cfg.t_ini_min+10,squeeze(HCST_X.data_lt_em(1,2:11)),'filled');
% % % % % plot_legend(2).Marker='s'; %square
% % % % % % s.MarkerEdgeColor=[0 128 0]/255; % dark green
% % % % % plot_legend(2).MarkerEdgeColor=[tmp.parula(1+100,:)]; 
% % % % % plot_legend(2).MarkerFaceColor = plot_legend(2).MarkerEdgeColor;
% % % % % plot_legend(2).MarkerFaceAlpha=1;
% % % % % plot_legend(2).MarkerEdgeAlpha=plot_legend(2).MarkerFaceAlpha;
% % % % % plot_legend(2).SizeData=200.*fig_weig;
% % % % % plot_legend(2).AlphaDataMode='manual';
% % % % % 
% % % % % 
% % % % % 
% % % % % %% normal scatter plot, lt=500 (tstart + 10 + pred_len ~ tstart + 10 + pred_len + 10 -1)
% % % % % for ensi=1:cfg.ensnum
% % % % %     s=scatter(cfg.t_ini_min+10+cfg.pred_len:cfg.t_ini_min+10+cfg.pred_len+10-1, ...
% % % % %         squeeze(HCST_X.data_lt(ensi,end,1+10:1+10+cfg.pred_len-1)),'filled');
% % % % %     s.Marker='p';
% % % % % %     s.MarkerEdgeColor=[255 208 0]/255; % yellow
% % % % %     s.MarkerEdgeColor=tmp.parula(cfg.pred_len*100+100,:); 
% % % % %     s.MarkerFaceColor = s.MarkerEdgeColor;
% % % % %     s.SizeData=100.*fig_weig;
% % % % %     s.MarkerFaceAlpha=1;  %0.8
% % % % %     s.MarkerEdgeAlpha=s.MarkerFaceAlpha;
% % % % %     s.AlphaDataMode='manual';
% % % % % end
% % % % % plot_legend(3)=scatter(cfg.t_ini_min+10+cfg.pred_len:cfg.t_ini_min+10+cfg.pred_len+10-1, ...
% % % % %     squeeze(HCST_X.data_lt_em(end,1+10:1+10+cfg.pred_len-1)),'filled');
% % % % % plot_legend(3).Marker='s';
% % % % % % s.MarkerEdgeColor=[255 140 0]/255; % dark orange
% % % % % plot_legend(3).MarkerEdgeColor=[tmp.parula(cfg.pred_len*100+100,:)]; 
% % % % % plot_legend(3).MarkerFaceColor = plot_legend(3).MarkerEdgeColor;
% % % % % plot_legend(3).MarkerFaceAlpha=1;
% % % % % plot_legend(3).MarkerEdgeAlpha=plot_legend(3).MarkerFaceAlpha;
% % % % % plot_legend(3).SizeData=200.*fig_weig;
% % % % % plot_legend(3).AlphaDataMode='manual';
% % % % % plot_legend(3).AlphaDataMapping='none';
% % % % % 
% % % % % hold off
% % % % % set(gca,'fontsize',15.*fig_weig);
% % % % % grid on;
% % % % % s=text(56,45,['Lorenz Model: ', 'X', ' versus ', '\tau'], 'fontsize', 22.*fig_weig, ...
% % % % %     'FontName', 'Arial', 'FontWeight', 'Bold', 'FontAngle', 'italic');
% % % % % % s=text(56,35,'Lorenz Model:', 'fontsize', 22, ...
% % % % % %     'FontName', 'Sans Serif', 'FontWeight', 'Bold');
% % % % % % s=text(59.9,34.9, '$$ x $$', 'Interpreter', 'latex', 'fontsize', 27, 'FontWeight', 'Bold' );
% % % % % % s=text(60.4,35,'versus', 'fontsize', 22, ...
% % % % % %     'FontName', 'Sans Serif');
% % % % % % s=text(62.0,34.9, '$$ t $$', 'Interpreter', 'latex', 'fontsize', 27, 'FontWeight', 'Bold' );
% % % % % 
% % % % % if ext_flag==1
% % % % %     ylabel('$$ x(\tau) + f(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% % % % % else
% % % % %     ylabel('$$ x(\tau) $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% % % % % end
% % % % % xlabel('$$ \tau $$', 'Interpreter', 'latex', 'fontsize', 25.*fig_weig)
% % % % % 
% % % % % pos=get(gca,'Position');
% % % % % outerpos=get(gca,'OuterPosition');
% % % % % 
% % % % % % legend([plot_legend], {'CTL  ', 'INI[LT=0]  ', ['INI[LT=', num2str(cfg.pred_len), ']  '], 'unINI  '}, ...
% % % % % %     'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 22.*fig_weig);
% % % % % % legend([plot_legend], ...
% % % % % %     {'$$ C  $$', '$$ I[lt=0] $$', ['$$ I[lt=', num2str(cfg.pred_len), ']  $$'], '$$ UI $$'}, ...
% % % % % %     'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 22.*fig_weig, ...
% % % % % %     'interpreter', 'latex');%% ver. 1
% % % % % legend([plot_legend], ...
% % % % %     {'$$ O  \hspace{5mm} $$ ', '$$ X_{I}[l\tau=0] \hspace{5mm} $$', ['$$ X_{I}[l\tau=', num2str(cfg.pred_len), '] \hspace{5mm} $$'], '$$ X_{UI} $$'}, ...
% % % % %     'location', 'northoutside', 'Orientation', 'horizontal', 'Fontsize', 35.*fig_weig, ...
% % % % %     'interpreter', 'latex'); %% fontsize 22*w = default
% % % % % 
% % % % % % print(gcf, 'abc.png', '-dpng');


%% corr coeficient calculation -------------------------------------------------------

% corrcoef(data_comb.X(901:100:3401),  HCST_X.data_lt(1,1,6:31))  % (lt=1, 60~85)
% corrcoef(data_comb.X(900:100:3400),  HCST_X.data_lt_em(500,1:26))  % (lt=500, 60~85)

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




%% additional read for trajectory plot
str.t_ini=num2str(55, '%04i');
for ensnum=1:cfg.ensnum
    str_ens=num2str(ensnum, '%02i');
%         disp('pred:', str.t_ini, 'ens:',str_ens);
    dir.hcst_tmp=[dir.hcst, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
    dir.hcst_lqini_tmp=[dir.hcst_lqini, '/', str.t_ini, '/', 'ens',str_ens, '/', 'out'];
    dir.hcst_unini_tmp=[dir.hcst_unini, '/', '0055', '/', 'ens',str_ens, '/', 'out'];

    for t_pred=1:50
        str.t_hcst=num2str(t_ini+t_pred-1, '%04i');
        fname_hcst=[dir.hcst_tmp, '/', 'out.', str.t_hcst, '.nc'];
        fname_hcst_lqini=[dir.hcst_lqini_tmp, '/', 'out.', str.t_hcst, '.nc'];
        fname_hcst_unini=[dir.hcst_unini_tmp, '/', 'out.', str.t_hcst, '.nc'];
        if ext_flag==0
        %% raw
            HCST_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'X'));
            HCST_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'Y'));
            HCST_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst, 'Z'));

        %% raw_unini
        HCST_unini_X.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'X'));
        HCST_unini_Y.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'Y'));
        HCST_unini_Z.(['ens',str_ens]).(['i',str.t_ini]).(['p',str.t_hcst])=single(ncread(fname_hcst_unini, 'Z')); 
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

% hold off
% % period=[1:100:50000];
% period=[1:1000];
% plot3(data_comb.X(period), data_comb.Y(period), data_comb.Z(period));
% hold on;
% 
% period2=[1:1000];
% plot3(data_comb_hcst.X(period2), data_comb_hcst.Y(period2), data_comb_hcst.Z(period2));



fig_cfg.fig_size = [0,0,13,14.5]; %% paper size (original)
fig_cfg.cb_size = [1, 13, 11, 0.3];
loc_column_first=1;
loc_row_first=7;
fig_h = figure('name', 'fig_lorenz','PaperUnits','inches', ...
        'PaperPosition',fig_cfg.fig_size,'position',fig_cfg.fig_size*get(groot,'ScreenPixelsPerInch')+[200,200,0,0],'visible','on');

%% plot using surf (CTL, OBS)
ax_m=subplot(2,2,1);
fig_cfg.ax_size = [loc_column_first, loc_row_first, 5, 5];
period=1:5000;
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
title1=title('(a) Observation', 'fontsize', 20);

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


%% plot using surf (ensemble mean)
ax_m=subplot(2,2,2);
fig_cfg.ax_size = [loc_column_first + 6, loc_row_first, 5, 5];
period=1:5000;
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
xlim([-20 20])
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















