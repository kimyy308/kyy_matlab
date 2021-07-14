clc;clear all;close all;

%Setting variable and path
result_path = '/home/ygkim/ipcc/data/zos'; %result path
root_path ='/data/ygkim/ipcc/'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%load ensemble data
load(strcat(result_path,'/processing_result_v3_new.mat'));

%find model name
for i = 1:length(root_path2)
    model_name(i) = {root_path2(i).name};
end

model_name_zos = model_name(zos_idx);

%----index plot(1D)---------%
[sort_idx, sort_num] = sort(total_RMSE);
sort_num_gl_zos = sort_num;
I = total_RMSE;
mean_I = mean(I);
I_std = se;
x_label = model_name_zos;

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',x_label(sort_num),...
    'XTick',[1:41],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');


%bar plot(blue to red)
bar(1:10,I(sort_num(1:10)),'BaseValue',mean_I,'FaceColor',rgb('indigo'));
bar(11:20,I(sort_num(11:20)),'BaseValue',mean_I,'FaceColor',rgb('royalblue'));
bar(21:30,I(sort_num(21:30)),'BaseValue',mean_I,'FaceColor',rgb('ForestGreen'));
bar(31:41,I(sort_num(31:41)),'BaseValue',mean_I,'FaceColor',rgb('darkorange'));
%bar(42:48,I(sort_num(42:48)),'BaseValue',1,'FaceColor',rgb('crimson'));

errorbar(I(sort_num),I_std(sort_num)*1.96,'LineStyle','none','Color','k','LineWidth',1);
set(gca,'FontSize',16,'xlim',[0,42]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('RMSE of sea level anomaly(m)','fontsize',25);
tt = strcat('RMSE of Sea level anomaly (Global)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat(tt,'.fig'),'fig');
close all;


%----index plot(1D,NWP)---------%
[sort_idx, sort_num] = sort(ocean_RMSE);
sort_num_nwp_zos = sort_num;
I = ocean_RMSE;
mean_I = mean(I);
I_std = se;
x_label = model_name_zos;

%----To plot bar color sync -------%
color_data = {};
for i = 1: length(sort_num)
    idx_nwp = sort_num_nwp_zos(i);
    idx_diff = find(idx_nwp == sort_num_gl_zos);
    if idx_diff < 11
        color_data{i} = 'indigo';
    elseif idx_diff>=11 && idx_diff<21
        color_data{i} = 'royalblue';
    elseif idx_diff>=21 && idx_diff<31
        color_data{i} = 'ForestGreen';
    elseif idx_diff>=31 && idx_diff<=41
        color_data{i} = 'darkorange';
    else
        color_data{i} = 'crimson';
    end
end



figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',x_label(sort_num),...
    'XTick',[1:41],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');


%bar plot(blue to red)
for i = 1:length(sort_num)
    bar(i,I(sort_num(i)),'BaseValue',mean_I,'FaceColor',rgb(color_data{i}));
end
errorbar(I(sort_num),I_std(sort_num)*1.96,'LineStyle','none','Color','k','LineWidth',1);
set(gca,'FontSize',16,'xlim',[0,42]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('RMSE of sea level anomaly(m)','fontsize',25);
tt = strcat('RMSE of Sea level anomaly (NWP)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat(tt,'.fig'),'fig');
close all;

cor_ssh = corr(total_RMSE,ocean_RMSE);


%% Heat contents
%clc;clear all;close all;

%Setting variable and path
result_path = '/home/ygkim/ipcc/data/zos'; %result path
root_path ='/data/ygkim/ipcc/'; %ipcc models root path
root_path2 = dir(strcat(root_path,'*')); 

%load ensemble data
load(strcat(result_path,'/processing_result_heatcontent_NWP.mat'));


%find model name
for i = 1:length(root_path2)
    model_name(i) = {root_path2(i).name};
end



%----index plot(1D)---------%
[sort_idx, sort_num] = sort(total_RMSE);
sort_num_gl_heat = sort_num;
I = total_RMSE;
mean_I = mean(I);
I_std = se;
x_label = model_name(4:end);

figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',x_label(sort_num),...
    'XTick',[1:48],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');


%bar plot(blue to red)
bar(1:10,I(sort_num(1:10)),'BaseValue',mean_I,'FaceColor',rgb('indigo'));
bar(11:20,I(sort_num(11:20)),'BaseValue',mean_I,'FaceColor',rgb('royalblue'));
bar(21:30,I(sort_num(21:30)),'BaseValue',mean_I,'FaceColor',rgb('ForestGreen'));
bar(31:41,I(sort_num(31:41)),'BaseValue',mean_I,'FaceColor',rgb('darkorange'));
bar(42:48,I(sort_num(42:48)),'BaseValue',mean_I,'FaceColor',rgb('crimson'));

errorbar(I(sort_num),I_std(sort_num)*1.96,'LineStyle','none','Color','k','LineWidth',1);
set(gca,'FontSize',16,'xlim',[0,49]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('RMSE of heat content(J/m^3)','fontsize',25);
tt = strcat('RMSE of upper ocean heat content (Global)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat(tt,'.fig'),'fig');
close all;

%----index plot(1D,NWP)---------%
[sort_idx, sort_num] = sort(ocean_RMSE);
sort_num_nwp_heat = sort_num;
I = ocean_RMSE;
mean_I = mean(I);
I_std = se;
x_label = model_name(4:end);


%----To plot bar color sync -------%
color_data = {};
for i = 1: length(sort_num)
    idx_nwp = sort_num_nwp_heat(i);
    idx_diff = find(idx_nwp == sort_num_gl_heat);
    if idx_diff < 11
        color_data{i} = 'indigo';
    elseif idx_diff>=11 && idx_diff<21
        color_data{i} = 'royalblue';
    elseif idx_diff>=21 && idx_diff<31
        color_data{i} = 'ForestGreen';
    elseif idx_diff>=31 && idx_diff<=41
        color_data{i} = 'darkorange';
    else
        color_data{i} = 'crimson';
    end
end


figure1 = figure;
set(figure1,'OuterPosition',[0,0,1000,400]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 3 1]*6);

% Create axes
axes1 = axes('Parent',figure1,...
    'XTickLabel',x_label(sort_num),...
    'XTick',[1:48],'XTickLabelRotation',45);
box(axes1,'on');
hold(axes1,'all');


%bar plot(blue to red)
for i =1:length(sort_num)
    bar(i,I(sort_num(i)),'BaseValue',mean_I,'FaceColor',rgb(color_data(i)));
end


errorbar(I(sort_num),I_std(sort_num)*1.96,'LineStyle','none','Color','k','LineWidth',1);
set(gca,'FontSize',16,'xlim',[0,49]);
% Create ylabel
xlabel('Model Name','fontsize',25);
ylabel('RMSE of heat content(J/m^3)','fontsize',25);
tt = strcat('RMSE of upper ocean heat content (NWP)');
title(tt,'fontsize',30,'fontweight','bold');
saveas(gca,strcat(tt,'.fig'),'fig');
close all;

cor_heat = corr(total_RMSE,ocean_RMSE);

