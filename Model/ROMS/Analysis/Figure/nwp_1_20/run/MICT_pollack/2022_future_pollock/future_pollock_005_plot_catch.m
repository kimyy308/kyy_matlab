% %  Updated 27-Apr-2021 by Yong-Yub Kim, 


close all; clear all;  clc;   


% % opts = spreadsheetImportOptions("NumVariables", 54);
% % 
% % opts.Sheet = "데이터";
% % opts.DataRange = "B1:BC26";
% % 
% % opts.VariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55"];
% % opts.VariableTypes = ["string", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical"];
% % 
% % opts = setvaropts(opts, "VarName2", "WhitespaceRule", "preserve");
% % opts = setvaropts(opts, ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName24", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54", "VarName55"], "EmptyFieldRule", "auto");
% % 
% % catch_data = readtable("D:\research\Ph_D_course\2022_pollock_future\catch\어업별_품종별_통계_20220425204305.xlsx", opts, "UseExcel", false);
% % 
% % clear opts

opts = spreadsheetImportOptions("NumVariables", 53);
opts.Sheet = "데이터";
opts.DataRange = "B1:BB22";
opts.VariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54"];
% opts.VariableTypes = ["string", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical", "categorical"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts = setvaropts(opts, "VarName2", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName24", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "VarName41", "VarName42", "VarName43", "VarName44", "VarName45", "VarName46", "VarName47", "VarName48", "VarName49", "VarName50", "VarName51", "VarName52", "VarName53", "VarName54"], "EmptyFieldRule", "auto");
catch_data = readtable("D:\Research\Ph_D_course\2022_pollock_future\catch\어업별_품종별_통계_20220425204305.xlsx", opts, "UseExcel", false);
clear opts





data_range=17:53;
year = table2array(catch_data(1,data_range));

pollock_juvenile.region = {'total', 'busan', 'kangwon', 's_chungchung', 'n_jeonra', 's_jeonra', 'n_kyeongsang', 's_kyeongsang', 'jeju'};
for regioni=1:length(pollock_juvenile.region)
    pollock_juvenile.(pollock_juvenile.region{regioni}) = table2array(catch_data(regioni+2,data_range));
    pollock_juvenile.(pollock_juvenile.region{regioni})(isnan(pollock_juvenile.(pollock_juvenile.region{regioni})))=0;
end
pollock_adult.region = {'total', 'busan', 'incheon', 'ulsan', 'kangwon', 's_chungchung', 'n_jeonra', 's_jeonra', 'n_kyeongsang', 's_kyeongsang', 'jeju'};
for regioni=1:length(pollock_adult.region)
    pollock_adult.(pollock_adult.region{regioni}) = table2array(catch_data(regioni+11,data_range));
    pollock_adult.(pollock_adult.region{regioni})(isnan(pollock_adult.(pollock_adult.region{regioni})))=0;
end

pollock_juvenile.eastsea = pollock_juvenile.kangwon + pollock_juvenile.n_kyeongsang + pollock_juvenile.s_kyeongsang + pollock_juvenile.busan;
% pollock_juvenile.eastsea(16:end)=NaN;
pollock_adult.eastsea = pollock_adult.kangwon + pollock_adult.n_kyeongsang + pollock_adult.s_kyeongsang + pollock_adult.busan;


jpgname = 'D:\Research\Ph_D_course\2022_pollock_future\figure\catch\Fig_catch_all.jpg';
% startind=1;
% endind=length(year);
% % endind=23;
% pollock_adult.plot=plot(year(startind:endind), pollock_adult.eastsea(startind:endind));
% hold on
% pollock_juvenile.plot=plot(year(startind:endind), pollock_juvenile.eastsea(startind:endind));
% hold off

pollock_adult.eastsea(end+1:end+2)=0;
pollock_juvenile.eastsea(end+1:end+2)=NaN;
pollock_juvenile.eastsea(16:end)=NaN;
year(end+1)=2020;
year(end+1)=2021;
pollock_adult.plot=plot(year, pollock_adult.eastsea);
hold on
pollock_juvenile.plot=plot(year, pollock_juvenile.eastsea);
hold off

xlabel('Year')
ylabel('Catch (metric tons)')
set(pollock_adult.plot,'Marker', 'o','LineWidth',2);
set(pollock_juvenile.plot,'Marker','^','LineWidth',2);
lgd=legend([pollock_adult.plot pollock_juvenile.plot], 'Adult', 'Juvenile');
axis tight
set(gca,'FontSize',22);
set(gcf,'PaperPosition', [0 0 36 12]) 
saveas(gcf,jpgname,'tif');
close all

[a, pval]=corrcoef(pollock_juvenile.eastsea(9:18), pollock_adult.eastsea(9:18))
[a, pval]=corrcoef(pollock_juvenile.eastsea(9:18), pollock_adult.eastsea(10:19))
[a, pval]=corrcoef(pollock_juvenile.eastsea(10:19), pollock_adult.eastsea(9:18))