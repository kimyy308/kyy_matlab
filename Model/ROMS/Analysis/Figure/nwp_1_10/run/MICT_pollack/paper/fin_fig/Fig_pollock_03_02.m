close all; clear all;  clc;   

opts = spreadsheetImportOptions("NumVariables", 40);

opts.Sheet = "데이터";
opts.DataRange = "B1:AO22";

opts.VariableNames = ["VarName2", "VarName3", "VarName4", "VarName5", "VarName6", "VarName7", "VarName8", "VarName9", "VarName10", "VarName11", "VarName12", "VarName13", "VarName14", "VarName15", "VarName16", "VarName17", "VarName18", "VarName19", "VarName20", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "p"];
opts.VariableTypes = ["string", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

opts = setvaropts(opts, "VarName2", "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["VarName2", "VarName3", "VarName4", "VarName11", "VarName21", "VarName22", "VarName23", "VarName24", "VarName25", "VarName26", "VarName27", "VarName28", "VarName29", "VarName30", "VarName31", "VarName32", "VarName33", "VarName34", "VarName35", "VarName36", "VarName37", "VarName38", "VarName39", "VarName40", "p"], "EmptyFieldRule", "auto");

catch_data = readtable("Z:\내 드라이브\research\Ph_D_course\Particle tracking of the walleye pollock\catch\어업별_품종별_통계_20210219_명태_노가리.xlsx", opts, "UseExcel", false);

clear opts

year = table2array(catch_data(1,4:39));

pollock_juvenile.region = {'total', 'busan', 'kangwon', 's_chungchung', 'n_jeonra', 's_jeonra', 'n_kyeongsang', 's_kyeongsang', 'jeju'};
for regioni=1:length(pollock_juvenile.region)
    pollock_juvenile.(pollock_juvenile.region{regioni}) = table2array(catch_data(regioni+2,4:39));
    pollock_juvenile.(pollock_juvenile.region{regioni})(isnan(pollock_juvenile.(pollock_juvenile.region{regioni})))=0;
end
pollock_adult.region = {'total', 'busan', 'incheon', 'ulsan', 'kangwon', 's_chungchung', 'n_jeonra', 's_jeonra', 'n_kyeongsang', 's_kyeongsang', 'jeju'};
for regioni=1:length(pollock_adult.region)
    pollock_adult.(pollock_adult.region{regioni}) = table2array(catch_data(regioni+11,4:39));
    pollock_adult.(pollock_adult.region{regioni})(isnan(pollock_adult.(pollock_adult.region{regioni})))=0;
end

pollock_juvenile.eastsea = pollock_juvenile.kangwon + pollock_juvenile.n_kyeongsang + pollock_juvenile.s_kyeongsang + pollock_juvenile.busan;
pollock_juvenile.eastsea(16:end)=NaN;
pollock_adult.eastsea = pollock_adult.kangwon + pollock_adult.n_kyeongsang + pollock_adult.s_kyeongsang + pollock_adult.busan;

pollock_adult.eastsea_half1(1:5)=mean(pollock_adult.eastsea(1:5));
pollock_adult.eastsea_half2(1:5)=mean(pollock_adult.eastsea(6:10));
pollock_juvenile.eastsea_half1(1:5)=mean(pollock_juvenile.eastsea(1:5));
pollock_juvenile.eastsea_half2(1:5)=mean(pollock_juvenile.eastsea(6:10));

jpgname = 'D:\research\Ph_D_course\2021_Particle tracking of the walleye pollock\figure\paper\Fig03_02.jpg';
startind=1;
% endind=length(year);
endind=10;



pollock_adult.plot=plot(year(startind:endind), pollock_adult.eastsea(startind:endind));
hold on
pollock_adult.halfplot1=plot(year(startind:5), pollock_adult.eastsea_half1(startind:5));
pollock_adult.halfplot2=plot(year(6:10), pollock_adult.eastsea_half2(startind:5));

pollock_juvenile.plot=plot(year(startind:endind), pollock_juvenile.eastsea(startind:endind));
pollock_juvenile.halfplot1=plot(year(startind:5), pollock_juvenile.eastsea_half1(startind:5));
pollock_juvenile.halfplot2=plot(year(6:10), pollock_juvenile.eastsea_half2(startind:5));
hold off
xlabel('Year')
ylabel('Catch (metric tons)')
set(pollock_adult.plot,'LineStyle','--', 'Marker', 'o','LineWidth',1, 'color', 'b');
set(pollock_adult.halfplot1,'LineWidth',3, 'color', 'b');
set(pollock_adult.halfplot2,'LineWidth',3, 'color', 'b');
set(pollock_juvenile.plot,'LineStyle','--','Marker','^','LineWidth',1, 'color', 'r');
set(pollock_juvenile.halfplot1,'LineWidth',3, 'color', 'r');
set(pollock_juvenile.halfplot2,'LineWidth',3, 'color', 'r');
lgd=legend([pollock_adult.plot pollock_juvenile.plot pollock_adult.halfplot1 pollock_juvenile.halfplot1], 'Adult', 'Juvenile', 'Adult (5-year mean)', 'Juvenile (5-year mean)');
axis tight
set(gca,'FontSize',22);
set(gcf,'PaperPosition', [0 0 36 12]) 
saveas(gcf,jpgname,'tif');
close all