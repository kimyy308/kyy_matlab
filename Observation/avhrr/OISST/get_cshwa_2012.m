
load ('D:\MEPL\project\SSH\data\OISST_monthly\AVHRR_2012\201201month')

pcolor(monthly_mean');
shading interp;
colormap(jet);
colorbar;


load ('D:\MEPL\project\SSH\data\OISST_monthly\AVHRR_climate(82-2016)\climate02')