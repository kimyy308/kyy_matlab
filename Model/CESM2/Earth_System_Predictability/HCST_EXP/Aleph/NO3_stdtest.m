% clear all; close all; clc;

addpath('/mnt/lustre/proj/kimyy/Dropbox/source/matlab/function')
years=[1960 1980 2000 2020];
varname='TEMP'; %TEMP SALT PD PO4 Fe WVEL photoC_TOT PAR_avg
xpoint=180;
ypoint=20;
levs=[1, 5, 10, 15];


head = ['/mnt/lustre/proj/earth.system.predictability/ASSM_EXP/archive_yearly_transfer/ocn/', varname];

a=dir([head, '/f09*']);

for zind = 1:length(levs)
    lev=levs(zind);
for yind=1:length(years)
    year=years(yind);
    yearstr=num2str(year);

for i=1:20
    expn=a(i).name;
    fname=[head,'/',expn,'/',varname,'_',expn,'.pop.h.', yearstr, '.nc'];
    grids.tlong=ncread(fname, 'TLONG');
    grids.tlat=ncread(fname, 'TLAT');

    varinfo=ncinfo(fname, varname);
    dim=length(varinfo.Dimensions);

    [grids.id_w, grids.id_e, grids.id_s, grids.id_n] = Func_0012_findind_Y(3, [xpoint, ypoint], ...
                    grids.tlong, ...
                    grids.tlat, 'CESM2'); % find valid lon, lat index near station
    
    if dim==3
        data(i)=ncread(fname, varname, [grids.id_w grids.id_s 1], [1 1 1]);
    elseif dim==4
        try
            grids.z_t=ncread(fname, 'z_t')./100;
        catch
            try
                grids.z_t=ncread(fname, 'z_t_150m')./100;
            catch
                grids.z_t=ncread(fname, 'z_w_top')./100;
            end
        end
        data(i)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);
    end
end
%     subplot(length(years)+1,length(levs),(zind-1)*length(years)+yind+1)
    subplot(length(years)+1,length(levs),length(levs)*yind+zind)

    plot(data)
    title(['memb ', num2str(year), 'Y, ', num2str(xpoint), 'E, ', num2str(ypoint), 'N, ', num2str(grids.z_t(lev)),'m' ])
    grid minor
end



for i=1960:2020
    fname = [ head, '/',  'ens_all/', varname, '_ensmean_', num2str(i), '.nc'];
    data_y(i-1959)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);
    fname = [ head, '/',  'ens_all/', varname, '_ensmax_', num2str(i), '.nc'];
    data_y_max(i-1959)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);
    fname = [ head, '/',  'ens_all/', varname, '_ensmin_', num2str(i), '.nc'];
    data_y_min(i-1959)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);

    fname = [ head, '/',  'ens_all/', varname, '_projdv7.3_ba_ensmean_', num2str(i), '.nc'];
    data_y_projd(i-1959)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);
    fname = [ head, '/',  'ens_all/', varname, '_en4.2_ba_ensmean_', num2str(i), '.nc'];
    data_y_en4(i-1959)=ncread(fname, varname, [grids.id_w grids.id_s lev 1], [1 1 1 1]);
end
% subplot(length(years)+1,length(levs),(zind-1)*length(years)+1)
subplot(length(years)+1,length(levs),zind)

plot(1960:2020, data_y, 'k');
hold on
plot(1960:2020, data_y_projd, 'r');
plot(1960:2020, data_y_en4, 'b');
plot(1960:2020, data_y_max, 'k--' );
plot(1960:2020, data_y_min, 'k--' );
hold off


title(['ensm t-s ', num2str(xpoint), 'E, ', num2str(ypoint), 'N, ', num2str(grids.z_t(lev)),'m' ])
grid minor


end

% kurtosis(squeeze(no3(200,255,:)))




