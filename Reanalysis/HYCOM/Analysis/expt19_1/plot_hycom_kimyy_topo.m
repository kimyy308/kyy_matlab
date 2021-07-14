
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_nwp'],[inputdir, 'input\'], [115 164 15 52], [0 5000], [0], testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_taiwan'],[inputdir, 'input\'], [117.5 123 22 27], [0 130], [30 50 80], testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_korea'],[inputdir, 'input\'], [128 132 32.5 36], [0 200],0:20:200, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_tsugaru'],[inputdir, 'input\'], [139 142 40 43], [0 200], 0:20:200, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_soya'],[inputdir, 'input\'], [141 143 45 47], [0 200],0:20:200, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_EJS'],[inputdir, 'input\'], [127 144 33 52], [0 2000],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_YS'],[inputdir, 'input\'], [116 128 32 42], [0 100],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_5000'],[inputdir, 'input\'], [132 143 27 35], [0 5000],0:500:5000, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_2000'],[inputdir, 'input\'], [132 143 27 35], [0 2000],0:200:2000, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_500'],[inputdir, 'input\'], [132 143 27 35], [0 500],0:50:500, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_200'],[inputdir, 'input\'], [132 143 27 35], [0 200],0:20:200, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_5000'],[inputdir, 'input\'], [117 124 16 23], [0 5000],0:500:5000, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_2000'],[inputdir, 'input\'], [117 124 16 23], [0 2000],0:200:2000, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_500'],[inputdir, 'input\'], [117 124 16 23], [0 500],0:50:500, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_200'],[inputdir, 'input\'], [117 124 16 23], [0 200],0:50:200, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_5000_nc'],[inputdir, 'input\'], [132 143 27 35], [0 5000],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_2000_nc'],[inputdir, 'input\'], [132 143 27 35], [0 2000],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_500_nc'],[inputdir, 'input\'], [132 143 27 35], [0 500],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_200_nc'],[inputdir, 'input\'], [132 143 27 35], [0 200],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_5000_nc'],[inputdir, 'input\'], [117 124 16 23], [0 5000],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_2000_nc'],[inputdir, 'input\'], [117 124 16 23], [0 2000],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_500_nc'],[inputdir, 'input\'], [117 124 16 23], [0 500],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_luzon_200_nc'],[inputdir, 'input\'], [117 124 16 23], [0 200],0, testname);
status=plot_HYCOM_bathy([figdir, 'Bathy\', 'bathy_kuro_2000'],[inputdir, 'input\'], [132 143 27 35], [0 2000],0:200:2000, testname);

% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_kuro_5000'],[inputdir, 'input\'], [132 143 27 35], [0 5000],0:500:5000, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_kuro_2000'],[inputdir, 'input\'], [132 143 27 35], [0 2000],0:200:2000, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_kuro_500'],[inputdir, 'input\'], [132 143 27 35], [0 500],0:50:500, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_kuro_200'],[inputdir, 'input\'], [132 143 27 35], [0 200],0:50:200, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_luzon_5000'],[inputdir, 'input\'], [117 124 16 23], [0 5000],0:500:5000, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_luzon_2000'],[inputdir, 'input\'], [117 124 16 23], [0 2000],0:200:2000, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_luzon_500'],[inputdir, 'input\'], [117 124 16 23], [0 500],0:50:500, testname);
% status=plot_etopo1([figdir, 'Bathy\', 'bathy_etopo1_luzon_200'],[inputdir, 'input\'], [117 124 16 23], [0 200],0:50:200, testname);
close all;

