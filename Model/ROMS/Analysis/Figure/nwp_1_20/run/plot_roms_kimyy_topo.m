% %  Updated 27-Apr-2018 by Yong-Yub Kim
% %  Updated 21-May-2018 by Yong-Yub Kim
% %  Updated 08-Jun-2018 by Yong-Yub Kim



if (windows ==1)
    % % for windows
    bathydir='Bathy\'; %% SNU_desktop
    % % set colorbar parameter
    load C:\Users\KYY\Dropbox\source\matlab\Common\Figure\jet_mod  % % set colormap (jet_modified)
elseif (linux==1)
    % % for linux
    bathydir='Bathy/'; %% Linux
end

status=plot_ROMS_bathy([figdir, bathydir, 'bathy_nwp'],inputdir, [115 164 15 52], [0 5000], [0], testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_taiwan'],inputdir, [117.5 123 22 27], [0 130], [30 50 80], testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_korea'],inputdir, [128 132 32.5 36], [0 200],0:20:200, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_tsugaru'],inputdir, [139 142 40 43], [0 200], 0:20:200, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_soya'],inputdir, [141 143 45 47], [0 200],0:20:200, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_EJS'],inputdir, [127 144 33 52], [0 200],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_EJS_4000'],inputdir, [127 144 33 52], [0 4000],0:500:4000, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_YS'],inputdir, [116 127 33 42], [0 2000],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_EKWC'],inputdir, [129 132 34 39], [0 300],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_YS'],inputdir, [116 128 32 42], [0 100],0:20:100, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_ECS'],inputdir, [116 128 25 42], [0 100],0:20:100, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_kuro_5000'],inputdir, [132 143 27 35], [0 5000],0:500:5000, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_kuro_2000'],inputdir, [132 143 27 35], [0 2000],0:200:2000, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_luzon_5000'],inputdir, [117 124 16 23], [0 5000],0:500:5000, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_luzon_2000'],inputdir, [117 124 16 23], [0 2000],0:200:2000, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_kuro_5000_nc'],inputdir, [132 143 27 35], [0 5000],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_kuro_2000_nc'],inputdir, [132 143 27 35], [0 2000],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_luzon_5000_nc'],inputdir, [117 124 16 23], [0 5000],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_luzon_2000_nc'],inputdir, [117 124 16 23], [0 2000],0, testname);
status=plot_ROMS_bathy([figdir, bathydir, 'bathy_kuro_2000'],inputdir, [132 143 27 35], [0 2000],0:200:2000, testname);
close all;

