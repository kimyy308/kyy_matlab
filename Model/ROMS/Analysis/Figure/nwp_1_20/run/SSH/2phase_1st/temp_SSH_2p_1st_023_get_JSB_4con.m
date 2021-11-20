close all; clc; clear all;


testname='test2102';
for year=1985:2014
    dirs.rootdir=['D:\Data\Model\ROMS\nwp_1_20\JSB\',testname '\harmonic_analysis'];

    tmp.filename=[dirs.rootdir, filesep, 'test2102_harmonic_analysis_zeta_JSB_', num2str(year), '.mat'];
    load(tmp.filename, 'tide_info', 'Model', 'cut_lon_rho', 'cut_lat_rho')
    Model_all(year-1984).Model=Model;
end

testname='test2107';
for year=2015:2050
    dirs.rootdir=['D:\Data\Model\ROMS\nwp_1_20\JSB\',testname '\harmonic_analysis'];
    tmp.filename=[dirs.rootdir, filesep, 'test2107_harmonic_analysis_zeta_JSB_', num2str(year), '.mat'];
    if (exist(tmp.filename)==2)
        load(tmp.filename, 'tide_info', 'Model', 'cut_lon_rho', 'cut_lat_rho')
        Model_all(year-1984).Model=Model;
    else
        Model_all(year-1984).Model=NaN;
    end
end


Model_all(47).Model=Model_all(46).Model;
Model_all(48).Model=Model_all(46).Model;

for i=1:66
    tmpval(i)=mean(mean(Model_all(i).Model.amp_M2(:,:), 'omitnan'), 'omitnan');
end
plot(1985:2050,tmpval)
% ylim([1.2 1.4])
% cut_lon_rho(13,6)
% cut_lat_rho(13,6)


% 
% 
% for i=1:23
%     for j=1:23
%         tmpval2(i,j)=Model_all(38).Model.amp_M2(i,j);
%     end
% end
% pcolor(cut_lon_rho', cut_lat_rho', tmpval2'); colorbar;

load('D:\Data\Model\ROMS\nwp_1_20\backup_surf\test2102\run\zeta\test2102_AKP4_RCM_ssh_1985_2014.mat')

[indw, inde, inds, indn] = Func_0012_findind_Y(1/20, ...
    [cut_lon_rho(1,1),cut_lon_rho(end,end),cut_lat_rho(1,1),cut_lat_rho(end,end)], ...
    RCM_grid.lon_rho, RCM_grid.lat_rho, 1)

Model_ssh=RCM_data.all(indw:inde, inds:indn,:);
figure; 
plot(squeeze(mean(mean(Model_ssh,1,'omitnan'),2,'omitnan')))