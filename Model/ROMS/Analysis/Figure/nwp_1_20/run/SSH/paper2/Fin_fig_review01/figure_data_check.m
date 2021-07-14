close all; clear all; clc;

% % % Data2
% pcolor(GCM_IPSL_L_lon, GCM_IPSL_L_lat, GCM_IPSL_L_topo);
% pcolor(GCM_Nor_lon, GCM_Nor_lat, GCM_Nor_topo);
% pcolor(GCM_MPI_lon, GCM_MPI_lat, GCM_MPI_topo);
% shading flat; 
% colorbar;

% % % Data3
% pcolor(AVHRR_lon, AVHRR_lat, AVHRR_SST);
% for i=1:4
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_SST{i})
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SST{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data4
% % pcolor(CMEMS_lon, CMEMS_lat, CMEMS_Svel.u);
% % pcolor(CMEMS_lon, CMEMS_lat, CMEMS_Svel.v);
% for i=1:4
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_Svel.u{i})
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_Svel.v{i})
% %     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_Svel.u{i})
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_Svel.v{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data5
% pcolor(CMEMS_lon, CMEMS_lat, CMEMS_SSH);
% for i=1:4
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_SSH{i})
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SSH{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end


% % % Data 6
% pcolor(CMEMS_lon, CMEMS_lat, CMEMS_SSH_std);
% for i=1:4
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_SSH_std{i})
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SSH_std{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data 7

% % % Data 8
% for i=1:4
% %     plot(GCM_time,GCM_MSL.rcp45(i,:))
% %     plot(GCM_time,GCM_MSL.rcp85(i,:))
%     plot(RCM_time,RCM_MSL.rcp45(i,:))
%     plot(RCM_time,RCM_MSL.rcp85(i,:))    
%     pause(2)
% end

% % % % Data 9
% for i=1:5
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_SLR{i})
% %     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SLR{i})
%     pcolor(DIFF_lon{i}, DIFF_lat{i}, DIFF_SLR{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data 10
% for i=1:5
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_SLR{i})
% %     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SLR{i})
%     pcolor(DIFF_lon{i}, DIFF_lat{i}, DIFF_SLR{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data 11
% for i=1:2
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_tidal_amp_change{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data 12
% for i=1:2
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_total_SLR{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % % % Data S01 S03
% for i=1:5
% %     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_steric_SLR{i})
% %     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_steric_SLR{i})
% %     pcolor(DIFF_lon{i}, DIFF_lat{i}, DIFF_steric_SLR{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data S02 S04
% for i=1:5
%     pcolor(GCM_lon{i}, GCM_lat{i}, GCM_manometric_SLR{i})
% %     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_manometric_SLR{i})
% %     pcolor(DIFF_lon{i}, DIFF_lat{i}, DIFF_manometric_SLR{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data S05 06 07 08 09
% for i=1:8
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_tidal_amp_change{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data S10
% for i=1:8
%     pcolor(RCM_lon{i}, RCM_lat{i}, RCM_SSH_std_future{i})
%     shading flat; 
%     colorbar;
%     pause(2)
% end

% % % Data S11
for i=1:8
    pcolor(RCM_lon{i}, RCM_lat{i}, RCM_total_SLR{i})
    shading flat; 
    colorbar;
    pause(2)
end