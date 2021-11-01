function [ssh_correction_for_fig] = Func_0017_SSH_correction_for_CMIP6_RMSE(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=Func_0017_SSH_correction_for_RMSE(testname);
%
% Getting the SSH correction value for RMSE (reference : 1995-2014 average)
%
%  input:
%  testname             Name of the dataset (string. RCM, GCM, Observation, etc.)
%
%  output:
%  ssh_correction_for_fig          correction value (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    07-Oct-2021 by Yong-Yub Kim
%  Updated    12-Oct-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
    case('CMEMS')
        ssh_correction_for_fig=0.5843; % 1995-2014 AKP4 mean       
    case{'test2102', 'test2107'}
        ssh_correction_for_fig=-0.3279; % 1995-2014 AKP4 mean
    case{'test2103', 'test2108'}
        ssh_correction_for_fig=-0.2331; % 1995-2014 AKP4 mean        
    case{'test2104', 'test2109'}
        ssh_correction_for_fig=-0.2493; % 1995-2014 AKP4 mean                
    case{'test2105', 'test2110'}
        ssh_correction_for_fig=-0.2783; % 1995-2014 AKP4 mean                
    case{'test2106', 'test2111'}
        ssh_correction_for_fig=-0.1649; % 1995-2014 AKP4 mean
    case{'RCM_ENSw_historical'}
        ssh_correction_for_fig=-0.2507; %mean of test2102~2106
        
    case{'test2117', 'test2127'}
        ssh_correction_for_fig=-0.3295; % 1995-2014 AKP4 mean
    case{'test2118', 'test2128'}
        ssh_correction_for_fig=-0.2384; % 1995-2014 AKP4 mean   
    case{'test2119', 'test2129'}
        ssh_correction_for_fig=-0.2835; % 1995-2014 AKP4 mean  
    case{'test2120', 'test2130'}
        ssh_correction_for_fig=-0.2898; % 1995-2014 AKP4 mean  
    case{'test2121', 'test2131'}
        ssh_correction_for_fig=-0.1346; % 1995-2014 AKP4 mean  
    case{'RCM_ENSg_historical'}
        ssh_correction_for_fig=-0.2552; %mean of test2117~2121
        
        
    case('CNRM-ESM2-1')
        ssh_correction_for_fig=-5.6232; % 1995-2014 AKP4 mean           
    case('EC-Earth3-Veg')
        ssh_correction_for_fig=-5.4209; % 1995-2014 AKP4 mean           
    case('ACCESS-CM2')
        ssh_correction_for_fig=0.8417; % 1995-2014 AKP4 mean                   
    case('CNRM-CM6-1-HR')
        ssh_correction_for_fig=-5.6777; % 1995-2014 AKP4 mean                           
    case('CMCC-ESM2')
        ssh_correction_for_fig=-5.6199; % 1995-2014 AKP4 mean
        
%     case('CMEMS')
%         ssh_correction_for_fig=0.6222; % 2014 AKP4 mean       
%     case{'test2102'}
%         ssh_correction_for_fig=-0.4135; % 2014 AKP4 mean
%     case{'test2103'}
%         ssh_correction_for_fig=-0.3453; % 2014 AKP4 mean        
%     case{'test2104'}
%         ssh_correction_for_fig=-0.3240; % 2014 AKP4 mean                
%     case{'test2105'}
%         ssh_correction_for_fig=-0.3655; % 2014 AKP4 mean                
%     case{'test2106'}
%         ssh_correction_for_fig=-0.2782; % 2014 AKP4 mean      
% 
%     case('CNRM-ESM2-1')
%         ssh_correction_for_fig=-5.6007; % 2014 AKP4 mean           
%     case('EC-Earth3-Veg')
%         ssh_correction_for_fig=-5.3777; % 2014 AKP4 mean           
%     case('ACCESS-CM2')
%         ssh_correction_for_fig=0.8866; % 2014 AKP4 mean                   
%     case('CNRM-CM6-1-HR')
%         ssh_correction_for_fig=-5.6533; % 2014 AKP4 mean                           
%     case('CMCC-ESM2')
%         ssh_correction_for_fig=-5.5835; % 2014 AKP4 mean    

%     case('CMEMS')
%         ssh_correction_for_fig=0.6163; % 2015 AKP4 mean       
%     case{'test2107'}
%         ssh_correction_for_fig=-0.4179; % 2015 AKP4 mean
%     case{'test2108'}
%         ssh_correction_for_fig=-0.3660; % 2015 AKP4 mean        
%     case{'test2109'}
%         ssh_correction_for_fig=-0.3160; % 2015 AKP4 mean                
%     case{'test2110'}
%         ssh_correction_for_fig=-0.2749; % 2015 AKP4 mean      (-0.3655; % 2014 AKP4 mean)           
%     case{'test2111'}
%         ssh_correction_for_fig=-0.3588; % 2015 AKP4 mean      (-0.2782; % 2014 AKP4 mean)

%     case('CNRM-ESM2-1')
%         ssh_correction_for_fig=-5.5736; % 2015 AKP4 mean           
%     case('EC-Earth3-Veg')
%         ssh_correction_for_fig=-5.3806; % 2015 AKP4 mean           
%     case('ACCESS-CM2')
%         ssh_correction_for_fig=0.8944; % 2015 AKP4 mean                   
%     case('CNRM-CM6-1-HR')
%         ssh_correction_for_fig=-5.6077; % 2015 AKP4 mean                           
%     case('CMCC-ESM2')
%         ssh_correction_for_fig=-5.5995; % 2015 AKP4 mean    

%     case{'test2107'}
%         ssh_correction_for_fig=-0.3316; % 2021 AKP4 mean
%     case{'test2108'}
%         ssh_correction_for_fig=-0.3475; % 2021 AKP4 mean        
%     case{'test2109'}
%         ssh_correction_for_fig=-0.2991; % 2021 AKP4 mean                
%     case{'test2110'}
%         ssh_correction_for_fig=-0.3485; % 2021 AKP4 mean      (-0.3655; % 2014 AKP4 mean)           
%     case{'test2111'}
%         ssh_correction_for_fig=-0.2894; % 2021 AKP4 mean      (-0.2782; % 2014 AKP4 mean)
%     
%     case('CNRM-ESM2-1')
%         ssh_correction_for_fig=-5.5476; % 2021 AKP4 mean           
%     case('EC-Earth3-Veg')
%         ssh_correction_for_fig=-5.3586; % 2021 AKP4 mean           
%     case('ACCESS-CM2')
%         ssh_correction_for_fig=0.9191; % 2021 AKP4 mean                   
%     case('CNRM-CM6-1-HR')
%         ssh_correction_for_fig=-5.6258; % 2021 AKP4 mean                           
%     case('CMCC-ESM2')
%         ssh_correction_for_fig=-5.5169; % 2021 AKP4 mean    
end

% % % % % get area_weighted_mean

% % % testnames={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
% % % for i=1:length(testnames)
% % %     testname=testnames{i}
% % %     load(['D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\', testname, '_AKP4_RCM_ssh_1993_2014.mat'])
% % %     for t=1:size(RCM_data.yearly_mean,3)
% % %         [ts_yearly_RCM(i,t), error_status] = Func_0011_get_area_weighted_mean(RCM_data.yearly_mean(:,:,t), RCM_grid.lon_rho, RCM_grid.lat_rho);
% % %     end
% % % end
% % % mean(ts_yearly_RCM(5,22))
% % % 
% % % testnames={'test2102', 'test2103', 'test2104', 'test2105', 'test2106'}
% % % for i=1:length(testnames)
% % %     testname=testnames{i}
% % %     load(['D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\', testname, '_AKP4_GCM_ssh_1993_2014.mat'])
% % %     for t=1:size(GCM_data.yearly_mean,3)
% % %         [ts_yearly_GCM(i,t), error_status] = Func_0011_get_area_weighted_mean(GCM_data.yearly_mean(:,:,t), GCM_grid.lon, GCM_grid.lat);
% % %     end
% % % end
% % % mean(ts_yearly_GCM(5,22))
% % % 
% % % 
% % % load(['D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\', testname, '_AKP4_CMEMS_ssh_1993_2014.mat'])
% % % for t=1:size(GCM_data.yearly_mean,3)
% % %     [ts_yearly_CMEMS(t), error_status] = Func_0011_get_area_weighted_mean(CMEMS_data.adt_yearly_mean(:,:,t), CMEMS_grid.lon2, CMEMS_grid.lat2);
% % % end


% % % % %  SSP 585
% % testnames={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'}
% % for i=1:length(testnames)
% %     testname=testnames{i}
% %     load(['D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\', testname, '_AKP4_RCM_ssh_2021_2050.mat'])
% %     for t=1:size(RCM_data.yearly_mean,3)
% %         [ts_yearly_RCM(i,t), error_status] = Func_0011_get_area_weighted_mean(RCM_data.yearly_mean(:,:,t), RCM_grid.lon_rho, RCM_grid.lat_rho);
% %     end
% % end
% % mean(ts_yearly_RCM(5,1))
% % 
% % testnames={'test2107', 'test2108', 'test2109', 'test2110', 'test2111'}
% % for i=1:length(testnames)
% %     testname=testnames{i}
% %     load(['D:\Data\Model\ROMS\nwp_1_20\backup_surf\', testname, '\run\zeta\', testname, '_AKP4_GCM_ssh_2021_2050.mat'])
% %     for t=1:size(GCM_data.yearly_mean,3)
% %         [ts_yearly_GCM(i,t), error_status] = Func_0011_get_area_weighted_mean(GCM_data.yearly_mean(:,:,t), GCM_grid.lon, GCM_grid.lat);
% %     end
% % end
% % mean(ts_yearly_GCM(5,1))




end


