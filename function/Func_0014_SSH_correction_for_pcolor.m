function [ssh_correction_for_fig] = Func_0014_SSH_correction_for_pcolor(testname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function var=Func_0014_SSH_correction_for_pcolor(testname);
%
% Getting the SSH correction value for pcolor
%
%  input:
%  testname             Name of the dataset (string. RCM, GCM, Observation, etc.)
%
%  output:
%  ssh_correction_for_fig          correction value (string)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    13-Jul-2021 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch(testname)
        case('CMEMS')
%         ssh_correction_for_fig=-0.1096; % 1993 NWP min
        ssh_correction_for_fig=0.7207; % 1993-2014 AKP4 mean       
    case('test2102')
%         ssh_correction_for_fig=-0.9915; % 1993 NWP min
        ssh_correction_for_fig=-0.3314; % 1993-2014 AKP4 mean
    case('test2103')
%         ssh_correction_for_fig=-1.0469; % 1993 NWP min
        ssh_correction_for_fig=-0.2376; % 1993-2014 AKP4 mean        
    case('test2104')
%         ssh_correction_for_fig=-1.0217; % 1993 NWP min
        ssh_correction_for_fig=-0.2518; % 1993-2014 AKP4 mean                
    case('test2105')
%         ssh_correction_for_fig=-0.9915; % 1993 NWP min
        ssh_correction_for_fig=-0.2822; % 1993-2014 AKP4 mean                
    case('test2106')
%         ssh_correction_for_fig=-1.0439; % 1993 NWP min
        ssh_correction_for_fig=-0.1680; % 1993-2014 AKP4 mean      
    case('RCM_ENS_historical')
        ssh_correction_for_fig=-0.2542; % 1993-2014 AKP4 mean
        
    case('CNRM-ESM2-1')
%         ssh_correction_for_fig=-5.4472; % 1985 NWP mean
%         ssh_correction_for_fig=-6.2406; % 1985 NWP min
%         ssh_correction_for_fig=-6.2427; % 1993 NWP min
        ssh_correction_for_fig=-5.4829; % 1993-2014 AKP4 mean           
    case('EC-Earth3-Veg')
%         ssh_correction_for_fig=-5.3736;
%         ssh_correction_for_fig=-6.3315; % 1985 NWP min
%         ssh_correction_for_fig=-6.3381; % 1993 NWP min
        ssh_correction_for_fig=-5.4008; % 1993-2014 AKP4 mean           
    case('ACCESS-CM2')
%         ssh_correction_for_fig=-5.3736;
%         ssh_correction_for_fig=-0.4512; % 1985 NWP min
%         ssh_correction_for_fig=-0.3847; % 1993 NWP min
        ssh_correction_for_fig=0.6046; % 1993-2014 AKP4 mean                   
    case('CNRM-CM6-1-HR')
%         ssh_correction_for_fig=-5.4724;
%         ssh_correction_for_fig=-6.3851; % 1985 NWP min
%         ssh_correction_for_fig=-6.2771; % 1993 NWP min
        ssh_correction_for_fig=-5.5167; % 1993-2014 AKP4 mean                           
    case('CMCC-ESM2')
%         ssh_correction_for_fig=-5.6812;
%         ssh_correction_for_fig=-6.6958; % 1985 NWP min
%         ssh_correction_for_fig=-6.7129; % 1993 NWP min
        ssh_correction_for_fig=-5.7181; % 1993-2014 AKP4 mean            
    case('GCM_ENS_historical')
        ssh_correction_for_fig=-4.1331; % 1993-2014 AKP4 mean
end

end

