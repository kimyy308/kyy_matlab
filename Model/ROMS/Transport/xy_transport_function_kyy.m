% %  Updated 22-Jun-2022 by Yong-Yub Kim, freshwater transport

% function [sum_segment,sum_segment_salt,sum_segment_fresh] = xy_transport_function(filename,g,selectdepth);
% ====================================================================
% [volume,heat,salinity] = ij_transport_onefile(filename,g,selectdepth);
% calculate transports of  volume(in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%
% across a line segment (slice along a constant I or J)
%
% 'filename'       = history or average file name
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
% 'selectdepth' = surface to which depth (m) such as 100, 500, 50000
%
% keep in mind that
% vertical coordinate changes in time = h + zeta(t) in ROMS
% ====================================================================
% version 1 for single   input file  and single segment.
% version 2 for multiple input files and single segment.
% onefile   for single   input file  and single segment.

% *******************************************************************
% FOR FIXED I (or longitude)
% *******************************************************************

% *******************************************************************
% FOR FIXED J (or latitude)
% *******************************************************************

% reference salt, 122.2E~123.4E; 23.8N~25N
switch testname
    case {'test2117', 'test2127'}
        salt_ref=34.7838; %1985~2014 mean
    case {'test2118', 'test2128'}
        salt_ref=33.6273; %1985~2014 mean
    case {'test2119', 'test2129'}
        salt_ref=33.8581; %1985~2014 mean
    case {'test2120', 'test2130'}
        salt_ref=34.8897; %1985~2014 mean
    case {'test2121', 'test2131'}
        salt_ref=33.7970; %1985~2014 mean
    case {'v04'} % SODA 3.4.2
        salt_ref=34.3742; %1985~2014 mean
    case {'v05'} % GLORYS
        salt_ref=34.3957; %1993~2014 mean
end

   if ( vpoint ==  1 )   

    % along a latitude
    %    _ __ 
    %     |
    %     | v  Hz(k+1,jj,i) 
    %    _|__
    %     |
    %     | v  Hz(k  ,jj,i)
    %    _|__
   
    jj=row_index(index);
    i1=col_index(index);
    i2=col_index(index);

    % actual calculation of transport

    sum_segment=0;
    sum_segment_salt=0;
    sum_segment_temp=0;
    sum_segment_fresh=0;

    for i=i1:i2
      Ip=find( z_v(:,jj,i) > selectdepth );
      lastlevel=min([N; min(Ip)]);
      
      sum=0;
      sum_salt=0;
      sum_fresh=0;
      sum_temp=0;
        seg_rho_fresh=1; %kg/m^3 

      for k=N:-1:lastlevel         % whole cell transport
         sum=sum+v(k,jj,i)*Hz_v(k,jj,i);
         sum_temp=sum_temp+v(k,jj,i)*Hz_v(k,jj,i)*temp_v(k,jj,i);
         temp_in_situ = sw_temp(salt_v(k,jj,i),temp_v(k,jj,i),-g.z_r(k,jj,i),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(k,jj,i));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
         seg_rho=sw_dens(salt_v(k,jj,i),temp_in_situ,-g.z_r(k,jj,i))/1000.0;
         sum_salt=sum_salt+v(k,jj,i)*Hz_v(k,jj,i)*salt_v(k,jj,i)/1000.0*seg_rho;
         sum_fresh=sum_fresh+v(k,jj,i)*Hz_v(k,jj,i)*(1-(salt_v(k,jj,i)/salt_ref))*seg_rho_fresh;
      end
      if (lastlevel == N)          % partial cell transport 
         partialH=zeta(jj,i)-selectdepth;
         sum=sum+v(N,jj,i)*partialH;
%          sum_salt=sum_salt+v(N,jj,i)*partialH*salt_v(N,jj,i);
         sum_temp=sum_temp+v(N,jj,i)*partialH*temp_v(N,jj,i);
         temp_in_situ = sw_temp(salt_v(N,jj,i),temp_v(N,jj,i),-g.z_r(N,jj,i),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(k,jj,i));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
         seg_rho=sw_dens(salt_v(N,jj,i),temp_in_situ,-g.z_r(N,jj,i))/1000.0;
         sum_salt=sum_salt+v(N,jj,i)*partialH*salt_v(N,jj,i)/1000.0*seg_rho;
         sum_fresh=sum_fresh+v(N,jj,i)*partialH*(1-(salt_v(N,jj,i)/salt_ref))*seg_rho_fresh;
      elseif (lastlevel ~= 1)         
         partialH=z_v(k,jj,i)-selectdepth;
         sum=sum+v(k-1,jj,i)*partialH;
         sum_salt=sum_salt+v(k-1,jj,i)*partialH*salt_v(k-1,jj,i);
         sum_temp=sum_temp+v(k-1,jj,i)*partialH*temp_v(k-1,jj,i);
         temp_in_situ = sw_temp(salt_v(k-1,jj,i),temp_v(k-1,jj,i),-g.z_r(k-1,jj,i),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(k,jj,i));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
         seg_rho=sw_dens(salt_v(k-1,jj,i),temp_in_situ,-g.z_r(k-1,jj,i))/1000.0;
         sum_salt=sum_salt+v(k-1,jj,i)*partialH*salt_v(k-1,jj,i)/1000.0*seg_rho;
         sum_fresh=sum_fresh+v(k-1,jj,i)*partialH*(1-(salt_v(k-1,jj,i)/salt_ref))*seg_rho_fresh;
      end
      sum_segment=sum_segment+sum*dx_v(jj,i);
      sum_segment_temp=sum_segment_temp+sum_temp*dx_v(jj,i);
      sum_segment_salt=sum_segment_salt+sum_salt*dx_v(jj,i);
      sum_segment_fresh=sum_segment_fresh+sum_fresh*dx_v(jj,i);
    end


  elseif  ( vpoint == 0 )  

    % along a longitude
    %    _ __ 
    %     |
    %     | u  Hz(k+1,j,ii) 
    %    _|__
    %     |
    %     | u  Hz(k  ,j,ii)
    %    _|__

    ii=col_index(index); 
    j1=row_index(index);
    j2=row_index(index);

    % actual calculation of transport

    sum_segment=0;
    sum_segment_salt=0;
    sum_segment_fresh=0;
    sum_segment_temp=0;
    for j=j1:j2
      Ip=find( z_u(:,j,ii) > selectdepth );
      lastlevel=min([N; min(Ip)]);
      sum=0;
      sum_salt=0;
      sum_temp=0;      
      sum_fresh=0;
      seg_rho_fresh=1; %kg/m^3 

      for k=N:-1:lastlevel         % whole cell transport
          sum=sum+u(k,j,ii)*Hz_u(k,j,ii);
%          sum_salt=sum_salt+u(k,j,ii)*Hz_u(k,j,ii)*salt_u(k,j,ii);
          sum_temp=sum_temp+u(k,j,ii)*Hz_u(k,j,ii)*temp_u(k,j,ii);
          temp_in_situ = sw_temp(salt_u(k,j,ii),temp_u(k,j,ii),-g.z_r(k,j,ii),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(k,j,ii));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
         seg_rho=sw_dens(salt_u(k,j,ii),temp_in_situ,-g.z_r(k,j,ii))/1000.0;
         sum_salt=sum_salt+u(k,j,ii)*Hz_u(k,j,ii)*salt_u(k,j,ii)/1000.0*seg_rho;
         sum_fresh=sum_fresh+u(k,j,ii)*Hz_u(k,j,ii)*(1-(salt_u(k,j,ii)/salt_ref))*seg_rho_fresh;
      end
      if (lastlevel == N)          % partial cell transport 
           partialH=zeta(j,ii)-selectdepth;
           sum=sum+u(N,j,ii)*partialH;
%           sum_salt =sum_salt +u(N,j,ii)*partialH*salt_u(N,j,ii);
           sum_temp=sum_temp+u(N,j,ii)*partialH*temp_u(N,j,ii);
           temp_in_situ = sw_temp(salt_u(N,j,ii),temp_u(N,j,ii),-g.z_r(N,j,ii),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(N,j,ii));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
             seg_rho=sw_dens(salt_u(N,j,ii),temp_in_situ,-g.z_r(N,j,ii))/1000.0;
             sum_salt=sum_salt+u(N,j,ii)*partialH*salt_u(N,j,ii)/1000.0*seg_rho;
             sum_fresh=sum_fresh+u(N,j,ii)*partialH*(1-(salt_u(N,j,ii)/salt_ref))*seg_rho_fresh;
      elseif (lastlevel ~= 1)     
          partialH=z_u(k,j,ii)-selectdepth;
          sum=sum+u(k-1,j,ii)*partialH;
%          sum_salt =sum_salt +u(k-1,j,ii)*partialH*salt_u(k-1,j,ii);
          sum_temp=sum_temp+u(k-1,j,ii)*partialH*temp_u(k-1,j,ii);
          temp_in_situ = sw_temp(salt_u(k-1,j,ii),temp_u(k-1,j,ii),-g.z_r(k-1,j,ii),0);
%          seg_rho_fresh=sw_dens(0,temp_in_situ,-g.z_r(k-1,j,ii));
%          (20degree,0m, 998.25 ~ 5000m, 1020)
             seg_rho=sw_dens(salt_u(k-1,j,ii),temp_in_situ,-g.z_r(k-1,j,ii))/1000.0;
             sum_salt=sum_salt+u(k-1,j,ii)*partialH*salt_u(k-1,j,ii)/1000.0*seg_rho;
             sum_fresh=sum_fresh+u(k-1,j,ii)*partialH*(1-(salt_u(k-1,j,ii)/salt_ref))*seg_rho_fresh;
      end
      
      sum_segment=sum_segment+sum*dy_u(j,ii);
      sum_segment_temp=sum_segment_temp+sum_temp*dy_u(j,ii);
      sum_segment_salt =sum_segment_salt +sum_salt*dy_u(j,ii);      
      sum_segment_fresh=sum_segment_fresh+sum_fresh*dy_u(j,ii);
    end


  end



