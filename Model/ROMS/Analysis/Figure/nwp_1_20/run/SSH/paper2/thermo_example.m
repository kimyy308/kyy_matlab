for yearij = 1:length(inputyear)
    if (exist('rho_0_va')==0)
        for monthij = 1:length(inputmonth)
            disp([num2str(yearij), 'y_',num2str(monthij),'m'])
            tempyear = inputyear(yearij);
            tempmonth = inputmonth(monthij);
            % ex : rootdir\test37\data\2001\test37_monthly_2001_01.nc
            filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');


           %% read model data
            if (exist('lon_min')==0)
                modelinfo=ncinfo(filename);
                lon = ncread(filename,'lon_rho',[1 1],[modelinfo.Dimensions(5).Length,1]);
                lat = ncread(filename,'lat_rho',[1 1],[1,modelinfo.Dimensions(6).Length]);

                [lon_min, lon_max, lat_min, lat_max] = ...
                    Func_0012_findind_Y(1, lonlat, lon, lat);

                lon = ncread(filename,'lon_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                lat = ncread(filename,'lat_rho', [lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);

                switch(regionname)
                    case('NWP') %% North western Pacific
                        mask_model(1:size(lon,1),1:size(lon,2))=1;
                    otherwise
                        mask_model = double(inpolygon(lon,lat,refpolygon(:,1),refpolygon(:,2)));
                        mask_model(mask_model==0)=NaN;
                end
            end
            data_info = ncinfo(filename, varname); 

            if (exist('h')==0)
                h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
                dA=1./pm .* 1./pn;  %% grid area
                Vtransform = ncread(filename, 'Vtransform');
                Vstretching = ncread(filename, 'Vstretching');
                theta_s = ncread(filename, 'theta_s');
                theta_b = ncread(filename, 'theta_b');
                hc = ncread(filename, 'hc');
                N = length(ncread(filename, 's_rho'));
            end
            data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]); % zeta (sea-level)
            PT_src(:,:,:,monthij) = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
            S_src(:,:,:,monthij) = ncread(filename,'salt',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
            S_src(S_src<=0)=0;

            vtype ='w';
            zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

            dz = diff(zw);
            dz = permute(dz, [2,3,1]);
%                         vol_0    = nansum( is_ocean.*dz.*dA, 'all' );       % global sum, 0-D
%                         area_0   = nansum( is_ocean(:,:,1).*dA, 'all' );    % global sum, 0-D
            vtype ='r';
            P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
            P = -permute(P, [2,3,1]);
            comb_P(:,:,:,monthij)=P;
        end
        %% get reference value (variable_0)
            PT_src = mean(PT_src, 4 );
            S_src = mean(S_src, 4 );
            P = mean(comb_P, 4);
            T_src = sw_temp(S_src, PT_src, P, zeros(size(P)) ); % potential temp -> in-situ temp
            rho_src = sw_dens(S_src, T_src, P ); % get density
            dep_0 = h + data;
            rho_0_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_0;       % vertical(local) average, 2-D
    end

    for monthij = 1:length(inputmonth)
        disp([num2str(yearij), 'y_',num2str(monthij),'m'])
        tic;
        tempyear = inputyear(yearij);
        tempmonth = inputmonth(monthij);
        filename = strcat(monfiledir, num2str(tempyear,'%04i'), '\', ...
                 testname, '_monthly_', num2str(tempyear,'%04i'), '_', num2str(tempmonth,'%02i'), '.nc');

        data_info = ncinfo(filename, varname);  %% [lon lat depth time] -> [1601 1201 33 1]

        if (exist('h')==0)
            h = ncread(filename,'h',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            pm = ncread(filename,'pm',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            pn = ncread(filename,'pn',[lon_min(1) lat_min(1)], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1]);
            dA=1./pm .* 1./pn;
            Vtransform = ncread(filename, 'Vtransform');
            Vstretching = ncread(filename, 'Vstretching');
            theta_s = ncread(filename, 'theta_s');
            theta_b = ncread(filename, 'theta_b');
            hc = ncread(filename, 'hc');
            N = length(ncread(filename, 's_rho'));
        end
        data = ncread(filename,varname,[lon_min(1) lat_min(1) 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 1]);
%%     thermosteric (S should be used as reference value for calculation of thermosteric change)
        PT_src = ncread(filename,'temp',[lon_min(1) lat_min(1) 1 1], [lon_max(1)-lon_min(1)+1 lat_max(1)-lat_min(1)+1 inf 1]);
        S_src(S_src<=0)=0;

        vtype ='w';
        zw=zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);

        dz = diff(zw);
        dz = permute(dz, [2,3,1]);

        vtype ='r';
        P = zlevs(Vtransform, Vstretching, h,data,theta_s,theta_b,hc,N,vtype);
        P = -permute(P, [2,3,1]);
        T_src = sw_temp (S_src, PT_src, P, zeros(size(P)) );
        rho_src = sw_dens (S_src, T_src, P );

        dep_n = h + data;
        rho_n_va = sum( rho_src.*dz, 3 , 'omitnan') ./ dep_n;      % vertical(local) average
%                         rho_n_va2 = sum( rho_src.*dz, 3 , 'omitnan') ./ (h+data);
%                         rho_n_ga = nansum( rho_src.*dz.*dA, 'all' ) ./ vol_0;  % global average
        zosto = dep_0       .* ( 1 - rho_n_va ./ rho_0_va );
        comb_zosto_halo(:,:,ind) = zosto;
       
    end
end