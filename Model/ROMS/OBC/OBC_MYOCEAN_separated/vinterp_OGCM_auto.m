function ROMS = vinterp_OGCM_auto(ROMS, OGCM_depth, ind_temp, ind_salt, ind_u, ind_v)

%
% Perform the vertical interpolations for the bry_file with OGCM for each boundary
%
% zr_bry = squeeze(zr_bry);
% zu_bry = squeeze(zu_bry);
% zv_bry = squeeze(zv_bry);
% dzr_bry = squeeze(dzr_bry);
% dzu_bry = squeeze(dzu_bry);
% dzv_bry = squeeze(dzv_bry);
%
% Add a level on top and bottom with no-gradient
%
for diri=1:4  % add first layer and last layer
    for vari=[ind_salt,ind_temp,ind_u,ind_v] % 3: salt, 4: temp, 7: u, 9: v
        ROMS.data(diri).var(vari).value = ...
            cat(3, ROMS.data(diri).var(vari).value(:,:,1), ...
            ROMS.data(diri).var(vari).value);
        ROMS.data(diri).var(vari).value = ...
            cat(3, ROMS.data(diri).var(vari).value, ...
            ROMS.data(diri).var(vari).value(:,:,end));
    end
end

% pcolor(squeeze(ROMS.data(diri).var(vari).value))
% shading flat;

%
% Perform the vertical interpolations
%
for diri=1:4
    for vari=[ind_salt, ind_temp] % 3: salt, 4: temp, 7: u, 9: v
        tempvar = flipdim(squeeze(ROMS.data(diri).var(vari).value)', 1);
        ROMS.data(diri).var(vari).value = ...
            ztosigma_1d(tempvar, ...
            ROMS.data(diri).zr, flipud(OGCM_depth));
    end
    tempu = flipdim(squeeze(ROMS.data(diri).var(ind_u).value)', 1);
    tempv = flipdim(squeeze(ROMS.data(diri).var(ind_v).value)', 1);
    ROMS.data(diri).var(ind_u).value = ...
        ztosigma_1d(tempu, ...
        ROMS.data(diri).zu, flipud(OGCM_depth));
    ROMS.data(diri).var(ind_v).value = ...
        ztosigma_1d(tempv, ...
        ROMS.data(diri).zv, flipud(OGCM_depth));
end

%
% Correct the horizontal transport
% i.e. remove the interpolated tranport and add
%      the OGCM transport
%
if ROMS.conserv == 1
    u_bry = u_bry-squeeze(tridim(squeeze(sum(u_bry.*dzu_bry)./sum(dzu_bry)),N));
    v_bry = v_bry-squeeze(tridim(squeeze(sum(v_bry.*dzv_bry)./sum(dzv_bry)),N));
    u_bry = u_bry + squeeze(tridim(ubar_bry,N));
    v_bry = v_bry + squeeze(tridim(vbar_bry,N));
end
%
% Barotropic velocities
%
% ubar_bry = squeeze(sum(u_bry.*dzu_bry)./sum(dzu_bry));
% vbar_bry = squeeze(sum(v_bry.*dzv_bry)./sum(dzv_bry));
% caculate barotropic u and v

for diri=1:4
    ROMS.data(diri).var(10).value = ...
        sum(ROMS.data(diri).var(ind_u).value .* ROMS.data(diri).dzu, 1) ...
        ./ sum(ROMS.data(diri).dzu, 1);
    ROMS.data(diri).var(11).value = ...
        sum(ROMS.data(diri).var(ind_v).value .* ROMS.data(diri).dzv, 1) ...
        ./ sum(ROMS.data(diri).dzv, 1);
end
%
return
