clear all;

base_dir = ['/Volumes/kyy_raid/kimyy/Research/Postdoc/03_IBS/', ...
            '2025_LIM_BGC_prediction/manuscripts/Revision_round1/To_Yong_Yub_n5/'];

for No = 1:50

    fprintf('=== Processing No_%d ===\n', No);

    %% ---------------------------------------------------------
    % 0. Load MAT file
    %% ---------------------------------------------------------
    matfile = sprintf('%sNPP_LIM_n_5_No_%d.mat', base_dir, No);
    assert(isfile(matfile), ['Missing file: ' matfile]);
    load(matfile);

    %% ---------------------------------------------------------
    % 1. Dimension inference & checks
    %% ---------------------------------------------------------
    n_state = 5;
    n_month = 12;
    n_lead  = 25;

    % --- xvec ---
    assert(size(xvec,1) == n_state, 'xvec: state dimension mismatch');
    n_time = size(xvec,2);


    %% ---------------------------------------------------------
    % 2. Create NetCDF file
    %% ---------------------------------------------------------
    ncfile = sprintf('%sNPP_LIM_n_5_No_%d.nc', base_dir, No);
    if isfile(ncfile)
        delete(ncfile);
    end

    ncid = netcdf.create(ncfile, 'NETCDF4');

    % ---- dimensions ----
    dim.state      = netcdf.defDim(ncid, 'state', n_state);
    dim.state_i    = netcdf.defDim(ncid, 'state_i', n_state);
    dim.state_j    = netcdf.defDim(ncid, 'state_j', n_state);
    dim.month      = netcdf.defDim(ncid, 'month', n_month);
    dim.lead       = netcdf.defDim(ncid, 'lead', n_lead);
    dim.time       = netcdf.defDim(ncid, 'time', n_time);
    % dim.init_year  = netcdf.defDim(ncid, 'init_year', n_init_year);

    %% ---------------------------------------------------------
    % 3. Define variables
    %% ---------------------------------------------------------
    varid = struct();

    % ===== A, D, Q =====
    mat_types = {'A','D','Q'};
    modes     = {'CS','ST'};
    colors    = {'CW','White','Colored'};

    for mt = 1:numel(mat_types)
        for c = 1:numel(colors)
            for s = 1:numel(modes)

                vname = sprintf('%s_%s_%s', mat_types{mt}, modes{s}, colors{c});
                if ~exist(vname,'var'); continue; end
                V = eval(vname);

                if ismember(mat_types{mt}, {'A','D'})
                    if strcmp(modes{s}, 'CS')
                        assert(isequal(size(V), [5 5 12]), [vname ' must be 5x5x12']);
                        varid.(vname) = netcdf.defVar(ncid, vname, 'double', ...
                            [dim.state_i, dim.state_j, dim.month]);
                    else
                        assert(isequal(size(V), [5 5]), [vname ' must be 5x5']);
                        varid.(vname) = netcdf.defVar(ncid, vname, 'double', ...
                            [dim.state_i, dim.state_j]);
                    end
                else % Q
                    assert(isequal(size(V), [5 5]), [vname ' must be 5x5']);
                    varid.(vname) = netcdf.defVar(ncid, vname, 'double', ...
                        [dim.state_i, dim.state_j]);
                end
            end
        end
    end

    % ===== K =====
    vname = sprintf('%s_%s_%s', 'K', 'CS', 'CW');
    if ~exist(vname,'var'); continue; end
    V = eval(vname);
    assert(isequal(size(V), [5 5 12 25]), [vname ' must be 5x5x12x25']);
    varid.(vname) = netcdf.defVar(ncid, vname, 'double', ...
        [dim.state_i, dim.state_j, dim.month, dim.lead]);

    % ===== Gamma =====
    gam_modes  = {'CS','ST'};
    gam_colors = {'CW','White','Colored'};

    for gm = 1:numel(gam_modes)
        for gc = 1:numel(gam_colors)

            vname = sprintf('Gam_%s_%s', gam_modes{gm}, gam_colors{gc});
            if ~exist(vname,'var'); continue; end
            V = eval(vname);

            assert(isvector(V) && numel(V) == n_state, ...
                [vname ' must be vector length 4']);

            varid.(vname) = netcdf.defVar(ncid, vname, 'double', dim.state);
        end
    end

    % ===== xvec =====
    varid.xvec = netcdf.defVar(ncid, 'xvec', 'double', ...
        [dim.state, dim.time]);

    % % ===== Y variants =====
    % Y_names = {'Y','Y_CS_CW','Y_CS_White','Y_CS_Colored', ...
    %            'Y_ST_CW','Y_ST_White','Y_ST_Colored'};
    % 
    % for i = 1:numel(Y_names)
    %     if exist(Y_names{i},'var')
    %         varid.(Y_names{i}) = netcdf.defVar(ncid, Y_names{i}, 'double', ...
    %             [dim.state, dim.init_year, dim.lead, dim.month]);
    %     end
    % end

    %% ---------------------------------------------------------
    % 4. End definition mode
    %% ---------------------------------------------------------
    netcdf.endDef(ncid);

    %% ---------------------------------------------------------
    % 5. Write variables
    %% ---------------------------------------------------------
    fnames = fieldnames(varid);

    for i = 1:numel(fnames)
        vname = fnames{i};
        v = eval(vname);

        % if iscell(v)
        %     Y_out = nan(n_state, n_init_year, n_lead, n_month);
        %     for m = 1:n_month
        %         for l = 1:n_lead
        %             Y_ml = v{m,l};
        %             if isempty(Y_ml); continue; end
        %             n_valid = size(Y_ml,2);
        %             Y_out(:,1:n_valid,l,m) = Y_ml;
        %         end
        %     end
        %     netcdf.putVar(ncid, varid.(vname), Y_out);
        % else
            netcdf.putVar(ncid, varid.(vname), v);
        % end
    end

    %% ---------------------------------------------------------
    % 6. Close
    %% ---------------------------------------------------------
    netcdf.close(ncid);

    fprintf('✔ Finished No_%d\n', No);
end

disp('🎉 All No_1 to No_50 NetCDF files successfully created.');

