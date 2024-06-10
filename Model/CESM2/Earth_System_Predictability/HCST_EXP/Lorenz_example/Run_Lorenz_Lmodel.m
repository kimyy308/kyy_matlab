close all; clear all; clc;


%% Model setting
rest = false; % -- True: Restart from file,  False: Start from the begining"
tstr = 0; % -- Start time
tend = 50; % -- End time
tstep = 1; % -- Restart time interval
Lroot= '/Volumes/kyy_raid/kimyy/Model/Lorenz/Lmodel';
rdir = [Lroot, '/rest/']; % -- Restart file directory
rfile = [rdir, 'restart.dddd.nc']; % -- Restart file name
odir = [Lroot, '/out/']; % -- Output file directory
outputf = [odir, 'out.dddd.nc']; % -- output file name
dt = 0.01;


%% run
nt = floor(tstep/dt);
if ~exist(odir, 'dir')
    mkdir(odir);
end
fprintf(' == Run the model ==\n');

for it = tstr:tstep:tend-tstep
    time = it + dt * (0:nt-1);
    fprintf(' time = %f\n', time(1));
    
    % Initial value
    if rest
        rin = strrep(rfile, 'dddd', sprintf('%04d', it));
        fin = load(rdir + rin);
        X = fin.X;
        Y = fin.Y;
        Z = fin.Z;
    else
        mkdir(rdir);
        X = zeros(1, 1);
        Y = 10;
        Z = zeros(1, 1);
    end

    % Model Run
    Xout = zeros(nt, 1);
    Yout = zeros(nt, 1);
    Zout = zeros(nt, 1);

    for i = 1:nt
        Xout(i) = X;
        Yout(i) = Y;
        Zout(i) = Z;
        [X, Y, Z] = Func_Lorenz_model(X, Y, Z, dt);
    end

    etime = time(end) + dt;
    
    

%     X = setfield(X, 'time', etime);
%     Y = setfield(Y, 'time', etime);
%     Z = setfield(Z, 'time', etime);

    % Model output
    fout = strrep(outputf, 'dddd', sprintf('%04d', it));

%%%    netcdf.create(fout)
Writing the code is not finished yet



    save(fullfile(odir, fout), 'Xout', 'Yout', 'Zout');

    % Restart file
    rout = strrep(rfile, 'dddd', sprintf('%04d', int32(etime)));
    save(fullfile(rdir, rout), 'X', 'Y', 'Z');
    rest = true;
end





    