function prog = Func_0030_AR1_prog(init, C_lambda, tsteps, noise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [psd, psd_sig, freq] = Func_0027_get_PSD_siglev(data, sig_level, DOF)
%
% Power Spectral Density using 1-d data [t] based on fft, with significant psd
%
%  input:
%  init         initial value
%  C_lambda     linear coefficient -> dT/dt = -lambda * T (dt = 1)
%  tstep        prognostic steps (lag)
%  noise        standard variation of noise
%
%  output:
%  prog         future prognostic value
%
%  e-mail:      kimyy308@pusan.ac.kr
%
%  Updated      11-Apr-2023 by Yong-Yub Kim using AR1 coef
%  Updated      17-Apr-2023 by Yong-Yub Kim using damped persistence

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% initialize
prog = init;

%% integration
for i=1:tsteps
    prog(i+1)= prog(i) + prog(i) .* C_lambda;
end




end


