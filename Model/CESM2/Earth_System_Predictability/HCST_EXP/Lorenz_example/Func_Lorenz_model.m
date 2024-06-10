function [Xout, Yout, Zout] = Func_Lorenz_model(Xin, Yin, Zin, dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% function [X, Y, Z] = Func_Lorenz_model(Xin, Yin, Zin, dt);
%
%  Lorenz model written by MATLAB based on Yoshi's python code(Lmodel-v2.py)
%
%  input:
%  Xin          X  (single or double)
%  Yin          Y  (single or double)
%  Zin          Z  (single or double)
%  dt           dt (single or double)  0.01
%
%  output:
%  Xout          X  (single or double)
%  Yout          Y  (single or double)
%  Zout          Z  (single or double)
%
%  e-mail:kimyy308@snu.ac.kr
%
%  Updated    29-Jan-2024 by Yong-Yub Kim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sigma = 10.;
r = 28.;
b = 8./3.

% -- Lorenz model --
Xout = Xin + dt*(-sigma*Xin + sigma *Yin);
Yout = Yin + dt*(-1.0*Xin*Zin + r *Xin - Yin);
Zout = Zin + dt*(Xin*Yin - b*Zin);


end



