function x=swp(p,t,s)
% SWP Sea Water Properties
%
% This function calculates some thermophysical properties of sea water 
% as a function of temperature in(c) and salinity in (ppm)
%
% Properties and their notations:
%    density in (kg/m^3) ---------------------------> 'rho'
%    enthalpy in (kJ/Kg) ---------------------------> 'h' 
%    constant pressure specific heat in(kJ/kg.K) ---> 'cp'
%    thermal conductivity in (kW/m.K) --------------> 'k'
%    dynamic viscosity in (Pa.s) -------------------> 'm'
%    boiling point elevation in (c)-----------------> 'bpe'
% 
% To call the function, type the following:
% swp(a property notation (as above),temperature in (c),salinity in (ppm)) 
% 
% Example: >> density=swp('rho',24,35000)
% 
%          density =
% 
%            1.0231e+003
s=s/10^6;lower(p);
switch p
    case 'rho'
        a=1002.4+754.8*s+236.3*s^2;
        b=(-0.1338-0.935*s-0.0976*s^2);
        c=(-0.003375+0.00996*s-0.439*s^2);
        d=(0.00000313-0.0000163*s+0.000244*s^2);
        x=a+b*t+c*t^2+d*t^3;
    case 'h'
        a=4206.8-6619.7*s+12288*s^2;
        b=-1.1262+54.17*s-227.19*s^2;
        c=0.012026-0.53566*s+1.8906*s^2;
        d=0.00000068774+0.001517*s-0.004268*s^2;
        x=4.1868*(2.3*s-103*s^2+0.000238846*(a*t+1/2*b*t^2+1/3*c*t^3+1/4*d*t^4));
    case 'cp'
        a=4206.8-6619.6*s+12288*s^2;
        b=-1.1262+54.178*s-227.19*s^2;
        c=-0.012026-0.53566*s+1.8906*s^2;
        d=0.00000068774+0.001517*s-0.0044268*s^2;
        x=(a+b*t+c*t^2+d*t^3)/1000;
    case 'k'
        m=28170*s/(1000-s*1000);
        a=576.6-34.64*m+7.286*m^2;
        b=0.001*(1526+466.2*m-226.8*m^2+28.67*m^3);
        c=-0.00001*(581+2055*m-991.6*m^2+146.4*m^3);
        x=(a+b*t+c*t^2)/1000000;
    case 'm'
        a=0.001474+0.000015*t-0.00000003927*t^2;
        b=0.000010734-0.000000085*t+0.000000000223*t^2;
        c=1+a*s*1000+(b*s*1000)^2;
        x=0.001*c*exp(-3.79418+(604.129/((t+273.15)-133.97)));
    case 'bpe'
        t=t+273.15;
        a=33.11452648-0.2008228*t+0.000402*t^2;
        b=-642.170248+3.74328*t-0.0052*t^2;
        c=3210.9472-19.392*t+0.03*t^2;
        x=a*s+b*s^2+c*s^3;
    otherwise
        'there is no such property name'
        help swp
end
% References:
% 1- Gruenberg, Properties of Seawater Concentrates, Dubrovnik 1970.
% 2- Hoemig, H.E, Seawater and Seawater Distillation, Vulkan-Verlag, 1978.