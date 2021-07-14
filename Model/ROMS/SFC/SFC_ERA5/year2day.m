function nd = yeardays(y)
% if mod(y,4)==0
if leapyear(y)==1
    nd = 366;
else
    nd=365;
end