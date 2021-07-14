function nd = year2day(y)
if mod(y,4)==0 && mod(y,100)~=0 || mod(y,400)==0
    nd=366;
else
    nd=365;
end
return