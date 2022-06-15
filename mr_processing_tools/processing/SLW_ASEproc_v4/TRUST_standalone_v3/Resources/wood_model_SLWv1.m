
% Author: Allison Griffin
% Date: 16 Oct 2019
% Purpose: calculate Y from T2 (ms), Hct (percent), Tcpmg (ms), model
% All values should be in seconds


function[Y_wood_model] = wood_model(T2, hct, Tcpmg)
%T2=T2/1000;
Tcpmg=Tcpmg/1000;

%% equation is 1/T2 = A1*hct(1-y)^2 + A2(1-y)^2 + A3*hct + A4

if Tcpmg == 0.01 %% Tcpmg in seconds
    A1= (77.5);
    A2= 27.8;
    A3= 6.95;
    A4= 2.34;
elseif Tcpmg == 0.02
A1= 167;
A2= 20.8;
A3= 10.3;
A4= 2.15;
else
    error('Invalid tau value')
end
disp('T2 is')
disp(T2)
disp('hct is')
disp(hct)

numer=-A3*hct-A4+1/T2;
denom=A1*hct+A2;

disp('numer is')
disp(numer)
disp('denom is')
disp(denom)

Y_wood_model=1-(numer/denom)^0.5;

end

