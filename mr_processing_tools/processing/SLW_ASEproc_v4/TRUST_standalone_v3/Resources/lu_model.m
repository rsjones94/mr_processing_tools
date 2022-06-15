
% Author: Allison Griffin
% Date: 11 Sep 2019
% Purpose: calculate Y from T2 (ms), Hct (percent), Tcpmg (ms), model
% All values should be in seconds

%% equation is T2 = A + B * (1-y) + C * (1-y)**2
function[Y_lu_model] = lu_model(T2, hct, Tcpmg, model)
%T2=T2/1000; % from T2 of 75 ms --> to .075 (seconds)
Tcpmg=Tcpmg/1000;%Tcpmg of 10 (ms) --> .01 (seconds)
 

if model == 1 %%% this is lu hba (bovine)
    if Tcpmg == 0.005 %% Tcpmg in seconds
        a1= (-4.4);
        a2= 39.1;
        a3= (-33.5);
        b1= (1.5);
        b2= 4.7;
        c1= 167.8;
    elseif Tcpmg == 0.01
        a1= (-13.5);
        a2= 80.2;
        a3= (-75.9);
        b1= (-0.5);
        b2= 3.4;
        c1=247.4;
    elseif Tcpmg == 0.015
        a1= (-12.0);
        a2= 77.7;
        a3= (-75.5);
        b1= (-6.6);
        b2= 31.4;
        c1= 249.4;
    elseif Tcpmg == 0.02
        a1= 7;
        a2= (-9.2);
        a3= 23.2;
        b1= (-4.5);
        b2= 5.3;
        c1= 310.8;
    else
        error('Invalid tau value')
    end 
end
    

if model == 2 %% This is hbf 
    a1= (-1.1);
    a2= 24.0;
    a3= (-21.4);
    b1= (-5.1);
    b2= 29.4;
    c1=247.9;
end 

A = a1 + a2 * hct + a3 * (hct^2);
B = b1 + hct + b2 * (hct^2);
C = c1 * hct * (1 - hct);

syms y
eqn = (1/T2) == (A + (B * (1-y)) + (C * (1-y)^2));
Y = solve(eqn,y);

% take the first root only
Y=Y(1);

Y_lu_model= double(Y);


end 

