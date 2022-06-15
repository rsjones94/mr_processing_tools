
% Author: Allison Griffin
% Date: 16 Oct 2019
% Purpose: calculate Y from T2 (ms), Hct (percent), Tcpmg (ms), model
% All values should be in seconds


function[Y_wood_model] = wood_model(T2, hct, Tcpmg, model)
%T2=T2/1000;
Tcpmg=Tcpmg/1000;
% hct=hct/100;

%% equation is 1/T2 = A1*hct(1-y)^2 + A2(1-y)^2 + A3*hct + A4
if model == 1 %%% this is wood hba 
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
syms y
eqn = (1/T2) == ((A1*hct*(1-y)^2)  + (A2*(1-y)^2)  + (A3*hct) + A4);
Y = solve(eqn,y);

% take the first root only
Y=Y(1);

Y_wood_model=double(Y);
end

%% equation is 1/T2 = A * (1-y)^2 + B
if model == 2 %%% this is wood hbs
    if Tcpmg == 0.01 %% Tcpmg in seconds
        A= 70.0;
        B= 5.75;
    elseif Tcpmg == 0.02
        A= 93.1;
        B= 7.16;
    else
        error('Invalid tau value')
    end
syms y
eqn = (1/T2) == ((A*(1-y)^2)  + B);
Y = solve(eqn,y);

% take the first root only
Y=Y(1);

Y_wood_model=double(Y);
end

end

