function [ logSs ] = ASE_logSignalFunc_shorttau_v3( x,xdata )
%%% SLW wrote this 12/1/2020
%%% x is unknowns, xdata is knowns
%%% rho is effective spin density
%%% lambda is vCBV
%%% dw is frequency shift between fully oxygenated and deoxygenated blood.
%%% dw= 4/3 .* pi .* dChi0 .* Hct .* B0 .* OEF
%%% dChi0 is the susceptibility difference between fully oxygenated and deoxygenated blood
%%% R2' = lambda .* dw
%%% B0 is scanner strength
%%% OEF is oxygen extraction fraction
%%% tau is refocusing pulse shift

if length(x)~=2 || size(xdata,2)~=1
    error('data was not entered correctly')
end

Cs=x(1);
R2primeSquareDivLambda=x(2);

tauvec=xdata(:,1);
tauTRANSFVEC=(tauvec*2).^2; %%% For short signal, we fit "delta TE" as the An lab calls it, which is 2*tau, squared
% TEvec=xdata(:,2)/1000;
% R2vec=xdata(:,3);
% % % R2primevec=xdata(:,4); %%% For short tau, only use the first echo.


logSs=zeros(size(tauTRANSFVEC));
for i = 1:length(logSs)
    tauTRANSF=tauTRANSFVEC(i);
    logSs(i)=Cs-0.3*R2primeSquareDivLambda*tauTRANSF;
end

end

