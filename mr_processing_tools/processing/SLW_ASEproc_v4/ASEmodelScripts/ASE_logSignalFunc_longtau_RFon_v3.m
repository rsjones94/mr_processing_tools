function [ logSl ] = ASE_logSignalFunc_longtau_RFon_v3( x,xdata )
%%% SLW Wrote this. 12/1/2020
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

if length(x)~=3 || size(xdata,2)~=2
    error('data was not entered correctly')
end

Cl=x(1);
R2=x(2);
R2prime=x(3);

tauvec=xdata(:,1);%%% leave in s because T2prime is in 1/s
tauTRANSFVEC=tauvec*2;
TEvec=xdata(:,2)/1000; %%%% put in s. R2 is in 1/s

logSl=zeros(size(tauvec));
for i = 1:length(logSl)
    tau=tauTRANSFVEC(i);
    TE=TEvec(i);
    logSl(i)=Cl-R2*TE-R2prime*tau;
end

end

