function [ xvec_smooth ] = line_filter_slw( xvec, filteropt )
%%% SLW wrote this to streamline the smoothing process I want. Matlabs
%%% smoothing is kindof a pain.
%xvec is a vector of values that should be smoothed. This function is built
%specifically for ASE data smoothing, and 'fft' option will likely not work for many
%other applications, as data is assumed to be a line.
%filteropt can be either 'smooth' or 'fft'. 

padval=5;

%%%%%%%% extrapolate values to account for ringing artifact and smoothing
%%%%%%%% ends
linefunc=@(x,xdata) xdata.*x(2) + x(1);
lenval=length(xvec);
fit1=lsqcurvefit(linefunc,[0 0],1:lenval,xvec,[],[],optimset('Display','off'));
xtemp=[linefunc(fit1,(-padval+1):0) xvec linefunc(fit1,(lenval+1):(lenval+padval))];

if strcmp(filteropt,'fft')
    passradius=6; %%% Higher value here means less smoothing
    freq1=fft(xtemp);
    freq1(passradius)=freq1(passradius)*0.5;
    freq1(end-passradius+1)=freq1(end-passradius+1)*0.5;
    freq1(passradius+1:end-passradius)=0;
    xvec_smooth=real(ifft(freq1));
end

if strcmp(filteropt,'smooth')
    smoothval=5; %%% Higher value here means more smoothing
    xvec_smooth=smooth(xtemp,smoothval);
end

%%%%%%%%% Truncate extrapolated values
xvec_smooth=xvec_smooth(padval+1:end-padval);

end

