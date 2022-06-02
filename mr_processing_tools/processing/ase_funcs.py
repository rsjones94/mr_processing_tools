#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Sky Jones
"""

import numpy as np
import scipy

def ase_logsignalfunc_shorttau(x, xdata):
    ### SLW wrote this 12/1/2020
    ### x is unknowns, xdata is knowns
    ### rho is effective spin density
    ### lambda is vCBV
    ### dw is frequency shift between fully oxygenated and deoxygenated blood.
    ### dw= 4/3 .* pi .* dChi0 .* Hct .* B0 .* OEF
    ### dChi0 is the susceptibility difference between fully oxygenated and deoxygenated blood
    ### R2' = lambda .* dw
    ### B0 is scanner strength
    ### OEF is oxygen extraction fraction
    ### tau is refocusing pulse shift

    Cs=x[1]
    R2primeSquareDivLambda=x[2]
    
    tauvec=xdata[:,1]
    tauTRANSFVEC = (tauvec * 2)^2 ### For short signal, we fit "delta TE" as the An lab calls it, which is 2*tau, squared
    # TEvec=xdata[:,2]/1000;
    # R2vec=xdata[:,3];
    # # # R2primevec=xdata[:,4]; ### For short tau, only use the first echo.
    
    
    logSs = np.zeros(tauTRANSFVEC.shape())
    for i,val in enumerate(logSs):
        tauTRANSF=tauTRANSFVEC(i)
        logSs[i]=Cs-0.3*R2primeSquareDivLambda*tauTRANSF
    
    
    



def ase_yablonskiy_fit(tauvec,TEvec,hct,RFon,sigData_all,B0):
    ### SLW wrote this.
    ### rho is effective spin density
    ### lambda is vCBV
    ### dw is frequency shift between fully oxygenated and deoxygenated blood.
    ### dw= 4/3 * pi * dChi0 * Hct * B0 * OEF
    ### dChi0 is the susceptibility difference between fully oxygenated and deoxygenated blood
    ### R2' = lambda * dw
    ### B0 is scanner strength
    ### OEF is oxygen extraction fraction
    ### tau is refocusing pulse shift
    ### x = fmincon(fun,x0,A,b)### starts at x0 to find minimizer of fun subject to A*x <= b
    ### x = lsqcurvefit(fun,x0,xdata,ydata)
    ### Hct is small vessel hematocrit
    
    # # #  In ASE_SignalFunction_v1
    # # lambda=x(1);
    # # dw=x(2);
    # # T1=x(3);
    # # T2=x(4);
    # # 
    # # tauvec=xdata(:,1);
    # # TR=xdata(:,2);
    # # TE=xdata(:,3);
    
    taucut=0.010;### How long is a "long tau"
    longtauinds=tauvec>=taucut;
    shorttauinds=tauvec<taucut;
    
    # B0=3;
    dChi0=0.18; ## ppm per unit Hct
    gamma=267.5; ## Mrad /s / T
    
    sigData=sigData_all(longtauinds);
    
    if RFon:
        ########################################
        ################## FIRST FIT, long tau
        # # x=[C, R2, R2prime]; 
        x0=[sigData(1)*1.2, 20, 20]; ### guesses.
        lb=[0,5,1];## lower bounds 
        ub=[sigData(1)*2.5,150,150];## upper bounds 
    
        xdata=zeros(sum(longtauinds),2);
        xdata(:,1)=tauvec(longtauinds);
        xdata(:,2)=TEvec(longtauinds);
        
        xOUT_long = lsqcurvefit(@ASE_logSignalFunc_longtau_RFon_v3,x0,xdata,log(sigData),lb,ub,optimset('Display','off'));
    #     Cl=xOUT_long(1)###intercept
        R2=xOUT_long(2);###R2 (RF on)
        otherR=R2;
        R2prime=xOUT_long(3); ### Used to calculaed OEF
        
        lnSl_extrap=ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,[0 min(TEvec)]);
        
    #     #### PLOT
    #     clf
    #     plot(xdata(:,1),log(sigData),'k.','MarkerSize',16)
    #     hold on
    #     plot(xdata(:,1),ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,xdata),'r.','MarkerSize',16)
    #     linevec=0:0.001:0.02;
    #     linexdata=zeros(length(linevec),2);
    #     linexdata(:,1)=linevec;
    #     linexdata(:,2)=min(TEvec);
    #     plot(linevec,ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,linexdata),'r-','LineWidth',3);
    #     #### END PLOT
    else
        ########################################
        ################## FIRST FIT, long tau
        # # x=[C, R2star, R2prime]; 
        x0=[sigData(1)*1.2, 50, 20]; ### guesses. 
        lb=[0,10,1];## lower bounds 
        ub=[sigData(1)*2.5,200,150];## upper bounds 
    
        dTEvec=TEvec-min(TEvec(longtauinds));# grab shorttauinds below
        xdata=zeros(sum(longtauinds),2);
        xdata(:,1)=tauvec(longtauinds);
        xdata(:,2)=dTEvec(longtauinds);
    
        xOUT_long = lsqcurvefit(@ASE_logSignalFunc_longtau_RFoff_v3,x0,xdata,log(sigData),lb,ub,optimset('Display','off'));
    #     Cl=xOUT_long(1)###intercept
        R2star=xOUT_long(2);###either R2 (RF on) or R2star (RF off)
        otherR=R2star;
        R2prime=xOUT_long(3); ### Used to calculaed OEF
        
        lnSl_extrap=xOUT_long(1);### Using TE1, because that's what Ss will use
    # ########### Plot
    #     clf
    #     plot(2*xdata(:,1),log(sigData),'k.','MarkerSize',16)
    #     hold on
    #     fitx=zeros(20,2);
    #     fitx(:,1)=0:max(xdata(:,1))/19:max(xdata(:,1));
    #     fitx(:,2)=0;
    #     plot(2*fitx(:,1),ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,fitx),'r-','LineWidth',2)
    #     fitx=zeros(10,2);
    #     fitx(:,1)=taucut:max(xdata(:,1))/18:max(xdata(:,1));
    #     fitx(:,2)=max(dTEvec);
    #     plot(2*fitx(:,1),ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,fitx),'r-','LineWidth',2,'HandleVisibility','off')
    # ############### Endplot
    end
    
    ########################################
    ################## SECOND FIT, short tau
    sigData=sigData_all(shorttauinds);
    # # x=[C, R2primeSquareDivLambda]; 
    x0=[sigData(1)*1.2, 8000]; ### guesses. 
    lb=[0,6.4];## lower bounds 
    ub=[sigData(1)*2.5,2e6];## upper bounds 
    
    xdata=zeros(sum(shorttauinds),1);
    xdata(:,1)=tauvec(shorttauinds);
    
    TE1=min(TEvec(shorttauinds));
    includeinds=TEvec(shorttauinds)==TE1; ## For short tau, only include first echo
    xdata=xdata(includeinds,:);
    sigData=sigData(includeinds);
    
    xOUT_short = lsqcurvefit(@ASE_logSignalFunc_shorttau_v3,x0,xdata,log(sigData),lb,ub,optimset('Display','off'));
    # Cs=xOUT_short(1)###intercept
    # R2primeSquareDivLambda=xOUT_short(2);### OEFcalc ;
    
    ######## Subtract fits to get lambda. Calculate dw from the results. Check,
    ######## second way to get lambda with slopes of two fits.
    
    lnSs0=xOUT_short(1);### Tau of 0
    lambda=lnSl_extrap-lnSs0;
    # lambda2=R2prime^2/R2primeSquareDivLambda;#### Second way to calculate
    # # # # lambda. Not how they do it at wash-u. Seems to give worse answers
    
    OEF=R2prime/(gamma*lambda*4/3*pi*dChi0*hct*B0);
    
    # ############### Plot
    # plot(xdata(:,1)*2,log(sigData),'k.','MarkerSize',16,'HandleVisibility','off')
    # hold on
    # plot(xdata(:,1)*2,ASE_logSignalFunc_shorttau_v3(xOUT_short,xdata),'b-','LineWidth',2)
    # pause
    # ########### Endplot
    
    
    #######################################
    ############ quantify fit quality
    
    calcData=zeros(1,length(tauvec));
    
    if RFon
        xdata=zeros(sum(longtauinds),2);
        xdata(:,1)=tauvec(longtauinds);
        xdata(:,2)=TEvec(longtauinds);
        calcData(longtauinds)=ASE_logSignalFunc_longtau_RFon_v3(xOUT_long,xdata);
    else
        dTEvec=TEvec-min(TEvec(longtauinds));# grab shorttauinds below
        xdata=zeros(sum(longtauinds),2);
        xdata(:,1)=tauvec(longtauinds);
        xdata(:,2)=dTEvec(longtauinds);
        calcData(longtauinds)=ASE_logSignalFunc_longtau_RFoff_v3(xOUT_long,xdata);
    end
    xdata=zeros(sum(shorttauinds),1);
    xdata(:,1)=tauvec(shorttauinds);
    
    TE1=min(TEvec(shorttauinds));
    includeinds=TEvec(shorttauinds)==TE1; ## For short tau, only include first echo
    xdata=xdata(includeinds,:);
    calcData(includeinds)=ASE_logSignalFunc_shorttau_v3(xOUT_short,xdata);
    
    shortecho2=TEvec~=TE1 & shorttauinds;
    ### TRuncating these at the very end to calculate Rsquared
    calcData(shortecho2)=[]; sigData_all(shortecho2)=[]; sigData=sigData_all';
    
    n=length(sigData);
    R2numer=n*sum(calcData.*sigData)-sum(calcData)*sum(sigData);
    R2denom=sqrt( (n*sum(calcData.^2)-(sum(calcData).^2))*(n*sum(sigData.^2)-(sum(sigData).^2)) );
    Rsquared=(R2numer/R2denom).^2; ### Calculate Rsquared
    if imag(Rsquared) ~= 0
        Rsquared = 0;
    end
    ####################
    
    end

