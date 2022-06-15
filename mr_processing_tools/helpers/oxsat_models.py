import numpy as np
import sympy as sym

def bovine_oxsat(T2, hct, hbs=None, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    
    trust_a1 = -13.5
    trust_a2 = 80.2
    trust_a3 = -75.9
    trust_b1 = -0.5
    trust_b2 = 3.4
    trust_c1 = 247.4

    # Compute Venous Oxygenation (Assuming Normal Curve)
    A = trust_a1 + trust_a2*hct + trust_a3*(hct)**2
    B = trust_b1*hct + trust_b2*(hct)**2
    C = trust_c1*hct*(1-hct)
   
    y = sym.symbols('y')
    solution = sym.solveset(A + B*(1-y) + C*(1-y)**2 - (1/T2))
    Y = list(solution)[0]
    
    return float(Y)


def aa_oxsat(T2, hct, hbs=None, tcpmg=10):
    #T2 in seconds
    #hct as float between 0 and 1
    #tcpmg in milliseconds
    
    # wood's AA model
    tcpmg=tcpmg/1000;
    
    # equation is 1/T2 = A1*hct(1-y)^2 + A2(1-y)^2 + A3*hct + A4
    if tcpmg == 0.01:
        A = (77.5)
        B = 27.8
        C = 6.95
        D = 2.34
    elif tcpmg == 0.02:
        A = 167
        B = 20.8
        C = 10.3
        D = 2.15
    else:
        Exception('Invalid tau value')
        
    #syms y
    #eqn = (1/T2) == ((A1*hct*(1-y)^2)  + (A2*(1-y)^2)  + (A3*hct) + A4)
    #Y = solve(eqn,y);
    
    y = sym.symbols('y')
    solution = sym.solveset(((A*hct*(1-y)**2)  + (B*(1-y)**2)  + (C*hct) + D) - (1/T2))
    Y = list(solution)[0]
    
    return float(Y)
    



def ss_oxsat(T2, hct=None, hbs=None, tcpmg=10):
    #T2 in seconds
    #hct as float between 0 and 1
    #tcpmg in milliseconds
    
    
    # wood's SS model
    tcpmg=tcpmg/1000;
    
    # equation is 1/T2 = A * (1-y)^2 + B
    
    if tcpmg == 0.01: # Tcpmg in seconds
        A= 70.0
        B= 5.75
    elif tcpmg == 0.02:
        A= 93.1
        B= 7.16
    else:
        Exception('Invalid tau value')
    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]

    #eqn = (1/T2) == ((A*(1-y)^2)  + B);
    #Y = solve(eqn,y);
    
    #numer = B*T2-1
    #denom = np.sqrt(-A*T2*(B*T2-1))
    
    #Y = numer / denom + 1
    
    return float(Y)
    



def f_oxsat(T2, hct, hbs=None, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    
    # lu model
    # equation is 1/T2 = A + B * (1-y) + C * (1-y)**2
    
    a1= -1.1
    a2= 24.0
    a3= -21.4
    b1= -5.1
    b2= 29.4
    c1= 247.9
    
    A = a1 + a2 * hct + a3 * (hct**2);
    B = b1 + hct + b2 * (hct**2);
    C = c1 * hct * (1 - hct);
    
    y = sym.symbols('y')
    solution = sym.solveset((A + (B * (1-y)) + (C * (1-y)**2)) - (1/T2))
    Y = list(solution)[0]
    
    
    #eqn = (1/T2) == (A + (B * (1-y)) + (C * (1-y)^2));
    #Y = solve(eqn,y);
    
    #Y = (-np.sqrt(T2*(-4*A*T2 + B**2*T2 + 4*C)) + B*T2 + 2*C*T2) / (2*C*T2)
    
    return float(Y)


def mixture_oxsat(T2, hct, hbs, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1

    # based off of DOI: 10.1002/mrm.28757, Bush et. al "Calibration of T2 oximetry MRI for subjects with sickle cell disease"
    # equation is 1/T2 = A * (1-y)^2 + B
    
    # note that this is not the MERGED model (which is just an SS model calibrated on a combined dataset)
    # this is the MIXTURE model, which essentially takes the weighted average of the merged model above and Wood's AA model, using HbS frac as the weighting (more HbS = weight merged model heavily)
    
    
    ay1hba = 77.5;
    ay2hba = 27.8;
    ay1hbs = 196.8;
    ay2hbs = 16.7;
    
    be1hba = 6.95;
    be2hba = 2.34;
    be1hbs = -6.6;
    be2hbs = 8.6;
    
    ay1 = (1 - (hbs/0.8)) * ay1hba + (hbs/0.8)*ay1hbs;
    ay2 = (1 - (hbs/0.8)) * ay2hba + (hbs/0.8)*ay2hbs;
    be1 = (1 - (hbs/0.8)) * be1hba + (hbs/0.8)*be1hbs;
    be2 = (1 - (hbs/0.8)) * be2hba + (hbs/0.8)*be2hbs;
    
    A = ay1*hct + ay2;
    B = be1*hct + be2;
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)


def lidonahue_mixture_oxsat(T2, hct, hbs, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1

    # based off of DOI: 10.1002/mrm.28757, Bush et. al "Calibration of T2 oximetry MRI for subjects with sickle cell disease"
    # but uses modified subparameters calculated from Li+Donahue data
    # note the scaling coef is 0.7, not 0.8
    # equation is 1/T2 = A * (1-y)^2 + B
    
    ay1hba = 77.5;
    ay2hba = 27.8;
    ay1hbs = 227.8;
    ay2hbs = 11.8;
    
    be1hba = 6.95;
    be2hba = 2.34;
    be1hbs = -23.9;
    be2hbs = 13.0;
    
    ay1 = (1 - (hbs/0.7)) * ay1hba + (hbs/0.7)*ay1hbs;
    ay2 = (1 - (hbs/0.7)) * ay2hba + (hbs/0.7)*ay2hbs;
    be1 = (1 - (hbs/0.7)) * be1hba + (hbs/0.7)*be1hbs;
    be2 = (1 - (hbs/0.7)) * be2hba + (hbs/0.7)*be2hbs;
    
    A = ay1*hct + ay2;
    B = be1*hct + be2;
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)


def libushdonahue_mixture_oxsat(T2, hct, hbs, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1

    # based off of DOI: 10.1002/mrm.28757, Bush et. al "Calibration of T2 oximetry MRI for subjects with sickle cell disease"
    # but uses modified subparameters calculated from Li+Bush+Donahue data
    # note the scaling coef is 0.7, not 0.8
    # equation is 1/T2 = A * (1-y)^2 + B
    
    ay1hba = 77.5;
    ay2hba = 27.8;
    ay1hbs = 170.8;
    ay2hbs = 25.5;
    
    be1hba = 6.95;
    be2hba = 2.34;
    be1hbs = -15.8;
    be2hbs = 10.9;
    
    ay1 = (1 - (hbs/0.7)) * ay1hba + (hbs/0.7)*ay1hbs;
    ay2 = (1 - (hbs/0.7)) * ay2hba + (hbs/0.7)*ay2hbs;
    be1 = (1 - (hbs/0.7)) * be1hba + (hbs/0.7)*be1hbs;
    be2 = (1 - (hbs/0.7)) * be2hba + (hbs/0.7)*be2hbs;
    
    A = ay1*hct + ay2;
    B = be1*hct + be2;
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)

def li_mixture_oxsat(T2, hct, hbs, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1

    # based off of DOI: 10.1002/mrm.28757, Bush et. al "Calibration of T2 oximetry MRI for subjects with sickle cell disease"
    # but uses modified subparameters calculated from Li only
    # equation is 1/T2 = A * (1-y)^2 + B
    
    ay1hba = 77.5;
    ay2hba = 27.8;
    ay1hbs = 231.2;
    ay2hbs = 11.9;
    
    be1hba = 6.95;
    be2hba = 2.34;
    be1hbs = -6.9;
    be2hbs = 8.8;
    
    ay1 = (1 - (hbs/0.8)) * ay1hba + (hbs/0.8)*ay1hbs;
    ay2 = (1 - (hbs/0.8)) * ay2hba + (hbs/0.8)*ay2hbs;
    be1 = (1 - (hbs/0.8)) * be1hba + (hbs/0.8)*be1hbs;
    be2 = (1 - (hbs/0.8)) * be2hba + (hbs/0.8)*be2hbs;
    
    A = ay1*hct + ay2;
    B = be1*hct + be2;
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)


def donahue_mixture_oxsat(T2, hct, hbs, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1

    # based off of DOI: 10.1002/mrm.28757, Bush et. al "Calibration of T2 oximetry MRI for subjects with sickle cell disease"
    # but uses modified subparameters calculated from Donahue only
    # note the scaling coef is 0.7, not 0.8
    # equation is 1/T2 = A * (1-y)^2 + B
    
    ay1hba = 77.5;
    ay2hba = 27.8;
    ay1hbs = 210.5;
    ay2hbs = 14.7;
    
    be1hba = 6.95;
    be2hba = 2.34;
    be1hbs = -65.5;
    be2hbs = 23.0;
    
    ay1 = (1 - (hbs/0.7)) * ay1hba + (hbs/0.7)*ay1hbs;
    ay2 = (1 - (hbs/0.7)) * ay2hba + (hbs/0.7)*ay2hbs;
    be1 = (1 - (hbs/0.7)) * be1hba + (hbs/0.7)*be1hbs;
    be2 = (1 - (hbs/0.7)) * be2hba + (hbs/0.7)*be2hbs;
    
    A = ay1*hct + ay2;
    B = be1*hct + be2;
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)


def generic_oxsat(ays, bees, T2, hct=None, tcpmg=None):
    #T2 in seconds
    #hct as float between 0 and 1
    #hbs as float between 0 and 1
    
    # equation is 1/T2 = A * (1-y)^2 + B
    # if ays is a float, A = ays
    # if ays is a list with len 2, A = ays[0]*hct+ays[1]
    # same for bees and B
    
    if type(ays) is list:
        A = ays[0]*hct + ays[1];
        B = bees[0]*hct + bees[1];
    else:
        A = ays
        B = bees
        
    #eqn = (1/T2) == ((A*(1-y)^2)  + B) -> solve for Y

    
    y = sym.symbols('y')
    solution = sym.solveset(((A*(1-y)**2)  + B) - (1/T2))
    Y = list(solution)[0]
    #eqn = eqn.subs(x, 2)
    
    #take the first root only
    #Y = 1 - np.sqrt((1/T2 - B)/A)
    
    return float(Y)




