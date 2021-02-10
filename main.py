from math import *
import numpy as np
from numpy.lib.function_base import kaiser
import matplotlib.pyplot as plt

class ProjectFailed(Exception):
    '''Custom Project Exception that can be thrown for error analysis'''
    pass

DtR = pi/180 # degrees to radians
RtD = 1/DtR # radians to degrees

def thetaTBM(M,beta,g):
    '''finds flow deflection angle in degrees for TBM shock given Mach number, shock angle, and ratio of specific heats'''
    numerator = (M*sin(DtR*beta))**2-1
    denominator = 0.5*tan(DtR*beta)*(M**2*(g+cos(2*DtR*beta))+2)
    theta = RtD*atan(numerator/denominator)
    return theta

def MTBM(theta,beta,g):
    '''finds Mach number for TBM shock given flow deflection angle, shock angle, and ratio of specific heats'''
    numerator = -2*(1+tan(DtR*theta)*tan(DtR*beta))
    denominator = tan(DtR*theta)*tan(DtR*beta)*(g+cos(2*DtR*beta))-2*(sin(DtR*beta))**2
    Mach = sqrt(numerator/denominator)
    return Mach

def betaMax(M,g):
    '''finds maximum shock angle given Mach number and ratio of specific heats'''
    gp1 = g+1 # gamma plus 1
    gm1 = g-1 # gamma minus 2
    # splitting up of beta_max equation
    isissq = gp1*(1+gm1*M*M/2+gp1/16*M**4)
    issq = 1/(g*M*M)*(gp1*M*M/4+sqrt(isissq)-1)
    bM = RtD*asin(sqrt(issq)) # beta max in deg
    return bM

#def betaMin(M,g):

def thetaMax(bM,M,g):
    '''finds maximum flow deflection angle given max shock angle, Mach number, and ratio of specific heats'''
    numerator = (M*M*(sin(DtR*bM))**2-1)/tan(DtR*bM)
    denominator = M*M*(g+1)/2-M*M*(sin(DtR*bM))**2+1
    tM = RtD*atan(numerator/denominator) # theta max in deg
    return tM
    

#def betaTBM(theta,M,g):


def TBMshockprop(theta,beta,M,g):
    '''finds all shock Mach numbers and thermo properties given T, B, and M; returns Mn1,Mn2,M2,p2/p1,T2/T1'''
    gp1 = g+1 # gamma plus 1
    gm1 = g-1 # gamma minus 2
    Mn1 = M*sin(DtR*beta)
    r2r1 = gp1*Mn1*Mn1/(gm1*Mn1*Mn1+2) # density ratio
    p2p1 = 1+2*g/gp1*(Mn1*Mn1-1) # pressure ratio
    T2T1 = p2p1/r2r1 # temperature ratio
    Mn2 = sqrt((Mn1*Mn1+2/gm1)/(2*g/gm1*Mn1*Mn1-1))
    M2 = Mn2/sin(DtR*(beta-theta))
    return Mn1,Mn2,M2,p2p1,T2T1




def main() -> int:
    gams = 1.4
    print("The flow deflection angle at Mach 6.7 that gives a shock angle of 20 deg is "+str(round(thetaTBM(6.7,20,gams),4))+" deg.")    
    print("The Mach number for a flow deflection angle of 6 deg and shock angle of 10 deg is "+str(round(MTBM(6,10,gams),4))+".")

    choice = input("What do you want to do?\n  (1) Find Theta (2) Find Beta (3) Find M? ")
    if choice == "1":
        M = float(input("Enter a value for Mach Number: "))
        bM = betaMax(M,gams)
        print("Maximum beta is "+str(round(bM,2))+"degrees\n")
        beta = float(input("Enter a value of Beta between these 2 values: "))
        theta = thetaTBM(M,beta,gams)
        Mn1,Mn2,M2,p2p1,T2T1 = TBMshockprop(theta,beta,M,gams)
        print("**************** Results ****************\n----- Given -----\n")
        print("Mach number: "+str(round(M,4))+"\n")
        print("Shock Angle: "+str(round(beta,4))+"\n\n")
        print("----- Results -----\n")
        print("Flow Deflection Angle: "+str(round(theta,4))+"\n")
        print("M1 normal: "+str(round(Mn1,4))+"\n")
        print("M2 normal: "+str(round(Mn2,4))+"\n")
        print("M2: "+str(round(M2,4))+"\n")
        print("P2/P1: "+str(round(p2p1,4))+"\n")
        print("T2/T1: "+str(round(T2T1,4))+"\n")
        print("**************** Done ****************\n")
    #elif choice == "2":

    elif choice == "3":
        beta = float(input("Enter a value of Beta: "))
        theta = float(input("Enter a value of Theta: "))
        M = MTBM(theta,beta,gams)
        Mn1,Mn2,M2,p2p1,T2T1 = TBMshockprop(theta,beta,M,gams)
        print("**************** Results ****************\n----- Given -----\n")
        print("Shock Angle: "+str(round(beta,4))+"\n\n")
        print("Flow Deflection Angle: "+str(round(theta,4))+"\n")
        print("----- Results -----\n")
        print("Mach number: "+str(round(M,4))+"\n")
        print("M1 normal: "+str(round(Mn1,4))+"\n")
        print("M2 normal: "+str(round(Mn2,4))+"\n")
        print("M2: "+str(round(M2,4))+"\n")
        print("P2/P1: "+str(round(p2p1,4))+"\n")
        print("T2/T1: "+str(round(T2T1,4))+"\n")
        print("**************** Done ****************\n")

    return 1




if __name__ == "__main__":
    if main() != 1:
        raise ProjectFailed("MAE 343 Project returned bad exit code!")



