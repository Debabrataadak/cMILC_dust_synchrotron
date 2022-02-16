import numpy as np
import scipy.constants as constants
import math

bita_s=-3.0
bita_d=1.53
Td=19.4
Tcmb=2.7255
xdd = constants.h * 353.0 * 1.e9 / constants.k / Td
###planck function ###
def B(nu):
    xd = constants.h * nu * 1.e9 / constants.k / Td
    return((nu/353.0)**(bita_d + 1.0)*(np.expm1(xdd)/np.expm1(xd)))
def planck(nu,T,bd):
    Temp_d=np.copy(T)
    spec_d=np.copy(bd)
    xd = constants.h * nu * 1.e9 / constants.k / Temp_d
    return((nu/353.0)**(spec_d + 1.0)*(np.expm1(xdd)/np.expm1(xd)))
############
frequency = np.array([23,30,33,41,44,61,70,94,100,143,217,353])
#frequency = np.array([23,33,41,61,94,30,44,70,100,143,217,353])
#frequency = np.array([28,44,70,100,143,217,353])
no_bands=np.size(frequency)

#####################



def Mixingmatrix(idn):
    if (idn=="01"):
                ncr=2
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        d=np.log(frequency[i]/353.0)
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        #f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-2.0)*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2

    if (idn=="02"):
                ncr=2
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        #F[i][1]=1.0
                        #F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][1] = (frequency[i]/30.0)**(bita_s)
    if (idn=="03"):
                ncr=3
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        #F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
    if (idn=="04"):
                ncr=3
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb 
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        #F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        #F[i][1] = (frequency[i]/30.0)**(bita_s)
                        F[i][2] = d*c
                        #F[i][5] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c

    if (idn=="05"):
                ncr=4
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        #F[i][4] = d*c
                        #F[i][5] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c

    if (idn=="06"):
                ncr=4
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        #F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][3] = d*c
                        #F[i][5] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c


    if (idn=="07"):
                ncr=5
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands): 
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        #F[i][5] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c
    if (idn=="10"):
                ncr=6
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c


    if (idn=="11"):
                ncr=7
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        #F[i][8] = d*e*c
    if (idn=="15"):
                ncr=8
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        #F[i][8] = d*e*c

    if (idn=="20"):
                ncr=9
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] =xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        F[i][8] = d*e*c


    if (idn=="08"):
                ncr=5
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        #F[i][4] = d*c
                        F[i][4] = e*c
                        #F[i][6] = f*e*c
                        #F[i][8] = d*e*c
    if (idn=="09"):
                ncr=5
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1]=xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        #F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][3] = d*c
                        F[i][4] = e*c
                        #F[i][3] = f*e*c
                        #F[i][8] = d*e*c

    if (idn=="12"):
                ncr=7
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        #F[i][8] = d*e*c
    if (idn=="13"):
                ncr=7
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*d*c
                        #F[i][8] = d*e*c

    if (idn=="14"):
                ncr=7
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        #F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*e*c
                        #F[i][8] = d*e*c
    if (idn=="16"):
                ncr=8
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*d*c
                        #F[i][8] = d*e*c
    if (idn=="17"):
                ncr=8
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*e*c
                        #F[i][8] = d*e*c
    if (idn=="18"):
                ncr=8
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = d*d*c
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        #F[i][8] = d*e*c


    if (idn=="19"):
                ncr=8
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] = xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = d*d*c
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*e*c
                        #F[i][8] = d*e*c
    if (idn=="21"):
                ncr=9
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] =xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = d*d*c
                        F[i][8] = d*e*c
    if (idn=="22"):
                ncr=9
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] =xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        F[i][8] = d*d*c
    if (idn=="23"):
                ncr=9
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] =xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = d*d*c
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        F[i][8] = d*e*c
    if (idn=="24"):
                ncr=10
                F=np.zeros((no_bands,ncr),dtype=np.float64)
                for i in range(no_bands):
                        xcmb=constants.h * frequency[i] * 1.e9 / constants.k / Tcmb
                        c=B(frequency[i])
                        d=np.log(frequency[i]/353.0)
                        xd = constants.h * frequency[i] * 1.e9 / constants.k / Td
                        e=(xd/Td)*np.exp(xd)/np.expm1(xd) - (xdd/Td)*np.exp(xdd)/np.expm1(xdd)
                        f=(xd*np.cosh(xd/2.0)/np.sinh(xd/2.0)-xdd*np.cosh(xdd/2.0)/np.sinh(xdd/2.0))*(1.0/Td)
                        F[i][0] = c
                        F[i][1] =xcmb**2*np.exp(xcmb)/np.expm1(xcmb)**2
                        F[i][3] = np.log(frequency[i]/30.0)*(frequency[i]/30.0)**(bita_s)
                        F[i][7] = (np.log(frequency[i]/30.0))**2*(frequency[i]/30.0)**(bita_s)
                        F[i][2] = (frequency[i]/30.0)**(bita_s)
                        F[i][4] = d*c
                        F[i][5] = e*c
                        F[i][6] = f*e*c
                        F[i][8] = d*e*c
                        F[i][9] = d*d*c
    return(F)


