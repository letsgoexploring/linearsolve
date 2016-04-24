from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import pandas as pd
import sys

''' This program defines a klein object class and provides functions for solving and creating impulse 
responses for linear dynamic models in the form of:
    
    a*Et[x(t+1)] = b*x(t) + c*z(t)

with x(t) = [k(t); u(t)] where k(t) is a vector of predetermined (state) variables and u(t) is a vector 
of nonpredetermined (forward-looking) variables. z(t) is a vector of exogenous driving variables with 
autocorrelation matrix rho. The solution to the model is a set of matrices f, n, p, l such that:

    u(t)   = f*k(t) + n*z(t)
    k(t+1) = p*k(t) + l*z(t).

The solution algorithm is based on Klein (2000) and is based on his solab.m Matlab program.'''

class klein:

    def __init__(self,a,b,c,rho,nk,endoNames=None,exoNames=None): #,shockNames=None):
        ''' Accepts as inputs the matrices of the linear system a, b, and c, the autocorrelation matrix of
        the exogenous driving processes, and the number of state variables nk. 

        Upon initialization, the model is solved. The solution matrices f, p, n, and l are returned as 
        attributes of the klein object. The solution program also gives the klein two other useful
        attributes:
            eig:  the eigenvalues from the genearlized Schur decomposition
            stab: equals 0 when the system has the correct number of stable eigenvalues, -1 when the system 
                    has too few stable eigenvalues, and 1 when the system has too many.
        '''
        
        f,n,p,l,stab,eig = self.solve(a,b,c,rho,nk)
        self.a   = a
        self.b   = b
        self.c   = c
        self.rho = rho
        self.f   = f
        self.n   = n
        self.p   = p
        self.l   = l
        self.stab= stab
        self.eig = eig
        self.nu = np.shape(self.n)[0]
        self.nk = np.shape(self.l)[0]
        self.nz = np.shape(self.l)[1]

        names = {}
        if endoNames == None or len(endoNames) != self.nu+self.nk:
            endoNames = []
            for i in range(nk):
                endoNames.append('state '+str(i+1))
            for i in range(self.nu):
                endoNames.append('forward '+str(i+1))

        if exoNames == None or len(exoNames) != self.nz:
            exoNames = []
            for i in range(self.nz):
                exoNames.append('exo '+str(i+1))

        # if shockNames == None or len(shockNames) != np.shape(c)[1]:
        #     shockNames = []
        #     for i in range(np.shape(c)[1]):
        #         shockNames.append('shock '+str(i+1))

        names['endo'] = endoNames
        names['exo'] = exoNames
        # names['shocks'] = shockNames

        self.names = names


    def solve(self,a,b,c,rho,nk):

        '''Method solves the linear model.'''

        s, t, alpha, beta, q, z = la.ordqz(A=a, B=b, sort='ouc')
        
        q=np.mat(q)
        z=np.mat(z)
        s=np.mat(s)
        t=np.mat(t)
        a=np.mat(a)
        b=np.mat(b)
        c=np.mat(c)
        rho=np.mat(rho)

        z11 = z[0:nk,0:nk]
        z12 = z[0:nk,nk:]
        z21 = z[nk:,0:nk]
        z22 = z[nk:,nk:]

        nu = np.shape(a)[0] - nk                # number of nonpredetermined variables
        nz = np.shape(c)[1]                     # number of forcing variables

        if nk>0:
            if np.linalg.matrix_rank(z11)<nk:
                sys.exit("Invertibility condition violated.")

        s11 = s[0:nk,0:nk];
        if nk>0:
            z11i = la.inv(z11)
            s11i = la.inv(s11)
        else:
            z11i = z11
            s11i = s11


        s12 = s[0:nk,nk:]
        s22 = s[nk:,nk:]
        t11 = t[0:nk,0:nk]
        t12 = t[0:nk,nk:]
        t22 = t[nk:,nk:]
        q1  = q[0:nk,:]
        q2  = q[nk:,:]

        # Verify that there are exactly nk stable (inside the unit circle) eigenvalues:

        stab = 0

        if nk>0:
            if np.abs(t[nk-1,nk-1])>np.abs(s[nk-1,nk-1]):
                print('Warning: Too few stable eigenvalues.')
                stab = -1

        if np.abs(t[nk,nk])<np.abs(s[nk,nk]):
            print('Warning: Too many stable eigenvalues.')
            stab = 1

        # Compute the generalized eigenvalues

        tii = np.diag(t)
        sii = np.diag(s)
        eig = np.zeros(np.shape(tii))

        for k in range(len(tii)):
            if np.abs(sii[k])>0:
                eig[k] = np.abs(tii[k])/np.abs(sii[k])
            else:
                eig[k] = np.inf

        mat1 = np.kron(np.transpose(rho),s22) - np.kron(np.identity(nz),t22)
        mat1i = la.inv(mat1)
        q2c = q2.dot(c)
        vecq2c = q2c.flatten(1).T
        vecm = mat1i.dot(vecq2c)
        m = np.transpose(np.reshape(np.transpose(vecm),(nz,nu)))
        if nk>0:
            dyn = np.linalg.solve(s11,t11)
        else:
            dyn = np.array([])


        n = np.real((z22 - z21.dot(z11i).dot(z12)).dot(m))
        l = np.real(-z11.dot(s11i).dot(t11).dot(z11i).dot(z12).dot(m) + z11.dot(s11i).dot(t12.dot(m) - s12.dot(m).dot(rho)+q1.dot(c)) + z12.dot(m).dot(rho))
        f = np.real(z21.dot(z11i))
        p = np.real(z11.dot(dyn).dot(z11i))

        return np.array(f),np.array(n),np.array(p),np.array(l),stab,eig


    def impulse(self,T=15,t0=1,shock=None):

        ''' Method for computing impulse responses for shocks arriving at date t0. t0 must be greater than 0.
        shock is an (nz x 1) list of shock values. If shock==None, then shock is set to a vector of ones.
        
        Returns a dictionary containing pandas dataframes. The dictionary has the form:

            self.irs['shock name']['endog var name']

        '''

        irsDict = {}

        nu = self.nu
        nk = self.nk
        nz = self.nz

        self.names['endo'] = self.names['endo']
        self.names['exo'] = self.names['exo']        

        # for j in range(nz):
        for j,name in enumerate(self.names['exo']):
            eps= np.zeros((T+1,nz))
            if shock==None:
                eps[t0-1][j] = 0.01
            else:
                eps[t0-1][j] = shock[j]
            k  = np.zeros((T+2,nk))
            z  = np.zeros((T+2,nz))
            u  = np.zeros((T+2,nu))

            x = ir(self.f,self.n,self.p,self.l,self.rho,z,k,u,eps)
            
            frameDict = {self.names['exo'][j]:x[j]}
            for i,endoName in enumerate(self.names['endo']):
                frameDict[endoName] = x[i+len(self.names['exo'])]

            irFrame = pd.DataFrame(frameDict, index = np.arange(T+1))
            irsDict[self.names['exo'][j]] = irFrame


        self.irs = irsDict

    def stochSim(self,T=50,dropFirst=100,covMat=None,seed=None):
        
        ''' Method for computing a stochastic simulation of the solved model. The simulation period is T 
        periods and the first 'dropFirst' periods aer excluded form the simulation. Shock realizations are
        drawn from the multivariate normal distribution with mean 0 and covariance matrix given. If the 
        covariance matrix is not given, it is set to the nz X nz identity matrix.

        'seed': integersets the seed 

        Returns a pandas dataframes with column names given by the names of the endogenous variables.

            self.sims['endog var name']

        '''

        simDict = {}

        nu = self.nu
        nk = self.nk
        nz = self.nz

        self.names['endo'] = self.names['endo']
        self.names['exo'] = self.names['exo']

        if covMat == None:
            covMat = 0.01*np.eye(nz)

        if seed!=None and type(seed)==int:

            np.random.seed(seed)

        epsSim = np.random.multivariate_normal(mean=np.zeros(nz),cov=covMat,size=[dropFirst+T+1])
        eps= epsSim
        if nz == 1:
            k  = np.zeros((dropFirst+T+2,nk))
            z  = np.zeros((dropFirst+T+2,nz))
            u  = np.zeros((dropFirst+T+2,nu))

            x = ir(self.f,self.n,self.p,self.l,self.rho,z,k,u,eps)

            frameDict = {self.names['exo'][0]:x[0][dropFirst:]}
            for i,endoName in enumerate(self.names['endo']):
                frameDict[endoName] = x[i+len(self.names['exo'])][dropFirst:]

            simFrame = pd.DataFrame(frameDict, index = np.arange(T+1))

        else:
            k  = np.zeros((dropFirst+T+2,nk))
            z  = np.zeros((dropFirst+T+2,nz))
            u  = np.zeros((dropFirst+T+2,nu))

            x = ir(self.f,self.n,self.p,self.l,self.rho,z,k,u,eps)

            frameDict={}
            
            for j,exoName in enumerate(self.names['exo']):
                frameDict[exoName] = x[j][dropFirst:]
            for i,endoName in enumerate(self.names['endo']):
                frameDict[endoName] = x[i+len(self.names['exo'])][dropFirst:]

            simFrame = pd.DataFrame(frameDict, index = np.arange(T+1))

        self.sims = simFrame


def ir(f,n,p,l,rho,z,k,u,eps):
    ''' Simulates a model given matrices and shock processes.

    u(t)   = f*k(t) + n*z(t)
    k(t+1) = p*k(t) + l*z(t)

    '''
    for i,e in enumerate(eps):
        rho = np.array(rho)
        z[i+1] = rho.dot(z[i]) + e
        k[i+1] = p.dot(k[i]) + l.dot(z[i])
        u[i]   = f.dot(k[i]) + n.dot(z[i])
    z = z[0:-1]
    k = k[0:-1]
    u = u[0:-1]
    return np.concatenate((z.T,k.T,u.T))

def shift(series,direction=None):
    values = series.values
    if direction=='forward':
        print('back')
        newValues = np.append(0,values[0:-1])
        print(newValues)
    elif direction=='backward':
        newValues = np.append(values[1:],np.nan)
    else:
        newSeries = pd.Series(newValues,index = series.index)
    return newValues

def trim(frame,lengthBegin=0,lengthEnd=0):
    if lengthBegin==0 and lengthEnd==0:
        return frame
    else:
        frame = frame.iloc[lengthBegin:-lengthEnd]
    return frame