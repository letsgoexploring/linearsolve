from __future__ import division,print_function
import numpy as np
import scipy.linalg as la
from statsmodels.tools.numdiff import approx_fprime_cs
from scipy.optimize import root,fsolve,broyden1,broyden2
import pandas as pd
import sys

class model:

    '''This program defines a class -- linearsolve.model -- with associated methods for solving and 
    simulating dynamic stochastic general equilibrium (DSGE) models.'''

    def __init__(self,equations=None,nstates=None,varNames=None,shockNames=None,parameters=None,parameterNames=None):
        
        '''Initializing an instance linearsolve.model requires values for the following variables:

        Args:
            equations:       (fun) A function that represents the equilibirum conditions for a DSGE model.
                                The function should accept three arguments:
                                    * vars_fwd:     endogenous variables dated t+1
                                    * vars_cur:     endogenous variables dated t
                                    * parameters:   the parameters of the model
                                The function should return an n-dimensional array. Each element of the returned
                                array is a list with the first element equaling the left-hand side of a single
                                equilibrium condition and the second elemtn equals the right-hand side.
            nstates:         (int) The number of state variables in the model.
            varNames:        (list) A list of strings with the names of the endogenous variables. The state
                                variables must be ordered first.
            shockNames:      (list) A list of strings with the names of the exogenous shocks to each state
                                variable. The order of names must agree with varNames.
            parameters:      (list or Pandas Series) Either a list of parameter values OR a Pandas Series object
                                with parameter name strings as the index.
            parameterNames:  (list) Optional. If parameters is given as a list,then this list of strings will
                                be used to save the parameters with names as a Pandas Series object.

        Returns:
            None

        Attributes:
            nvars:              (int) The number of variables in the model.
            nstates:            (int) The number of state variables in the model.
            ncostates:          (int) The number of costate or control variables in the model.
            parameters:         (Pandas Seris) A Pandas Series with parameter name strings as the 
                                    index. If parameterNames wasn't supplied,then parameters are labeled 
                                    'parameter 1',parameter2',etc.
            names:              (dict) A dictionary with keys 'variables','shocks',and 'param' that
                                    stores the names of the model's variables and parameters.
            equilibrium_fun:    (fun) Function that returns the equilibrium comditions of the model
            ss:                 (None)

        '''
        
        self.equilibrium_fun= equations
        self.nvars = len(varNames)
        self.nstates = nstates
        self.costates=self.nvars-nstates

        names = {}

        names['variables'] = varNames

        if shockNames is not None:
            if len(shockNames)<self.nstates:
                shockNames_temp = []
                for i in range(self.nstates):
                    try:
                        shockNames_temp.append(shockNames[i])
                    except:
                        shockNames_temp.append('e_'+varNames[i])
                shockNames = shockNames_temp

        
        else:
            shockNames = []
            for i in range(self.nstates):
                shockNames.append('e_'+varNames[i])

        names['shocks'] = shockNames

        if isinstance(parameters,pd.Series):
            self.parameters = parameters

        else:
            if parameterNames is None:
                parameterNames = ['parameter '+str(i+1) for i in range(len(parameters))]

            self.parameters = pd.Series(parameters,index=parameterNames)

        names['param'] = parameterNames
        self.names = names
        self.steady_state = None



    def compute_ss(self,guess=None,method='fsolve',options={}):

        '''This method attempts to solve for the steady state of the model.

        Args:
            guess:      (Pandas Series,Numpy array,or list) An initial guess for the 
                            steady state solution. The result is highly sensisitve to the intial 
                            guess chosen,so be careful. If the guess is a Numpy ndarray or a list
                            then the elements must be ordered to conform with self.names['variables'].
            method:     (str) The function from the Scipy library to use. Your choices are:
                        a. root
                        b. fsolve (default)
                        c. broyden1
                        d. broyden2
            options:    (dict) A dictionary of optional arguments to pass to the numerical solver.
                            Check out the Scipy documentation to see the options available for each routine:
                                http://docs.scipy.org/doc/scipy/reference/optimize.html

        Returns:
            None

        Attributes:
            ss: (Pandas Series) Steady state values of endogenous variables

            '''

        if guess is None:
            guess = np.ones(self.nvars)
        else:
            try:
                guess = guess[self.names['variables']]
            except:
                pass


        # Create function for nonlinear solver
        def ss_fun(variables):

            variables = pd.Series(variables,index = self.names['variables'])

            return self.equilibrium_fun(variables,variables,self.parameters)

        if method == 'fsolve':
            steady_state =fsolve(ss_fun,guess,**options)

        elif method == 'root':
            steady_state =root(ss_fun,guess,**options)['x']

        elif method == 'broyden1':
            steady_state =broyden1(ss_fun,guess,**options)

        elif method == 'broyden2':
            steady_state =broyden2(ss_fun,guess,**options)

        # Add ss attribute
        self.ss = pd.Series(steady_state,index=self.names['variables'])



    def set_ss(self,steady_state):

        '''This method allows the user to directly set the steady state of the model Uses Numpy.isclose()
        to print whether each steady state equilibrium condition evaluates to something close to zero.

        Args:
            steady_state:   (Pandas Series,Numpy array,or list)

        Returns:
            None

        Attributes:
            ss: (Pandas Series) Steady state values of endogenous variables

        '''

        try:
            self.ss = steady_state[self.names['variables']]
        except:

            self.ss = pd.Series(steady_state,index=self.names['variables'])


    def check_ss(self):

        '''This method allows the user to directly set the steady state of the model.

        Args:
            None

        Returns:
            None

        Attributes:
            None

        '''

        try:
            print(np.isclose(self.equilibrium_fun(self.ss,self.ss,self.parameters),0))
        except:
            print('Set the steady state first.')

    def linear_approximation(self,steady_state=None):

        ''' Given a nonlinear rational expectations model in the form:

                    psi_1[x(t+1),x(t)] = psi_2[x(t+1),x(t)]

            this method returns the linear approximation of the model with matrices a and b such that:

                    a * y(t+1) = b * y(t)

            where y(t) = x(t) - x is the log deviation of the vector x from its steady state value.

        Args:
            steady_state:   (Pandas Series or numpy array)

        Returns:
            None

        Attributes:
            loglinear:      (bool) Whether the model is loglinear. Sets to False.
            a:              (Numpy ndarray)
            b:              (Numpy ndarray)

        '''

        # Set loglinear attribute
        self.loglinear=False

        # Warn if steady state attribute ss has not been assigned
        if steady_state is None:

            try:
                steady_state = self.ss
            except :
                raise ValueError('You must specify a steady state for the model before attempting to linearize.')

        # Compute approximation
        def equilibrium(vars_fwd,vars_cur):

            vars_fwd = pd.Series(vars_fwd,index = self.names['variables'])
            vars_cur = pd.Series(vars_cur,index = self.names['variables'])

            equilibrium_left = self.equilibrium_fun(vars_fwd,vars_cur,self.parameters)
            equilibrium_right = np.ones(len(self.names['variables']))

            return equilibrium_left - equilibrium_right

        equilibrium_fwd = lambda fwd: equilibrium(fwd,steady_state)
        equilibrium_cur = lambda cur: equilibrium(steady_state,cur)

        # Assign attributes
        self.a= approx_fprime_cs(steady_state.ravel(),equilibrium_fwd)
        self.b= -approx_fprime_cs(steady_state.ravel(),equilibrium_cur)


    def log_linear_approximation(self,steady_state=None):

        ''' Given a nonlinear rational expectations model in the form:

                    psi_1[x(t+1),x(t)] = psi_2[x(t+1),x(t)]

            this method returns the log-linear approximation of the model with matrices a and b such that:

                    a * y(t+1) = b * y(t)

            where y(t) = log x(t) - log x is the log deviation of the vector x from its steady state value.

        Args:
            steady_state:   (Pandas Series or numpy array)

        Returns:
            None

        Attributes:
            loglinear:      (bool) Whether the model is loglinear. Sets to True.
            a:              (Numpy ndarray)
            b:              (Numpy ndarray)

        '''

        # Set loglinear attribute
        self.loglinear=True

        # Warn if steady state attribute ss has not been assigned
        if steady_state is None:

            try:
                steady_state = self.ss
            except :
                raise ValueError('You must specify a steady state for the model before attempting to linearize.')

        # Compute approximation
        def log_equilibrium(log_vars_fwd,log_vars_cur):

            log_vars_fwd = pd.Series(log_vars_fwd,index = self.names['variables'])
            log_vars_cur = pd.Series(log_vars_cur,index = self.names['variables'])

            equilibrium_left = self.equilibrium_fun(np.exp(log_vars_fwd),np.exp(log_vars_cur),self.parameters)+1
            equilibrium_right = np.ones(len(self.names['variables']))

            return np.log(equilibrium_left) - np.log(equilibrium_right)

        log_equilibrium_fwd = lambda log_fwd: log_equilibrium(log_fwd,np.log(steady_state))
        log_equilibrium_cur = lambda log_cur: log_equilibrium(np.log(steady_state),log_cur)

        # Assign attributes
        self.a= approx_fprime_cs(np.log(steady_state).ravel(),log_equilibrium_fwd)
        self.b= -approx_fprime_cs(np.log(steady_state).ravel(),log_equilibrium_cur)


    def solve_klein(self,a=None,b=None):

        '''Method solves a linear rational expectations model of the form:

                a * x(t+1) = b * x(t) + e(t)

        where z(t) is a VAR(1) exogenous forcing process with autocorrelation matrix rho.

        The method returns the solution to the law of motion:

                u(t)   = f*s(t) + e(t)
                s(t+1) = p*s(t)

        Args:
            a:  (Numpy ndarray) coefficient matrix
            b:  (Numpy ndarray) coefficient matrix

        Returns:
            None

        Attributes:
            f:      (Numpy ndarray) Solution matrix coeffients on s(t)
            p:      (Numpy ndarray) Solution matrix coeffients on s(t)
            stab:   (int) Indicates solution stability and uniqueness
                        stab == 1: too many stable eigenvalues
                        stab == -1: too few stable eigenvalues
                        stab == 0: just enoughstable eigenvalues
            eig:    The generalized eigenvalues from the Schur decomposition

         '''

        if a is None and b is None:
            
            a = self.a
            b = self.b

        self.f,n,self.p,l,self.stab,self.eig = klein(a=a,b=b,c=None,rho=None,nstates=self.nstates)


    def impulse(self,T=51,t0=1,shocks=None,percent=False,diff=True):

        ''' Method for computing impulse responses for shocks to each state variable.

        Arguments:
                T:          (int) Number of periods to simulate. Default: 51
                t0:         (int) Period in which the shocks are realized. May be equal to 0. Default: 1
                shocks:     (list or Numpy array) An (ns x 1) list of shock values. If shocks==None,shocks is set to a vector of 0.01s. Default = None
                percent:    (bool) Whether to multiply simulated values by 100. Only works for log-linear approximations. Default: False
                diff:       (bool) Subtract steady state for linear approximations (or log steady state for log-linear approximations). Default: True
        
        Returns
            None

        Attributes:
            irs:    (dict) A dictionary containing Pandas DataFrames. Has the form: self.irs['shock name']['endog var name']

        '''

        # Initialize dictionary
        irsDict = {}

        # Set numbers of costate and state variables
        ncostates = self.costates
        nstates = self.nstates
        
        # iterate over all shocks, compute impulse responses, and add results to dictionary
        for j,name in enumerate(self.names['shocks']):

            s0 = np.zeros([1,nstates])
            eps= np.zeros([T,nstates])
            if shocks is not None:
                try:
                    eps[t0][j] = shocks[name]
                except:
                    try:
                        eps[t0][j] = shocks[j]
                    except:
                        if self.loglinear:
                            eps[t0][j] = 0.01
                        else:
                            eps[t0][j] = 1

            else:
                if self.loglinear:
                    eps[t0][j] = 0.01
                else:
                    eps[t0][j] = 1

            x = ir(self.f,self.p,eps,s0)

            frameDict = {self.names['shocks'][j]:eps.T[j]}
            for i,endoName in enumerate(self.names['variables']):
                if diff:
                    frameDict[endoName] = x[i]
                else:
                    if not self.loglinear:
                        frameDict[endoName] = x[i] + self.ss[endoName]
                    else:
                        frameDict[endoName] = x[i] + np.log(self.ss[endoName])

            irFrame = pd.DataFrame(frameDict,index = np.arange(T))

            if percent==True and self.loglinear:
                irFrame = 100*irFrame

            if shocks is None or len(shocks)>j:
                irsDict[self.names['shocks'][j]] = irFrame

        # Set attribute
        self.irs = irsDict

    def stoch_sim(self,T=51,dropFirst=300,covMat=None,seed=None,percent=False,diff=True):
        
        ''' Method for computing impulse responses for shocks to each state variable.

        Arguments:
                T:          (int) Number of periods to simulate. Default: 51
                dropFirst:  (int) Number of periods to simulate before saving output. Default: 300
                covMat:     (list or Numpy.ndarray) Covariance matrix shocks. If covMat is None,it's set to eye(nstates). Default: None
                seed:       (int) Sets the seed for the Numpy random number generator. Default: None
                percent:    (bool) Whether to multiply simulated values by 100. Only works for log-linear approximations. Default: False
                diff:       (bool) Subtract steady state for linear approximations (or log steady state for log-linear approximations). Default: True
        
        Returns
            None

        Attributes:
            simulated:    (Pandas DataFrame)

        '''

        # Set numbers of costate and state variables
        ncostates = self.costates
        nstates = self.nstates

        # Initialize states
        s0 = np.zeros([1,nstates])

        # Set covMat if not given
        if covMat is None:
            covMat = np.eye(nstates)

        # Set seed for the Numpy random number generator
        if seed is not None and type(seed)==int:

            np.random.seed(seed)

        # Simulate shocks
        eps = np.zeros([dropFirst+T,nstates])
        eps[:,:len(covMat)] = np.random.multivariate_normal(mean=np.zeros(len(covMat)),cov=covMat,size=[dropFirst+T])

        # Compute impulse responses given shocks
        x = ir(self.f,self.p,eps,s0)

        # Construct DataFrame
        frameDict = {}
        for j,exoName in enumerate(self.names['shocks']):
            frameDict[exoName] = eps.T[j][dropFirst:]
        for i,endoName in enumerate(self.names['variables']):
            if diff:
                frameDict[endoName] = x[i][dropFirst:]
            else:
                frameDict[endoName] = x[i][dropFirst:] + self.ss[endoName]

        simFrame = pd.DataFrame(frameDict,index = np.arange(T))
        if percent==True:
            simFrame = 100*simFrame

        # Assign attribute
        self.simulated = simFrame

    def approximate_and_solve(self,loglinear=True):

        '''Method approximates and solves a dynamic stochastic general equilibrium (DSGE) model by 
        constructing the log-linear approximation (if the model isn't log-linear) and solving the model
        using Klein's (2000) method. Arguments:

        Args:
            loglinear:  (bool) whether to compute log-linear or linear approximation. Default: True

        Returns:
            None

        Attributes:
            a:          (Numpy ndarray)
            b:          (Numpy ndarray)
            f:          (Numpy ndarray) Solution matrix coeffients on s(t)
            p:          (Numpy ndarray) Solution matrix coeffients on s(t)
            stab:       (int) Indicates solution stability and uniqueness
                            stab == 1: too many stable eigenvalues
                            stab == -1: too few stable eigenvalues
                            stab == 0: just enoughstable eigenvalues
            eig:        The generalized eigenvalues from the Schur decomposition
            loglinear:  (bool) Whether the model is loglinear. Sets to loglinear.

        '''

        # Set attribute
        self.loglinear = loglinear

        # Approximate
        if loglinear == True:
            self.log_linear_approximation()
        else:
            self.linear_approximation()

        # Solve the model
        self.solve_klein(self.a,self.b)

    def approximated(self,round=True,precision=4):

        '''Returns a string containing the log-linear approximation to the equilibrium conditions

        Args:
            round:       (bool) Whether to round the coefficents in the linear equations. Default: True
            precisions:  (int) Number of decimals to round the coefficients. Default: 4

        Returns:
            String with the log-linear approximation to the equilibrium conditions.

        Attributes:
            None

        '''

        if round is True:
            a = np.round(self.a,precision)
            b = np.round(self.b,precision)

            
        leftsides =  []
        rightsides = []
        if self.loglinear==True:
            lines ='Log-linear equilibrium conditions:\n\n'
        else:
            lines ='Linear equilibrium conditions:\n\n'


        left_length = 1


        for i in range(self.nvars):
            
            left = ''
            right = ''
            
            left_plus_flag = 0
            right_plus_flag = 0
            if all(np.isclose(0,a[i])):
                left += '0'
            else:
                for j in range(self.nvars):

                    if not np.isclose(0,a[i][j]):

                        name = self.names['variables'][j]

                        if j >self.nstates-1:
                            name +='[t+1|t]'

                        else:
                            name +='[t+1]'

                        if np.isclose(1,a[i][j]):
                            coeff = ''

                        elif np.isclose(-1,a[i][j]):
                            coeff = '-'

                        else:
                            coeff = str(a[i][j])+'·'

                        if left_plus_flag == 0:
                            left += coeff+name
                            left_plus_flag+=1

                        else:
                            if a[i][j] > 0:
                                left += '+'+coeff+name
                            else:
                                left += coeff+name
            
            if all(np.isclose(0,b[i])):
                right += '0'
            else:
                for j in range(self.nvars):

                    if not np.isclose(0,b[i][j]):

                        name = self.names['variables'][j]+'[t]'

                        if np.isclose(1,b[i][j]):
                            coeff = ''

                        elif np.isclose(-1,b[i][j]):
                            coeff = '-'

                        else:
                            coeff = str(b[i][j])+'·'

                        if right_plus_flag == 0:
                            right += coeff+name
                            right_plus_flag+=1

                        else:
                            if b[i][j] > 0:
                                right += '+'+coeff+name
                            else:
                                right += coeff+name
                                
                        
            leftsides.append(left)
            rightsides.append(right)
            
            if len(left)>left_length:
                left_length = len(left)
                

        for i in range(self.nvars):
            leftsides[i] = leftsides[i].rjust(left_length)
            lines+=leftsides[i]+' = '+rightsides[i]+'\n\n'

        lines = lines[:-2]

        return lines

    def solved(self,round=True,precision=4):

        '''Returns a string containing the solution to the log-linear system

        Args:
            round:       (bool) Whether to round the coefficents in the solution equations. Default: True
            precisions:  (int) Number of decimals to round the coefficients. Default: 4

        Returns:
            String with the log-linear approximation to the equilibrium conditions.

        Attributes:
            None

        '''

        if round is True:
            f = np.round(self.f,precision)
            p = np.round(self.p,precision)

            
        leftsides =  []
        rightsides = []
        if self.loglinear==True:
            lines ='Solution to the log-linear system:\n\n'
        else:
            lines ='Solution to the linear system:\n\n'

        left_length = 1


        for i in range(self.nstates):
            
            left = ''
            right = ''
            right_plus_flag = 0

            left+= self.names['variables'][i]+'[t+1]'
            
            if all(np.isclose(0,p[i])):
                right += self.names['shocks'][i]+'[t+1]'
            else:
                for j in range(self.nstates):

                    if not np.isclose(0,p[i][j]):

                        if right_plus_flag == 0:
                            right += str(p[i][j])+'·'+self.names['variables'][j]+'[t]'
                            right_plus_flag+=1

                        else:
                            if p[i][j] > 0:
                                right += '+'+str(p[i][j])+'·'+self.names['variables'][j]+'[t]'

                            else:
                                right += str(p[i][j])+'·'+self.names['variables'][j]+'[t]'
                right+='+'+self.names['shocks'][i]+'[t+1]'  
            leftsides.append(left)
            rightsides.append(right)
                            
                        
            if len(left)>left_length:
                left_length = len(left)

                
        for i in range(self.nvars-self.nstates):
            
            left = ''
            right = ''
            right_plus_flag = 0
            
            left+= self.names['variables'][self.nstates+i]+'[t]'
            
            if all(np.isclose(0,f[i])):
                right += '0'
            else:
                for j in range(self.nstates):
                    if not np.isclose(0,f[i][j]):

                        name = self.names['variables'][j]+'[t]'

                        if np.isclose(1,f[i][j]):
                            coeff = ''

                        elif np.isclose(-1,f[i][j]):
                            coeff = '-'

                        else:
                            coeff = str(f[i][j])+'·'

                        if right_plus_flag == 0:
                            right += coeff+name
                            right_plus_flag+=1

                        else:
                            if f[i][j] > 0:
                                right += '+'+coeff+name
                            else:
                                right += coeff+name
                                
            leftsides.append(left)
            rightsides.append(right)
                            
                        
            if len(left)>left_length:
                left_length = len(left)
                
                
        for i in range(self.nvars):
            leftsides[i] = leftsides[i].rjust(left_length)
            lines+=leftsides[i]+' = '+rightsides[i]+'\n\n'
        lines = lines[:-2]
                    

        return lines


def ir(f,p,eps,s0=None):

    ''' Function for simulating a model in the following form:

            u(t)   = f*s(t) + e(t)
            s(t+1) = p*s(t)

    where s(t) is an (nstates x 1) vector of state variables,u(t) is an (ncostate x 1) vector of costate
    variables,and e(t) is an (nstates x 1) vector of structural shocks.

    Args:
        f:      (Numpy ndarray) matrix of appropriate size 
        p:      (Numpy ndarray) matrix of appropriate size
        eps:    (Numpy ndarray) T x nstates array of exogenous strucural shocks. 
        s0:     (Numpy ndarray) 1 x nstates array of zeros of initial state value. Optional; Default: 0.

    Returns
        s:   (Numpy ndarray) states simulated from t = 0,1,...,T-1
        u:   (Numpy ndarray) costates simulated from t = 0,1,...,T-1

    '''

    T = np.max(eps.shape)
    nstates = np.shape(p)[0]
    ncostates = np.shape(f)[0]

    if s0 is None:

        s0 = np.zeros([1,nstates])

    s = np.array(np.zeros([T+1,nstates]))
    u = np.array(np.zeros([T,ncostates]))

    s[0]=s0

    for i,e in enumerate(eps):
        s[i+1] = p.dot(s[i]) + e
        u[i]   = f.dot(s[i+1])

    s = s[1:]

    return np.concatenate((s.T,u.T))


def klein(a=None,b=None,c=None,rho=None,nstates=None):

    '''This function solves linear dynamic models with the form of:
    
                a*Et[x(t+1)] = b*x(t) + c*z(t)

        with x(t) = [s(t); u(t)] where s(t) is a vector of predetermined (state) variables and u(t) is
        a vector of nonpredetermined costate variables. z(t) is a vector of exogenous driving variables with 
        autocorrelation matrix rho. The solution to the model is a set of matrices f,n,p,l such that:

                u(t)   = f*s(t) + n*z(t)
                s(t+1) = p*s(t) + l*z(t).

        The solution algorithm is based on Klein (2000) and is based on his solab.m Matlab program.

    Args:
        a:          (Numpy ndarray) coefficient matrix
        b:          (Numpy ndarray) coefficient matrix
        c:          (Numpy ndarray) coefficient matrix
        rho:        (Numpy ndarray) autocorrelation of exogenous forcing variables
        nstates:    (int) number of state variables

    Returns:
        f:          (Numpy ndarray) Solution matrix coeffients on s(t)
        p:          (Numpy ndarray) Solution matrix coeffients on s(t)
        n:          (Numpy ndarray) Solution matrix coeffients on z(t)
        l:          (Numpy ndarray) Solution matrix coeffients on z(t)
        stab:       (int) Indicates solution stability and uniqueness
                        stab == 1: too many stable eigenvalues
                        stab == -1: too few stable eigenvalues
                        stab == 0: just enoughstable eigenvalues
        eig:        The generalized eigenvalues from the Schur decomposition

    '''

    s,t,alpha,beta,q,z = la.ordqz(A=a,B=b,sort='ouc')
    
    q=np.mat(q)
    z=np.mat(z)
    s=np.mat(s)
    t=np.mat(t)
    a=np.mat(a)
    b=np.mat(b)
    

    forcingVars = False
    if len(np.shape(c))== 0:
        nz = 0
        rho = np.empty([0,0])
    else:
        forcingVars = True
        nz = np.shape(c)[1]
        

    # Components of the z matrix
    z11 = z[0:nstates,0:nstates]
    z12 = z[0:nstates,nstates:]
    z21 = z[nstates:,0:nstates]
    z22 = z[nstates:,nstates:]

    # number of nonpredetermined variables
    ncostates = np.shape(a)[0] - nstates
    
    if nstates>0:
        if np.linalg.matrix_rank(z11)<nstates:
            sys.exit("Invertibility condition violated.")

    s11 = s[0:nstates,0:nstates];
    if nstates>0:
        z11i = la.inv(z11)
        s11i = la.inv(s11)
    else:
        z11i = z11
        s11i = s11

    # Components of the s,t,and q matrices
    s12 = s[0:nstates,nstates:]
    s22 = s[nstates:,nstates:]
    t11 = t[0:nstates,0:nstates]
    t12 = t[0:nstates,nstates:]
    t22 = t[nstates:,nstates:]
    q1  = q[0:nstates,:]
    q2  = q[nstates:,:]

    # Verify that there are exactly nstates stable (inside the unit circle) eigenvalues:
    stab = 0

    if nstates>0:
        if np.abs(t[nstates-1,nstates-1])>np.abs(s[nstates-1,nstates-1]):
            print('Warning: Too few stable eigenvalues.')
            stab = -1

    if nstates<nstates+ncostates:
        if np.abs(t[nstates,nstates])<np.abs(s[nstates,nstates]):
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


    # Solution matrix coefficients on the endogenous state
    if nstates>0:
            dyn = np.linalg.solve(s11,t11)
    else:
        dyn = np.array([])

    f = np.real(z21.dot(z11i))
    p = np.real(z11.dot(dyn).dot(z11i))

    # Solution matrix coefficients on the exogenous state
    if not forcingVars:
        n = np.empty([ncostates,0])
        l = np.empty([nstates,0])
    else:
        mat1 = np.kron(np.transpose(rho),s22) - np.kron(np.identity(nz),t22)
        mat1i = la.inv(mat1)
        q2c = q2.dot(c)
        vecq2c = q2c.flatten(1).T
        vecm = mat1i.dot(vecq2c)
        m = np.transpose(np.reshape(np.transpose(vecm),(nz,ncostates)))
        
        n = np.real((z22 - z21.dot(z11i).dot(z12)).dot(m))
        l = np.real(-z11.dot(s11i).dot(t11).dot(z11i).dot(z12).dot(m) + z11.dot(s11i).dot(t12.dot(m) - s12.dot(m).dot(rho)+q1.dot(c)) + z12.dot(m).dot(rho))


    return f,n,p,l,stab,eig