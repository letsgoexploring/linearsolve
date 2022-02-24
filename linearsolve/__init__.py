from __future__ import division,print_function
import numpy as np
import scipy.linalg as la
from statsmodels.tools.numdiff import approx_fprime_cs, approx_fprime
from scipy.optimize import root,fsolve,broyden1,broyden2
import pandas as pd
import sys

class model:

    '''Defines a class -- linearsolve.model -- with associated methods for solving and simulating dynamic 
    stochastic general equilibrium (DSGE) models.'''

    def __init__(self,equations=None,n_states=None,n_exo_states=None,var_names=None,shock_names=None,parameters=None):
        
        '''Initializing an instance linearsolve.model requires values for the following variables:

        Args:
            equations:          (fun) A function that represents the equilibirum conditions for a DSGE model.
                                    The function should accept three arguments:
                                        * vars_fwd:     endogenous variables dated t+1
                                        * vars_cur:     endogenous variables dated t
                                        * parameters:   the parameters of the model
                                    The function should return an n-dimensional array with each element of 
                                    the returned array being equaling an equilibrium condition of the model 
                                    solved for zero. 
            n_states:           (int) The number of state variables in the model.
            n_exo_states:       (int) The number of state variables with exogenous shocks. If None, then it's assumed
                                    that all state variables have exogenous shocks. Default: None
            var_names:          (list) A list of strings with the names of the endogenous variables. The 
                                    state variables with exogenous shocks must be ordered first, followed by state
                                    variables without exogenous shocks, followed by control variables. E.g., for a 
                                    3-variables RBC model, var_names = ['a','k','c']
            shock_names:        (list) A list of strings with the names of the exogenous shocks to each state
                                    variable. The order of names must agree with the relevant elements of var_names.
            parameters:         (Pandas Series) Pandas Series object with parameter name strings as the index.

        Returns:
            None

        Attributes:
            equilibrium_fun:    (fun) Function that returns the equilibrium comditions of the model.
            n_vars:             (int) The number of variables in the model.
            n_states:           (int) The number of state variables in the model.
            n_exo_states:       (int) The number of exogenous state variables.
            n_endo_states:      (int) The number of endogenous state variables.
            n_costates:         (int) The number of costate or control variables in the model.
            names:              (dict) A dictionary with keys 'variables', 'shocks', and 'param' that
                                    stores the names of the model's variables, shocks, and parameters.
            parameters:         (Pandas Series) A Pandas Series with parameter name strings as the 
                                    index.
        '''
        
        self.equilibrium_fun= equations
        self.n_vars = len(var_names)
        self.n_states = n_states
        self.n_exo_states = n_exo_states or n_states
        self.n_endo_states = n_states - n_exo_states
        self.n_costates=self.n_vars-n_states
        self.parameters = parameters

        names = {}

        names['variables'] = var_names
        
        names['shocks'] = shock_names or ['e_'+var_names[i] for i in range(self.n_exo_states)]
        
        if len(names['shocks']) != self.n_exo_states:
            raise Exception('Length of shock_names doesn\'t match number of exogenous states')
        
        names['param'] = list(self.parameters.index)
        
        self.names = names

    # Methods

    def approximate_and_solve(self,log_linear=True,eigenvalue_warnings=True):

        '''Method approximates and solves a dynamic stochastic general equilibrium (DSGE) model by 
        constructing the log-linear approximation (if the model isn't log-linear) and solving the model
        using Klein's (2000) method.

        Args:
            log_linear:             (bool) Whether to compute log-linear or linear approximation. Default: True
            eigenvalue_warnings:    (bool) Whether to print warnings that there are too many or few eigenvalues. Default: True

        Returns:
            None

        Attributes:
            a:          (Numpy ndarray) Coefficient matrix on forward-dated variables.
            b:          (Numpy ndarray) Coefficient matrix on current-dated variables.
            f:          (Numpy ndarray) Solution matrix coeffients on s(t) in control equation.
            p:          (Numpy ndarray) Solution matrix coeffients on s(t) in state equation.
            stab:       (int) Indicates solution stability and uniqueness
                            stab == 1: too many stable eigenvalues
                            stab == -1: too few stable eigenvalues
                            stab == 0: just enoughstable eigenvalues
            eig:        The generalized eigenvalues from the Schur decomposition
            log_linear: (bool) Whether the model is log-linear. Sets to log-linear.

        '''

        # Set attribute
        self.log_linear = log_linear

        # Approximate
        if log_linear == True:
            self.log_linear_approximation()
        else:
            self.linear_approximation()

        # Solve the model
        self.solve_klein(self.a,self.b,eigenvalue_warnings=eigenvalue_warnings)

    def approximated(self,round=True,precision=4):

        '''Returns a string containing the log-linear approximation to the equilibrium conditions

        Args:
            round:      (bool) Whether to round the coefficents in the linear equations. Default: True
            precision:  (int) Number of decimals to round the coefficients. Default: 4

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
        if self.log_linear==True:
            lines ='Log-linear equilibrium conditions:\n\n'
        else:
            lines ='Linear equilibrium conditions:\n\n'


        left_length = 1


        for i in range(self.n_vars):
            
            left = ''
            right = ''
            
            left_plus_flag = 0
            right_plus_flag = 0
            if all(np.isclose(0,a[i])):
                left += '0'
            else:
                for j in range(self.n_vars):

                    if not np.isclose(0,a[i][j]):

                        name = self.names['variables'][j]

                        if j >self.n_states-1:
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
                for j in range(self.n_vars):

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
                

        for i in range(self.n_vars):
            leftsides[i] = leftsides[i].rjust(left_length)
            lines+=leftsides[i]+' = '+rightsides[i]+'\n\n'

        lines = lines[:-2]

        return lines


    def check_ss(self):

        '''Uses Numpy.isclose() to print whether each steady state equilibrium condition evaluates to 
        something close to zero.

        Args:
            None

        Returns:
            Numpy ndarray

        Attributes:
            None

        '''

        try:
            print(np.isclose(self.equilibrium_fun(self.ss,self.ss,self.parameters),0))
        except:
            print('Set the steady state first.')

    
    def compute_ss(self,guess=None,method='fsolve',options={}):

        '''Attempts to solve for the steady state of the model.

        Args:
            guess:      (Pandas Series, Numpy array, or list) An initial guess for the 
                            steady state solution. The result is highly sensisitve to the intial 
                            guess chosen, so be careful. If the guess is a Numpy ndarray or a list
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
            guess = np.ones(self.n_vars)
        else:
            if isinstance(guess, list):
                guess = np.array(guess)


        # Create function for nonlinear solver
        def ss_fun(variables):

            variables = pd.Series(variables,index = self.names['variables'])

            return self.equilibrium_fun(variables,variables,self.parameters)
        
        def real_ss_fun(variables_transformed):

            z_vals = [a + b*1j for a,b in zip(variables_transformed[:int(len(variables_transformed)/2)],variables_transformed[int(len(variables_transformed)/2):])]
            actual_ss_fun = ss_fun(z_vals)
            retval = np.r_[np.real(actual_ss_fun),np.imag(actual_ss_fun)]
            return retval

        if not np.iscomplexobj(self.parameters) and not np.iscomplexobj(guess):
            
            if method == 'fsolve':
                steady_state =fsolve(ss_fun,guess,**options)

            elif method == 'root':
                steady_state =root(ss_fun,guess,**options)['x']

            elif method == 'broyden1':
                steady_state =broyden1(ss_fun,guess,**options)

            elif method == 'broyden2':
                steady_state =broyden2(ss_fun,guess,**options)
            
            
        else:
            
            if method == 'fsolve':
                steady_state =fsolve(real_ss_fun,np.r_[np.real(guess),np.imag(guess)],**options)

            elif method == 'root':
                steady_state =root(real_ss_fun,np.r_[np.real(guess),np.imag(guess)],**options)['x']

            elif method == 'broyden1':
                steady_state =broyden1(real_ss_fun,np.r_[np.real(guess),np.imag(guess)],**options)

            elif method == 'broyden2':
                steady_state =broyden2(real_ss_fun,np.r_[np.real(guess),np.imag(guess)],**options)
                
            steady_state = [a + b*1j for a,b in zip(steady_state[:int(len(steady_state)/2)],steady_state[int(len(steady_state)/2):])]


        # Add ss attribute
        self.ss = pd.Series(steady_state,index=self.names['variables'])


    def impulse(self,T=51,t0=1,shocks=None,percent=False,diff=True):

        '''Computes impulse responses for shocks to each state variable.

        Arguments:
                T:          (int) Number of periods to simulate. Default: 51
                t0:         (int) Period in which the shocks are realized. Must be greater than or equal to 
                                0. Default: 1
                shocks:     (list or Numpy array) An (ns x 1) array of shock values. If shocks==None, and 
                                log_linear==True, shocks is set to a vector of 0.01s. If shocks==None and
                                log_linear==False, shocks is set to a vector of 1s. Default = None
                percent:    (bool) Whether to multiply simulated values by 100. Only works for log-linear 
                                approximations. Default: False
                diff:       (bool) Subtract steady state for linear approximations (or log steady state for 
                                log-linear approximations). Default: True
        
        Returns
            None

        Attributes:
            irs:    (dict) A dictionary containing Pandas DataFrames. Has the form: 
                        self.irs['shock name']['endog var name']

        '''

        # Initialize dictionary
        irsDict = {}

        # Set numbers of costate and state variables
        n_costates = self.n_costates
        n_states = self.n_states
        
        # iterate over all shocks, compute impulse responses, and add results to dictionary
        for j,name in enumerate(self.names['shocks']):

            s0 = np.zeros([1,n_states])
            eps= np.zeros([T,n_states])
            if shocks is not None:
                try:
                    eps[t0][j] = shocks[name]
                except:
                    try:
                        eps[t0][j] = shocks[j]
                    except:
                        if self.log_linear:
                            eps[t0][j] = 0.01
                        else:
                            eps[t0][j] = 1

            else:
                if self.log_linear:
                    eps[t0][j] = 0.01
                else:
                    eps[t0][j] = 1

            x = ir(self.f,self.p,eps,s0)

            frameDict = {self.names['shocks'][j]:eps.T[j]}
            for i,endoName in enumerate(self.names['variables']):
                if diff:
                    frameDict[endoName] = x[i]
                else:
                    if not self.log_linear:
                        frameDict[endoName] = x[i] + self.ss[endoName]
                    else:
                        frameDict[endoName] = x[i] + np.log(self.ss[endoName])

            irFrame = pd.DataFrame(frameDict,index = np.arange(T))

            if percent==True and self.log_linear:
                irFrame = 100*irFrame

            if shocks is None or len(shocks)>j:
                irsDict[self.names['shocks'][j]] = irFrame

        # Set attribute
        self.irs = irsDict


    def linear_approximation(self,steady_state=None):

        ''' Given a nonlinear rational expectations model in the form:

                    psi_1[x(t+1),x(t)] = psi_2[x(t+1),x(t)]

            this method returns the linear approximation of the model with matrices a and b such that:

                    a * y(t+1) = b * y(t)

            where y(t) = x(t) - x is the log deviation of the vector x from its steady state value.

        Args:
            steady_state:   (Pandas Series or numpy array or list)

        Returns:
            None

        Attributes:
            log_linear:     (bool) Whether the model is log-linear. Sets to False.
            a:              (Numpy ndarray)
            b:              (Numpy ndarray)

        '''

        # Set log_linear attribute
        self.log_linear=False

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

        if not np.iscomplexobj(self.parameters):
            
            # Assign attributes
            self.a= approx_fprime_cs(steady_state.ravel(),equilibrium_fwd)
            self.b= -approx_fprime_cs(steady_state.ravel(),equilibrium_cur)

        else:

            # Assign attributes
            self.a= approx_fprime(steady_state.ravel(),equilibrium_fwd)
            self.b= -approx_fprime(steady_state.ravel(),equilibrium_cur)


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
            log_linear:     (bool) Whether the model is log_linear. Sets to True.
            a:              (Numpy ndarray)
            b:              (Numpy ndarray)

        '''

        # Set log_linear attribute
        self.log_linear=True

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

        log_equilibrium_fwd = lambda log_fwd: log_equilibrium(log_fwd,np.log(self.ss))
        log_equilibrium_cur = lambda log_cur: log_equilibrium(np.log(self.ss),log_cur)
        

        if not np.iscomplexobj(self.parameters):

            # Assign attributes
            self.a= approx_fprime_cs(np.log(self.ss).ravel(),log_equilibrium_fwd)
            self.b= -approx_fprime_cs(np.log(self.ss).ravel(),log_equilibrium_cur)

        else:

            # Assign attributes
            self.a= approx_fprime(np.log(self.ss).ravel(),log_equilibrium_fwd,centered=True)
            self.b= -approx_fprime(np.log(self.ss).ravel(),log_equilibrium_cur,centered=True)
            
            


    def set_ss(self,steady_state):

        '''Directly set the steady state of the model. 

        Args:
            steady_state:   (Pandas Series, Numpy array, or list)

        Returns:
            None

        Attributes:
            ss: (Pandas Series) Steady state values of endogenous variables

        '''

        try:
            self.ss = steady_state[self.names['variables']]
            self.ss = self.ss.astype(np.promote_types(float,np.promote_types(self.parameters.dtype,steady_state.dtype)))

        except:
            self.ss = pd.Series(steady_state,index=self.names['variables'],dtype=np.promote_types(float,np.promote_types(self.parameters.dtype,np.array(steady_state).dtype)))
            
        self.parameters = self.parameters.astype(np.promote_types(self.ss.dtype,self.parameters.dtype))


    def solve_klein(self,a=None,b=None,eigenvalue_warnings=True):

        '''Solves a linear rational expectations model of the form:

                a * x(t+1) = b * x(t) + e(t)

        The method returns the solution to the law of motion:

                u(t)   = f*s(t) + e(t)
                s(t+1) = p*s(t)

        Args:
            a:                      (Numpy ndarray) coefficient matrix
            b:                      (Numpy ndarray) coefficient matrix
            eigenvalue_warnings:    (bool) Whether to print warnings that there are too many or few eigenvalues. Default: True

        Returns:
            None

        Attributes:
            f:      (Numpy ndarray) Solution matrix coeffients on s(t)
            p:      (Numpy ndarray) Solution matrix coeffients on s(t)
            stab:   (int) Indicates solution stability and uniqueness
                        stab == 1: too many stable eigenvalues
                        stab == -1: too few stable eigenvalues
                        stab == 0: just enough stable eigenvalues
            eig:    The generalized eigenvalues from the Schur decomposition

         '''

        if a is None and b is None:
            
            a = self.a
            b = self.b

        self.f,n,self.p,l,self.stab,self.eig = klein(a=a,b=b,c=None,phi=None,n_states=self.n_states,eigenvalue_warnings=eigenvalue_warnings)

        if not np.iscomplexobj(self.parameters):

            self.f = np.real(self.f)
            self.p = np.real(self.p)
            l = np.real(l)
            n = np.real(n)

            # self.f = np.abs(self.f)
            # self.p = np.abs(self.p)
            # l = np.abs(l)
            # n = np.abs(n)

            



    def solved(self,round=True,precision=4):

        '''Returns a string containing the solution to the linear system

        Args:
            round:       (bool) Whether to round the coefficents in the solution equations. Default: True
            precisions:  (int) Number of decimals to round the coefficients. Default: 4

        Returns:
            String with the linear approximation to the equilibrium conditions.

        Attributes:
            None

        '''

        if round is True:
            f = np.round(self.f,precision)
            p = np.round(self.p,precision)

            
        leftsides =  []
        rightsides = []
        if self.log_linear==True:
            lines ='Solution to the log-linear system:\n\n'
        else:
            lines ='Solution to the linear system:\n\n'

        left_length = 1


        for i in range(self.n_states):
            
            left = ''
            right = ''
            right_plus_flag = 0

            left+= self.names['variables'][i]+'[t+1]'
            
            if all(np.isclose(0,p[i])):
                right += self.names['shocks'][i]+'[t+1]'
            else:
                for j in range(self.n_states):

                    if not np.isclose(0,p[i][j]):

                        if right_plus_flag == 0:
                            right += str(p[i][j])+'·'+self.names['variables'][j]+'[t]'
                            right_plus_flag+=1

                        else:
                            if p[i][j] > 0:
                                right += '+'+str(p[i][j])+'·'+self.names['variables'][j]+'[t]'

                            else:
                                right += str(p[i][j])+'·'+self.names['variables'][j]+'[t]'
                if i <self.n_exo_states:
                    right+='+'+self.names['shocks'][i]+'[t+1]'
            leftsides.append(left)
            rightsides.append(right)
                            
                        
            if len(left)>left_length:
                left_length = len(left)

                
        for i in range(self.n_vars-self.n_states):
            
            left = ''
            right = ''
            right_plus_flag = 0
            
            left+= self.names['variables'][self.n_states+i]+'[t]'
            
            if all(np.isclose(0,f[i])):
                right += '0'
            else:
                for j in range(self.n_states):
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
                
                
        for i in range(self.n_vars):
            leftsides[i] = leftsides[i].rjust(left_length)
            lines+=leftsides[i]+' = '+rightsides[i]+'\n\n'
        lines = lines[:-2]
                    

        return lines
    

    def stoch_sim(self,T=51,drop_first=300,cov_mat=None,seed=None,percent=False,diff=True):
        
        '''Computes a stohcastic simulation of the model.

        Arguments:
                T:          (int) Number of periods to simulate. Default: 51
                drop_first: (int) Number of periods to simulate before generating the simulated periods. 
                                Default: 300
                cov_mat:    (list or Numpy.ndarray) Covariance matrix shocks.
                seed:       (int) Sets the seed for the Numpy random number generator. Default: None
                percent:    (bool) Whether to multiply simulated values by 100. Only works for log-linear 
                                approximations. Default: False
                diff:       (bool) Subtract steady state for linear approximations (or log steady state for 
                                log-linear approximations). Default: True
        
        Returns
            None

        Attributes:
            simulated:    (Pandas DataFrame)

        '''

        # Set numbers of costate and state variables
        n_costates = self.n_costates
        n_states = self.n_states

        # Initialize states
        s0 = np.zeros([1,n_states])

        # Set cov_mat if not given
        cov_mat = np.array(cov_mat)
        if cov_mat.ndim==1:
            cov_mat = np.array([cov_mat])
            
            if len(cov_mat) != self.n_exo_states:
                raise Exception('Length of cov_mat doesn\'t match number of exogenous states')

        # Set seed for the Numpy random number generator
        if seed is not None and type(seed)==int:

            np.random.seed(seed)

        # Simulate shocks
        eps = np.zeros([drop_first+T,n_states])
        eps[:,:len(cov_mat)] = np.random.multivariate_normal(mean=np.zeros(len(cov_mat)),cov=cov_mat,size=[drop_first+T])

        # Compute impulse responses given shocks
        x = ir(self.f,self.p,eps,s0)

        # Construct DataFrame
        frameDict = {}
        for j,exoName in enumerate(self.names['shocks']):
            frameDict[exoName] = eps.T[j][drop_first:]
        for i,endoName in enumerate(self.names['variables']):
            if diff:
                frameDict[endoName] = x[i][drop_first:]
            else:
                frameDict[endoName] = x[i][drop_first:] + self.ss[endoName]

        simFrame = pd.DataFrame(frameDict,index = np.arange(T))
        if percent==True:
            simFrame = 100*simFrame

        # Assign attribute
        self.simulated = simFrame

    


### End of model class ####################################################################################


def ir(f,p,eps,s0=None):

    '''Simulates a model in the following form:

            u(t)   = f*s(t) + e(t)
            s(t+1) = p*s(t)

    where s(t) is an (n_states x 1) vector of state variables, u(t) is an (n_costates x 1) vector of costate
    variables, and e(t) is an (n_states x 1) vector of exogenous shocks.

    Args:
        f:      (Numpy ndarray) Coefficnent matrix of appropriate size 
        p:      (Numpy ndarray) Coefficnent matrix of appropriate size
        eps:    (Numpy ndarray) T x n_states array of exogenous shocks. 
        s0:     (Numpy ndarray) 1 x n_states array of zeros of initial state value. Optional; Default: 0.

    Returns
        s:   (Numpy ndarray) states simulated from t = 0,1,...,T-1
        u:   (Numpy ndarray) costates simulated from t = 0,1,...,T-1

    '''

    T = np.max(eps.shape)
    n_states = np.shape(p)[0]
    n_costates = np.shape(f)[0]

    if s0 is None:

        s0 = np.zeros([1,n_states])
    
    s = np.zeros([T+1,n_states],dtype=np.promote_types(f.dtype,p.dtype))
    u = np.zeros([T,n_costates],dtype=np.promote_types(f.dtype,p.dtype))

    s[0]=s0

    for i,e in enumerate(eps):
        s[i+1] = p.dot(s[i]) + e
        u[i]   = f.dot(s[i+1])

    s = s[1:]

    return np.concatenate((s.T,u.T))


def klein(a=None,b=None,c=None,phi=None,n_states=None,eigenvalue_warnings=True):

    '''Solves linear dynamic models with the form of:
    
                a*Et[x(t+1)] = b*x(t) + c*z(t)

        with x(t) = [s(t); u(t)] where s(t) is a vector of predetermined (state) variables and u(t) is
        a vector of nonpredetermined costate variables. z(t) is a vector of exogenous forcing variables with 
        autocorrelation matrix phi. The solution to the model is a set of matrices f, n, p, l such that:

                u(t)   = f*s(t) + n*z(t)
                s(t+1) = p*s(t) + l*z(t).

        The solution algorithm is based on Klein (2000) and his solab.m Matlab program.

    Args:
        a:                      (Numpy ndarray) Coefficient matrix on future-dated variables
        b:                      (Numpy ndarray) Coefficient matrix on current-dated variables
        c:                      (Numpy ndarray) Coefficient matrix on exogenous forcing variables
        phi:                    (Numpy ndarray) Autocorrelation of exogenous forcing variables
        n_states:               (int) Number of state variables
        eigenvalue_warnings:    (bool) Whether to print warnings that there are too many or few eigenvalues. Default: True

    Returns:
        f:          (Numpy ndarray) Solution matrix coeffients on s(t) for u(t)
        n:          (Numpy ndarray) Solution matrix coeffients on z(t) for u(t)
        p:          (Numpy ndarray) Solution matrix coeffients on s(t) for s(t+1)
        l:          (Numpy ndarray) Solution matrix coeffients on z(t) for s(t+1)
        stab:       (int) Indicates solution stability and uniqueness
                        stab == 1: too many stable eigenvalues
                        stab == -1: too few stable eigenvalues
                        stab == 0: just enoughstable eigenvalues
        eig:        The generalized eigenvalues from the Schur decomposition

    '''

    s,t,alpha,beta,q,z = la.ordqz(A=a,B=b,sort='ouc',output='complex')

    # print('type of s,',s.dtype)
    # print('type of t,',t.dtype)

    forcingVars = False
    if len(np.shape(c))== 0:
        nz = 0
        phi = np.empty([0,0])
    else:
        forcingVars = True
        nz = np.shape(c)[1]
        

    # Components of the z matrix
    z11 = z[0:n_states,0:n_states]
    z12 = z[0:n_states,n_states:]
    z21 = z[n_states:,0:n_states]
    z22 = z[n_states:,n_states:]

    # number of nonpredetermined variables
    n_costates = np.shape(a)[0] - n_states
    
    if n_states>0:
        if np.linalg.matrix_rank(z11)<n_states:
            sys.exit("Invertibility condition violated.")

    s11 = s[0:n_states,0:n_states];
    if n_states>0:
        z11i = la.inv(z11)
        s11i = la.inv(s11)
    else:
        z11i = z11
        s11i = s11

    # Components of the s,t,and q matrices
    s12 = s[0:n_states,n_states:]
    s22 = s[n_states:,n_states:]
    t11 = t[0:n_states,0:n_states]
    t12 = t[0:n_states,n_states:]
    t22 = t[n_states:,n_states:]
    q1  = q[0:n_states,:]
    q2  = q[n_states:,:]

    # Verify that there are exactly n_states stable (inside the unit circle) eigenvalues:
    stab = 0

    if n_states>0:
        if np.abs(t[n_states-1,n_states-1])>np.abs(s[n_states-1,n_states-1]):
            if eigenvalue_warnings:
                print('Warning: Too few stable eigenvalues.')
            stab = -1

    if n_states<n_states+n_costates:
        if np.abs(t[n_states,n_states])<np.abs(s[n_states,n_states]):
            if eigenvalue_warnings:
                print('Warning: Too many stable eigenvalues.')
            stab = 1

    # Compute the generalized eigenvalues
    tii = np.diag(t)
    sii = np.diag(s)
    eig = np.zeros(np.shape(tii),dtype=np.complex128)
    # eig = np.zeros(np.shape(tii))

    for k in range(len(tii)):
        if np.abs(sii[k])>0:
            # eig[k] = np.abs(tii[k])/np.abs(sii[k])
            eig[k] = tii[k]/sii[k]    
        else:
            eig[k] = np.inf

    # Solution matrix coefficients on the endogenous state
    if n_states>0:
            dyn = np.linalg.solve(s11,t11)
    else:
        dyn = np.array([])


    f = z21.dot(z11i)
    p = z11.dot(dyn).dot(z11i)

    # Solution matrix coefficients on the exogenous state
    if not forcingVars:
        n = np.empty([n_costates,0])
        l = np.empty([n_states,0])
    else:
        mat1 = np.kron(np.transpose(phi),s22) - np.kron(np.identity(nz),t22)
        mat1i = la.inv(mat1)
        q2c = q2.dot(c)
        vecq2c = q2c.flatten(1).T
        vecm = mat1i.dot(vecq2c)
        m = np.transpose(np.reshape(np.transpose(vecm),(nz,n_costates)))
        
        n = (z22 - z21.dot(z11i).dot(z12)).dot(m)
        l = -z11.dot(s11i).dot(t11).dot(z11i).dot(z12).dot(m) + z11.dot(s11i).dot(t12.dot(m) - s12.dot(m).dot(phi)+q1.dot(c)) + z12.dot(m).dot(phi)

    return f,n,p,l,stab,eig