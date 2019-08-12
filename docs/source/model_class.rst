.. fredpy documentation master file, created by
   sphinx-quickstart on Fri Aug 19 15:23:34 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

``linearsolve.model``
==================================




.. py:class:: linearsolve.model(equations=None,n_states=None,var_names=None,shock_names=None,parameters=None,parameter_names=None)
	
	Creates an instance of :py:class:`linearsolve.model` that stores equilibrium conditions for a DSGE model.

	:param function equations: A function that represents the equilibirum conditions for a DSGE model. The function should return an n-dimensional array with each element of the returned array being equaling an equilibrium condition of the model solved for zero. The function should accept three arguments:

		* **variables_forward**: Endogenous variables dated t+1.
		* **variables_current**: Endogenous variables dated t.
		* **parameters**: The parameters of the model.


	:param int n_states: The number of state variables in the model.
	:param list var_names: A list of strings with the names of the endogenous variables. The state variables must be ordered first.
	:param list shock_names: A list of strings with the names of the exogenous shocks to each state variable. The order of names must agree with var_names.
	:param parameters: Either a Pandas Series object with parameter name strings as the index OR a list or an array of parameter values.
	:type parameters: Pandas.Series or Numpy.ndarray or list
	:param list parameter_names: (Optional) If parameters is given as a list, then this list of strings will be used to save the parameters with names as a Pandas Series object.


	**Attributes:**
    

		:equilibrium_fun: (:py:obj:`function`) -- Function input with the **equations** parameter.
		:n_vars: (:py:obj:`int`) -- Number of endogenous variables.
		:n_states: (:py:obj:`Numpy.ndarray`) -- Number of state variables.
		:n_costates: (:py:obj:`int`) -- Number of control or costate variables.
		:names: (:py:obj:`dict`) -- A dictionary with keys 'variables', 'shocks', and 'param' that stores the names of the model's variables, shocks, and parameters.
		:parameters: (:py:obj:`Pandas.Series`) -- A Pandas Series with parameter name strings as the index. If **parameter_names** wasn't supplied, then parameters are labeled 'parameter 1', 'parameter2', etc.

		
	**Methods:**


		.. py:function:: approximate_and_solve(log_linear=True)

			Method approximates and solves a dynamic stochastic general equilibrium (DSGE) model by constructing the log-linear approximation (if the model isn't log-linear) and solving the model using Klein's (2000) method.


			:param log_linear: Whether to compute log-linear or linear approximation. Default: True
   			:type log_linear: bool

   			:var a: Coefficient matrix on forward-dated variables.
			:vartype a: Numpy.ndarray

			:var b: Coefficient matrix on current-dated variables.
			:vartype b: Numpy.ndarray

			:var f: Solution matrix coeffients on s(t) in control equation.
			:vartype f: Numpy.ndarray

			:var p: Solution matrix coeffients on u(t) in state equation.
			:vartype p: Numpy.ndarray

			:var stab: Indicates solution stability and uniqueness. stab =1: too many stable eigenvalues, stab = -1: too few stable eigenvalues, stab = 0: just enough stable eigenvalues

			:vartype stab: int

			:var eig: Generalized eigenvalues from the Schur decomposition

			:vartype eig: Numpy.ndarray



		.. py:function:: approximated(round=True,precision=4)

			Returns a string containing the log-linear approximation to the equilibrium conditions.

			:param round: Whether to round the coefficents in the linear equations. Default: True
   			:type round: bool

   			:param precision: Number of decimals to round the coefficients. Default: 4
   			:type precision: int

   			:returns: string containing the log-linear approximation to the equilibrium conditions

   			:rtype: string



   		.. py:function:: check_ss()

   			Uses Numpy.isclose() to print whether each steady state equilibrium condition evaluates to something close to zero.

   			:returns: Numpy.ndarry of booleans indicating whether the stored staeady state satisfies each equilibrium condition

   			:rtype: Numpy.ndarray



		.. py:function:: compute_ss(guess=None,method='fsolve',options={})

			Attempts to solve for the steady state of the model. Stores results as **ss** attribute.

			:param guess: An initial guess for the steady state solution. The result is highly sensisitve to the intial guess chosen, so be careful. If the guess is a Numpy ndarray or a list then the elements must be ordered to conform with self.names['variables'].
   			:type guess: Pandas.Series or Numpy.ndarray or list

			:param string method: The function from the Scipy library to use. Your choices are: 'root', 'fsolve' (default), broyden1, broyden2.
			:param dictionary options: A dictionary of optional arguments to pass to the numerical solver. Check out the Scipy documentation to see the options available for each routine: http://docs.scipy.org/doc/scipy/reference/optimize.html

			:var ss: The steady state of the model.
			:vartype ss: Pandas.Series



		.. py:function:: impulse(T=51,t0=1,shocks=None,percent=False,diff=True)

			Computes impulse responses for shocks to each state variable.

			:param T: Number of periods to simulate. Default: 1
   			:type T: int

   			:param t0: Period in which the shocks are to be realized. Must be greater than or equal to 0. default: 1
   			:type t0: int

   			:param shocks: An array of shock values with length equal to the number of shocks. If shocks=None and log_linear=True, shocks is set to a vector of 0.01s. If shocks=None and log_linear=False, shocks is set to a vector of 1s. Default: None
   			:type shocks: list or Numpy.ndaray

			:param percent: Whether to multiply simulated values by 100. Only works for log-linear approximations. Default: False
   			:type percent: bool

   			:param diff: Subtract steady state for linear approximations (or log steady state for log-linear approximations). Default: True
   			:type diff: bool


   			:var irs: A dictionary containing Pandas DataFrames. Has the form: self.irs['shock name']['endog var name']
			:vartype irs: dict



		.. py:function:: linear_approximation(steady_state=None)

			Given a nonlinear rational expectations model in the form:
			
			.. math::
				\psi_1[x_{t+1},x_t] = \psi_2[x_{t+1},x_t]

			this method returns the linear approximation of the model with matrices :math:`A` and :math:`B` such that:
			
			.. math::
				A y_{t+1} = B y_t

			where :math:`y_t = x_t - x` is the log deviation of the vector :math:`x` from its steady state value.

			:param steady_state: Coefficient matrix on forward-dated variables.
   			:type steady_state: Pandas.Series

			:var a: Coefficient matrix on forward-dated variables.
   			:vartype a: Numpy.ndaray

			:var b: Coefficient matrix on current-dated variables.
   			:vartype b: Numpy.ndaray

			:var log_linear: Whether the model is log-linear. Sets to False.
			:vartype ss: bool



		.. py:function:: log_linear_approximation(steady_state=None)

			Given a nonlinear rational expectations model in the form:
			
			.. math::
				\psi_1[x_{t+1},x_t] = \psi_2[x_{t+1},x_t]

			this method returns the log-linear approximation of the model with matrices :math:`A` and :math:`B` such that:
			
			.. math::
				A y_{t+1} = B y_t

			where :math:`y_t = \log x_t - \log x` is the log-deviation of the vector :math:`x` from its steady state value.

			:param steady_state: Coefficient matrix on forward-dated variables.
   			:type steady_state: Pandas.Series

   			:var a: Coefficient matrix on forward-dated variables.
   			:vartype a: Numpy.ndaray

			:var b: Coefficient matrix on current-dated variables.
   			:vartype b: Numpy.ndaray

			:var log_linear: Whether the model is log-linear. Sets to True.
			:vartype ss: bool



		.. py:function:: set_ss(steady_state)

			Directly set the steady state of the model. Stores results as **ss** attribute.

			:param steady_state: The steady state of the model.
   			:type steady_state: Pandas.Series or Numpy.ndarray or list

			:var ss: The steady state of the model.
			:vartype ss: Pandas.Series

		


		.. py:function:: solve_klein(a=None,b=None)

			Solves a linear rational expectations model of the form:
			
			.. math::
				A x_{t+1} = B x_t

			this method computes the log-linear approximation of the model with matrices :math:`A` and :math:`B` such that:

			.. math::

   				u_t  &= fs_t + \epsilon_t\\

   			and:

   			.. math::
   				s_{t+1}  &= ps_t\\

   			where :math:`s_t` denotes the vector of state variables and :math:`f_t` denotes the vector of forward-looking variables.

	        
			:param a: Coefficient matrix on forward-dated variables.
   			:type a: Numpy.ndaray

			:param b: Coefficient matrix on current-dated variables.
   			:type b: Numpy.ndaray

   			:var f: coeficient matrix.
			:vartype f: Numpy.ndarray

			:var p: coeficient matrix.
			:vartype p: Numpy.ndarray

			:var stab: Indicates solution stability and uniqueness. stab =1: too many stable eigenvalues, stab = -1: too few stable eigenvalues, stab = 0: just enough stable eigenvalues

			:vartype stab: int

			:var eig: Generalized eigenvalues from the Schur decomposition

			:vartype eig: Numpy.ndarray

		

		.. py:function:: solved(round=True,precision=4)

			Returns a string containing the solution to the linear system

			:param round: Whether to round the coefficents in the linear equations. Default: True
   			:type round: bool

   			:param precision: Number of decimals to round the coefficients. Default: 4
   			:type precision: int

   			:returns: string containing the solution to linear system
   			
   			:rtype: string



		.. py:function:: stoch_sim(T=51,drop_first=300,cov_mat=None,seed=None,percent=False,diff=True)

			Computes a stohcastic simulation of the model.

			:param T: Number of periods to simulate. Default: 1
   			:type T: int

   			:param drop_first: Number of periods to simulate before generating the simulated periods. Default: 300
   			:type drop_first: int

   			:param cov_mat: Covariance matrix shocks. If cov_mat is None, it's set to Numpy.eye(n_states). Default: None
   			:type cov_mat: list or Numpy.ndarray

   			:param seed: Sets the seed for the Numpy random number generator. Default: None
   			:type seed: int

			:param percent: Whether to multiply simulated values by 100. Only works for log-linear approximations. Default: False
   			:type percent: bool

   			:param diff: Subtract steady state for linear approximations (or log steady state for log-linear approximations). Default: True
   			:type diff: bool

   			:var simulated: A DataFrame with a column for each variable.
			:vartype simulated: Pandas.DataFrame