
Additional ``linearsolve`` Functions
====================================


.. py:function:: linearsolve.ar1(mu=0,rho=0.5)

            Initiates an AR(1) model in linearsolve. Model has the form:

            .. math::
                  x_{t+1} &= (1-\rho)\mu \rho x_t + \epsilon_{t+1}\\

            Computes the steady state, and solves the model.

            :param mu: Mean of the process
            :type mu: float
            :param rho: Autoregressive coeficcient of the process
            :type rho: float
            :return: :py:class:`linearsolve.model`

.. py:function:: linearsolve.ir(f,p,eps,s0=None)

            Simulates a model in the following form:

            .. math:: 
                  u_t   &= fs_t + \epsilon_{t}  \\

            .. math:: 
                  s_{t+1} &= ps_t \\

            where :math:`s_t` is an (n_states x 1) vector of state variables, :math:`u_t` is an (n_costates x 1) vector of costate variables, and :math:`\epsilon_t` is an (n_states x 1) vector of exogenous shocks.

            :param f: Coefficnent matrix of appropriate size 
            :type f: Numpy.ndarray
            :param p: Coefficnent matrix of appropriate size
            :type p: Numpy.ndarray
            :param eps: T x n_states array of exogenous shocks. 
            :type eps: Numpy.ndarray or list
            :param s0: 1 x n_states array of zeros of initial state value. Optional; Default: 0.
            :type s0: Numpy.ndarray or list
            :returns: Array containing two arrays: one corresponding to simulated :math:`u_t` values and the other to simulated :math:`s_t` values.
            :rtype: Numpy.ndarray

.. py:function:: linearsolve.klein(a=None,b=None,c=None,phi=None,n_states=None)

            Solves linear dynamic models with the form of:
    
            .. math::
                A E_tx_{t+1} = Bx_t + Cz_t

            with :math:`x_t = [s_t; u_t]` where :math:`s_t` is a vector of predetermined (state) variables and :math:`u_t` is a vector of nonpredetermined costate variables. :math:`z_t` is a vector of exogenous forcing variables with autocorrelation matrix :math:`\Phi`. The solution to the model is a set of matrices :math:`f,n,p,l` such that:

            .. math::

                u_t   = fs_t + nz_t

            .. math::

                s_{t+1} = ps_t + lz_t.

            The solution algorithm is based on Klein (2000) and his solab.m Matlab program.

            :param a: Coefficient matrix on future-dated variables
            :type a: Numpy.ndarray
            :param b: Coefficient matrix on current-dated variables
            :type b: Numpy.ndarray
            :param p: Coefficient matrix on exogenous forcing variables
            :type p: Numpy.ndarray
            :param phi: Autocorrelation of exogenous forcing variables
            :type phi: Numpy.ndarray or list
            :param n_states: number of state variables
            :type n_states: int


            :returns:
                * **f** (:py:obj:`Numpy.ndarray`) - Coefficient matrix on state vector in control equation
                * **p** (:py:obj:`Numpy.ndarray`) - Coefficient matrix on state vector in state equation
                * **n** (:py:obj:`Numpy.ndarray`) - Coefficient matrix on forcing vector in control equation
                * **l** (:py:obj:`Numpy.ndarray`) - Coefficient matrix on forcing vector in state equation
                * **stab** (:py:class:`int`) - Indicates solution stability and uniqueness. stab =1: too many stable eigenvalues, stab = -1: too few stable eigenvalues, stab = 0: just enough stable eigenvalues
                * **eig** (:py:obj:`Numpy.ndarray`) - Generalized eigenvalues from the Schur decomposition