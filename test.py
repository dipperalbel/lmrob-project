class psi_func:


    """
    The class "psi_func" is used to store ψ (psi) functions for M-estimation. In particular, an object of the class contains ρ(x) (rho), its derivative ψ(x) (psi), the weight function ψ(x)/x, and ﬁrst derivative of ψ, Dpsi = ψ0(x)
    """

    def __init__(self, rho, psi,wgt, Dpsi, Dwgt, tDefs, Erho, Epsi2, EDpsi, name, xtras):
        """Example of docstring on the __init__ method.

        The __init__ method may be documented in either the class level
        docstring, or as a docstring on the __init__ method itself.

        Either form is acceptable, but the two should not be mixed. Choose one
        convention to document the __init__ method and be consistent with it.

        Note:
            Do not include the `self` parameter in the ``Args`` section.

        Args:
            rho (function): the ρ() function. This is used to formulate the objective function; ρ() can be regarded as generalized negative log-likelihood.
            psi (function): ψ() is the derivative of ρ.
            wgt (function): The weight function ψ(x)/x.
            Dpsi (function): the derivative of ψ, Dpsi(x) = psi0(x).
            Dwgt (function): the derivative of the weight function.
            tDefs (list): named numeric vector oftuning parameterDefault values.
            Epsi2 (function): The function for computing E[ψ2(X)] when X is standard normal.
            EDpsi (function): The function for computing E[ψ0(X)] when X is standard normal.
            name (string): Name of ψ-function used for printing.
            xtras (string): Potentially further information.
        """
        self.rho = rho
        self.psi = psi
        self.wgt = wgt
        self.Dpsi = Dpsi
        self.Dwgt = Dwgt
        self.tDefs = tDefs
        self.Erho = Erho
        self.Epsi2 = Epsi2
        self.EDpsi = EDpsi
        self.name = name
        self.xtras = xtras
    def __str__(self):
        """ 
        The function returns the name of the ψ-function class when you want to print the class.

        Returns:
            The name of the ψ-function class.        
        """
        return str(self.name) + ' psi function'