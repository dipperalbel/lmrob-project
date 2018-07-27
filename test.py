class psi_func:
    def __init__(self, rho, psi,wgt, Dpsi, Dwgt, tDefs, Erho, Epsi2, EDpsi, name, xtras):
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
        return str(self.name) + ' psi function'
    

    