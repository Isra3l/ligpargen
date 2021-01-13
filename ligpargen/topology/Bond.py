"""

Bond Class definition

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""
class Bond(object):
    """
    docstring
    """
    def __init__(self, atomA, atomB, R0, K0):

        self.atomA = atomA
        self.atomB = atomB
        self.R0 = float(R0)
        self.K0 = float(K0)

        self.R0_B = float(R0)
        self.K0_B = float(K0)

    def getAtoms(self):
        return [self.atomA, self.atomB]


    def __str__(self) -> str: 
        return ' '.join([str(item) for item in [self.atomA.serial,self.atomB.serial,
        self.R0,self.K0,self.R0_B,self.K0_B]])

