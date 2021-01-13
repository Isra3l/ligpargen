"""

Angle Class definition

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""
class Angle(object):
    """
    docstring
    """
    def __init__(self, atomA, atomB, atomC, angle0, K0):

        self.atomA = atomA
        self.atomB = atomB
        self.atomC = atomC

        self.angle0 = float(angle0)
        self.K0 = float(K0)

        self.angle0_B = float(angle0)
        self.K0_B = float(K0)

    
    def getAtoms(self):
        return [self.atomA, self.atomB, self.atomC]


    def __eq__(self,other):

        if (self.atomA == other.atomA and self.atomB == other.atomB and self.atomC == other.atomC) or (self.atomA == other.atomC 
        and self.atomB == other.atomB and self.atomC == other.atomA):
            return True
        else: return False

    def __str__(self) -> str: 
        return ' '.join([str(item) for item in [self.atomA.serial,self.atomB.serial,self.atomC.serial,
        self.angle0,self.K0,self.angle0_B,self.K0_B]])

