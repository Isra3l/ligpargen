"""

Torsion Class definition

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""
class Torsion(object):
    """
    Torsion Class
    """
    def __init__(self, atomA, atomB, atomC, atomD, typeInitial, typeFinal, V1, V2, V3, V4, phase1, phase2, phase3, displacement, improper=False):

        self.atomA = atomA
        self.atomB = atomB
        self.atomC = atomC
        self.atomD = atomD

        self.typeInitial = int(typeInitial)
        self.typeFinal = int(typeFinal)


        self.V1 = float(V1)
        self.V2 = float(V2)
        self.V3 = float(V3)
        self.V4 = float(V4)


        self.V1_B = float(V1)
        self.V2_B = float(V2)
        self.V3_B = float(V3)
        self.V4_B = float(V4)


        self.phase1 = float(phase1)
        self.phase2 = float(phase2)
        self.phase3 = float(phase3)

        self.displacement = float(displacement)

        self.improper = improper


    def __str__(self) -> str: 
        return ' '.join([str(item) for item in [self.atomA.serial,self.atomB.serial,self.atomC.serial,self.atomD.serial,
        self.V1,self.V2,self.V3,self.V4,self.V1_B,self.V2_B,self.V3_B,self.V4_B, self.typeInitial, self.typeFinal]])
    
    def getAtoms(self):
        return [self.atomA, self.atomB, self.atomC, self.atomD]

