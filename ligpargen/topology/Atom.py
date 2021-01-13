"""

Atom Class definition

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""
class Atom(object):
    """
    docstring
    """

    def __init__(self, index, serial, name, typeA, typeB, parent, parentParent, parentParentParent, r, angle, dihedral, resname, 
        resnum, atomicNumber, atomTypeOPLS, charge, sigma, epsilon, x, y, z, mass, element):

        self.serial = int(serial)
        self.index = int(index)

        self.serialOriginal = int(serial) # In case the order of the atoms is changed
        self.indexOriginal = int(index) # In case the order  of the atoms is changed

        self.name = name
        self.nameOriginal = name

        self.typeA = int(typeA)
        self.typeB = int(typeB)
        self.parent = int(parent)
        self.parentParent = int(parentParent)
        self.parentParentParent = int(parentParentParent)
        self.r = float(r)
        self.angle = float(angle)
        self.dihedral = float(dihedral)

        self.resname = resname
        self.resnum = int(resnum)
        
        self.atomicNumber = int(atomicNumber)
        self.atomTypeOPLS = atomTypeOPLS
        self.charge = float(charge)
        self.sigma = float(sigma)
        self.epsilon = float(epsilon)

        # B corresponds to the final state in casa of alchemical transformation

        self.atomicNumber_B = int(atomicNumber)
        self.atomTypeOPLS_B = atomTypeOPLS
        self.charge_B = float(charge)
        self.sigma_B = float(sigma)
        self.epsilon_B = float(epsilon)

        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

        self.mass = float(mass)
        self.element = element

        self.mass_B = float(mass)
        self.element_B = element

        self.dual_just_stateB = False


        if self.typeA !=-1: 

            if len(self.element)==1:
                self.type_q = element+str(typeA)
            else:
                self.type_q = element+str(int(typeA)-800)

            self.type_gmx = 'opls_'+str(typeA)
            self.type_charmm = element+str(typeA)

            self.type_q_B = element+str(typeA)
            self.type_gmx_B = 'opls_'+str(typeA)
        else:
            self.type_q = element
            self.type_gmx = 'opls_'+str(typeA)
            self.type_charmm = element+str(typeA)

            self.type_q_B = element
            self.type_gmx_B = 'opls_'+str(typeA)
            

        self.bondedAtoms = []

    @classmethod
    def fromZmat(cls, zmatAtomLine,  atomicNumber, atomTypeOPLS, charge, sigma, epsilon,x,y,z, mass, element):

        serial = zmatAtomLine[:4].strip()
        index = int(serial)-1
        name = zmatAtomLine[4:8].strip()
        typeA = zmatAtomLine[9:13].strip()
        typeB = zmatAtomLine[14:18].strip()
        parent = zmatAtomLine[19:23].strip()
        parentParent = zmatAtomLine[35:39].strip()
        parentParentParent = zmatAtomLine[51:55].strip()
        r = zmatAtomLine[23:35].strip()
        angle = zmatAtomLine[39:51].strip()
        dihedral = zmatAtomLine[55:67].strip()

        resname = zmatAtomLine[67:71].strip()
        resnum = zmatAtomLine[71:].strip()

        return cls(index, serial, name, typeA, typeB, parent, parentParent, parentParentParent, r, angle, dihedral, resname, 
        resnum, atomicNumber, atomTypeOPLS, charge, sigma, epsilon, x, y, z, mass, element)        


    @classmethod
    def fromRDkit(cls, indexa, a, b, c, d, r0, angle0, dihedral0, element, atomicNumber, mass):
        """
        Alternative constructor using RDkit information

        Parameters
        ----------
        indexa : int
            Atom index
        a : int
            Atom position
        b : int
            Parent atom
        c : int
            Parent of the parent atom
        d : int
            Parent of the parent of the parent atom
        r0 : Float
            Distance between a and b
        angle0 : Float
            Angle between a-b-c
        dihedral0 : Float
            Dihedral angle between a-b-c-d
        element : int
            Atomic element
        atomicNumber : int
            Atomic number
        mass : Float
            Atomic mass

        Returns
        -------
        Atom object
            An atom object generated using the atom class
        """

        serial = a
        name = element
        typeA = atomicNumber
        typeB = atomicNumber
        parent = b 
        parentParent = c
        parentParentParent = d
        r = r0
        angle = angle0
        dihedral = dihedral0

        resname = 'MOL'
        resnum = 1
        
        atomicNumber = 0
        atomTypeOPLS = 0
        charge = 0.0
        sigma = 0.0
        epsilon = 0.0

        # B corresponds to the final state in casa of alchemical transformation

        x = 0.0
        y = 0.0
        z = 0.0

        return cls(indexa, serial, name, typeA, typeB, parent, parentParent, parentParentParent, r, angle, dihedral, resname, 
            resnum, atomicNumber, atomTypeOPLS, charge, sigma, epsilon, x, y, z, mass, element)


    def __str__(self):

        data = [self.serial,self.name,self.typeA,self.typeB,self.parent,self.parentParent,
        self.parentParentParent,self.r,self.angle,self.dihedral,self.resname, self.resnum,
        self.atomicNumber,self.atomicNumber_B,self.atomTypeOPLS,self.atomTypeOPLS_B,self.charge,self.sigma,
        self.epsilon,self.x,self.y,self.z, self.mass, self.element, self.type_q, self.type_q_B]

        data = [str(ele) for ele in data]

        return '  '.join(data)

