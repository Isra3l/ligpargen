"""

Module with functions to generate TINKER software inputs (XYZ, KEY)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os

headerKEY = '''!---------------------------------------------
!        Generated with LigParGen 
!        William L. Jorgensen Lab
!     Author: israel.cabezadevaca@yale.edu
!    OPLS Force Field with CM1A derived Atomic Charges
!---------------------------------------------
#Force Field Definition

forcefield              OPLS-AA

vdwindex                TYPE
vdwtype                 LENNARD-JONES
radiusrule              GEOMETRIC
radiustype              SIGMA
radiussize              DIAMETER
epsilonrule             GEOMETRIC
torsionunit             1.0
imptorunit              1.0
vdw-14-scale            2.0
chg-14-scale            2.0
electric                332.06
dielectric              1.0

'''


def printBond(bond, alchemical = False):
    """Generate bond line

    Parameters
    ----------
    bond : Bond Object
        Bond Object

    Returns
    -------
    bondLine : str
        Bond line data
    """

    label = 'BOND %10s %5s %8.1f %8.4f\n' % (bond.atomA.typeA, bond.atomB.typeA, bond.K0, bond.R0)

    if alchemical: label = 'BOND %10s %5s %8.1f %8.4f\n' % (bond.atomA.typeB, bond.atomB.typeB, bond.K0_B, bond.R0_B)

    return label

def printAngle(angle, alchemical = False):
    """Generate angle line

    Parameters
    ----------
    angle : Angle Object
        Angle Object

    Returns
    -------
    angleLine : str
        Angle line data
    """
    
    label = 'angle %9s %5s %5s %8.1f %8.2f\n' % (angle.atomA.typeA, angle.atomB.typeA, angle.atomC.typeA, angle.K0, angle.angle0)

    if alchemical: label = 'angle %9s %5s %5s %8.1f %8.2f\n' % (angle.atomA.typeB, angle.atomB.typeB, angle.atomC.typeB, angle.K0_B, angle.angle0_B)

    return label

def printDihedral(dihedral, alchemical = False):
    """Generate dihedral line

    Parameters
    ----------
    dihedral : dihedral Object
        dihedral Object

    Returns
    -------
    dihedralLine : str
        dihedral line data
    """

    V1 = dihedral.V1*0.5
    V2 = dihedral.V2*0.5
    V3 = dihedral.V3*0.5

    V1_B = dihedral.V1_B*0.5
    V2_B = dihedral.V2_B*0.5
    V3_B = dihedral.V3_B*0.5

    g = [0.0, 180.0, 0.0]
    n = [1, 2, 3]

    label = 'torsion %7s %5s %5s %5s %6.3f %4.1f %2d %6.3f %4.1f %2d %6.3f %4.1f %2d\n' % \
        (dihedral.atomA.typeA, dihedral.atomB.typeA, dihedral.atomC.typeA, dihedral.atomD.typeA, V1, g[0], n[0], V2, g[1], n[1], V3, g[2], n[2])

    if alchemical: label = 'torsion %7s %5s %5s %5s %6.3f %4.1f %2d %6.3f %4.1f %2d %6.3f %4.1f %2d\n' % \
        (dihedral.atomA.typeB, dihedral.atomB.typeB, dihedral.atomC.typeB, dihedral.atomD.typeB, V1_B, g[0], n[0], V2_B, g[1], n[1], V3_B, g[2], n[2])


    return label


def printImproperDihedral(dihedral, alchemical = False):
    """Generate improper dihedral line

    Parameters
    ----------
    dihedral : dihedral Object
        dihedral Object

    Returns
    -------
    dihedralLine : str
        Improper dihedral line data
    """

    V2 = dihedral.V2*0.5
    V2_B = dihedral.V2_B*0.5

    label = 'imptors %7s %5s %5s %5s %8.3f %4.1f %2d\n' % \
        (dihedral.atomA.typeA, dihedral.atomB.typeA, dihedral.atomC.typeA, dihedral.atomD.typeA, V2, 180.0, 2)


    if alchemical: label = 'imptors %7s %5s %5s %5s %8.3f %4.1f %2d\n' % \
        (dihedral.atomA.typeB, dihedral.atomB.typeB, dihedral.atomC.typeB, dihedral.atomD.typeB, V2_B, 180.0, 2)

    return label


def writeKEY(molecule, keyFile):
    """Generate KEY file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    prmFile : str
        KEY file name
    """

    with open(keyFile,'w') as ofile:

        ofile.write(headerKEY)

        ofile.write('''
        #############################
        ##                         ##
        ##  Atom Type Definitions  ##
        ##                         ##
        #############################
''')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write('atom %10s %6d %5s  \"%s\" %10d %10.3f %5d\n' %(atom.typeA, atom.typeB, atom.element, atom.nameOriginal, atom.atomicNumber, atom.mass, len(atom.bondedAtoms)))

        ofile.write('''
        ################################
        ##                            ##
        ##  Van der Waals Parameters  ##
        ##                            ##
        ################################
''')

        for i, atom in enumerate(atomsToWrite, start =1): ofile.write('vdw %11s %11.4f %8.4f\n' % (atom.typeA, atom.sigma, atom.epsilon))

        if molecule.alchemicalTransformation: 
            for i, atom in enumerate(atomsToWrite, start =1): ofile.write('vdw %11s %11.4f %8.4f\n' % (atom.typeB, atom.sigma_B, atom.epsilon_B))


        ofile.write('''
        ########################################
        ##                                    ##
        ##  Atomic Partial Charge Parameters  ##
        ##                                    ##
        ########################################
''')

        for i, atom in enumerate(atomsToWrite, start =1): ofile.write('charge %8s %11.4f\n' % (atom.typeA, atom.charge))

        if molecule.alchemicalTransformation: 
            for i, atom in enumerate(atomsToWrite, start =1): ofile.write('charge %8s %11.4f\n' % (atom.typeB, atom.charge_B))


        ofile.write('''
        ##################################
        ##                              ##
        ##  Bond Stretching Parameters  ##
        ##                              ##
        ##################################
''')
              
        for bond in molecule.bondsVariable: ofile.write(printBond(bond))
        for bond in molecule.bondsAdditional: ofile.write(printBond(bond))

        if molecule.alchemicalTransformation: 

            for bond in molecule.bondsVariable: ofile.write(printBond(bond, True))
            for bond in molecule.bondsAdditional: ofile.write(printBond(bond, True))


        ofile.write('''
        ################################
        ##                            ##
        ##  Angle Bending Parameters  ##
        ##                            ##
        ################################
''')

        for angle in molecule.anglesVariable: ofile.write(printAngle(angle))
        for angle in molecule.anglesAdditional: ofile.write(printAngle(angle))

        if molecule.alchemicalTransformation: 

            for angle in molecule.anglesVariable: ofile.write(printAngle(angle, True))
            for angle in molecule.anglesAdditional: ofile.write(printAngle(angle, True))


        ofile.write('''
        ############################
        ##                        ##
        ##  Torsional Parameters  ##
        ##                        ##
        ############################
''')

        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))

        if molecule.alchemicalTransformation: 

            for dihedral in molecule.torsionsVariable: 
                if dihedral.improper==False: ofile.write(printDihedral(dihedral, True))
            for dihedral in molecule.torsionsAdditional: 
                if dihedral.improper==False: ofile.write(printDihedral(dihedral, True))


        ofile.write('''
        #####################################
        ##                                 ##
        ##  Improper Torsional Parameters  ##
        ##                                 ##
        #####################################
''')

        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))

        if molecule.alchemicalTransformation: 

            for dihedral in molecule.torsionsVariable: 
                if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral, True))
            for dihedral in molecule.torsionsAdditional: 
                if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral, True))


def writeXYZ(molecule, xyzFile):
    """Generate XYZ file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    rtfFile : str
        XYZ file name
    """

    with open(xyzFile,'w') as ofile:

        ofile.write('%6d\n' % len(molecule.atoms[molecule.numberOfStructuralDummyAtoms:]))

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            bondedAtoms = []

            for atomBonded in atom.bondedAtoms:

                atomIndex = atomBonded -1 - molecule.numberOfStructuralDummyAtoms
                atomSerial = molecule.atoms[atomIndex].serial

                bondedAtoms.append(atomSerial)

            bondedAtoms = ''.join(['%5d' % atomBonded for atomBonded in bondedAtoms])

            ofile.write('%6d %4s %14.6f %14.6f %14.6f %6s%s\n' %(i, atom.nameOriginal, atom.x, atom.y, atom.z, atom.typeA, bondedAtoms))



def getFileNames(molname, workdir):
    """Return output file names

    Parameters
    ----------
    molname : str
        Molecule name
    workdir : str
        Working folder path

    Returns
    -------
    xyzFile : str
        XYZ file name
    keyFile : str
        KEY file name
    """

    xyzFile = os.path.join(workdir,molname+'.tinker.xyz')
    keyFile = os.path.join(workdir,molname+'.tinker.key')

    return xyzFile, keyFile


def write(molecule, molName, workdir):

    xyzFile, keyFile = getFileNames(molName, workdir)

    writeXYZ(molecule, xyzFile)
    
    writeKEY(molecule, keyFile)

    #TODO: Add exclusion list if it is needed!!!!