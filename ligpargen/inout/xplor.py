"""

Module with functions to generate CNS/X-PLOR software inputs (TOP, PARAM)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os
import math

headerTOPPARAM = 'Remarks generated with LigParGen (israel.cabezadevaca@yale.edu)\n\n'

kcalToKj = 4.184


def printBondedAtoms(bond, shiftdummySerial):
    """Generate bond atoms line

    Parameters
    ----------
    bond : Bond Object
        Bond Object
    shiftdummySerial : int
        Shift to remove structural dummy atoms and get index (should be 4)

    Returns
    -------
    bondLine : str
        Bond line data
    """

    return '<Bond from=\"%d\" to=\"%d\"/>\n' % (bond.atomA.serialOriginal- shiftdummySerial, bond.atomB.serialOriginal - shiftdummySerial)


def printBond(bond):
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

    R0 = 0.1*bond.R0
    K0 = bond.K0*836.80

    return '<Bond class1=\"%s\" class2=\"%s\" length=\"%0.6f\" k=\"%0.6f\"/>\n' % (bond.atomA.type_q, bond.atomB.type_q, R0, K0)



def printAngle(angle):
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
    
    k0 = angle.K0*8.3680
    angle0 = math.radians(angle.angle0)

    return '<Angle class1=\"%s\" class2=\"%s\" class3=\"%s\" angle=\"%0.6f\" k=\"%0.6f\"/>\n' % \
        (angle.atomA.type_q, angle.atomB.type_q, angle.atomC.type_q, angle0, k0)


def printDihedral(dihedral):
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

    V1 = dihedral.V1*0.5*kcalToKj
    V2 = dihedral.V2*0.5*kcalToKj
    V3 = dihedral.V3*0.5*kcalToKj
    V4 = dihedral.V4*0.5*kcalToKj

    label = '<Proper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\" ' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)
    label += 'k1=\"%0.6f\" k2=\"%0.6f\" k3=\"%0.6f\" k4=\"%0.6f\" ' % (V1, V2, V3, V4)
    label += 'periodicity1=\"1\" periodicity2=\"2\" periodicity3=\"3\" periodicity4=\"4\" '
    label += 'phase1=\"0.00\" phase2=\"3.141592653589793\" phase3=\"0.00\" phase4=\"3.141592653589793\"/>\n'

    return label


def printImproperDihedral(dihedral):
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

    V1 = dihedral.V1*0.5*kcalToKj
    V2 = dihedral.V2*0.5*kcalToKj
    V3 = dihedral.V3*0.5*kcalToKj
    V4 = dihedral.V4*0.5*kcalToKj

    label = '<Improper class1=\"%s\" class2=\"%s\" class3=\"%s\" class4=\"%s\" ' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)
    label += 'k1=\"%0.6f\" k2=\"%0.6f\" k3=\"%0.6f\" k4=\"%0.6f\" ' % (V1, V2, V3, V4)
    label += 'periodicity1=\"1\" periodicity2=\"2\" periodicity3=\"3\" periodicity4=\"4\" '
    label += 'phase1=\"0.00\" phase2=\"3.141592653589793\" phase3=\"0.00\" phase4=\"3.141592653589793\"/>\n'

    return label


def writeTOP(molecule, topFile):
    """Generate TOP file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    topFile : str
        TOP file name
    """

    with open(topFile,'w') as ofile:

        ofile.write(headerTOPPARAM)

        ofile.write('set echo=false end\n\nautogenerate angles=True dihedrals=True end\n\n{ atomType mass }\n')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for atom in atomsToWrite:

            ofile.write('MASS %3s %5.4f\n' %(atom.type_q, atom.mass))

        ofile.write('\nRESIdue %5s\n' % molecule.residueName)

        ofile.write('\nGROUP\n')

        ofile.write('\n{ atomName atomType Charge }\n')

        for atom in atomsToWrite:

            ofile.write('ATOM %6s TYPE= %6s CHARGE= %8.4f END\n' %(atom.nameOriginal, atom.type_q, atom.charge))

        ofile.write('\n{ Bonds: atomName1 atomName2 }\n')

        bonds = molecule.bondsVariable + molecule.bondsAdditional

        for bond in bonds: 
            ofile.write('BOND %s %s\n' %( bond.atomA.nameOriginal, bond.atomB.nameOriginal))

        ofile.write('\n{ Improper Dihedrals: aName1 aName2 aName3 aName4 }\n')

        torsions = molecule.torsionsVariable + molecule.torsionsAdditional
        impropers = [torsion for torsion in torsions if torsion.improper]

        for torsion in impropers: 
            ofile.write('IMPRoper %s %s %s %s\n' %( torsion.atomA.nameOriginal, torsion.atomB.nameOriginal, torsion.atomC.nameOriginal, torsion.atomD.nameOriginal)) 

        ofile.write('\nEND {RESIdue %s}\n' % molecule.atoms[0].resname)

        ofile.write('\nset echo=true end\n\n')


def writePARAM(molecule, paramFile):
    """Generate PARAM file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    paramFile : str
        PARAM file name
    """
    
    with open(paramFile, 'w') as ofile:

        ofile.write(headerTOPPARAM)

        ofile.write('set echo=false end\n\n')

        ofile.write('{ BOND: atomType1 atomType2 kb r0 }\n')

        bonds = molecule.bondsVariable + molecule.bondsAdditional

        for bond in bonds: 
            ofile.write('BOND   %s %5s %8.1f %8.4f\n' %( bond.atomA.type_q, bond.atomB.type_q, bond.K0, bond.R0))

        ofile.write('\n{ ANGLE: aType1 aType2 aType3 kt t0 }\n')

        angles = molecule.anglesVariable + molecule.anglesAdditional

        for angle in angles: 
            ofile.write('ANGLE %5s %5s %5s %8.1f %8.2f\n' %( angle.atomA.type_q, angle.atomB.type_q, angle.atomC.type_q, angle.K0, angle.angle0))


        ofile.write('\n{ Proper Dihedrals: aType1 aType2 aType3 aType4 kt period phase }\n')

        torsions = molecule.torsionsVariable + molecule.torsionsAdditional
        propers = [torsion for torsion in torsions if not torsion.improper]

        for torsion in propers: 

            ofile.write('DIHEDRAL %5s %5s %5s %5s MULT 4' %(torsion.atomA.type_q, torsion.atomB.type_q, torsion.atomC.type_q, torsion.atomD.type_q))

            ofile.write('%10.3f %4d %8.2f\n' % (torsion.V1/2.0, 1, 0.0))
            ofile.write('%49.3f %4d %8.2f\n' % (torsion.V2/2.0, 2, 180.0))
            ofile.write('%49.3f %4d %8.2f\n' % (torsion.V3/2.0, 3, 0.0))
            ofile.write('%49.3f %4d %8.2f\n' % (torsion.V4/2.0, 4, 180.0))


        ofile.write('\n{ Improper Dihedrals: aType1 aType2 aType3 aType4 kt period phase }\n')


        impropers = [torsion for torsion in torsions if torsion.improper]

        for torsion in impropers: 

            ofile.write('IMPROPER %5s %5s %5s %5s %8.5f %d %4.4f \n' %(torsion.atomA.type_q, torsion.atomB.type_q, torsion.atomC.type_q, torsion.atomD.type_q,\
                torsion.V2/2.0, 2, 180.0))


        ofile.write('\n{ Nonbonded: Type Emin sigma; (1-4): Emin/2 sigma }\n')


        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for atom in atomsToWrite:

            ofile.write('NONBONDED %5s %11.6f %11.6f %11.6f %11.6f\n' %(atom.type_q, atom.epsilon, atom.sigma, 0.5*atom.epsilon, atom.sigma))
 
        ofile.write('\nset echo=true end\n\n\n')


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
    topFile : str
        TOP file name
    paramFile : str
        PARAM file name
    """

    topFile = os.path.join(workdir,molname+'.xplor.top')
    paramFile = os.path.join(workdir,molname+'.xplor.param')

    return topFile, paramFile


def write(molecule, molName, workdir):

    topFile, paramFile = getFileNames(molName, workdir)

    writeTOP(molecule, topFile)

    writePARAM(molecule, paramFile)

