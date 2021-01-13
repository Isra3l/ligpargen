"""

Module with functions to generate Q software inputs (PRM, LIB, PDB, FEP)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os
import math

headerPRM = '''*---------------------------------------------
*        Generated with LigParGen 
*        William L. Jorgensen Lab
*     Author: israel.cabezadevaca@yale.edu
*    OPLS Force Field with CM1A derived Atomic Charges
*---------------------------------------------
[options]
vdw_rule     geometric !vwd combination rule
scale_14     0.5       !electrostatics 1-4 scaling factor
switch_atoms on
improper_potential      periodic
force_field AMBER
improper_definition explicit

[atom_types]
*tac--Avdw1--Avdw2--Bvdw1--Avdw3--Bvdw2&3--mass--comment
'''

headerLIB = '*Generated with LigParGen (israel.cabezadevaca@yale.edu)\n\n'


headerFEP = '''!info: %s
[FEP]
states 2
!softcore_use_max_potential on
offset_residue 1

'''

def printAtom(atom, alchemical=False):
    """Generate atom line

    Parameters
    ----------
    atom : Atom class
        Atom object
    alchemical : bool, optional
        If alchemical data, by default False

    Returns
    -------
    atomLine : str
        Atom line data
    """

    if alchemical:

        lj_A = 2.0*(atom.sigma_B**6)*math.sqrt(atom.epsilon_B)
        lj_B = 2.0*(atom.sigma_B**3)*math.sqrt(atom.epsilon_B)

        type_q = atom.type_q_B
        mass = atom.mass_B
        
    else:

        lj_A = 2.0*(atom.sigma**6)*math.sqrt(atom.epsilon)
        lj_B = 2.0*(atom.sigma**3)*math.sqrt(atom.epsilon)
        

        type_q = atom.type_q
        mass = atom.mass


    return '%-5s%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f ! %5s\n' % (type_q, 
            lj_A, lj_A, lj_B, lj_A/math.sqrt(2.0), lj_B/math.sqrt(2.0), mass, atom.resname)


def printAtomFEP(atom, alchemical=False):
    """Generate atom line if FEP

    Parameters
    ----------
    atom : Atom class
        Atom object
    alchemical : bool, optional
        If alchemical data, by default False

    Returns
    -------
    atomLine : str
        Atom line data
    """

    if alchemical:

        lj_A = 2.0*(atom.sigma_B**6)*math.sqrt(atom.epsilon_B)
        lj_B = 2.0*(atom.sigma_B**3)*math.sqrt(atom.epsilon_B)

        type_q = atom.type_q_B
        mass = atom.mass_B
        
    else:

        lj_A = 2.0*(atom.sigma**6)*math.sqrt(atom.epsilon)
        lj_B = 2.0*(atom.sigma**3)*math.sqrt(atom.epsilon)
        

        type_q = atom.type_q
        mass = atom.mass


    return '%-5s%12.4f%12.4f    0.0    0.0%12.4f%12.4f %12.4f ! %5s\n' % (type_q, 
            lj_A, lj_B, lj_A/math.sqrt(2.0), lj_B/math.sqrt(2.0), mass, atom.resname)


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

    k0 = bond.K0*2.0

    return '%-11s%-9s%12.1f%10.3f\n' % (bond.atomA.type_q, bond.atomB.type_q, k0, bond.R0)

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

    k0 = angle.K0*2.0

    return '%-11s%-11s%-8s%11.2f%12.3f\n' % (angle.atomA.type_q, angle.atomB.type_q, angle.atomC.type_q, k0, angle.angle0)


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

    V1 = dihedral.V1*0.5
    V2 = dihedral.V2*0.5
    V3 = dihedral.V3*0.5
    V4 = dihedral.V4*0.5

    label = '%-11s%-11s%-11s%-8s' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)

    torsion1 = '%s%8.3f     -1.000     0.000     1.000\n' % (label, V1)
    torsion2 = '%s%8.3f     -2.000   180.000     1.000\n' % (label, V2)
    torsion3 = '%s%8.3f     -3.000     0.000     1.000\n' % (label, V3)
    torsion4 = '%s%8.3f      4.000   180.000     1.000\n' % (label, V4)  

    return torsion1+torsion2+torsion3+torsion4


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

    V2 = dihedral.V2*0.5

    label = '%-11s%-11s%-11s%-8s' % (dihedral.atomA.type_q, dihedral.atomB.type_q, dihedral.atomC.type_q, dihedral.atomD.type_q)

    torsion2 = '%s%13.3f   180.000\n' % (label, V2)

    return torsion2



def writePRM(molecule, prmFile):
    """Generate PRM file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    prmFile : str
        PRM file name
    """

    with open(prmFile,'w') as ofile:

        ofile.write(headerPRM)

        qatomtypeInMoleculeLst = []

        for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]: 
            
            if atom.type_q in qatomtypeInMoleculeLst: continue

            ofile.write(printAtom(atom))
            qatomtypeInMoleculeLst.append(atom.type_q)


        ofile.write('\n[bonds]\n*iaci--iacj--force.c--dist.\n')

        for bond in molecule.bondsVariable: ofile.write(printBond(bond))
        for bond in molecule.bondsAdditional: ofile.write(printBond(bond))

        ofile.write('\n[angles]\n*iaci--iacj--iack--force.c--angle.\n')

        for angle in molecule.anglesVariable: ofile.write(printAngle(angle))
        for angle in molecule.anglesAdditional: ofile.write(printAngle(angle))

        ofile.write('\n[torsions]\n*iaci--iacj--iack--iackl--force.c--minima--phase--path.\n')

        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==False: ofile.write(printDihedral(dihedral))

        ofile.write('!  X    X    X    X    0.00000 1 0.000000 1 ! WILD CARD FOR MISSING TORSION PARAMETERS\n')

        ofile.write('\n[impropers]\n*iaci--iacj--iack--iackl--force.c--phase.\n')

        for dihedral in molecule.torsionsVariable: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))
        for dihedral in molecule.torsionsAdditional: 
            if dihedral.improper==True: ofile.write(printImproperDihedral(dihedral))

        ofile.write('!  X    X    X    X    0.00000 2 0.000000 1 ! WILD CARD FOR MISSING IMPROPER PARAMETERS\n')


def writeLIB(molecule, libFile):
    """Generate LIB file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    libFile : str
        LIB file name
    """


    with open(libFile,'w') as ofile:

        ofile.write(headerLIB)
        ofile.write('{'+molecule.atoms[0].resname+'}\n')
        ofile.write('[atoms]\n')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atomOriginalOrder in enumerate(atomsToWrite, start = 1):

            ofile.write('%4d   %-10s%-10s%9.4f\n' %(i, atomOriginalOrder.nameOriginal.upper(), atomOriginalOrder.type_q, atomOriginalOrder.charge))

        ofile.write('\n[bonds]\n')

        for bond in molecule.bondsVariable: ofile.write('%-10s%-10s\n' % (bond.atomA.nameOriginal.upper(), bond.atomB.nameOriginal.upper()))
        for bond in molecule.bondsAdditional: ofile.write('%-10s%-10s\n' % (bond.atomA.nameOriginal.upper(), bond.atomB.nameOriginal.upper()))

        ofile.write('\n[impropers]\n')

        for dihedral in molecule.torsionsVariable:
            if dihedral.improper==True:
                ofile.write('%-10s%-10s%-10s%-10s\n' % (dihedral.atomA.nameOriginal.upper(), dihedral.atomB.nameOriginal.upper(),
                            dihedral.atomC.nameOriginal.upper(),dihedral.atomD.nameOriginal.upper()))

        for dihedral in molecule.torsionsAdditional:
            if dihedral.improper==True:
                ofile.write('%-10s%-10s%-10s%-10s\n' % (dihedral.atomA.nameOriginal.upper(), dihedral.atomB.nameOriginal.upper(),
                            dihedral.atomC.nameOriginal.upper(),dihedral.atomD.nameOriginal.upper()))

        ofile.write('\n[charge_groups]\n')

        for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]: ofile.write('%-4s' % atom.nameOriginal.upper())

        ofile.write('\n\n')

def writeFEP(molecule, fepFile):
    """Generate FEP file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    fepFile : str
        FEP file name
    """

    with open(fepFile,'w') as ofile:

        ofile.write(headerFEP % molecule.atoms[0].resname.upper())

        ofile.write('[atoms]\n')

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atomOriginalOrder in enumerate(atomsToWrite, start = 1):

            ofile.write('%s %5d\n' % (str(i).ljust(3),atomOriginalOrder.serial-molecule.numberOfStructuralDummyAtoms))


        ofile.write('\n[PBC]\n')

        for i, atom in enumerate(molecule.atoms[molecule.numberOfStructuralDummyAtoms:], start=1): 
            ofile.write('switching_atom  %d\n' % (atom.serial-molecule.numberOfStructuralDummyAtoms))


        ofile.write('\n[change_charges]\n')

        for i, atom in enumerate(molecule.atoms[molecule.numberOfStructuralDummyAtoms:], start=1): 
            ofile.write('%s %10.4f %10.4f\n' % (str(i).ljust(3),atom.charge, atom.charge_B))


        ofile.write('\n[atom types]\n')


        qatomtypeInMoleculeLst = []

        for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]: 
            
            if atom.type_q in qatomtypeInMoleculeLst: continue

            ofile.write(printAtomFEP(atom))
            qatomtypeInMoleculeLst.append(atom.type_q)

        for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]: 

            if atom.type_q_B in qatomtypeInMoleculeLst: continue

            ofile.write(printAtomFEP(atom, True))

            qatomtypeInMoleculeLst.append(atom.type_q_B)

            
        ofile.write('\n[change atoms]\n')

        for i, atom in enumerate(molecule.atoms[molecule.numberOfStructuralDummyAtoms:], start=1): 
            ofile.write('%s %10s %10s\n' % (str(i).ljust(3),atom.type_q, atom.type_q_B))


        if len(molecule.excludedListDict)!=0: ofile.write('\n[excluded pairs]\n')

        for atomA, excludedAtoms in molecule.excludedListDict.items():
            for atomB in excludedAtoms:
                ofile.write('%5d%5d%5d%5d\n' % (atomA, atomB, 1, 1))


        totalBonds = molecule.bondsVariable + molecule.bondsAdditional

        bondTypes, bondChange = generateChangeList(totalBonds, molecule.numberOfStructuralDummyAtoms)

        if len(bondTypes)!=0:

            ofile.write('\n[bond_types]\n')

            for bond in bondTypes: ofile.write('%d %5.1f %5.1f\n' % (bond[0], bond[1], bond[2]))

            ofile.write('\n[change_bonds]\n')

            for bond in bondChange: ofile.write('%d %d %d %d\n' % (bond[0], bond[1], bond[2], bond[3]))


        totalAngles = molecule.anglesVariable + molecule.anglesAdditional

        angleTypes, angleChange = generateChangeList(totalAngles, molecule.numberOfStructuralDummyAtoms)

        if len(angleTypes)!=0:

            ofile.write('\n[angle_types]\n')

            for angle in angleTypes: ofile.write('%d %5.1f %5.1f\n' % (angle[0], angle[1], angle[2]))

            ofile.write('\n[change_angles]\n')

            for angle in angleChange: ofile.write('%d %d %d %d %d\n' % (angle[0], angle[1], angle[2], angle[3], angle[4]))

        totalTorsions = [torsion for torsion in molecule.torsionsVariable + molecule.torsionsAdditional if not torsion.improper] 

        torsionTypes, torsionChange = generateChangeList(totalTorsions, molecule.numberOfStructuralDummyAtoms)

        if len(torsionTypes)!=0:

            ofile.write('\n[torsion_types]\n')

            for torsion in torsionTypes: ofile.write('%d %5.1f %5.1f %5.1f %5.1f\n' % (torsion[0], torsion[1], torsion[2], torsion[3], torsion[4]))

            ofile.write('\n[change_torsions]\n')

            for torsion in torsionChange: ofile.write('%d %d %d %d %d\n' % (torsion[0], torsion[1], torsion[2], torsion[3], torsion[4]))


        totalTorsionsImproper = [torsion for torsion in molecule.torsionsVariable + molecule.torsionsAdditional if torsion.improper] 

        torsionTypes, torsionChange = generateChangeList(totalTorsionsImproper, molecule.numberOfStructuralDummyAtoms)

        if len(torsionTypes)!=0:

            ofile.write('\n[improper_types]\n')

            for torsion in torsionTypes: ofile.write('%d %5.1f %5.1f %5.1f %5.1f\n' % (torsion[0], torsion[1], torsion[2]))

            ofile.write('\n[change_impropers]\n')

            for torsion in torsionChange: ofile.write('%d %d %d %d %d\n' % (torsion[0], torsion[1], torsion[2], torsion[3], torsion[4]))



def generateChangeList(bondedTerms, numberOfStructuralDummyAtoms):

    termTypes = []
    termChange = []

    tmpType = []
    tmpChanges = []

    if len(bondedTerms[0].getAtoms())==2:

        for term in bondedTerms: 

            if term.R0!=term.R0_B or term.K0!=term.K0_B:

                tmpType.append((term.K0, term.R0))
                tmpType.append((term.K0_B, term.R0_B))

                tmpChanges.append([term.atomA.serial-numberOfStructuralDummyAtoms, term.atomB.serial-numberOfStructuralDummyAtoms, term.K0, term.R0, term.K0_B, term.R0_B])

        tmpType = list(dict.fromkeys(tmpType))

        termTypes = [[i]+ list(term) for i, term in enumerate(tmpType, start =1)]

        for item in tmpChanges:

            indexA = tmpType.index((item[2], item[3]))+1
            indexB = tmpType.index((item[4], item[5]))+1

            termChange.append([item[0], item[1], indexA, indexB])

    if len(bondedTerms[0].getAtoms())==3:
    
        for term in bondedTerms: 

            if term.angle0!=term.angle0_B or term.K0!=term.K0_B:

                tmpType.append((term.K0, term.angle0))
                tmpType.append((term.K0_B, term.angle0_B))

                tmpChanges.append([term.atomA.serial-numberOfStructuralDummyAtoms, term.atomB.serial-numberOfStructuralDummyAtoms, 
                    term.atomC.serial-numberOfStructuralDummyAtoms, term.K0, term.angle0, term.K0_B, term.angle0_B])

        tmpType = list(dict.fromkeys(tmpType))

        termTypes = [[i]+ list(term) for i, term in enumerate(tmpType, start =1)]

        for item in tmpChanges:

            indexA = tmpType.index((item[3], item[4]))+1
            indexB = tmpType.index((item[5], item[6]))+1

            termChange.append([item[0], item[1], item[2], indexA, indexB])


    if len(bondedTerms[0].getAtoms())==4:

        if bondedTerms[0].improper:

            for term in bondedTerms: 

                if term.V2!=term.V2_B:

                    tmpType.append((term.V2))
                    tmpType.append((term.V2_B))

                    tmpChanges.append([term.atomA.serial-numberOfStructuralDummyAtoms, term.atomB.serial-numberOfStructuralDummyAtoms,
                        term.atomC.serial-numberOfStructuralDummyAtoms,  term.atomD.serial-numberOfStructuralDummyAtoms, term.V2, term.V2_B])


            tmpType = list(dict.fromkeys(tmpType))

            termTypes = [[i]+ list(term) for i, term in enumerate(tmpType, start =1)]

            for item in tmpChanges:

                indexA = tmpType.index((item[4]))+1
                indexB = tmpType.index((item[5]))+1

                termChange.append([item[0], item[1], item[2], item[3], indexA, indexB])


        else:

            for term in bondedTerms: 

                if term.V1!=term.V1_B or term.V2!=term.V2_B or term.V3!=term.V3_B or term.V4!=term.V4_B:

                    tmpType.append((term.V1, term.V2, term.V3, term.V4))
                    tmpType.append((term.V1_B, term.V2_B, term.V3_B, term.V4_B))

                    tmpChanges.append([term.atomA.serial-numberOfStructuralDummyAtoms, term.atomB.serial-numberOfStructuralDummyAtoms, 
                    term.atomC.serial-numberOfStructuralDummyAtoms,  term.atomD.serial-numberOfStructuralDummyAtoms, term.V1, 
                    term.V2, term.V3, term.V4, term.V1_B, term.V2_B, term.V3_B, term.V4_B])


            tmpType = list(dict.fromkeys(tmpType))

            termTypes = [[i]+ list(term) for i, term in enumerate(tmpType, start =1)]

            for item in tmpChanges:

                indexA = tmpType.index((item[4], item[5], item[6], item[7]))+1
                indexB = tmpType.index((item[8], item[9], item[10], item[11]))+1

                termChange.append([item[0], item[1], item[2], item[3], indexA, indexB])



    return termTypes, termChange

def writePDB(molecule, pdbFile):
    """Generate PDB file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    pdbFile : str
        PDB file name
    """

    with open(pdbFile, 'w') as ofile:

        for i, atom in enumerate(molecule.atoms[molecule.numberOfStructuralDummyAtoms:], start=1):

            atomOriginalOrder = molecule.atoms[atom.serialOriginal -1]

            ofile.write('ATOM%7d%5s%4s%6d%12.3f%8.3f%8.3f  1.00  0.00%12s  \n' % (i, atomOriginalOrder.nameOriginal.upper(), atomOriginalOrder.resname,
                        1, atomOriginalOrder.x + molecule.shiftX, atomOriginalOrder.y + molecule.shiftY, atomOriginalOrder.z + molecule.shiftZ,
                        atomOriginalOrder.element))


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
    prmFile : str
        PRM file name
    libFile : str
        LIB file name
    pdbFile : str
        PDB file name
    fepFile : str
        FEP file name
    """

    prmFile = os.path.join(workdir,molname+'.q.prm')
    libFile = os.path.join(workdir,molname+'.q.lib')
    pdbFile = os.path.join(workdir,molname+'.q.pdb')
    fepFile = os.path.join(workdir,molname+'.q.fep')

    return prmFile, libFile, pdbFile, fepFile


def write(molecule, molName, workdir):
    """Write Q output formats

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    molName : str
        Molecule name
    workdir : str
        Working folder
    """

    prmFile, libFile, pdbFile, fepFile = getFileNames(molName, workdir)

    writePRM(molecule, prmFile)

    writeLIB(molecule, libFile)

    writePDB(molecule, pdbFile)

    if molecule.alchemicalTransformation: 
        
        writeFEP(molecule, fepFile)

