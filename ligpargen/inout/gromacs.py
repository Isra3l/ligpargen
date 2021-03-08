"""

Module with functions to generate GROMACS software inputs (GRO, ITP)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os

headerITP = ''';---------------------------------------------
;        Generated with LigParGen 
;        William L. Jorgensen Lab
;     Author: israel.cabezadevaca@yale.edu
;    OPLS Force Field with CM1A derived Atomic Charges
;---------------------------------------------
[ defaults ]
; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
    1               3              yes            0.5     0.5

[ atomtypes ]
'''

headerGRO = 'Generated with LigParGen (israel.cabezadevaca@yale.edu)\n'

kcalToKj = 4.184


def printAtom( serial, molecule, atom, alchemicalTransformation):
    """Generate atom line   

    Parameters
    ----------
    serial : int
        Atom serial
    alchemicalTransformation : bool
        True if alchemical transformation

    Returns
    -------
    atomLine : str
        Atom line 
    """

    line = ''

    if alchemicalTransformation:
        line = ' %5d %10s %6d %6s %5s %6d %10.4f %10.4f%11s%11.4f%11.4f \n' % (serial, atom.type_gmx, 1, 
        molecule.residueName, atom.nameOriginal, 1, atom.charge, atom.mass, atom.type_gmx_B, atom.charge_B, atom.mass_B)
    else:
        line = ' %5d %10s %6d %6s %5s %6d %10.4f %10.4f \n' % (serial, atom.type_gmx, 1, 
        molecule.residueName, atom.nameOriginal, 1, atom.charge, atom.mass)

    return line

def printAtomType(atom, alchemicalTransformation):
    """Generate atom type line

    Parameters
    ----------
    atom : atom class
        atom class
    alchemicalTransformation : bool
        True if there is alchemical

    Returns
    -------
    atomLine : str 
        atomline
    """
    
    sigma = 0.1*atom.sigma
    epsilon = 4.184*atom.epsilon

    line = ''

    if alchemicalTransformation:

        sigma_B = 0.1*atom.sigma_B
        epsilon_B = 4.184*atom.epsilon_B

        line = '%10s %5s %10.4f     0.000    A    %10.5E   %10.5E\n' % (atom.type_gmx_B, atom.atomTypeOPLS_B, atom.mass_B, sigma_B, epsilon_B)
    else:

        line = '%10s %5s %10.4f     0.000    A    %10.5E   %10.5E\n' % (atom.type_gmx, atom.atomTypeOPLS, atom.mass, sigma, epsilon)

    return line

def printBond(bond, shift, molecule, alchemicalTransformation):
    """Generate bond line

    Parameters
    ----------
    bond : Bond Object
        Bond Object
    shift : int
        Shift produced by structural dummy atoms
    molecule : molecule object
        Molecule object
    alchemicalTransformation : bool
        True if alchemical transformation

    Returns
    -------
    bondLine : str
        Bond line data
    """


    k0 = bond.K0*836.80
    R0 = bond.R0*0.1
    ftype = 1

    line = ''

    atomAbond = molecule.atoms[bond.atomA.serialOriginal -1].serial-shift
    atomBbond = molecule.atoms[bond.atomB.serialOriginal -1].serial-shift


    if alchemicalTransformation:

        k0_B = bond.K0_B*836.80
        R0_B = bond.R0_B*0.1

        line = '%5d%5d%5d %11.4f %13.3f%12.4f%11.3f\n' % (atomAbond, atomBbond, ftype, R0, k0, R0_B, k0_B)
        # line = '%5d%5d%5d %11.4f %13.3f%12.4f%11.3f\n' % (bond.atomA.serial-shift, bond.atomB.serial-shift, ftype, R0, k0, R0_B, k0_B)

    else: line = '%5d%5d%5d %11.4f %13.3f\n' % (atomAbond, atomBbond, ftype, R0, k0)

    return line

def printAngle(angle, shift, molecule, alchemicalTransformation):
    """Generate angle line

    Parameters
    ----------
    angle : Angle Object
        Angle Object
    shift : int
        Shift produced by structural dummy atoms
    molecule : molecule object
        Molecule object
    alchemicalTransformation : bool
        True if alchemical transformation

    Returns
    -------
    angleLine : str
        Angle line data
    """

    
    k0 = angle.K0*8.3680
    angle0 = angle.angle0
    ftype = 1

    line = ''

    atomAangle = molecule.atoms[angle.atomA.serialOriginal -1].serial-shift
    atomBangle = molecule.atoms[angle.atomB.serialOriginal -1].serial-shift
    atomCangle = molecule.atoms[angle.atomC.serialOriginal -1].serial-shift


    if alchemicalTransformation: 

        k0_B = angle.K0_B*8.3680
        angle0_B = angle.angle0_B

        line = '%5d%5d%5d%5d    %10.3f %10.3f%11.3f%11.3f\n' % (atomAangle, atomBangle, 
            atomCangle, ftype, angle0, k0, angle0_B, k0_B)

    else: 
        
        line = '%5d%5d%5d%5d    %10.3f %10.3f\n' % (atomAangle, atomBangle, 
            atomCangle, ftype, angle0, k0)

    return line

def printDihedral(dihedral, shift, molecule, alchemicalTransformation):
    """Generate dihedral line

    Parameters
    ----------
    dihedral : Angle Object
        Angle Object
    shift : int
        Shift produced by structural dummy atoms
    molecule : molecule object
        Molecule object
    alchemicalTransformation : bool
        True if alchemical transformation

    Returns
    -------
    dihedralLine : str
        Dihedral line data
    """

    a = (dihedral.V2 + (dihedral.V1 + dihedral.V3) * 0.5)*kcalToKj

    V1 = (1.5*dihedral.V3 - 0.5*dihedral.V1)*kcalToKj
    V2 = (4.0*dihedral.V4 - dihedral.V2)*kcalToKj
    V3 = (-2.0*dihedral.V3)*kcalToKj
    V4 = (-4.0*dihedral.V4)*kcalToKj

    ftype = 3.0

    line = ''

    atomAdihedral = molecule.atoms[dihedral.atomA.serialOriginal -1].serial-shift
    atomBdihedral = molecule.atoms[dihedral.atomB.serialOriginal -1].serial-shift
    atomCdihedral = molecule.atoms[dihedral.atomC.serialOriginal -1].serial-shift
    atomDdihedral = molecule.atoms[dihedral.atomD.serialOriginal -1].serial-shift


    if alchemicalTransformation:

        a_B = (dihedral.V2_B + (dihedral.V1_B + dihedral.V3_B) * 0.5)*kcalToKj

        V1_B = (1.5*dihedral.V3_B - 0.5*dihedral.V1_B)*kcalToKj
        V2_B = (4.0*dihedral.V4_B - dihedral.V2_B)*kcalToKj
        V3_B = (-2.0*dihedral.V3_B)*kcalToKj
        V4_B = (-4.0*dihedral.V4_B)*kcalToKj

        line = '%5d%5d%5d%5d%5d%12.3f%8.3f%9.3f%8.3f%8.3f%8.3f%12.3f%8.3f%9.3f%8.3f%8.3f%8.3f\n' \
        % (atomAdihedral, atomBdihedral, atomCdihedral, atomDdihedral, ftype, a, V1, V2, V3, V4, 0.000, a_B, V1_B, V2_B, V3_B, V4_B, 0.000)

    else: line = '%5d%5d%5d%5d%5d%12.3f%8.3f%9.3f%8.3f%8.3f%8.3f\n' \
        % (atomAdihedral, atomBdihedral, atomCdihedral, atomDdihedral, ftype, a, V1, V2, V3, V4, 0.000)


    return line

def printImproperDihedral(dihedral, shift, molecule, alchemicalTransformation):
    """Generate improper dihedral line

    Parameters
    ----------
    dihedral : Angle Object
        Angle Object
    shift : int
        Shift produced by structural dummy atoms
    molecule : molecule object
        Molecule object
    alchemicalTransformation : bool
        True if alchemical transformation

    Returns
    -------
    dihedralLine : str
        Dihedral line data
    """


    a = 180.0
    V1 = (dihedral.V2*0.5)*kcalToKj
    V2 = 2.0
    ftype = 4

    line = ''

    atomAdihedral = molecule.atoms[dihedral.atomA.serialOriginal -1].serial-shift
    atomBdihedral = molecule.atoms[dihedral.atomB.serialOriginal -1].serial-shift
    atomCdihedral = molecule.atoms[dihedral.atomC.serialOriginal -1].serial-shift
    atomDdihedral = molecule.atoms[dihedral.atomD.serialOriginal -1].serial-shift


    if alchemicalTransformation:

        V1_B = (dihedral.V2_B*0.5)*kcalToKj

        line = '%5d%5d%5d%5d     %d     %6.3f   %6.3f  %6d    %6.3f   %6.3f  %6d  \n' % \
        (atomAdihedral, atomBdihedral, atomCdihedral, atomDdihedral, ftype, a, V1, V2, a, V1_B, V2)

    else: line = '%5d%5d%5d%5d     %d     %6.3f   %6.3f  %6d  \n' % \
        (atomAdihedral, atomBdihedral, atomCdihedral, atomDdihedral, ftype, a, V1, V2)

    return line


def print14list(pair14, shift, molecule):
    """Generate pair 14 list line

    Parameters
    ----------
    pair14 : Angle Object
        Angle Object
    shift : int
        Shift produced by structural dummy atoms
    molecule : molecule object
        Molecule object

    Returns
    -------
    pair14 : str
        Pair14 line data
    """

    ftype = 1

    atomApair = molecule.atoms[pair14[0].serialOriginal -1].serial-shift
    atomBdihedral = molecule.atoms[pair14[1].serialOriginal -1].serial-shift


    return '%5d%5d%5d\n' % (atomApair, atomBdihedral, ftype)

def writeITP(molecule, itpFile):
    """Generate ITP file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    itpFile : str
        ITP file name
    """

    serialAtomShift = molecule.numberOfStructuralDummyAtoms
    alchemicalTransformation = molecule.alchemicalTransformation

    remarks = {'molType':'\n[ moleculetype ]\n; Name               nrexcl\n%5s                  3\n\n[ atoms ]\n',
                'atoms':';   nr       type  resnr residue  atom   cgnr     charge       mass\n',
                'bonds':'\n[ bonds ]\n;  ai    aj   funct      c0       c1\n',
                'angles':'\n[angles]\n;  ai    aj   ak  funct      c0       c1 \n',
                'dihedrals':'[ dihedrals ]\n;%8s DIHEDRAL ANGLES\n;  ai    aj   ak   al funct      c0       c1      c2      c3       c4     c5\n',
                'pairs':'[ pairs ]\n;  ai    aj funct\n'}

    if alchemicalTransformation:

        remarks['atoms'] = ';   nr       type  resnr    res   atom   cgnr     charge       mass      typeB    chargeB      massB comments\n'
        remarks['bonds'] = '\n[ bonds ]\n;  ai    aj   funct      r0       k0            r0_B     k0_B\n'
        remarks['angles'] = '\n[angles]\n;  ai    aj   ak  funct     angle0       k0      angle0_B    k0_B \n'
        remarks['dihedrals'] = '[ dihedrals ]\n;%8s DIHEDRAL ANGLES\n;  ai    aj   ak   al funct      c0       c1     \
c2      c3       c4       c5      c0_B       c1_B      c2_B      c3_B       c4_B     c5_B\n'


    with open(itpFile,'w') as ofile:

        ofile.write(headerITP)

        for atom in molecule.atoms: 
            if atom.typeA!=-1: ofile.write(printAtomType(atom, False))

        for atom in molecule.atoms: 
            if atom.typeA!=-1 and atom.typeA!=atom.typeB: 
                ofile.write(printAtomType(atom, alchemicalTransformation))


        ofile.write(remarks['molType'] % molecule.residueName)

        ofile.write(remarks['atoms'])

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for serial, atomOriginalOrder in enumerate(atomsToWrite, start = 1):

            ofile.write(printAtom(serial, molecule, atomOriginalOrder, alchemicalTransformation))

        ofile.write(remarks['bonds'])

        for bond in molecule.bondsVariable: ofile.write(printBond(bond, serialAtomShift, molecule, alchemicalTransformation))
        for bond in molecule.bondsAdditional: ofile.write(printBond(bond, serialAtomShift, molecule, alchemicalTransformation))

        ofile.write(remarks['angles'])

        for angle in molecule.anglesVariable: ofile.write(printAngle(angle, serialAtomShift, molecule, alchemicalTransformation))
        for angle in molecule.anglesAdditional: ofile.write(printAngle(angle, serialAtomShift, molecule, alchemicalTransformation))

        ofile.write(remarks['dihedrals'] % 'PROPER')

        dihedrals = [dihedral for dihedral in molecule.torsionsVariable + molecule.torsionsAdditional if dihedral.improper==False and dihedral.typeInitial!=0]

        for dihedral in dihedrals: ofile.write(printDihedral(dihedral, serialAtomShift, molecule, alchemicalTransformation))

        ofile.write(remarks['dihedrals'] % 'IMPROPER')

        dihedralsImproper = [dihedral for dihedral in molecule.torsionsVariable + molecule.torsionsAdditional if dihedral.improper==True]

        for dihedral in dihedralsImproper: ofile.write(printImproperDihedral(dihedral, serialAtomShift, molecule, alchemicalTransformation))

        if len(molecule.excludedListDict)!=0: ofile.write('[ exclusions ]\n')

        for atomA, excludedAtoms in molecule.excludedListDict.items():
            ofile.write('%5d' % atomA +''.join(['%5d' % atom for atom in excludedAtoms])+'\n')

        ofile.write(remarks['pairs'])

        for pair14 in molecule.lst14pairs: ofile.write(print14list(pair14, serialAtomShift, molecule))


def writeGRO(molecule, groFile):
    """Generate GRO file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    groFile : str
        GRO file name
    """

    with open(groFile,'w') as ofile:

        ofile.write(headerGRO)
        ofile.write('%5d\n' % (len(molecule.atoms)-molecule.numberOfStructuralDummyAtoms))

        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atomOriginalOrder in enumerate(atomsToWrite, start = 1):

            ofile.write('%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n' % (1, molecule.residueName, atomOriginalOrder.nameOriginal, i, 0.1*(atomOriginalOrder.x + molecule.shiftX), 
                0.1*(atomOriginalOrder.y + molecule.shiftY), 0.1*(atomOriginalOrder.z + molecule.shiftZ)))

        ofile.write('%10.5f%10.5f%10.5f\n' % (10.0000, 10.0000, 10.0000))



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
    groFile : str
        GRO file name
    itpFile : str
        ITP file name
    """

    groFile = os.path.join(workdir,molname+'.gmx.gro')
    itpFile = os.path.join(workdir,molname+'.gmx.itp')

    return groFile, itpFile


def write(molecule, molName, workdir):
    """Write GROMACS output formats

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    molName : str
        Molecule name
    workdir : str
        Working folder
    """

    groFile, itpFile = getFileNames(molName, workdir)

    writeITP(molecule, itpFile)
    
    writeGRO(molecule, groFile)

