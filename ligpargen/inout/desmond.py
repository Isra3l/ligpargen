"""

Module with functions to generate DESMOND software inputs (CMS)

Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""

import os

part_0 = '''{
  s_m_m2io_version
  :::
  2.0.0
}
f_m_ct {
  s_m_title
  r_chorus_box_ax
  r_chorus_box_ay
  r_chorus_box_az
  r_chorus_box_bx 
  r_chorus_box_by
  r_chorus_box_bz
  r_chorus_box_cx 
  r_chorus_box_cy
  r_chorus_box_cz
  s_ffio_ct_type
  :::
  "Generated with LigParGen (israel.cabezadevaca@yale.edu)" 
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
                   0.0
                   0.0
                   0.0
                  70.0
  full_system
  m_atom[NUMBEROFATOMS] {
    i_m_mmod_type
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord 
    i_m_residue_number
    s_m_pdb_residue_name
    i_m_atomic_number 
    s_m_atom_name
    r_ffio_x_vel
    r_ffio_y_vel
    r_ffio_z_vel
    :::

'''

part_1 = '''  m_bond[NUMBEROFBONDS] {
    i_m_from
    i_m_to
    i_m_order
    i_m_from_rep 
    i_m_to_rep
    :::
'''


part_2 = '''      
f_m_ct {
  s_m_title
  r_chorus_box_ax
  r_chorus_box_ay
  r_chorus_box_az
  r_chorus_box_bx 
  r_chorus_box_by
  r_chorus_box_bz
  r_chorus_box_cx 
  r_chorus_box_cy
  r_chorus_box_cz
  s_ffio_ct_type
  :::
  RESNAME 
                70.0
                 0.0
                 0.0
                 0.0
                70.0
                 0.0
                 0.0
                 0.0
                70.0
  solute
  m_atom[NUMBEROFATOMS] {
    i_m_mmod_type
    r_m_x_coord
    r_m_y_coord
    r_m_z_coord 
    i_m_residue_number
    s_m_pdb_residue_name
    i_m_atomic_number 
    s_m_atom_name
    r_ffio_x_vel
    r_ffio_y_vel
    r_ffio_z_vel
    :::
'''


part_3 = '''  ffio_ff {
    s_ffio_name
    s_ffio_comb_rule
    i_ffio_version
    :::
    "RESNAME"
    GEOMETRIC
    1.0.0
    ffio_vdwtypes[VDWTYPENUMBER] {
      s_ffio_name
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''

part_4 = '''    ffio_sites[NUMBEROFATOMS] {
      s_ffio_type
      r_ffio_charge
      r_ffio_mass
      s_ffio_vdwtype
      i_ffio_resnr
      s_ffio_residue
      :::
'''

part_5 = '''     ffio_bonds[NUMBEROFBONDS] {
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''

part_6 = '''    ffio_angles[NUMBEROFANGLES] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''

part_7 = '''    ffio_dihedrals[NUMBEROFTORSIONS] {
      i_ffio_ai
      i_ffio_aj
      i_ffio_ak
      i_ffio_al
      s_ffio_funct
      r_ffio_c0
      r_ffio_c1
      r_ffio_c2
      r_ffio_c3
      r_ffio_c4
      :::
'''


part_7_2 = '''    ffio_torsion_torsion[0] {
      :::
      :::
    }
'''


part_8 = '''    ffio_exclusions[NUMBEROFEXCLUDED] {
      i_ffio_ai
      i_ffio_aj
      :::
'''


part_9 = '''    ffio_pairs[NUMBEROF14PAIRS] {
      i_ffio_ai
      i_ffio_aj
      s_ffio_funct
      r_ffio_c1
      r_ffio_c2
      :::
'''

part_10 = '''    ffio_constraints[0] {
      :::
      :::
    }
  }
}


'''


def writeCMS(molecule, cmsFile):
    """Generate CMS file

    Parameters
    ----------
    molecule : molecule class
        Molecule class
    cmsFile : str
        CMS file name
    """

    numberOfAtoms = len(molecule.atoms[molecule.numberOfStructuralDummyAtoms:])

    with open(cmsFile,'w') as ofile:

        ofile.write(part_0.replace('NUMBEROFATOMS', str(len(molecule.atoms[molecule.numberOfStructuralDummyAtoms:]))))


        atomsToWrite = sorted([atom for atom in molecule.atoms[molecule.numberOfStructuralDummyAtoms:]], key = lambda x: x.serialOriginal)

        for i, atom in enumerate(atomsToWrite, start =1):

            ofile.write('%7d%3d%8.4f%8.4f%8.4f 1 \"%s\"%3d \"%s\" %8.4f%8.4f%8.4f\n' %(i, 1, atom.x, atom.y, atom.z, atom.resname, atom.atomicNumber, atom.nameOriginal, 0.0, 0.0, 0.0))

        ofile.write('    :::\n  }\n')

        numberofbonds = len(molecule.bondsVariable) + len(molecule.bondsAdditional)

        ofile.write(part_1.replace('NUMBEROFBONDS',str(numberofbonds)))

        bonds = molecule.bondsVariable + molecule.bondsAdditional

        shift = molecule.numberOfStructuralDummyAtoms

        for i, bond in enumerate(bonds):

            ofile.write('%7d%4d%4d%4d%4d%4d\n' % (i+1, bond.atomA.serialOriginal - shift, bond.atomB.serialOriginal - shift, 1, 1, 1))

        ofile.write('    :::\n  }\n}\n')

        ofile.write(part_2.replace('NUMBEROFATOMS', str(numberOfAtoms)).replace('RESNAME', molecule.atoms[0].resname))

        for i, atom in enumerate(atomsToWrite, start =1):
    
            ofile.write('%7d%4d%8.4f%8.4f%8.4f 1 \"%s\"%4d \"%s\" %8.4f%8.4f%8.4f\n' %(i, 1, atom.x, atom.y, atom.z, atom.resname, atom.atomicNumber, atom.nameOriginal, 0.0, 0.0, 0.0))

        ofile.write('    :::\n  }\n')

        ofile.write(part_1.replace('NUMBEROFBONDS',str(numberofbonds)))

        bonds = molecule.bondsVariable + molecule.bondsAdditional

        shift = molecule.numberOfStructuralDummyAtoms

        for i, bond in enumerate(bonds):

            ofile.write('%7d%4d%4d%4d%4d%4d\n' % (i+1, bond.atomA.serialOriginal - shift, bond.atomB.serialOriginal - shift, 1, 1, 1))

        ofile.write('    :::\n  }\n}\n')

        ofile.write(part_3.replace('VDWTYPENUMBER', str(numberOfAtoms)).replace('RESNAME', molecule.atoms[0].resname))

        for i, atom in enumerate(atomsToWrite, start =1):
        
            ofile.write('%9d%10s LJ12_6_sig_epsilon %14.8f%14.8f\n' %(i, atom.type_gmx, atom.sigma, atom.epsilon))

        ofile.write('    :::\n  }\n')

        ofile.write(part_4.replace('NUMBEROFATOMS', str(numberOfAtoms)))

        for i, atom in enumerate(atomsToWrite, start =1):
        
            ofile.write('%9d  atom %14.8f%14.8f%10s 1 %4s\n' % (i, atom.charge, atom.mass, atom.type_gmx, atom.resname))

        ofile.write('      :::\n     }\n')

        ofile.write(part_5.replace('NUMBEROFBONDS', str(len(bonds))))

        for i, bond in enumerate(bonds):
    
            ofile.write('%9d%4d%4d HARM  %14.8f%14.8f\n' % (i+1, bond.atomA.serialOriginal - shift, bond.atomB.serialOriginal - shift, bond.R0, bond.K0))

        ofile.write('    :::\n  }\n')

        angles = molecule.anglesVariable + molecule.anglesAdditional

        ofile.write(part_6.replace('NUMBEROFANGLES', str(len(angles))))

        for i, angle in enumerate(angles):
        
            ofile.write('%9d%4d%4d%4d HARM  %14.8f%14.8f\n' % (i+1, angle.atomA.serialOriginal - shift, angle.atomB.serialOriginal - shift, angle.atomC.serialOriginal - shift, \
                 angle.angle0, angle.K0))

        ofile.write('    :::\n  }\n')


        torsions = molecule.torsionsVariable + molecule.torsionsAdditional

        ofile.write(part_7.replace('NUMBEROFTORSIONS', str(len(torsions))))

        for i, torsion in enumerate(torsions):

            if torsion.improper==False:
        
                ofile.write('%9d%4d%4d%4d%4d Opls_proper    0   %14.8f%14.8f%14.8f%14.8f\n' % (i+1, torsion.atomA.serialOriginal - shift, \
                    torsion.atomB.serialOriginal - shift, torsion.atomC.serialOriginal - shift, torsion.atomD.serialOriginal - shift, \
                    torsion.V1, torsion.V2, torsion.V3, torsion.V4))

            if torsion.improper==True:
            
                ofile.write('%9d%4d%4d%4d%4d Opls_improper    0   %14.8f%14.8f%14.8f%14.8f\n' % (i+1, torsion.atomA.serialOriginal - shift, \
                    torsion.atomB.serialOriginal - shift, torsion.atomC.serialOriginal - shift, torsion.atomD.serialOriginal - shift, \
                    torsion.V1, torsion.V2, torsion.V3, torsion.V4))

        ofile.write('    :::\n  }\n')

        ofile.write(part_7_2)

        exclusionList = getExclusionList(molecule)

        ofile.write(part_8.replace('NUMBEROFEXCLUDED', str(len(exclusionList))))

        for i, excludedPair in enumerate(exclusionList, start=1):

            ofile.write('%9d%4d%4d\n' % (i, excludedPair[0], excludedPair[1]))
        
        ofile.write('    :::\n  }\n')


        ofile.write(part_9.replace('NUMBEROF14PAIRS', str(2*len(molecule.lst14pairs))))

        i = 1
        for pair14 in molecule.lst14pairs:

            ofile.write('%9d%4d%4d Coulomb 0.5 <>\n' % (i, pair14[0].serialOriginal - shift, pair14[1].serialOriginal - shift))
            ofile.write('%9d%4d%4d LJ      0.5 <>\n' % (i + 1, pair14[0].serialOriginal - shift, pair14[1].serialOriginal - shift))

            i +=2

        ofile.write('    :::\n   }\n')

        ofile.write(part_10)


def getExclusionList(molecule):

    serialShift = molecule.numberOfStructuralDummyAtoms

    pairSet = set()

    bonds = molecule.bondsVariable + molecule.bondsAdditional

    for bond in bonds: pairSet.add(tuple(sorted([bond.atomA.serialOriginal-serialShift, bond.atomB.serialOriginal-serialShift])))

    angles = molecule.anglesVariable + molecule.anglesAdditional

    for angle in angles: pairSet.add(tuple(sorted([angle.atomA.serialOriginal-serialShift, angle.atomC.serialOriginal-serialShift])))


    torsions = molecule.torsionsVariable + molecule.torsionsAdditional

    for torsion in torsions: pairSet.add(tuple(sorted([torsion.atomA.serialOriginal-serialShift, torsion.atomD.serialOriginal-serialShift])))

    return list(pairSet)

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
    cmsFile : str
        CMS file name
    """

    return os.path.join(workdir,molname+'.desmond.cms')


def write(molecule, molName, workdir):

    cmsFile = getFileNames(molName, workdir)

    writeCMS(molecule, cmsFile)

