import unittest
import os
import logging
import shutil

import copy

import ligpargen.tools.utilities as utilities
import ligpargen.topology.Molecule as molecule
import ligpargen.tools.boss as boss 
import ligpargen.inout.zmat as zmat
import ligpargen.inout.gromacs as gromacs
import ligpargen.inout.Q as Q
from ligpargen.inout import charmm
import ligpargen.tools.alchemify as alchemify

class test_ligpargen(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.workdir = os.path.join(os.getcwd(), 'test_ligpargen')
        cls.molname = 'MOL'

        logging.disable(logging.NOTSET)

        cls.benzene_RDkit, cls.newIndexToOriginalIndex, cls.atomsNameOriginal, cls.residueNameOriginal = utilities.generateRDkitMolecule(None, 'c1ccccc1', cls.workdir, cls.molname)

        cls.molecule_benzene_fromRDkit = molecule.fromRDkitMolecule(cls.benzene_RDkit)

        cls.goldenDataPath = os.path.join(os.getcwd(),'test_ligpargen','goldenData')

        zmatFile_benzene = os.path.join(cls.goldenDataPath, 'benzene_after_BOSS.z')
        outFile_benzene = os.path.join(cls.goldenDataPath, 'benzene_after_BOSS_out')
        pdbFile_benzene = os.path.join(cls.goldenDataPath, 'benzene_after_BOSS.pdb')

        cls.molecule_benzene = molecule.fromBOSS(zmatFile_benzene, outFile_benzene, pdbFile_benzene)

        cls.phenol_RDkit, cls.newIndexToOriginalIndex_phenol, cls.atomsNameOriginal_phenol, cls.residueNameOriginal_phenol = utilities.generateRDkitMolecule(None, 'c1ccc(cc1)O', cls.workdir, cls.molname)

        cls.molecule_phenol_fromRDkit = molecule.fromRDkitMolecule(cls.phenol_RDkit)

        zmatFile_phenol = os.path.join(cls.goldenDataPath, 'phenol_after_BOSS.z')
        outFile_phenol = os.path.join(cls.goldenDataPath, 'phenol_after_BOSS_out')
        pdbFile_phenol = os.path.join(cls.goldenDataPath, 'phenol_after_BOSS.pdb')

        cls.molecule_phenol = molecule.fromBOSS(zmatFile_phenol, outFile_phenol, pdbFile_phenol)

    @classmethod
    def tearDownClass(cls):

        os.remove(os.path.join(cls.workdir, 'MOL-debug.pdb'))

    def test_RDkitMolecule(self):

        atomNames_mol1, atomAtomicnumber_mol1, atomPositions = getAtomProperties(self.benzene_RDkit)

        golden_mol1_atomNames = ['C', 'C', 'C', 'C', 'C', 'C', 'H', 'H', 'H', 'H', 'H', 'H']
        golden_mol1_atomAtomicnumber = [6, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 1]
        golden_mol1_atomPositions = [[1.298, 0.393, -0.042], [0.304, 1.343, -0.014], [1.021, -0.963, -0.029], [-0.299, -1.351, 0.014], [-1.306, -0.403, 0.043], [-1.031, 0.95, 0.029], [-0.503, -2.402, 0.024], [-2.338, -0.724, 0.077], [-1.83, 1.679, 0.052], [0.506, 2.412, -0.024], [2.329, 0.722, -0.076], [1.849, -1.656, -0.053]]

        golden_newIndexToOriginalIndex = {0: 0, 1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10, 11: 11}
        golden_atomsNameOriginal = []
        golden_residueNameOriginal = []


        self.assertEqual(atomNames_mol1, golden_mol1_atomNames)
        self.assertEqual(atomAtomicnumber_mol1, golden_mol1_atomAtomicnumber)
        self.assertEqual(atomPositions, golden_mol1_atomPositions)
        self.assertEqual(self.newIndexToOriginalIndex, golden_newIndexToOriginalIndex)
        self.assertEqual(self.atomsNameOriginal, golden_atomsNameOriginal)
        self.assertEqual(self.residueNameOriginal, golden_residueNameOriginal)

    def test_bonds(self):

        bondsVariables = [[bond.atomA.serial, bond.atomB.serial] for bond in self.molecule_benzene_fromRDkit.bondsVariable]
        golden_mol1_bondsVariables = [[5, 4], [6, 4], [7, 6], [8, 7], [9, 8], [10, 7], [11, 8], [12, 9], [13, 5], [14, 4], [15, 6]]

        golden_mol1_bondsAdditionals = [[9, 5]]
        bondsAdditionals = [[bond.atomA.serial, bond.atomB.serial] for bond in self.molecule_benzene_fromRDkit.bondsAdditional]


        self.assertEqual(bondsVariables, golden_mol1_bondsVariables)
        self.assertEqual(bondsAdditionals, golden_mol1_bondsAdditionals)


        ethane_RDkit,  _ , _ , _ = utilities.generateRDkitMolecule(None, 'CC', self.workdir, self.molname)

        molecule_ethane = molecule.fromRDkitMolecule(ethane_RDkit)

        bondsVariables = [[bond.atomA.serial, bond.atomB.serial] for bond in molecule_ethane.bondsVariable]

        golden_ethane_bondsVariables = [[5, 4], [6, 4], [7, 4], [8, 4], [9, 5], [10, 5], [11, 5]]

        bondsAdditionals = [[bond.atomA.serial, bond.atomB.serial] for bond in molecule_ethane.bondsAdditional]

        golden_ethane_bondsAdditionals = []


        self.assertEqual(bondsVariables, golden_ethane_bondsVariables)
        self.assertEqual(bondsAdditionals, golden_ethane_bondsAdditionals)


    def test_angles(self):

        anglesVariables = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial] for angle in self.molecule_benzene_fromRDkit.anglesVariable]
        golden_mol1_anglesVariables = [[6, 4, 5], [7, 6, 4], [8, 7, 6], [9, 8, 7], [10, 7, 6], [11, 8, 7], [12, 9, 8], [13, 5, 4], [14, 4, 5], [15, 6, 4]]

        anglesAdditionals = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial] for angle in self.molecule_benzene_fromRDkit.anglesAdditional]
        golden_mol1_anglesAdditionals = [[6, 4, 14], [7, 6, 15], [8, 7, 10], [9, 8, 11], [9, 5, 4], [9, 5, 13], [12, 9, 5], [5, 9, 8]]

        self.assertEqual(anglesVariables, golden_mol1_anglesVariables)
        self.assertEqual(anglesAdditionals, golden_mol1_anglesAdditionals)

    def test_dihedrals(self):
        
        dihedralsVariables = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial, angle.atomD.serial] for angle in self.molecule_benzene_fromRDkit.torsionsVariable]
        golden_mol1_dihedralsVariables = [[7, 6, 4, 5], [8, 7, 6, 4], [9, 8, 7, 6], [10, 7, 6, 4], [11, 8, 7, 6], [12, 9, 8, 7], [13, 5, 4, 6], [14, 4, 5, 9], [15, 6, 4, 5]]

        dihedraslAdditionals = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial, angle.atomD.serial] for angle in self.molecule_benzene_fromRDkit.torsionsAdditional]
        golden_mol1_dihedralsAdditionals = [[6, 4, 5, 9], [7, 6, 4, 14], [8, 7, 6, 15], [9, 8, 7, 10], [10, 7, 6, 15], [10, 7, 8, 11], [11, 8, 9, 5], [11, 8, 9, 12], [12, 9, 5, 4], [12, 9, 5, 13], [13, 5, 4, 14], [13, 5, 9, 8], [14, 4, 6, 15], [5, 9, 8, 7], [8, 9, 5, 4], [14, 6, 4, 5], [13, 9, 5, 4], [15, 7, 6, 4], [10, 8, 7, 6], [11, 9, 8, 7], [12, 5, 9, 8]]

        self.assertEqual(dihedralsVariables, golden_mol1_dihedralsVariables)
        self.assertEqual(dihedraslAdditionals, golden_mol1_dihedralsAdditionals)

    def test_atoms(self):

        atoms_benzene = [str(atom) for atom in self.molecule_benzene_fromRDkit.atoms]


        with open(os.path.join(self.goldenDataPath,'benzene_atoms.txt')) as f: 
            golden_atoms_benzene = f.read().splitlines()

        self.assertEqual(atoms_benzene, golden_atoms_benzene)


    def test_molecule_writeZmatInitial(self):

        zmat.write(self.molecule_benzene_fromRDkit, self.molname, self.workdir, False)

        moleculeFile = os.path.join(self.workdir, self.molname+'.z')

        with open(os.path.join(self.goldenDataPath, 'benzene_before_BOSS.z')) as f: 
            golden_benzene = f.read().splitlines()

        with open(moleculeFile) as f: benzene = f.read().splitlines()

        # for line in benzene: print(line)

        self.assertEqual(golden_benzene, benzene)

        os.remove(moleculeFile)


    # def test_runBOSS(self):

    #     bnz_charge = 0.0
    #     bnz_opt = 0
    #     bnz_cgen = 'CM1A-LBCC'

        
    #     zmatNameGolden = os.path.join(self.goldenDataPath,'benzene_before_BOSS.z')
    #     runNameZmat = os.path.join(self.workdir, 'benzene_after_BOSS.z')

    #     shutil.copyfile(zmatNameGolden, runNameZmat)

    #     boss.run(runNameZmat, bnz_cgen, bnz_opt, bnz_charge, self.molname, self.workdir, False, False)

    #     with open(os.path.join(self.goldenDataPath, 'benzene_after_BOSS.z')) as f: golden_benzene = f.read().splitlines()
    #     with open(runNameZmat) as f: benzene = f.read().splitlines()

    #     self.assertEqual(golden_benzene, benzene)

    #     os.remove(runNameZmat)

    #     for fileName in ['out', 'MOL.pdb']:
    #         os.remove(os.path.join(self.workdir, fileName))


    #     bnz_opt = 3
    #     bnz_cgen = 'CM1A-LBCC'

        
    #     zmatNameGolden = os.path.join(self.goldenDataPath, 'benzene_before_BOSS.z')
    #     runNameZmat = os.path.join(self.workdir,'benzene_after_BOSS.z')
    #     shutil.copyfile(zmatNameGolden, runNameZmat)

    #     boss.run(runNameZmat, bnz_cgen, bnz_opt, bnz_charge, self.molname, self.workdir, False)

    #     with open(os.path.join(self.goldenDataPath, 'benzene_after_BOSS_opt3.z')) as f: golden_benzene = f.read().splitlines()
    #     with open(runNameZmat) as f: benzene = f.read().splitlines()

    #     self.assertEqual(golden_benzene, benzene)

    #     os.remove(runNameZmat)

    #     for fileName in ['out', 'MOL.pdb']:
    #         os.remove(os.path.join(self.workdir, fileName))

    def test_generateMoleculeFromZmat(self):

        bonds_benzeneVariables = [[bond.atomA.serial, bond.atomB.serial, bond.R0, bond.K0] for bond in self.molecule_benzene.bondsVariable]
    
        bonds_benzene_Variables_golden = [[5, 4, 1.4, 469.0], [6, 4, 1.4, 469.0], [7, 6, 1.4, 469.0], [8, 7, 1.4, 469.0], [9, 8, 1.4, 469.0], [10, 7, 1.08, 367.0], [11, 8, 1.08, 367.0], [12, 9, 1.08, 367.0], [13, 5, 1.08, 367.0], [14, 4, 1.08, 367.0], [15, 6, 1.08, 367.0]]

        self.assertEqual(bonds_benzeneVariables, bonds_benzene_Variables_golden)

        bonds_benzeneAdditionals = [[bond.atomA.serial, bond.atomB.serial, bond.R0, bond.K0] for bond in self.molecule_benzene.bondsAdditional]
    
        bonds_benzene_Additionals_golden = [[9, 5, 1.4, 469.0]]

        self.assertEqual(bonds_benzeneAdditionals, bonds_benzene_Additionals_golden)


        angles_benzeneVariables = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial, angle.angle0, angle.K0] for angle in self.molecule_benzene.anglesVariable]

        angles_benzene_Variables_golden = [[6, 4, 5, 120.0, 63.0], [7, 6, 4, 120.0, 63.0], [8, 7, 6, 120.0, 63.0], [9, 8, 7, 120.0, 63.0], [10, 7, 6, 120.0, 35.0], [11, 8, 7, 120.0, 35.0], [12, 9, 8, 120.0, 35.0], [13, 5, 4, 120.0, 35.0], [14, 4, 5, 120.0, 35.0], [15, 6, 4, 120.0, 35.0]]

        self.assertEqual(angles_benzeneVariables, angles_benzene_Variables_golden)

        angles_benzeneAdditionals = [[angle.atomA.serial, angle.atomB.serial, angle.atomC.serial, angle.angle0, angle.K0] for angle in self.molecule_benzene.anglesAdditional]

        angles_benzene_Additionals_golden = [[14, 4, 6, 120.0, 35.0], [15, 6, 7, 120.0, 35.0], [10, 7, 8, 120.0, 35.0], [11, 8, 9, 120.0, 35.0], [4, 5, 9, 120.0, 63.0], [13, 5, 9, 120.0, 35.0], [5, 9, 12, 120.0, 35.0], [8, 9, 5, 120.0, 63.0]]

        self.assertEqual(angles_benzeneAdditionals, angles_benzene_Additionals_golden)


        dihedrals_benzeneVariables = [[dihedral.atomA.serial, dihedral.atomB.serial, dihedral.atomC.serial, dihedral.atomD.serial,
            dihedral.V1, dihedral.V2, dihedral.V3, dihedral.V4] for dihedral in self.molecule_benzene.torsionsVariable]

        dihedrals_benzene_Variables_golden = [[7, 6, 4, 5, 0.0, 7.25, 0.0, 0.0], [8, 7, 6, 4, 0.0, 7.25, 0.0, 0.0], [9, 8, 7, 6, 0.0, 7.25, 0.0, 0.0], [10, 7, 6, 4, 0.0, 7.25, 0.0, 0.0], [11, 8, 7, 6, 0.0, 7.25, 0.0, 0.0], [12, 9, 8, 7, 0.0, 7.25, 0.0, 0.0], [13, 5, 4, 6, 0.0, 7.25, 0.0, 0.0], [14, 4, 5, 9, 0.0, 7.25, 0.0, 0.0], [15, 6, 4, 5, 0.0, 7.25, 0.0, 0.0]]

        self.assertEqual(dihedrals_benzeneVariables, dihedrals_benzene_Variables_golden)


        dihedrals_benzeneAdditionals = [[dihedral.atomA.serial, dihedral.atomB.serial, dihedral.atomC.serial, dihedral.atomD.serial,
            dihedral.V1, dihedral.V2, dihedral.V3, dihedral.V4] for dihedral in self.molecule_benzene.torsionsAdditional]

        dihedrals_benzene_Additionals_golden = [[9, 5, 4, 6, 0.0, 7.25, 0.0, 0.0], [14, 4, 6, 7, 0.0, 7.25, 0.0, 0.0], [15, 6, 7, 8, 0.0, 7.25, 0.0, 0.0], [10, 7, 8, 9, 0.0, 7.25, 0.0, 0.0], [15, 6, 7, 10, 0.0, 7.25, 0.0, 0.0], [11, 8, 7, 10, 0.0, 7.25, 0.0, 0.0], [5, 9, 8, 11, 0.0, 7.25, 0.0, 0.0], [12, 9, 8, 11, 0.0, 7.25, 0.0, 0.0], [4, 5, 9, 12, 0.0, 7.25, 0.0, 0.0], [13, 5, 9, 12, 0.0, 7.25, 0.0, 0.0], [14, 4, 5, 13, 0.0, 7.25, 0.0, 0.0], [8, 9, 5, 13, 0.0, 7.25, 0.0, 0.0], [15, 6, 4, 14, 0.0, 7.25, 0.0, 0.0], [7, 8, 9, 5, 0.0, 7.25, 0.0, 0.0], [4, 5, 9, 8, 0.0, 7.25, 0.0, 0.0], [5, 4, 6, 14, 0.0, 5.0, 0.0, 0.0], [4, 5, 9, 13, 0.0, 5.0, 0.0, 0.0], [4, 6, 7, 15, 0.0, 5.0, 0.0, 0.0], [6, 7, 8, 10, 0.0, 5.0, 0.0, 0.0], [7, 8, 9, 11, 0.0, 5.0, 0.0, 0.0], [8, 9, 5, 12, 0.0, 5.0, 0.0, 0.0]]

        self.assertEqual(dihedrals_benzeneAdditionals, dihedrals_benzene_Additionals_golden)


        atoms_benzene = [str(atom) for atom in self.molecule_benzene.atoms]

        with open(os.path.join(self.goldenDataPath, 'benzene_atoms_fromZMAT.txt')) as f: golden_atoms_benzene = f.read().splitlines()

        self.assertEqual(atoms_benzene, golden_atoms_benzene)


    def test_gromacs_outputs(self):

        gromacs.write(self.molecule_benzene, self.molname, self.workdir)

        groFileName, itpFileName = gromacs.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, itpFileName)) as f: benzene_itp = f.read().splitlines()
        with open(os.path.join(self.workdir, groFileName)) as f: benzene_gro = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.gmx.itp')) as f: benzene_itp_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.gmx.gro')) as f: benzene_gro_golden = f.read().splitlines()

        self.assertEqual(benzene_itp, benzene_itp_golden)
        self.assertEqual(benzene_gro, benzene_gro_golden)

        os.remove(itpFileName)
        os.remove(groFileName)

    def test_gromacs_FEP_singleTopoloogy_outputs(self):

        phenol_RDkit = copy.deepcopy(self.phenol_RDkit)
        molecule_phenol = copy.deepcopy(self.molecule_phenol)

        BtoAserialCorrespondency, umatchB = alchemify.alignMolecules(self.benzene_RDkit, phenol_RDkit, self.workdir, self.debug)

        moleculeAB_singleTopology = alchemify.generateMoleculeAB_singleTopology(self.molecule_benzene, molecule_phenol, BtoAserialCorrespondency, umatchB)

        gromacs.write(moleculeAB_singleTopology, self.molname, self.workdir)

        groFileName, itpFileName = gromacs.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, itpFileName)) as f: benzene_itp = f.read().splitlines()
        with open(os.path.join(self.workdir, groFileName)) as f: benzene_gro = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.gmx.itp')) as f: benzene_itp_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.gmx.gro')) as f: benzene_gro_golden = f.read().splitlines()

        self.assertEqual(benzene_itp, benzene_itp_golden)
        self.assertEqual(benzene_gro, benzene_gro_golden)

        os.remove(itpFileName)
        os.remove(groFileName)

        os.remove(os.path.join(self.workdir, 'molecule_A.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B_superposedToA.pdb'))



    def test_Q_outputs(self):

        Q.write(self.molecule_benzene, self.molname, self.workdir)

        prmFileName, libFileName, pdbFileName , _ = Q.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, libFileName)) as f: benzene_lib = f.read().splitlines()
        with open(os.path.join(self.workdir, prmFileName)) as f: benzene_prm = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene.q.lib')) as f: benzene_lib_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.q.prm')) as f: benzene_prm_golden = f.read().splitlines()


        self.assertEqual(benzene_lib, benzene_lib_golden)
        self.assertEqual(benzene_prm, benzene_prm_golden)

        os.remove(libFileName)
        os.remove(prmFileName)
        os.remove(pdbFileName)


    def test_Q_FEP_singleTopology_outputs(self):

        phenol_RDkit = copy.deepcopy(self.phenol_RDkit)
        molecule_phenol = copy.deepcopy(self.molecule_phenol)

        BtoAserialCorrespondency, umatchB = alchemify.alignMolecules(self.benzene_RDkit, phenol_RDkit, self.workdir, self.debug)

        moleculeAB_singleTopology = alchemify.generateMoleculeAB_singleTopology(self.molecule_benzene, molecule_phenol, BtoAserialCorrespondency, umatchB)
    
        Q.write(moleculeAB_singleTopology, self.molname, self.workdir)

        prmFileName, libFileName, pdbFileName , fepFileName = Q.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, libFileName)) as f: benzene_lib = f.read().splitlines()
        with open(os.path.join(self.workdir, prmFileName)) as f: benzene_prm = f.read().splitlines()
        with open(os.path.join(self.workdir, fepFileName)) as f: benzene_fep = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.q.lib')) as f: benzene_lib_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.q.prm')) as f: benzene_prm_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.q.fep')) as f: benzene_fep_golden = f.read().splitlines()


        self.assertEqual(benzene_lib, benzene_lib_golden)
        self.assertEqual(benzene_prm, benzene_prm_golden)
        self.assertEqual(benzene_fep, benzene_fep_golden)

        os.remove(libFileName)
        os.remove(prmFileName)
        os.remove(pdbFileName)
        os.remove(fepFileName)

        os.remove(os.path.join(self.workdir, 'molecule_A.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B_superposedToA.pdb'))


    def test_CHARMM_FEP_dualTopology_outputs(self):

        molecule_phenol = copy.deepcopy(self.molecule_phenol)
        moleculeAB_dualTopology = alchemify.generateMoleculeAB_dualTopology(self.molecule_benzene, molecule_phenol)

        charmm.write(moleculeAB_dualTopology, self.molname, self.workdir)

        prmFileName, rtfFileName, pdbFileName , fepFileName = charmm.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, rtfFileName)) as f: benzene_rtf = f.read().splitlines()
        with open(os.path.join(self.workdir, prmFileName)) as f: benzene_prm = f.read().splitlines()
        with open(os.path.join(self.workdir, fepFileName)) as f: benzene_fep = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.charmm.rtf')) as f: benzene_rtf_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.charmm.prm')) as f: benzene_prm_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene_phenol.charmm.fep')) as f: benzene_fep_golden = f.read().splitlines()


        self.assertEqual(benzene_rtf, benzene_rtf_golden)
        self.assertEqual(benzene_prm, benzene_prm_golden)
        self.assertEqual(benzene_fep, benzene_fep_golden)

        os.remove(rtfFileName)
        os.remove(prmFileName)
        os.remove(pdbFileName)
        os.remove(fepFileName)

    def test_CHARMM_outputs(self):
        
        charmm.write(self.molecule_benzene, self.molname, self.workdir)

        prmFileName, rtfFileName, pdbFileName , _ = charmm.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, rtfFileName)) as f: benzene_rtf = f.read().splitlines()
        with open(os.path.join(self.workdir, prmFileName)) as f: benzene_prm = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene.charmm.rtf')) as f: benzene_rtf_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.charmm.prm')) as f: benzene_prm_golden = f.read().splitlines()


        self.assertEqual(benzene_rtf, benzene_rtf_golden)
        self.assertEqual(benzene_prm, benzene_prm_golden)

        os.remove(rtfFileName)
        os.remove(prmFileName)
        os.remove(pdbFileName)


    def test_alignMolecules(self):

        phenol_RDkit = copy.deepcopy(self.phenol_RDkit)

        BtoAserialCorrespondency, umatchB = alchemify.alignMolecules(self.benzene_RDkit, phenol_RDkit, self.workdir, self.debug)
        
        golden_BtoAserialCorrespondency = {0: 0, 1: 1, 6: 5, 5: 4, 3: 3, 2: 2, 10: 11, 11: 6, 7: 7, 8: 8, 9: 9, 4: 10}
        golden_umatch = [12]

        self.assertEqual(BtoAserialCorrespondency, golden_BtoAserialCorrespondency)
        self.assertEqual(umatchB, golden_umatch)

        os.remove(os.path.join(self.workdir, 'molecule_A.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B_superposedToA.pdb'))


    def test_ZMAT_FEP_singleTopology(self):

        phenol_RDkit = copy.deepcopy(self.phenol_RDkit)
        molecule_phenol = copy.deepcopy(self.molecule_phenol)

        BtoAserialCorrespondency, umatchB = alchemify.alignMolecules(self.benzene_RDkit, phenol_RDkit, self.workdir, self.debug)

        moleculeAB_singleTopology = alchemify.generateMoleculeAB_singleTopology(self.molecule_benzene, molecule_phenol, BtoAserialCorrespondency, umatchB)

        zmat.write(moleculeAB_singleTopology, self.molname, self.workdir, writeAtomParameters = True)

        zmatFileName = zmat.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, zmatFileName)) as f: benzene_phenol_single = f.read().splitlines()  

        with open(os.path.join(self.goldenDataPath, 'benzene_phenol_singleTopology.z')) as f: benzene_phenol_single_golden = f.read().splitlines()

        self.assertEqual(benzene_phenol_single, benzene_phenol_single_golden)

        os.remove(os.path.join(self.workdir, zmatFileName))
        os.remove(os.path.join(self.workdir, 'molecule_A.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B.pdb'))
        os.remove(os.path.join(self.workdir, 'molecule_B_superposedToA.pdb'))

    def test_ZMAT_FEP_dualTopology(self):

        molecule_phenol = copy.deepcopy(self.molecule_phenol)

        moleculeAB_dualTopology = alchemify.generateMoleculeAB_dualTopology(self.molecule_benzene, molecule_phenol)

        zmat.write(moleculeAB_dualTopology, self.molname, self.workdir, writeAtomParameters = True)

        zmatFileName = zmat.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, zmatFileName)) as f: benzene_phenol_dual = f.read().splitlines()  

        with open(os.path.join(self.goldenDataPath, 'benzene_phenol_dualTopology.z')) as f: benzene_phenol_dual_golden = f.read().splitlines()

        self.assertEqual(benzene_phenol_dual, benzene_phenol_dual_golden)

        os.remove(os.path.join(self.workdir, zmatFileName))

    def test_ZMAT_outputs(self):

        zmat.write(self.molecule_benzene, self.molname, self.workdir, writeAtomParameters = True)
        
        zmatFileName = zmat.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, zmatFileName)) as f: benzene_zmat = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.boss.z')) as f: benzene_zmat_golden = f.read().splitlines()

        self.assertEqual(benzene_zmat, benzene_zmat_golden)

        os.remove(zmatFileName)

    def test_Tinker_outputs(self):

        import ligpargen.inout.tinker as tinker
        
        tinker.write(self.molecule_benzene, self.molname, self.workdir)

        xyzFile, keyFile = tinker.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, xyzFile)) as f: benzene_xyz = f.read().splitlines()
        with open(os.path.join(self.workdir, keyFile)) as f: benzene_key = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene.tinker.xyz')) as f: benzene_xyz_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.tinker.key')) as f: benzene_key_golden = f.read().splitlines()


        self.assertEqual(benzene_xyz, benzene_xyz_golden)
        self.assertEqual(benzene_key, benzene_key_golden)

        os.remove(xyzFile)
        os.remove(keyFile)

    def test_Lammps_outputs(self):
    
        import ligpargen.inout.lammps as lammps
        
        lammps.write(self.molecule_benzene, self.molname, self.workdir)

        lmpFile = lammps.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, lmpFile)) as f: benzene_lmp = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.lammps.lmp')) as f: benzene_lmp_golden = f.read().splitlines()

        self.assertEqual(benzene_lmp, benzene_lmp_golden)

        os.remove(lmpFile)

    def test_pqr_outputs(self):
        
        import ligpargen.inout.pqr as pqr
        
        pqr.write(self.molecule_benzene, self.molname, self.workdir)

        pqrFile = pqr.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, pqrFile)) as f: benzene_pqr = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.pqr')) as f: benzene_pqr_golden = f.read().splitlines()

        self.assertEqual(benzene_pqr, benzene_pqr_golden)

        os.remove(pqrFile)


    def test_Desmond_outputs(self):
    
        import ligpargen.inout.desmond as desmond
        
        desmond.write(self.molecule_benzene, self.molname, self.workdir)

        cmsFile = desmond.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, cmsFile)) as f: benzene_cms = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.desmond.cms')) as f: benzene_cms_golden = f.read().splitlines()

        self.assertEqual(benzene_cms, benzene_cms_golden)

        os.remove(cmsFile)

    def test_openMM_outputs(self):

        import ligpargen.inout.openmm as openmm
        
        openmm.write(self.molecule_benzene, self.molname, self.workdir)

        pdbFileName, xmlFileName = openmm.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, pdbFileName)) as f: benzene_pdb = f.read().splitlines()
        with open(os.path.join(self.workdir, xmlFileName)) as f: benzene_xml = f.read().splitlines()


        with open(os.path.join(self.goldenDataPath, 'benzene.openmm.pdb')) as f: benzene_pdb_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.openmm.xml')) as f: benzene_xml_golden = f.read().splitlines()


        self.assertEqual(benzene_pdb, benzene_pdb_golden)
        self.assertEqual(benzene_xml, benzene_xml_golden)

        os.remove(pdbFileName)
        os.remove(xmlFileName)

    def test_xplor_outputs(self):
    
        import ligpargen.inout.xplor as xplor
        
        xplor.write(self.molecule_benzene, self.molname, self.workdir)

        topFileName, paramFileName = xplor.getFileNames(self.molname, self.workdir)

        with open(os.path.join(self.workdir, topFileName)) as f: benzene_top = f.read().splitlines()
        with open(os.path.join(self.workdir, paramFileName)) as f: benzene_param = f.read().splitlines()

        with open(os.path.join(self.goldenDataPath, 'benzene.xplor.top')) as f: benzene_top_golden = f.read().splitlines()
        with open(os.path.join(self.goldenDataPath, 'benzene.xplor.param')) as f: benzene_param_golden = f.read().splitlines()


        self.assertEqual(benzene_top, benzene_top_golden)
        self.assertEqual(benzene_param, benzene_param_golden)

        os.remove(topFileName)
        os.remove(paramFileName)


def getAtomProperties(molecule1):

    conf = molecule1.GetConformer()

    atomNames = [atom.GetSymbol() for atom in molecule1.GetAtoms()]
    atomAtomicnumber = [atom.GetAtomicNum() for atom in molecule1.GetAtoms()]
    atomPositions = [[conf.GetAtomPosition(atom.GetIdx()).x, conf.GetAtomPosition(atom.GetIdx()).y, \
        conf.GetAtomPosition(atom.GetIdx()).z] for atom in molecule1.GetAtoms()]

    return atomNames, atomAtomicnumber, atomPositions


if __name__ == '__main__'  and __package__ is None:

    unittest.main()