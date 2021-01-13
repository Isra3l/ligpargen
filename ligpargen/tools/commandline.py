"""

Module with command line parsing related functions (compute angles, distances, ...)
Author: Israel Cabeza de Vaca Lopez
Email:  israel.cabezadevaca@icm.uu.se

"""


from argparse import ArgumentParser, RawDescriptionHelpFormatter
from os import getcwd
from sys import argv

def commandline_options():
   """
   Parse argumments from command line

   Returns
   -------
   argsparse
       Arguments
   """

   parser = ArgumentParser(
         prog='LigParGen2.1',
         formatter_class=RawDescriptionHelpFormatter,
         description='''
      OPLSAA with CM1A/CM1A-LBCC charge models parameter Generator
      This version includes Alchemical templates generation 

      Author: Israel Cabeza de Vaca Lopez
      Email: israel.cabezadevaca@yale.edu

      Output Files Generated : OpenMM, CHARMM/NAMD, GROMACS, CNS/X-PLOR, Q, DESMOND, BOSS/MCPRO, PDB2PQR

      Input Files Supported: SMILES, PDB, MOL and MOL2 files + any input format supported by Open Babel

      Examples:

      Template generation:

            ligpargen -i phenol.pdb 
            ligpargen -i phenol.pdb    -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A
            ligpargen -s 'c1ccc(cc1)O' -n phenol -p phenol -r MOL -c 0 -o 0 -cgen CM1A-LBCC

      Alchemical Template generation:
      
            ligpargen -s 'c1ccc(cc1)O' -sb 'c1ccccc1' -n phenol_benzene -p phenol_benzene -r MOL -c 0 -o 0 -cgen CM1A-LBCC 
      
      Requirements:
      
      BOSS software  (It can be downloaded for free from: http://zarbi.chem.yale.edu/software.html)


      Please cite the following references:

      1. LigParGen web server: an automatic OPLS-AA parameter generator for organic ligands  
         Leela S. Dodda  Israel Cabeza de Vaca  Julian Tirado-Rives William L. Jorgensen 
         Nucleic Acids Research, Volume 45, Issue W1, 3 July 2017, Pages W331–W336
      2. 1.14*CM1A-LBCC: Localized Bond-Charge Corrected CM1A Charges for Condensed-Phase Simulations
   '''
   )
   parser.add_argument("-p",       "--path",     help="Path to the working directory",                         type=str, default=getcwd())
   parser.add_argument("-r",       "--resname",  help="A three-letter abbriviation of molecule name",          type=str, default=None)
   parser.add_argument("-n",       "--molname",  help="Molecule name",                                         type=str, default=None)
   parser.add_argument("-s",       "--smile",    help="SMILES code from CHEMSPIDER or PubChem",                type=str, default=None)
   parser.add_argument("-i",       "--ifile",    help="Molecule file from CHEMSPIDER or PubChem",              type=str, default=None)
   parser.add_argument("-cgen",    "--cgen",     help="Charge generating algorithm: choices[CM1A, CM1A-LBCC]", type=str, default=None, choices=['CM1A', 'CM1A-LBCC'])
   parser.add_argument("-c",       "--charge",   help="0: Neutral <0: Anion >0: Cation ",                      type=int, choices=[charge for charge in range(-10,11)], default=None) 
   parser.add_argument("-o",       "--opt",      help="Optimization or Single Point Calculation",              type=int, choices=[0, 1, 2, 3], default=0)
   parser.add_argument("-debug",   "--debug",    help="Turn on debug mode",                                    action="store_true")
   parser.add_argument("-verbose", "--verbose",  help="Print mor info",                                        action="store_true")

   # Alchemical transformations

   parser.add_argument("-sb",       "--smileB",    help="SMILES code from CHEMSPIDER or PubChem",                type=str, default=None)
   parser.add_argument("-ib",       "--ifileB",    help="Molecule file from CHEMSPIDER or PubChem",              type=str, default=None)

   # Alchemical optional arguments

   parser.add_argument("-cgenb",    "--cgenB",     help="Charge generating algorithm: choices[CM1A, CM1A-LBCC]", type=str, default=None, choices=['CM1A', 'CM1A-LBCC'])
   parser.add_argument("-cb",       "--chargeB",   help="0: Neutral <0: Anion >0: Cation ",                      type=int, choices=[charge for charge in range(-10,11)], default=None) 
   parser.add_argument("-ob",       "--optB",      help="Optimization or Single Point Calculation",              type=int, choices=[0, 1, 2, 3], default=None)

   # Checker 

   parser.add_argument("-check",    "--checker",  help="Flag to compare GAS phase energies between softwares to verify  consistency in the outputs", action="store_true", default=False)


   if len(argv)==1:
      parser.print_help()
      exit(1)

   args = parser.parse_args()

   return args


def checkArguments(args):
   """
   Check argumments 

   Parameters
   ----------
   args : argparse
       Arguments by command line

   Returns
   -------
   bool
       True if there is alchemical transformations
   """
       
   if args.smile==None and args.ifile==None:

      print("\n       ERROR!    Missing input molecule.\n")
      exit()

   alchemicalTransformation = False

   if args.smileB!=None or args.ifileB!=None: alchemicalTransformation = True

   if alchemicalTransformation:
    
      if args.chargeB==None: args.chargeB = args.charge
      if args.cgenB==None: args.cgenB = args.cgen
      if args.optB==None: args.optB = args.opt

   return alchemicalTransformation


