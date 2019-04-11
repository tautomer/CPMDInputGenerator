"""Generate CPMD input files."""
import sys
import json
import re
import numpy as np


class GenerateCPMDInput:

    def __init__(self):
        self.inp = ''
        # general
        self.pimd = False
        # not implementing geometry optimization by now
        # convergence and minimizer hard-coded, too
        self.opt = False
        self.md = False
        self.cp = False
        self.bo = False
        self.rst = False
        self.rstopt = ''
        self.symm = int(-1)
        self.solver = ''
        self.cellparmtype = ''
        self.cellparm = ''
        self.cutoff = 0
        self.functional = ''
        self.psptype = ''
        # implementing in a general way, though only D works now
        self.isotope = False
        self.isoidx = {}
        self.cnstr = {}
        self.restr = {}
        self.ncnstr = int(0)
        self.nrestr = int(0)
        # pimd
        self.nb = int(0)
        self.nm = False
        self.wmass = 0
        self.facmass = 0
        self.debrogile = False
        self.inittemp = 0
        # since generate is not necessary, only read is implemented
        self.init = False
        self.repfname = ''
        self.prtlvl = int(1)
        self.procgrp = int(1)
        # md
        self.keepwfn = False
        self.traj = False
        self.sample = False
        self.step = int(0)
        self.temp = 0
        self.maxstep = int(0)
        self.dt = 0
        self.nose = False
        self.ion = {}
        self.elec = {}
        self.center = False
        self.quench = False
        # flux side
        self.fs = False
        self.fsidx = ''
        self.ranvel = False


    def ReadConfig(self, conf):
        """read 
        
        Arguments:
            conf {str} -- configuration json
        """
        print("This script is not able to handle typos, etc.\nDouble-check you"
              "r json file and resulting input file.")
        with open(conf) as conf_file:    
            confdata = json.load(conf_file)
        config = {}
        # convert all strings to caps, since CPMD input file only takes caps
        for key, value in confdata.items():
            if isinstance(key, str):
                key = key.upper()
            # TODO: better way to handle this?
            if isinstance(value, str) and key != 'INIT' and 'inp' not in value:
                value = value.upper()
            elif isinstance(value, dict):
                tmp = dict(value)
                value = {}
                for k, v in tmp.items():
                    if isinstance(k, str):
                        k = k.upper()
                    if isinstance(v, str):
                        v = v.upper()
                    value.update({k:v})
            config.update({key:value})
        config = self.ParseConfig(config)
        self.CheckParm()
        if len(config):
            print("Some keys in the json file does not match the records:")
            for key, value in config.items():
                print("key: {}, value: {}".format(key, value))
            print("Double-check these keys!")


    def ParseConfig(self, config):
        """parse json file
        
        Arguments:
            config {dictionary} -- loaded from json
        
        Returns:
            dictionary -- reduced config for future uses
        """
        # parse the dictionary
        if 'INP' in config:
            self.inp = config['INP']
            config.pop('INP')
        # pimd or classical md
        if 'PI' in config:
            self.pimd = True
            config.pop('PI')
            # pimd section
            if 'NB' in config:
                self.nb = config['NB']
                config.pop('NB')
            if 'NM' in config:
                self.nm = True
                self.wmass = config['NM']
                config.pop('NM')
            if 'FACMASS' in config:
                self.facmass = config['FACMASS']
                config.pop('FACMASS')
            if 'DEBROGILE' in config:
                self.debrogile = True
                self.inittemp = config['DEBROGILE']
                config.pop('DEBROGILE')
            if 'INIT' in config:
                self.init = True
                self.repfname = config['INIT']
                config.pop('INIT')
            if 'PRINT' in config:
                self.prtlvl = config['PRINT']
                config.pop('PRINT')
            if 'PROCGRP' in config:
                self.procgrp = config['PROCGRP']
                config.pop('PROCGRP')
        else:
            self.nb = 1
        # job type
        if 'OPT' in config:
            self.opt = True
            config.pop('OPT')
        elif 'MD' in config:
            self.md = True
            if config['MD'] == 'BO':
                self.bo = True
            elif config['MD'] == 'CP' or not config['MD']:
                self.cp = True
            else:
                raise ValueError("MD should have an value of BO, "
                                 "CP or omitted.")
            config.pop('MD')
            # general md
            if 'KEEPWFN' in config:
                self.keepwfn = True
                config.pop('KEEPWFN')
            if 'QUENCH' in config:
                self.quench = True
                config.pop('QUENCH')
            if 'TRAJ' in config:
                self.traj = True
                if 'SAMPLE' in config['TRAJ']:
                    self.sample = True
                    self.step = config['TRAJ']['SAMPLE']
                config.pop('TRAJ')
            if 'TEMP' in config:
                self.temp = config['TEMP']
                config.pop('TEMP')
            if 'STEP' in config:
                self.maxstep = config['STEP']
                config.pop('STEP')
            if 'DT' in config:
                self.dt = config['DT']
                config.pop('DT')
            if 'NOSE' in config:
                self.nose = True
                if 'ION' in config['NOSE']:
                    # use massive by default
                    self.ion = config['NOSE']['ION']
                if 'ELEC' in config['NOSE']:
                    # use massive by default
                    self.elec = config['NOSE']['ELEC']
                config.pop('NOSE')
            if 'CENTER' in config:
                self.center = True
                config.pop('CENTER')
            # flux side
            if 'FS' in config:
                self.fs = True
                self.fsidx = config['FS']
                config.pop('FS')
            if 'RANVEL' in config:
                self.ranvel = True
                config.pop('RANVEL')
        # restart or not
        if 'RESTART' in config:
            self.rst = True
            self.rstopt = config['RESTART']
            config.pop('RESTART')
        # cell type
        if 'SYMM' in config:
            self.symm = config['SYMM']
            config.pop('SYMM')
        # poisson solver
        if 'SOLVER' in config:
            self.solver = config['SOLVER']
            config.pop('SOLVER')
        # cell parameters
        if 'CELL' in config:
            self.cellparmtype = config['CELL']['TYPE']
            self.cellparm = config['CELL']['SIZE']
            config.pop('CELL')
        # cutoff
        if 'CUTOFF' in config:
            self.cutoff = config['CUTOFF']
            config.pop('CUTOFF')
        # pesudopotential and functional
        if 'PSP' in config:
            self.psptype = config['PSP']['TYPE']
            self.functional = config['PSP']['FUNC']
            config.pop('PSP')
        # constraints and restraints
        if 'ISOTOPE' in config:
            self.isotope = True
            self.isoidx = config['ISOTOPE']
            config.pop('ISOTOPE')
        if 'CNSTR' in config:
            for key, value in config['CNSTR'].items():
                if 'CONS' in key:
                    self.ncnstr += 1
                    self.cnstr[self.ncnstr] = value
                elif 'REST' in key:
                    self.nrestr += 1
                    self.restr[self.nrestr] = value
                else:
                    raise ValueError('Error in CONSTRAINTS section')
            config.pop('CNSTR')
        return config


    def CheckParm(self):
        """check conflicts in the configurations
        
        """
        # TODO: a better way to handle errors and warnings.
        # job type
        if not self.opt and not self.md:
            raise ValueError('Either OPT or MD has to be specified')
        # general md
        if self.symm < 0:
            self.symm = 0
            raise Warning('Symmetry not specified. Reset to isolated box.')
        if not self.cellparm:
            raise ValueError('Cell parameters must be specified.')
        if not self.solver:
            self.solver = 'TUCKERMAN'
            raise Warning('Poisson solver not specified, use Tuckerman')
        if self.md:
            if self.temp == 0:
                raise ValueError('Temperature must be specified for MD')
            if self.maxstep <= 0:
                raise ValueError('Number of steps missing or incorrect.')
            if self.dt == 0:
                raise Warning('dt not specified. CPMD default dt is very '
                              'dangerous. Better re-write the json file.')
            elif self.dt < 0:
                raise ValueError('dt value incorrect.')
            if self.nose and self.fs:
                raise Warning('Using thermostat with fs is meaningless.')
            if self.ranvel and self.fs:
                raise Warning('ranvel is redundant when fs is enabled')
        # bomd and nose
            if self.bo and self.elec:
                raise Warning('Running BOMD, electron thermostat ignored.')
            if self.quench:
                if self.sample and self.step != 0:
                    self.step = 0
                    raise Warning('No tracjectory saved with quench enabled.')
                else:
                    self.traj = True
                    self.sample = True
                    self.step = 0
        # pimd
        if self.pimd:
            if self.nb == 0:
                raise ValueError('nb must be specified.')
            if self.md and not self.rst:
                raise ValueError('PI MD has to restart.')
            if self.debrogile and self.init:
                raise ValueError('DEBROGILE cannot work with READ REPLICAS.')


    def ReadInput(self, fname):
        """read input file and sort atom types and indices
        
        Arguments:
            fname {list} -- geo and traj files
        
        Returns:
            [list, np.array] -- atomic labels and coordinates
        """
        fac = 0.52917721092
        atomlabel = np.genfromtxt(fname[0], dtype="a2", skip_header=2)
        if len(fname) == 2:
            mol = np.genfromtxt(fname[1], usecols=(1, 2, 3))
        else:
            mol = np.genfromtxt(fname[0], usecols=(1, 2, 3))
        mol1 = [i[0].decode('ascii') for i in atomlabel]
        if self.isotope:
            for k, v in self.isoidx.items():
                mol1[int(k)-1] = v
        atomtype = sorted(list(set(mol1)))
        indices = []
        length = []
        for i in atomtype:
            tmp = [j for j, x in enumerate(mol1) if x == i]
            indices.append(tmp)
            length.append(len(tmp))
        atomtype = list(zip(atomtype, indices, length))
        natom = len(atomlabel)
        if self.nb > 1:
            mol2 = np.zeros((natom, 3))
            if self.init:
                crd = open(self.repfname, 'w')
            for i in range(self.nb):
                ind = natom * i
                mol2 += mol[ind:ind+natom]
                if self.init:
                    print("{}\n{}".format(' ', re.sub('[\[\]]', ' ',
                          np.array_str(mol[ind:ind+natom]))), file=crd)
            mol2 /= self.nb
        else:
            mol2 = mol
        mol2 *= fac
        return atomtype, mol2


    def WriteInput(self, atomtype, xyz):
        """write CPMD input file
        
        Arguments:
            atomtype {list} -- labels
            xyz {np.array} -- geometry
        """
        # elements dictionary
        # guess only Deuterium will be used in the near future!
        #      label type   amu   lmax
        elem = {"H": ["H", 1.008, "S"],
                "D": ["H", 2.014, "S"],
                "C": ["C", 12.011, "P"],
                "N": ["N", 14.007, "P"],
                "O": ["O", 15.999, "P"]}
        out = open(self.inp, 'w')
        # &CPMD section
        print("&CPMD", file=out)
        if self.pimd:
            print("  PATH INTEGRAL", file=out)
        if self.opt:
            if self.rst:
                print("  RESTART {}".format(self.rstopt), file=out)
            print("  OPTIMIZE WAVEFUNCTION\n  CONVERGENCE ORBITALS\n    1d-7\n"
                  "  PCG MINIMIZE\n  TIMESTEP\n    20\n   MEMORY BIG\n  CENTER"
                  " MOLECULE ON", file=out)
        else:
            if self.cp:
                print("  MOLECULAR DYNAMICS CP", file=out)
            else:
                print("  MOLECULAR DYNAMICS BO", file=out)
            if self.rst:
                print("  RESTART {}".format(self.rstopt), file=out)
            # TODO: f.write may work better in this case
            if self.keepwfn:
                print("  REAL SPACE WFN KEEP", file=out)
            if self.quench:
                print("  QUENCH IONS ELECTRONS", file=out)
            if self.traj:
                if self.sample:
                    print("  TRAJECTORY SAMPLE\n    {}".format(self.step),
                          file=out)
                else:
                    print("  TRAJECTORY", file=out)
            print("  TEMPERATURE\n    {}\n  MAXSTEP\n    {}".format(self.temp,
                  self.maxstep), file=out)
            if self.dt != 0:
                print("  TIMESTEP\n    {}".format(self.dt), file=out)
            if self.nose:
                if self.ion:
                    print("  NOSE IONS MASSIVE\n    {}".format(self.ion),
                          file=out)
                if self.elec and self.cp:
                    print("  NOSE ELECTRONS\n    {}".format(self.elec),
                          file=out)
            if self.center:
                print("  CENTER MOLECULE ON", file=out)
            if self.fs:
                print("  FLUX SIDE\n    {}".format(self.fsidx), file=out)
            if self.ranvel:
                print("  RAN INIT VEL", file=out)
        print("&END\n\n", file=out)
        # &SYSTEM section
        print("&SYSTEM\n  SYMMETRY\n    {}\n  POISSON SOLVER {}\n   ANGSTROM"
              "\n  CELL {}\n    {}\n  CUTOFF\n    {}\n&END\n\n\n&ATOMS".format(
              self.symm, self.solver, self.cellparmtype, self.cellparm,
              self.cutoff), file=out)
        # &ATOMS section
        # TODO: ability to handle other isotopes (not really necessary though)
        mass = ''
        for i in range(len(atomtype)):
            label = atomtype[i][0]
            if self.isotope:
                mass += '\n' + '    ' + str(elem[label][1]) 
            pspfname = '*' + elem[label][0] + '_' + self.psptype + '_' + \
                       self.functional
            print("{}\nLMAX={}\n{}".format(pspfname, elem[label][2],
                  atomtype[i][2]), file=out)
            for j in atomtype[i][1]:
                print("    {0:11.6f}{1:11.6f}{2:11.6f}".format(xyz[j][0],
                      xyz[j][1], xyz[j][2]), file=out)
        if self.isotope:
            print("  ISOTOPES{}".format(mass), file=out)
        if self.ncnstr or self.nrestr:
            print("  CONSTRAINTS", file=out)
            if self.ncnstr:
                print("    FIX STRUCTURE\n      {}".format(self.ncnstr),
                      file=out)
                for i in range(self.ncnstr):
                    print("    {}".format(self.cnstr[i+1]), file=out)
            if self.nrestr:
                print("    RESTRAINTS\n      {}".format(self.nrestr), file=out)
                for i in range(self.nrestr):
                    print("    {}".format(self.restr[i+1]), file=out)
            print("  END CONSTRAINTS", file=out)
        # &DFT section, pretty much hardcoded
        print("&END\n\n\n&DFT\n  FUNCTIONAL {}\n  GC-CUTOFF\n    1d-8\n&END".
             format(self.functional), file=out)
        if self.pimd:
            print("\n\n&PIMD\n  TROTTER DIMENSION\n    {}".format(self.nb),
                  file=out)
            if self.nm: 
                print("  NORMAL MODES\n    {}".format(self.wmass), file=out)
            print("  FACMASS\n    {}".format(self.facmass), file=out)
            if self.debrogile:
                print("  DEBROGLIE CENTROID\n    {}".format(self.inittemp),
                      file=out)
            elif self.init:
                print("  INITIALIZATION\n  READ REPLICAS\n    {}".format(
                      self.repfname), file=out)
            print("  PRINT LEVEL\n    {}\n  PROCESSOR GROUPS\n    {}\n&END".
                  format(self.prtlvl, self.procgrp), file=out)
        out.close()
        print("Successfully Created CPMD Input File {}.\nPlease check the "
              "resulting file carefully before running MD".format(self.inp))


gi = GenerateCPMDInput()
conf = sys.argv[1]
fname = sys.argv[2:4]
gi.ReadConfig(conf)
atomtype, xyz = gi.ReadInput(fname)
gi.WriteInput(atomtype, xyz)
