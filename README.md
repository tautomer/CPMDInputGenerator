# CPMDInputGenerator

I wrote this script to generate CPMD input files by reading configurations from
a json file for my daily usage. Geometry is read from a xyz file and CPMD 
TRAJECTORY file.

Its functionality is very limited, since lots of keywords are hardcoded, 
especially for optimization and DFT.

Since CPMD only takes capitals (correct me if it is wrong), almost everything 
in json will be converted to caps after being loaded, excpet stuff like 
filenames. So it does not matter whether you use upper or lower cases in the 
file, but input file name and replica coordinate file name must be right.

List of keywords:

<pre>
* pi        -- PIMD or classical MD
* nb        -- number of beads for PIMD
* nm        -- use normal modes
* facmss    -- factor to manipulate masses in nm
* debrogile -- generate bead postions by initial guess
* init      -- hardcoded for 'INITILIZATION' & 'READ REPLICAS'
* print     -- print level (trivial)
* procgrp   -- processor groups, often limited to 1 due to CPMD bugs
* opt       -- hardcoded for 'OPTIMIZATION WAVEFUNCTION'
* MD        -- molecular dynamics, both BO and CP are supported
* keepwfn   -- store wavefunction in real spcae in memory
* traj      -- write TRAJECTORY (support 'SAMPLE')
* temp      -- temperature
* step      -- total steps
* dt        -- timestep
* nose      -- use Nose-Hoover thermostat, only massive ion and electron are supported
* center    -- center molecule
* fs        -- flux-side (only work for my modified CPMD)
* ranvel    -- randomize initial velocities from Boltzmann distribution (only work for my modified CPMD)
* restart   -- use RESTART files
* symm      -- simulation box type
* solver    -- Poisson solver to use
* cell      -- how cell parameters are specified and what the parameters are
* cutoff    -- cutoff for the plane wave basis
* psp       -- pseduopotential type and functional type (also used in DFT part)
* isotope   -- specify isotopes with atomic masses
* cnstr     -- specify any number of constraints and restraints
</pre>
