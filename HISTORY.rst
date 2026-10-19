=======
History
=======
2026.10.19 -- Bugfix: a convergence preset without the dialog
    * The dialog sets the convergence formula and criteria when a preset, such as
      'Gaussian (tight)', is chosen, but a flowchart made without it -- by
      seamm-flowchart or the MCP server -- had only the preset's name: the
      optimization stopped with KeyError '', or with the formula set, ran with the
      default criteria. They now come from the preset when the step runs.

2026.3.1 -- Internal: switching from deprecated library pkg_resources to importlib

2025.5.14 -- Added search for stereoisomers
    * Added stereoisomers as a target, This will take a structure and generate
      structures for stereoisomers, including that caused by double bonds.
      
2025.5.8 -- Bugfixes: using 'natoms' in steps, and continuing when not converged.
    * Corrected errors handling the use of 'natoms' in the maximum number of steps
      allowed, e.g. '6 * natoms'.
    * Correctly handle convergence failure in  geomeTRIC
      
2024.10.20 -- Bugfixes & internal changes
    * Reorganized the code, moving the ASE specific parts to seamm-ase.
    * Fixed formatting of output
    * Fixed typos in the references

2024.10.15 -- Enhancement: improved output and aded to results
    * Added energy, etc. to results
    * Added more detail to printed output
    * Fixed bugs with the GUI and with labeling substeps.

2024.10.13 -- Added geomeTRIC minimizer
    * Added transition states as a target
    * Added the geomeTRIC minimizer
    * Rationalized the convergence criteria, allowing them to be used across minimizers.
      
2024.7.31 -- Plug-in created using the SEAMM plug-in cookiecutter.

