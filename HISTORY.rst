=======
History
=======
2026.11.25.1 -- Every optimizer works from the Model Chemistry
    * A new choice, 'Energies and gradients from:', with the Thermochemistry step's
      choices: the flowchart's 'Model Chemistry' (the default) or the step's
      'subflowchart'. The geomeTRIC and ASE optimizers can now use the Model
      Chemistry: over a warm MDI engine when the program has one here, otherwise as
      one separate calculation per step. Before, only the SEAMM optimizer could,
      and the others stopped at the first gradient with "There are no energy and
      gradients in properties.json". The SEAMM optimizer always uses the Model
      Chemistry, and now also runs over separate calculations when the program has
      no MDI engine here, rather than stopping.
    * Flowcharts saved before this release keep their source: the subflowchart for
      the geomeTRIC and ASE optimizers, the Model Chemistry for the SEAMM optimizer.
    * When the chosen source is missing and the other is there -- no Model
      Chemistry before the step, or an empty subflowchart -- the other is used, with
      a note, and the flowchart check warns, offering to change the choice.
    * A Structure step with neither is reported by the flowchart check, so the
      editor's Run, ``seamm-flowchart validate`` and the MCP server stop it before
      it is queued; at run time it stops with an error that says what to add. Sub-
      steps that save no energy and gradients to JSON also get an error saying so.
    * The dialog shows the Flowchart tab and the step-subdirectory options only for
      the subflowchart, and says in red what is missing for the choice.

2026.11.25 -- The SEAMM optimizer's minimum check is reported
    * With the SEAMM optimizer and a cheap level, the results now say whether the end
      point is a minimum, its lowest frequency, the saddle escapes and the
      cheap-level modes confirmed, and how many cheap-level Hessians were computed;
      before, the check ran but its verdict was not recorded.
    * Results no longer carry over from an earlier run of the step, e.g. in a loop.

2026.11.24.1 -- Bugfix: the presets' gradient criteria are in E_h/bohr
    * The convergence presets (Gaussian, NWChem, QChem, Turbomole, Inter-fragment) gave
      the programs' gradient criteria, which are in E_h/bohr, labeled E_h/Å. With
      geomeTRIC the two errors cancelled, so its presets are unchanged; ASE used them
      1.89 times too tight (Gaussian's maximum force 0.0123 eV/Å instead of 0.0231).
      The default custom gradient criteria are also in E_h/bohr now. Needs
      seamm-geometric 2026.11.24.1 or later.
    * structure-step now requires seamm-optimizer and seamm-exec 2026.11.24, so the
      SEAMM optimizer works in any installation; before, an installation without
      seamm-optimizer stopped with an error when it was chosen.

2026.11.24 -- The SEAMM optimizer
    * A new optimizer, SEAMM, run over the flowchart's Model Chemistry (a warm MDI
      engine): internal coordinates (TRIC by default), Lindh's or geomeTRIC's model
      Hessian or the Hessian of a cheap level such as GFN2-xTB, and a check that the
      end point is a minimum. Needs seamm-exec 2026.11.24 and seamm-optimizer.
    * The dialog shows the coordinate system for the optimizers that have one (also
      geomeTRIC, where it was not shown) and only the options that apply.

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

