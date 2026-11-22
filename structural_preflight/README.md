# PABPC1–RNA preflight — incomplete research tool

This is a **preflight**, not a complete structural preparation, mutation, minimization, MD, or MM/GBSA pipeline. No output should be treated as a validated or clash-repaired model. It was made because the requested revised 31-nt PDB is not available, and the preceding soft-repair procedure was found capable of changing stereochemistry. Neither an output PDB nor claims of scientific validation can be provided on that basis.

## Notebook

Open `PABPC1_RNA_structural_preflight.ipynb` in Google Colab. Run the cells in order. Upload `revised_template(1).pdb` and optionally the two revised Excel workbooks. The notebook installs pinned libraries, runs limited tests, compares the input sequence against UniProt P11940, performs a preliminary structural audit and produces a report ZIP. The notebook includes its own code; the separate `.py` files are provided for inspection and local testing.

The workbook report follows the stated convention: take 16 observed UTR bases ending at the most 3′ non-A, **retain** the original terminal A run and then append A15. `candidate_sequence` is only a planning sequence. It must not be used as a fixed-length 31-nt variant without a validated extension workflow. It is not an experimentally validated construct. Do not add terminal A residues twice in downstream preparations.

The raw PDB is preserved. A 1-Å CRYST1 box is reported; this tool does not remove it. The application must remove or override any placeholder box before simulation. Contact count is the number of protein–RNA heavy-atom pairs <2 Å, not the number of residues or a force-field clash score. Protein/RNA cannot be deemed stereochemically sound by this audit. Atom completeness checks cover canonical RNA heavy atoms (with special handling of the 5′ phosphate) and protein N/CA/C/O, not all protein side chains. Residue names and chain IDs must still be reviewed. Bio.PDB strict parsing can stop on malformed duplicate records.

## What is and is not implemented

| Feature | Status |
| --- | --- |
| Workbook row counts, UTR consistency, native terminal-A and candidate-length report | Implemented; tested on supplied revised workbooks |
| Input chain, sequence, limited missing-atom, RNA-linkage, peptide-linkage, CRYST1 and protein–RNA contact checks | Implemented; synthetic fixture tested, plus real run on an earlier, different PDB |
| CSV checks for 31-nt canonical sequences and explicit modification annotation | Implemented; test fixtures only; rejects unsupported longer inputs |
| OpenMM ff19SB/OL3/TIP3P XML loading and harmonic-restraint conversion | Tested in OpenMM 8.6.1; tiny mathematical minimization only |
| Full side-chain completeness, absolute stereochemistry, force-field topology matching or terminal chemistry | Not implemented |
| Revised reference audit and stereochemistry-preserving minimization | Blocked: revised PDB missing; no validated algorithm here |
| Canonical A→G/C/U builder, ChimeraX execution/scripts, mutation smoke tests | Not implemented or tested |
| Extension of the RNA beyond 31 nt | Not implemented; required for 47 of 61 planning constructs |
| Modified nucleotide geometry and parameters | Not implemented; unsupported modifications are marked `requires_parameterization` in CSV checks |
| Full-batch structural generation, per-variant minimization and restart | Not implemented |

`tests/test_physics.py` uses a one-particle harmonic model to confirm force constant and OpenMM minimization, **not** a test of a molecular complex. Example reports are not included in the repository because they contain user-provided sequences or model information; generated reports should be handled as private research data.

## Locally

With Python, OpenMM, NumPy, SciPy, Biopython and openpyxl installed:

```bash
python -m unittest discover -s tests -v
python workbook_report.py /path/Muller_benchmark.xlsx /path/Muller_Position_Top_Two_Chemical_Modification_Library.xlsx workbook_length_planning.csv
```

The next critical step is to obtain the revised PDB and independently audit the complete atom and stereochemical geometry, then validate a conservative clash-handling strategy and a real mutation smoke test. The 47 longer constructs additionally need explicit direction, positioning, and validation for RNA extension. It is scientifically incorrect to report all clashes resolved or any binding affinity from these artifacts.

Sources: [OpenMM force fields](https://docs.openmm.org/latest/userguide/application/02_running_sims.html); [ChimeraX swap command source](https://github.com/RBVI/ChimeraX/blob/develop/src/bundles/swap_res/src/cmd.py); [UniProt P11940](https://www.uniprot.org/uniprotkb/P11940/entry).
