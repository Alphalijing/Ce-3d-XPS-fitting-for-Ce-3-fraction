# English GUI Guide

Open the English Windows application by double-clicking `Ce3d_XPS_Fitter_v2.0.1_EN.exe` inside its complete application folder. Keep `_internal/` and every other bundled file beside the executable. Copy or move the **whole folder**, not the `.exe` alone. The English copy is separate from the original Chinese application.

If you are using the source folder, double-click `Launch_GUI.bat`. It starts the packaged English application when available, or uses the local Python installation. Before using the source version, install the dependencies from the same folder:

```powershell
python -m pip install -r requirements.txt
```

## Choose the fit settings

1. **Peak model**: LA is selected by default; Voigt is also available. LA uses a Gaussian-convolved LA-type line shape with fixed alpha=0.9 and beta=2. This implementation differs from the CasaXPS `LA(50)` definition; the name does not imply identical parameters or line shapes.
2. **Satellite peaks**: Exclude is selected by default. Select Include to add the satellite peaks.
3. **Sampling mode**: Quick fit runs 8000 steps by default. Convergence validation keeps sampling until the engineering convergence criteria are met, with no preset maximum step count. You can request **Stop and save** during either mode.

## Add files and select an output folder

Use **Add files...** to select one or more `.xls`, `.xlsx`, `.vms`, `.vamas`, `.csv`, `.txt`, `.dat`, or `.xy` files. The format is detected automatically. Files are fitted sequentially, with separate results folders for each input and each run. **Remove selected** removes the highlighted entries; **Clear** empties the input list. Drag-and-drop is available when the optional `tkinterdnd2` dependency is installed.

Use **Browse...** under **Output folder** to choose a writable location. Each run creates a separate folder so earlier results are preserved. GUI activity and complete error details are also saved under `gui_logs/` in the selected output folder.

The **Run log** shows feedback from the sampler at the start and every 2000 steps. The progress line shows actual sampled steps, elapsed time, and available Ce3+ quantiles. Convergence validation has no percentage progress because it has no preset step limit. Batch feedback identifies each file as it runs.

## Advanced settings

**Advanced settings** provides the following controls:

- **Walkers**: 80 by default; must be an even integer of at least 22 without satellites, or at least 30 with satellites.
- **Quick steps**: 8000 by default; must be a multiple of 2000.
- **Burn-in steps**: 4000 by default; must be a multiple of 2000 and at least 2000. For a quick fit, it must be less than the total steps.
- **Seed**: 42 by default; an integer from 0 through 4294967295.
- **VAMAS block** and **Signal index**: zero-based indices; leave blank for automatic selection.
- **KE offset / eV**: required when converting a kinetic-energy axis; confirm the instrument's energy convention.
- **Satellite gamma / eV**: 0.8 eV by default; must be positive.
- **Background**: Fixed Shirley is selected by default. Use **Use background from file** only when the input actually contains the required background column.

Convergence validation disables **Quick steps**. Its runtime depends on the data and diagnostic results. Opening Advanced settings increases the minimum window height so the inputs, output folder, progress, and log remain visible.

## Select a VAMAS block and intensity channel

If a VAMAS file contains multiple blocks or signal channels, select that file in the input list, then click **Select VAMAS block...**. The dialog lists the file's actual XPS blocks, Ce 3d candidates, sample labels, energy axes, and channels. Applying a selection fills the corresponding advanced settings.

When several Ce 3d candidates exist, the first is not silently selected. Choose the measured intensity channel; transmission or TXFN functions cannot serve as spectral intensity. A kinetic-energy axis requires an explicit KE offset based on the instrument convention. A binding-energy axis does not use this offset.

The selected indices and KE offset apply to all VAMAS files in the current run. Run files separately if their block or channel layouts differ. File paths, sample labels, and other imported metadata retain their original language.

## Start, stop, and review results

Click **Start fit** to begin. **Stop and save** requests a safe exit from the sampler and saves available results. Wait for saving to finish. Closing an active window asks whether to stop sampling and exits after saving; choosing No continues the fit.

A stopped run and a completed quick fit do not establish convergence. The interface reports engineering convergence in green only when the saved results show that every requested file has status `converged`. Even after these criteria are met, assess the spectra, residuals, diagnostics, and model assumptions.

When the run ends, the final result line shows **Ce3+ percentage ± one posterior standard deviation**. This display uses the posterior mean; the median and 95% credible interval remain available separately. ± is in percentage points and does not denote the 95% interval.

Click **Open report** to open `Results.html` directly. The report starts with the Ce3+ result, then the component plot and concise selected-scan conditions, followed by component parameters and posterior distributions. Other parameter distributions and complete records appear farther down or in expandable sections. Batch reports list the sample results and link to each report.

The run folder contains `Results.html` plus `figures/`, `tables/` and `details/`. The report embeds its figures, so the HTML alone can be copied for viewing. Keep the folders when you need editable figure/table exports or complete records. Experimental conditions JSON, raw retained posterior samples, configuration, checkpoints and chain blocks are in `details/`.

If a run fails, the interface reports the error and writes complete details to the GUI log. For a batch with failed files, review the batch report and the affected file logs. A stopped run with insufficient retained samples does not display a fabricated Ce3+ estimate.

## Build the English Windows application

```powershell
python -m pip install pyinstaller
python build_windows_exe.py
```

The build script uses PyInstaller from the current Python environment to create a windowed, **one-folder** application. If `tkinterdnd2` is installed, its required data are collected for drag-and-drop. Otherwise, use **Add files...**.

The application is generated under `dist/Ce3d_XPS_Fitter_v2.0.1_EN/`. Distribute that entire folder, preserving `_internal/` and all other contents. Open `Ce3d_XPS_Fitter_v2.0.1_EN.exe` inside it.

## Hidden checks without fitting

Create and destroy a hidden interface without sampling:

```powershell
python ce3d_xps_fitter_gui.py --self-test
```

Write a detailed hidden runtime report, including English interface wording, defaults, four window layouts, and basic model predictions:

```powershell
python ce3d_xps_fitter_gui.py --self-check english_gui_self_check.json
```

These checks do not run MCMC. Window layout checks use a fully transparent, borderless window positioned off-screen. Optional VAMAS inspection is available with `--self-check-input PATH_TO_FILE.vms`.
