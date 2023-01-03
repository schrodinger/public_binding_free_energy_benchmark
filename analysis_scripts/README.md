## Analysis scripts
The following scripts were used to generate the main results and plots for accompanying maunscript. 
* `process_fep_benchmark.py`: calculate and print the overall statistics for the FEP+ benchmark. FMP files or CSV files 
can be used as inputs. Using FMP files produces the most accurate results owing to the proper handling of ligands with 
multple protomers or tautomers.
* `process_experimental_survey`: calculate and print the overall statistics of the experimental reproduceability survey.
* `generate_boxplots_and_histograms.py`: generate boxplots and histograms that analyze the distribution of errors in 
both the experimental survey and FEP+ benchmark.
* `generate_scatter_plots.py`: generate example scatter plots from the experimental survey and FEP+ benchmark.
* `scatterplot_data/`: the directory that contains the files used in the scatter plot.
* `write_group_summary_tables.py`: Write tables that summarize the error of each data set in the FEP+ benchmark.
* `print_latex_tables.py`: Print out latex formatted tables of each groups results. Requires a Schrodinger installation.
Please see the doc-strings in each script to see how run each script. The above scripts use functions in the following
 files:
* `helper_functions.py`
* `analysis_functions.py`

## Python dependencies
* `numpy`
* `scipy`
* `matplotlib`
* `pandas`
* `schrodinger` (optional)

If `schrodinger` is not installed (or accessible in python), then all the scrips listed above can use `csv` files as 
input (as opposed to `fmp` files in `schrodinger` was available). The correctly formatted `csv` files for all the 
analysis can be found in `../21_4_results/ligand_predictions` of this repository. However, using `csv` files instead of
`fmp` files as inputs produces statistics that miscount ligands that have multiple tautomers or protomers in 
perturbation graph. This miscounting does not occur with schrodinger is intalled and `fmp` files are used as input. 

## FEP data structure
The FEP+ benchmark data is expected to be contained in subdirectories of a single directory. If the uppermost directory 
is called 'upper_dir', then the FMP (or CSV) files are expected to be in subdirectories. For instance, with 2 
groups 'group_1' and 'group_2' - with 3 output files each - the following directory structure should be used: 
```
    upper_dir
        ├── group_1 
              ├── system_1_out.fmp 
              ├── system_2_out.fmp 
              ├── system_3_out.fmp 
        ├── group_2 
              ├── system_A_out.fmp 
              ├── system_B_out.fmp 
              ├── system_C_out.fmp 
```
This structure is present in `../21_4_results/`.