# Unrelated Machine Scheduling with Speed Scaling

This archive is distributed under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper *Unrelated Machine Scheduling with Speed Scaling* by the authors.

## Description
This project implements the experiments for the unrelated machine scheduling problem with speed scaling, written in C. It requires **Gurobi 13.01** to be configured. All programs are compiled with Microsoft Visual Studio 2022 and executed on a machine equipped with an Intel Core i5-1135G7 CPU (2.4 GHz) and 16 GB RAM, running Windows 11 (64-bit).

### Code Structure
- `Unrelated/src/`: Contains C source code files.
  - `main.c`: Program entry point.
  - `globals.c`, `MKP.c`, `PBD.c`, `SOCP.c`, `utils.c`: Core logic and algorithmic implementations.
  - `def.h`: Definitions and configurations.
- `Unrelated/data/`: Data required/used for experiments.
- `Unrelated/results/`: Outputs of the experiments.
- `x64/Release/`: Directory containing the compiled executable (`Unrelated.exe`), run script (`run.bat`), data and results folder containing `analysis.py`.

### Benchmark Instances Data
The benchmark instances used in this study are adapted from the classical GAP instances, available at [http://www.al.cm.is.nagoya-u.ac.jp/~yagiura/gap/](http://www.al.cm.is.nagoya-u.ac.jp/~yagiura/gap/). Specifically, this dataset consists of five series, denoted by A--E, with the number of machines $m$ ranging from 5 to 80 and the number of jobs $n$ ranging from 100 to 1600. For example, instance `a05100` corresponds to assigning 100 jobs to 5 machines in series A, while `e801600` represents assigning 1600 jobs to 80 machines in series E. In total, the dataset contains 57 instances of varying sizes.

The data format of the instances located in `data/unrelated` is as follows:

```text
number of agents (m), number of jobs (n)
 for each agent i (i=1,...,m) in turn:
     cost of allocating job j to agent i (j=1,...,n)
 for each agent i (i=1,...,m) in turn:
     resource consumed in allocating job j to agent i (j=1,...,n)
 resource capacity of agent i (i=1,...,m)
 for each agent i (i=1,...,m) in turn:
     lower bound of resource consumed in allocating job j to agent i (j=1,...,n)
 for each agent i (i=1,...,m) in turn:
     upper bound of resource consumed in allocating job j to agent i (j=1,...,n)
```

### Using the Compiled Executable
The compiled program `Unrelated.exe` can be found in the `x64\Release` folder. It can be executed with the following arguments:

```cmd
Usage: Unrelated.exe <input_file> <P> <Method> [mode]

Arguments:
  <input_file> : Path to the file containing instance names.
                 Optional parameters for different instance sets: A_.txt, B_.txt, C_.txt, D_.txt, E_.txt.
  <P>          : The P parameter (double value). It controls the function $\phi_i(x)$.
                 In the first four settings, all machines share the same function $\phi_i(x) = x^{1-P}$:
                 - 1.5 : $\phi_i(x) = x^{-0.5}$
                 - 2.0 : $\phi_i(x) = x^{-1}$
                 - 2.5 : $\phi_i(x) = x^{-1.5}$
                 - 3.0 : $\phi_i(x) = x^{-2}$
                 - 0   : Mixed configuration, $\phi_i(x)$ depends on machine index $i$ and cyclically takes values in $\{x^{-0.5}, x^{-1}, x^{-1.5}, x^{-2}\}$.
  <Method>     : 'SOCP' or 'PBD'.
                 - 'SOCP' uses Gurobi to solve the MISOCP model.
                 - 'PBD' uses the Benders Decomposition method.
  [mode]       : (Optional) 'single' or 'multi'. Defaults to 'single' if omitted.
                 For PBD, this selects between the single-cut version and the multi-cut version.

Examples:
  Unrelated.exe list.txt 2.0 SOCP
  Unrelated.exe list.txt 2.0 PBD          (Runs as single)
  Unrelated.exe list.txt 2.0 PBD multi    (Runs as multi)
```

It is generally recommended to execute the experiments via the batch script provided. 

### Reproducing Experiments via `run.bat`
To extract the raw data required for the experimental section of the paper, follow these steps:
1. Open a command prompt.
2. Navigate to the x64 Release directory:
   ```cmd
   cd C:\Users\Young\Desktop\ORArticle\Codes\x64\Release
   ```
3. Run the batch file to execute the experiments:
   ```cmd
   run.bat
   ```
   This script automatically runs `Unrelated.exe` with the configurations used in our experiments.

### Generating the Report Table
Once the experiments have finished generating data via `run.bat`, follow these steps to process the raw output into a final table formatted for the paper (`results_final_summary.xlsx`):

1. Switch into the results folder:
   ```cmd
   cd C:\Users\Young\Desktop\ORArticle\Codes\x64\Release\results\unrelated
   ```
2. Run the data processing script using Python:
   ```cmd
   python analysis.py
   ```
   This will summarize the raw metrics and generate `results_final_summary.xlsx`, replicating the experiment tables shown in the paper.


