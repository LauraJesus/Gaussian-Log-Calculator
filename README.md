# scriptLogs

Used to extract ESP and NBO charges from Gaussian output files.

## Description

The `scriptLog.py` script reads `.log` and `.txt` files within a folder and generates HTML reports containing the found ESP and NBO charges.

Currently, it only works for ESP and NBO charges under the specific Gaussian calculation parameters used by the script.

## Features

- Processes `.log` and `.txt` files.
- Extracts ESP charges from the Gaussian `ESP charges:` text block.
- Extracts NBO charges from the Gaussian `Summary of Natural Population Analysis:` text block.
- Generates HTML tables displaying the mean and standard deviation per solute atom.

## Usage

```bash
python scriptLog.py <folder> --solute <n> --output <output_folder>
# or 
python scriptLog.py # (The script will manually prompt for the folder name and number of solute atoms)