scVarID
=============
**scVarID** is a streamlined pipeline designed for processing and classifying single-cell variant data. By integrating variant files (VCF/BCF/TXT) with BAM files, scVarID efficiently extracts read information, identifies overlaps with variants, and categorizes reads into reference (ref), alternative (alt), missing, or unknown classifications. Leveraging parallel processing, scVarID ensures rapid and scalable analysis suitable for large genomic datasets.

## Table of Contents
- [Features](#features)
- [Installation](#installation)
- [Usage](#usage)
  - [Command-Line Arguments](#command-line-arguments)
  - [Example](#example)
- [Project Structure](#project-structure)
- [Dependencies](#dependencies)
- [Logging](#logging)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Features
* **Variant Processing** Handles multiple variant file formats (VCF, BCF, TXT) and filters variants based on specific criteria.
* **Read Extraction**: Extracts detailed read information from BAM files, including cell barcodes.
* **Overlap Analysis**: Utilizes PyRanges for efficient overlap detection between reads and variants.
* **Classification**: Categorizes reads into ref, alt, missing, or unknown using parallel processing for enhanced performance.
* **Logging**: Comprehensive logging to monitor progress and debug issues.
* **Scalable**: Supports parallel processing to handle large datasets efficiently.
* **Output**: Generates sparse matrices in HDF5 format and stores variant-barcode mappings for downstream analysis.

## Installation
1. **Clone the Repository**
```bash
git clone https://github.com/NCC-sdklab/scVarID.git
cd scVarID
```
<pre data-lang="bash">
<code>git clone https://github.com/NCC-sdklab/scVarID.git</code>
</pre>
2. **Install Dependencies**
```bash
pip install -r requirements.txt
```

## Usage
scVarID is a command-line tool. Below are instructions on how to use it effectively.
### Command-Line Arguments
* `--variant-files`: **(Required)** List of sample and variant file path pairs. Each origin must be paried with a file path.
    * **Format**: `SAMPLE1 FILE_PATH1 SAMPLE2 FILE_PATH2 ...`
    * **Example Samples**: `normal`, `tumor`, `HG002`
* `--bam-path`: **(Required)** Path to the BAM file.
* `--barcode-path`: **(Optional)** Path to the barcode file. If not provided, barcodes will be generated from the data.
* `--save-dir`: **(Required)** Directory where the results will be saved.
* `--num-cores`: **(Optional)** Number of CPU cores to use for parallel processing. Default is `4`.
* `--ref-alt-only`: **(Optional)** If set, only ref and alt classifications will be computed, skipping missing and unknown.
