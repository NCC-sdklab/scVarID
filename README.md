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
2. **Install Dependencies**
```bash
pip install -r requirements.txt
```

## Usage
scVarID is a command-line tool. Below are instructions on how to use it effectively.
### Command-Line Arguments
* `--variant-files`: **(Required)** List of sample and variant file path pairs. Each origin must be paired with a file path.
    * **Format**: `SAMPLE1 FILE_PATH1 SAMPLE2 FILE_PATH2 ...`
    * **Example Samples**: `normal`, `tumor`, `HG002`
* `--bam-path`: **(Required)** Path to the BAM file.
* `--barcode-path`: **(Optional)** Path to the barcode file. If not provided, barcodes will be generated from the data.
* `--save-dir`: **(Required)** Directory where the results will be saved.
* `--num-cores`: **(Optional)** Number of CPU cores to use for parallel processing. Default is `4`.
* `--ref-alt-only`: **(Optional)** If set, only ref and alt classifications will be computed, skipping missing and unknown.
### Example
```bash
python main.py \
  --variant-files normal variants_normal.vcf tumor variants_tumor.vcf \
  --bam-path sample.bam \
  --barcode-path barcodes.txt \
  --save-dir results \
  --num-cores 8
```
#### Explanation:
* Processes two variants files:
    * `variant_normal.vcf` with sample `normal`
    * `variant_tumor.vcf` with sample `tumor`
* Uses `sample.bam` as the BAM file.
* Reads barcodes from `barcodes.txt`
* Saves all outputs to the `results` directory.
* Utilizes `8` CPU cores for parallel processing.
#### Running with one variant file
(Edit)
#### Running without barcode file
If you do not have a barcode file, you can omit the `--barcode-path` argument:
```bash
python main.py \
  --variant-files normal variants_normal.vcf tumor variants_tumor.vcf \
  --bam-path sample.bam \
  --save-dir results \
  --num-cores 8
```

## Project Structure
```
scVarID/
├── README.md
├── requirements.txt
├── src/
│   ├── classification.py
│   ├── main.py
│   ├── read_processing.py
│   ├── variant_processing.py
│   └── utils.py
└── results/
    ├── processing.log
    ├── ref_matrix.h5
    ├── alt_matrix.h5
    ├── missing_matrix.h5
    ├── unknown_matrix.h5
    └── variant_barcode_mappings.pkl
```

## Dependencies
All dependencies are listed in `requirements.txt`. Below is a summary of key packages used:
* **h5py**: Handling HDF5 file formats.
* **joblib**: For parallel processing.
* **logging**: Comprehensive logging of the pipeline's progress.
* **pandas**: Data manipulation and analysis.
* **pyranges**: Efficient genomic interval operations.
* **pysam**: For reading and manipulating BAM files.
* **scipy**: Scientific computations, especially for sparse matirces.

To install all dependencies, run:
```bash
pip install -r requirements.txt
```

## Logging
scVarID utilizes Python's `logging` module to provide detailed logs of the pipeline's execution. Logs are saved both to the console and to a log file named `processing.log` within the specified `--save-dir`.
**Log Levels:**
* `INFO`: General progress and import milestones.
* `WARNING`: Issues that do not halt the pipeline but may require attention.
* `ERROR`: Critical issues that cause the pipeline to terminate.
**Example Log Entries:**
```yaml
2024-11-28 10:00:00,000:INFO:=== [scVarID] Program Started ===

2024-11-28 10:00:00,001:INFO:=== Step Start: Process variant files ===
2024-11-28 10:05:30,500:INFO:Extracted 1500 variants from variants_normal.vcf
2024-11-28 10:05:30,501:INFO:Extracted 2000 variants from variants_tumor.vcf
2024-11-28 10:05:30,502:INFO:Filtered variants count: 3000
2024-11-28 10:05:30,503:INFO:=== Step End: Process variant files | Elapsed Time: 00h 00m 03s ===
...
2024-11-28 10:20:30,003:INFO:=== All steps completed successfully ===
```

## Contributing
(Edit)

## Contact
For any questions, issues, or suggestions, please open an issue on [GitHub Issues](https://github.com/NCC-sdklab/scVarID/issues) or contact [dshin@ncc.re.kr](mailto:dshin@ncc.re.kr)
