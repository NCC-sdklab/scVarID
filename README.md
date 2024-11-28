scVarID
=============
**scVarID** is a streamlined pipeline designed for processing and classifying single-cell variant data. By integrating variant files (VCF/BCF/TXT) with BAM files, scVarID efficiently extracts read information, identifies overlaps with variants, and categorizes reads into reference (ref), alternative (alt), missing, or unknown classifications. Leveraging parallel processing, scVarID ensures rapid and scalable analysis suitable for large genomic datasets.

## Table of Contents
* Features
* Installation
* Usage
    * Command-Line Arguments
    * Example
* Project Structure
* Dependencies
* Logging
* Contributing
* License
* Contact

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
git clone https://gihub.com/NCC-sdklab/scVarID.git
cd scVarID
```
2. **Install Dependencies**
```bash
pip install -r requirements.txt
```

## Usage
scVarID is a command-line tool. Below are instructions on how to use it effectively.
### Command-Line Arguments
* `--variant-files`: **(Required)** List of sample and variant file path pairs. Each origin must be paried with a file path.
    * **Format**: `SAMPLE1 FILE_PATH1 SAMPLE2 FILE_PATH2 ...`
    * **Eample Samples: `normal`, `tumor`, `HG002`
* `--bam-path`: **(Required)** 
