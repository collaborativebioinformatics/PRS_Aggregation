# PRS Locus Viewer

**PRS Locus Viewer** is an interactive Dash application designed for visualizing Polygenic Risk Score (PRS) variants in their genomic context. It allows researchers to inspect SNP effect sizes across multiple scoring files, map variants to nearby genes, and explore specific loci dynamically.

## Key Features
- **Multi-Score Comparison:** Overlay effect weights from multiple PRS files side-by-side.
- **Interactive Visualization:** Clickable SNP tracks with heatmap-style coloring based on effect size.
- **Gene Mapping:** Automatically identifies and visualizes genes within a configurable window (e.g., ±25kb) of reported SNPs.
- **Locus Zoom:** Search by rsID or Gene Symbol to zoom into specific genomic regions.

---

## Installation

### 1. Clone the Repository
```bash
git clone https://github.com/humengying0907/prs-locus-viewer.git
cd prs-locus-viewer
```

### 2. Create the Environment

You can set up the required dependencies using Conda. Create a file named `environment.yml` with the content below, or install manually.

**Option A: Using Conda (Recommended)**

```bash
conda env create -f environment.yml
conda activate prs_viewer
```

**Option B: Using pip**

```bash
pip install pandas dash requests pyranges
```

---
## Data Preparation

This tool runs locally and requires **Harmonized Scoring Files** (containing `hm_chr` and `hm_pos` columns) as input. A future web-based version is planned to eliminate the need for manual file uploads.

### 1. Downloading Harmonized Scores

You can download harmonized files directly from the [PGS Catalog FTP site](https://www.pgscatalog.org/):


- **Navigate to:** `PGS_ID → ScoringFiles → Harmonized`
- **Example URL:**  
  https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/PGS003402/ScoringFiles/Harmonized/
- **Example file:**
  PGS003402_hmPOS_GRCh37.txt.gz

Ensure the downloaded files include at least the following columns:

- `hm_chr`
- `hm_pos`
- `effect_weight`
- `rsID`

### 2. Using Example Data

For quick testing, sample harmonized scoring files are provided in the `example_input/` directory of this repository.

---

## Usage

To launch the viewer, run the main script with your input PRS text files and the target chromosome.

```bash
python prs_locus_viewer.py [input_file_1] [input_file_2] ... --chr [chromosome]
```

- `input_file`: Path to PRS scoring files (must contain columns: `hm_chr`, `hm_pos`, `effect_weight`, `rsID`).
- `--chr`: The specific chromosome to load (e.g., `chr5`).

### Example Run
```bash
python prs_locus_viewer.py example_input/PGS001818_hmPOS_GRCh37.txt example_input/PGS000020_hmPOS_GRCh37.txt example_input/PGS003402_hmPOS_GRCh37.txt --chr chr5
```

---

## Startup Output

Upon running the command, you should see output like:

```text
============================================================
  PRS LOCUS VIEWER | Chromosome: chr5
============================================================

[1/3] Loading data for chr5...
[2/3] Calculating gene coverage (Window: ±25kb)...
      Targeting 2,838 unique SNPs

[info] Calculating gene coverage (±25kb) for 2838 unique SNPs...
[==============================] 100% (processing chr5)
[info] Gene coverage calculation finished.
[3/3] Calculation complete.

------------------------------------------------------------
 [READY] Viewer active at: http://127.0.0.1:8050/
         (Press CTRL+C to stop)
------------------------------------------------------------

Dash is running on http://127.0.0.1:8050/

 * Serving Flask app 'prs_locus_viewer'
 * Debug mode: on
```

Once the server is ready, open your web browser and navigate to `http://127.0.0.1:8050/` to start using the tool.

## Interactive Preview
