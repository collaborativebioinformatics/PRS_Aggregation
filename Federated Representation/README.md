#  Federated Representation Learning for Polygenic Risk Scores (PRS / PGS)

This repository implements a **modular, interpretable representation learning framework for Polygenic Risk Scores (PGS)** that integrates heterogeneous PGS Catalog scoring files across ancestries, harmonizes variants at the locus level, enriches them with biological annotations, and learns **PGS-level and variant-level embeddings**.

The framework is explicitly designed to align with **federated learning principles**, enabling analysis of **shared vs ancestry-specific genetic architecture** without requiring individual-level genotype data.

---

##  Motivation

Polygenic risk scores are:

- ancestry-specific  
- heterogeneous in format  
- difficult to compare across studies  
- often biologically opaque  

At the same time, federated learning in genomics requires:

- **client-level representations**
- **shared vs private signal decomposition**
- **privacy-preserving aggregation**

This project bridges these gaps by:

- treating each PGS as a *client*
- learning structured embeddings at multiple biological levels
- explicitly analyzing cross-ancestry sharing of variants

---

##  Core Ideas

- **PGS as representations**  
  Each PGS is embedded as a vector capturing how it distributes genetic risk across loci.

- **Variants as representations**  
  Each locus is embedded based on its effect profile across multiple PGSs.

- **Cross-ancestry structure**  
  Variants are classified as ancestry-specific or shared across populations.

- **Federated alignment**  
  The entire pipeline mirrors client-local learning + global aggregation.

---

##  End-to-End Workflow

### **Step 1 — Parse & Harmonize PGS Files**
**Script:** `01_parse_and_harmonize_pgs.py`

- Parses heterogeneous PGS Catalog scoring files  
- Extracts metadata (ancestry, trait, genome build, weight type)  
- Harmonizes coordinates using `hm_chr` / `hm_pos`  
- Defines canonical loci as: locus_id = hm_chr : hm_pos

- Normalizes effect sizes into a common numeric space  

**Output**
- `pgs_variants_master.parquet`
- `pgs_metadata.parquet`

---

### **Step 2 — Variant Annotation (Biological Context)**
**Script:** `01b_annotate_variants_ensembl_vep.py`

Annotates each locus using Ensembl VEP:

- associated gene(s)
- gene region (exonic / intronic / intergenic)
- mutation type (missense, synonymous, splice, etc.)
- most severe consequence
- number of genes affected

**Output**
- `pgs_variants_master_annotated.parquet`

---

### **Step 3 — Sparse Feature Matrix Construction**
**Script:** `02_build_sparse_matrices.py`

Constructs a **PGS × locus sparse matrix**:

- rows = PGSs (clients)
- columns = loci
- values = normalized effect size (`weight_scaled`)

Options include:

- top-K variant selection per PGS
- annotation-aware feature schemas
- memory-safe sparse storage (CSR)

**Output**
- `X_pgs_feat_csr.npz`
- `locus_index.json`
- `feature_schema.json`

---

### **Step 4 — PGS-Level Representation Learning**
**Script:** `03_train_eval_embed.py`

Learns **PGS embeddings** using:

- SVD or autoencoder
- train/test splits
- masked reconstruction loss
- cosine distance evaluation

Each PGS → a compact latent vector representing its genetic architecture.

**Output**
- `pgs_embeddings.pt`
- `pgs_input_cosine_distance.csv`

---

### **Step 5 — Variant-Level Representation Learning**
**Script:** `03_train_eval_embed.py`

Defines: Variant = [effect in PGS₁, effect in PGS₂, …]

Trains a **variant autoencoder** to learn:

- variant embeddings capturing cross-PGS behavior
- consistency vs heterogeneity across ancestries

Variants can be filtered by minimum PGS support.

**Output**
- `variant_embeddings.csv`

---

### **Step 6 — Visualization & Interpretation**

**Scripts**
- `03b_plot_embeddings.py`
- `03c_plot_variant_embeddings*.py`

#### PGS-level plots
- PCA / UMAP
- colored by ancestry or trait
- labeled by ancestry + disease
- cosine distance heatmaps

#### Variant-level plots
- PCA / UMAP
- colored by:
  - primary ancestry
  - primary trait
  - shared vs ancestry-specific
  - number of ancestries observed
- top variants labeled for interpretation

---

##  Variant Sharing Definitions

| Category | Meaning |
|--------|--------|
| ancestry_specific | Variant appears in only 1 ancestry |
| shared_2 | Appears in exactly 2 ancestries |
| shared_3 | Appears in exactly 3 ancestries |
| shared_4 | Appears in exactly 4 ancestries |
| shared_5plus | Appears in ≥5 ancestries |

These categories are central to interpreting **cross-population genetic architecture**.

---

##  Federated Learning Perspective

This framework naturally aligns with federated learning:

| Genomics Concept | Federated Learning Analogy |
|------------------|----------------------------|
| PGS | Client |
| Ancestry | Client domain |
| Variant | Model parameter |
| Shared variants | Global parameters |
| Ancestry-specific variants | Client-private parameters |
| PGS embedding | Client representation |
| Variant embedding | Parameter representation |

No individual-level genotypes are required.

---

##  Example Use Cases

- Cross-ancestry PRS comparison  
- Identifying ancestry-robust disease loci  
- Federated PRS modeling  
- Variant prioritization across populations  
- Methods development for privacy-preserving genomics  

---
