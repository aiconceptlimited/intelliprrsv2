# IntelliPRRSV2

**IntelliPRRSV2: An Integrated Computational Framework for Comparative and Evolutionary Analysis of PRRSV Genomes**

IntelliPRRSV2 is a computational framework developed for comparative and evolutionary analysis of porcine reproductive and respiratory syndrome virus (PRRSV) genomes.

The platform integrates sequence acquisition, quality control, multiple-sequence analysis, nucleotide-discordance profiling, phylogenetic analysis, lineage classification, geographic metadata extraction, exploratory miRNA interaction analysis, and vaccine-reference sequence profiling.

## Scientific scope

The documented manuscript analysis comprised:

- **1,866** isolates in the validated phylogenetic/lineage cohort
- **1,627** sequences in the multiple-sequence alignment
- **17,741** alignment positions
- **27,164** observed nucleotide-discordance records
- **12,965** unique alignment positions containing observed discordance
- **1,953** geographic metadata records covering **29 countries**
- **681,525** exploratory miRNA-interaction records involving **388 miRNAs**
- vaccine-reference profiling involving **11 reference sequences**

These values describe the documented manuscript analysis and should not be interpreted as live database counts.

## Main computational workflow

The authoritative production workflow is:

    scripts/run_pipeline_ultrafast_real.sh

The workflow comprises the following major stages:

1. PRRSV sequence acquisition from NCBI
2. Sequence quality control and cleaning
3. Nucleotide-discordance profiling
4. ORF annotation
5. Phylogenetic reconstruction
6. Lineage classification
7. Vaccine-reference sequence profiling
8. Exploratory miRNA interaction analysis
9. Geographic metadata extraction
10. Dashboard data refresh

## Repository structure

    app/                  Application package
    dashboard/            Streamlit dashboard
    scripts/              Production and supporting computational scripts
    tests/                Test package
    utils/                Utility modules
    reference/            Curated reference sequences and resources
    paper/                Manuscript-supporting figures, data and scripts
    figures/              Submitted manuscript figures
    data/                 Runtime data directory
    results/              Generated analysis results

Generated datasets, runtime outputs, logs, local environments, credentials, and other deployment artifacts are excluded from version control.

## Reproducibility

The repository separates reproducible computational source code from generated runtime data.

The production pipeline is implemented through the scripts under `scripts/`, with project-relative paths defined in:

    scripts/paths.py

Database configuration is supplied through environment variables rather than committed credentials.

Create a local environment using:

    python -m venv .venv
    source .venv/bin/activate
    pip install -r requirements.txt

Copy the example configuration when database-backed execution is required:

    cp .env.example .env

Populate the required environment variables with credentials appropriate to the local deployment.

## External computational dependencies

The workflow uses established computational software including:

- Python
- Biopython
- MAFFT
- FastTree
- NumPy
- pandas
- SciPy
- scikit-learn
- SQLAlchemy
- MySQL Connector/Python
- Streamlit
- Plotly

Exact Python package versions are specified in `requirements.txt`.

## Manuscript-supporting materials

The `paper/` directory contains supporting figure data and figure-generation scripts associated with the documented manuscript analysis.

The submitted manuscript figures preserved under:

    figures/submitted_manuscript/

are retained as the authoritative submitted figure assets.

The figure-generation scripts and supporting data are provided as computational documentation and should not be assumed to reproduce every submitted figure without the original analysis environment and data.

## Scientific limitations

The IntelliPRRSV2 framework includes exploratory computational modules whose outputs require biological interpretation and validation.

In particular:

- the nucleotide-discordance metric is descriptive and is not a substitution-rate estimate;
- geographic counts describe available metadata records and are not prevalence estimates;
- the lineage classifier is implementation-specific and is not presented as a contemporary genotype-specific classification framework;
- the miRNA module produces hypothesis-generating computational predictions requiring biological validation;
- vaccine-reference profiling is not experimentally calibrated as a probability of vaccine escape;
- the documented study did not perform formal recombination, codon-selection, or temporal/phylodynamic inference.

## Data provenance

The genomic analyses use publicly available PRRSV sequence data obtained from the NCBI sequence database.

Derived computational outputs are generated from the analysis workflow. Large raw datasets and generated runtime outputs are intentionally not committed to this repository.

## Manuscript

**IntelliPRRSV2: An Integrated Computational Framework for Comparative and Evolutionary Analysis of PRRSV Genomes**

Author:

**Abubakar Garba**

The manuscript associated with this repository is a scientific work describing the computational framework and its documented analysis.

## Contact

**Abubakar Garba**

AI Concept Limited  
Computational Research Laboratory

Email: `abubakar.garba@aiconceptlimited.com.ng`

## License

A repository license will be added after the publication and intellectual-property requirements for the associated scientific work and software have been finalized.
