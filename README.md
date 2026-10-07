# Obesity Target Prioritization Platform

Interactive Streamlit platform for exploring the obesity gene candidate pool,
prioritized Top 50, pathway enrichment, score breakdown, Ensembl records,
hallmarks, and PubMed search summaries. The app uses a target-and-molecule logo
and a matching target icon in the browser tab.

The Overview tab includes an interactive Established/Emerging/Novel gene
classification pie chart. Click a slice to filter its genes, then click a gene
row to display matching Biological Process, Cellular Component, and Molecular
Function enrichment terms. GO links open the matching term in QuickGO.

## Run on Windows

From PowerShell, change to the dashboard directory and install the listed
dependencies:

```powershell
cd C:\Users\Curad\clawbio-lab\obesity_dashboard
python -m pip install -r requirements.txt
python -m streamlit run app.py
```

To install and run with the workspace virtual environment instead:

```powershell
cd C:\Users\Curad\clawbio-lab\obesity_dashboard
& ..\.venv\Scripts\python.exe -m pip install -r requirements.txt
& ..\.venv\Scripts\python.exe -m streamlit run app.py
```

Streamlit prints a local URL (normally `http://localhost:8501`). Open it in a
browser. Press `Ctrl+C` in the terminal to stop the server.

## Data loading

The app first loads files from an included `data` directory or the repository
root, so it works when deployed to Streamlit Community Cloud and does not
depend on a local Windows path. When neither location contains the candidate
CSV, it falls back to the original OneDrive data directory on the developer's
machine. Supporting both layouts lets the dashboard run from either a
`data/` subfolder or a flat repository containing the files alongside `app.py`.

Use the sidebar uploader to upload one or more CSVs. Uploaded files replace the
matching default dataset by filename. Recognized dataset names include
`candidate_genes_union`, `top50`, `go_bp`, `go_cc`, `go_mf`, `kegg`, `reactome`,
`ensembl_gene_records`, `gene_specific_pubmed`, `hallmark_mapping`, and
`score_breakdown`.

The included Excel workbook is also loaded automatically when present. Its `TOP_50_GENES`,
`HALLMARK_MAPPING`, `EVIDENCE_TABLE`, and `SCORE_BREAKDOWN` sheets supply the
prioritized-gene, hallmark, evidence, and score views. Those sheets take
precedence over matching CSV defaults; an uploaded XLSX workbook takes
precedence over the attached workbook for any matching sheets. Other supplied
CSV files, including the candidate pool and pathway enrichments, continue to
load alongside it.

To upload a workbook, select an `.xlsx` file in the sidebar uploader. The
dashboard reads the workbook locally and recognizes the sheet names above; no
Excel-reading package is required.

Gene and trait filters apply to gene tables and pathway gene overlaps. Trait
matching uses the phenotype/evidence text in the supplied Top 50 file; it does
not infer phenotype associations for candidate-only genes. The GO ID filter
applies to GO identifiers parsed from pathway terms. Tables show direct Ensembl,
PubMed, QuickGO, Reactome, and KEGG links where identifiers are present.

The dashboard uses Streamlit charts plus Altair for the interactive evidence
classification pie chart.

## Create a shareable Streamlit Community Cloud URL

1. Create a **public GitHub repository** for this dashboard.
2. Upload `app.py`, `requirements.txt`, `obesity_target_logo.svg`, and the
   supplied data files. The files may be placed in the repository root or in
   its `data` folder. Only publish these files if you are comfortable
   making the included gene evidence and analysis outputs public.
3. Sign in to [Streamlit Community Cloud](https://share.streamlit.io/) with the
   GitHub account that owns the repository.
4. Choose **Create app**, select the repository and branch, and set the main
   file path to `app.py`.
5. Deploy. Streamlit will provide a shareable `https://<app-name>.streamlit.app`
   URL. Anyone with the link can access the public app.

The workbook is optional for the deployed app because its four displayed
tables are also available as CSVs in `data`. The dashboard's built-in Excel
reader needs no additional requirements.
