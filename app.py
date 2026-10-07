from __future__ import annotations

import re
import math
import posixpath
import zipfile
from pathlib import Path
from typing import Any
from urllib.parse import quote
import xml.etree.ElementTree as ET

import altair as alt
import pandas as pd
import streamlit as st


LOGO_PATH = Path(__file__).resolve().parent / "obesity_target_logo.svg"

st.set_page_config(
    page_title="Obesity Target Prioritization Platform",
    page_icon="🎯",
    layout="wide",
    initial_sidebar_state="expanded",
)

APP_DIR = Path(__file__).resolve().parent
BUNDLED_DATA_ROOT = APP_DIR / "data"
LOCAL_DATA_ROOT = Path(
    r"C:\Users\Curad\OneDrive\Desktop\Obesity\Trial results"
)
DATA_ROOT = (
    BUNDLED_DATA_ROOT
    if (BUNDLED_DATA_ROOT / "candidate_genes_union_434.csv").is_file()
    else APP_DIR
    if (APP_DIR / "candidate_genes_union_434.csv").is_file()
    else LOCAL_DATA_ROOT
)
WORKBOOK_NAME = "Obesity_Top50_Gene_Target_Discovery_ClawBio.xlsx"
LOCAL_WORKBOOK_PATH = (
    APP_DIR.parents[0]
    / "output"
    / "obesity_target_discovery"
    / "run_20260930_1045"
    / "final_results"
    / WORKBOOK_NAME
)
WORKBOOK_PATH = (
    BUNDLED_DATA_ROOT / WORKBOOK_NAME
    if (BUNDLED_DATA_ROOT / WORKBOOK_NAME).is_file()
    else APP_DIR / WORKBOOK_NAME
    if (APP_DIR / WORKBOOK_NAME).is_file()
    else LOCAL_WORKBOOK_PATH
)
GO_TABLES_DIR = (
    Path(__file__).resolve().parents[1]
    / "output"
    / "obesity_target_discovery"
    / "run_20260930_1045"
    / "evidence"
    / "pathway_enrichment"
    / "tables"
)


def go_table_path(name: str, legacy_name: str) -> Path:
    bundled_path = BUNDLED_DATA_ROOT / name
    source_path = GO_TABLES_DIR / name
    if bundled_path.is_file():
        return bundled_path
    if source_path.is_file():
        return source_path
    root_path = APP_DIR / name
    if root_path.is_file():
        return root_path
    return DATA_ROOT / legacy_name

DATASETS: dict[str, tuple[str, Path]] = {
    "candidates": (
        "Candidate pool",
        DATA_ROOT / "candidate_genes_union_434.csv",
    ),
    "top50": (
        "Prioritized Top 50",
        DATA_ROOT / "Obesity_Top50_Gene_Target_Discovery_TOP50.csv",
    ),
    "go_bp": (
        "GO Biological Process",
        go_table_path("go_bp_enrichment.csv", "go_bp_enrichment_MC_done.csv"),
    ),
    "go_cc": (
        "GO Cellular Component",
        go_table_path("go_cc_enrichment.csv", "go_cc_enrichment_MC_done.csv"),
    ),
    "go_mf": (
        "GO Molecular Function",
        go_table_path("go_mf_enrichment.csv", "go_mf_enrichment_MC_done.csv"),
    ),
    "kegg": (
        "KEGG",
        DATA_ROOT / "kegg_enrichment_MC_done.csv",
    ),
    "reactome": (
        "Reactome",
        DATA_ROOT / "reactome_enrichment_MC_done.csv",
    ),
    "ensembl": (
        "Ensembl records",
        DATA_ROOT / "ensembl_gene_records.csv",
    ),
    "pubmed": (
        "PubMed searches",
        DATA_ROOT / "gene_specific_pubmed_searches.csv",
    ),
    "hallmark": (
        "Hallmark mapping",
        DATA_ROOT / "hallmark_mapping.csv",
    ),
    "scores": (
        "Score breakdown",
        DATA_ROOT / "score_breakdown.csv",
    ),
}

PATHWAY_KEYS = ("go_bp", "go_cc", "go_mf", "kegg", "reactome")
WORKBOOK_SHEETS = {
    "TOP_50_GENES": "top50",
    "HALLMARK_MAPPING": "hallmark",
    "EVIDENCE_TABLE": "evidence",
    "SCORE_BREAKDOWN": "scores",
}
DATASET_LABELS = {
    **{key: label for key, (label, _) in DATASETS.items()},
    "evidence": "Evidence table",
}
XLSX_NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "relationships": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
TRAIT_PATTERNS: dict[str, tuple[str, ...]] = {
    "BMI / body mass index": (r"\bbmi\b", r"\bbody mass index\b"),
    "Body weight": (r"\bbody weight\b", r"\bweight\b"),
    "Body-fat percentage": (r"\bbody[- ]fat percentage\b", r"\bbody fat\b"),
    "Waist circumference": (r"\bwaist circumference\b", r"\bwaist\b"),
    "Obesity / adiposity": (r"\bobesity\b", r"\badiposity\b"),
}
GO_ID_PATTERN = re.compile(r"\bGO:\d{7}\b", re.IGNORECASE)
PMID_PATTERN = re.compile(r"\bPMID\s*:?\s*(\d{5,10})\b", re.IGNORECASE)
URL_PATTERN = re.compile(r"https?://[^\s;|,]+", re.IGNORECASE)


def normalise_symbol(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper()


def find_column(frame: pd.DataFrame, *names: str) -> str | None:
    by_lower = {str(column).strip().lower(): column for column in frame.columns}
    for name in names:
        if name.lower() in by_lower:
            return str(by_lower[name.lower()])
    return None


def read_csv_source(source: Any, label: str) -> pd.DataFrame | None:
    try:
        return pd.read_csv(source, encoding="utf-8-sig")
    except (OSError, UnicodeDecodeError, pd.errors.ParserError, ValueError) as exc:
        st.error(f"Could not read **{label}**: {exc}")
        return None


def column_number(cell_reference: str) -> int:
    letters = re.match(r"[A-Z]+", cell_reference.upper())
    if not letters:
        return 0
    number = 0
    for letter in letters.group(0):
        number = number * 26 + ord(letter) - ord("A") + 1
    return number - 1


def read_xlsx_workbook(source: Any) -> dict[str, pd.DataFrame]:
    """Read workbook sheets with the standard library, avoiding an Excel engine dependency."""
    with zipfile.ZipFile(source) as archive:
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(
            archive.read("xl/_rels/workbook.xml.rels")
        )
        targets = {
            item.attrib["Id"]: item.attrib["Target"]
            for item in relationships
        }
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in archive.namelist():
            shared_root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            shared_strings = [
                "".join(text.text or "" for text in item.findall(".//main:t", XLSX_NS))
                for item in shared_root.findall("main:si", XLSX_NS)
            ]

        loaded: dict[str, pd.DataFrame] = {}
        sheets = workbook.findall("main:sheets/main:sheet", XLSX_NS)
        for sheet in sheets:
            dataset_key = WORKBOOK_SHEETS.get(sheet.attrib["name"])
            if dataset_key is None:
                continue
            relationship_id = sheet.attrib[f"{{{XLSX_NS['relationships']}}}id"]
            target = targets[relationship_id]
            if target.startswith("/"):
                sheet_path = target.lstrip("/")
            else:
                sheet_path = posixpath.normpath(
                    posixpath.join("xl", target)
                )
            worksheet = ET.fromstring(archive.read(sheet_path))
            table: list[list[Any]] = []
            for row in worksheet.findall("main:sheetData/main:row", XLSX_NS):
                values: dict[int, Any] = {}
                for cell in row.findall("main:c", XLSX_NS):
                    index = column_number(cell.attrib.get("r", ""))
                    cell_type = cell.attrib.get("t")
                    if cell_type == "inlineStr":
                        value = "".join(
                            text.text or ""
                            for text in cell.findall(".//main:t", XLSX_NS)
                        )
                    else:
                        value_node = cell.find("main:v", XLSX_NS)
                        value = value_node.text if value_node is not None else None
                        if cell_type == "s" and value is not None:
                            value = shared_strings[int(value)]
                        elif cell_type == "b" and value is not None:
                            value = value == "1"
                    values[index] = value
                if values:
                    row_values = [None] * (max(values) + 1)
                    for index, value in values.items():
                        row_values[index] = value
                    table.append(row_values)

            if not table:
                continue
            width = max(map(len, table))
            header = [
                str(value).strip() if value is not None else f"Unnamed: {index}"
                for index, value in enumerate(table[0] + [None] * (width - len(table[0])))
            ]
            normalized_rows = [
                row + [None] * (width - len(row))
                for row in table[1:]
                if any(value is not None for value in row)
            ]
            frame = pd.DataFrame(normalized_rows, columns=header)
            frame = frame.dropna(how="all").reset_index(drop=True)
            loaded[dataset_key] = frame
        return loaded


def detect_uploaded_dataset(filename: str) -> str | None:
    name = filename.casefold().replace("-", "_").replace(" ", "_")
    if "candidate_genes_union" in name:
        return "candidates"
    if "go_bp" in name or "go_biological_process" in name:
        return "go_bp"
    if "go_cc" in name or "go_cellular_component" in name:
        return "go_cc"
    if "go_mf" in name or "go_molecular_function" in name:
        return "go_mf"
    if "kegg" in name:
        return "kegg"
    if "reactome" in name:
        return "reactome"
    if "ensembl_gene_records" in name:
        return "ensembl"
    if "top50" in name or "top_50" in name:
        return "top50"
    if "gene_specific_pubmed" in name or "pubmed" in name:
        return "pubmed"
    if "hallmark_mapping" in name:
        return "hallmark"
    if "score_breakdown" in name:
        return "scores"
    return None


def load_datasets(
    uploaded_files: list[Any],
) -> tuple[dict[str, pd.DataFrame], dict[str, str]]:
    sources: dict[str, Any] = {
        key: path for key, (_, path) in DATASETS.items() if path.is_file()
    }
    loaded: dict[str, pd.DataFrame] = {}
    origins: dict[str, str] = {}
    for key, source in sources.items():
        frame = read_csv_source(source, DATASETS[key][0])
        if frame is not None:
            loaded[key] = frame
            origins[key] = str(source)

    if WORKBOOK_PATH.is_file():
        try:
            workbook_data = read_xlsx_workbook(WORKBOOK_PATH)
            loaded.update(workbook_data)
            origins.update(
                {key: str(WORKBOOK_PATH) for key in workbook_data}
            )
        except (
            OSError,
            KeyError,
            ValueError,
            zipfile.BadZipFile,
            ET.ParseError,
        ) as exc:
            st.error(f"Could not read attached Excel workbook: {exc}")

    uploaded_keys: set[str] = set()
    for uploaded_file in uploaded_files:
        if uploaded_file.name.casefold().endswith(".xlsx"):
            try:
                workbook_data = read_xlsx_workbook(uploaded_file)
                loaded.update(workbook_data)
                origins.update(
                    {key: uploaded_file.name for key in workbook_data}
                )
                if not workbook_data:
                    st.warning(
                        f"No recognized data sheets were found in "
                        f"**{uploaded_file.name}**."
                    )
            except (
                OSError,
                KeyError,
                ValueError,
                zipfile.BadZipFile,
                ET.ParseError,
            ) as exc:
                st.error(f"Could not read **{uploaded_file.name}**: {exc}")
            continue

        key = detect_uploaded_dataset(uploaded_file.name)
        if key is None:
            st.warning(
                f"Unrecognized CSV filename **{uploaded_file.name}**; "
                "it was not loaded. Use a supplied filename or include the "
                "dataset name (for example, `go_bp` or `top50`)."
            )
            continue
        if key in uploaded_keys:
            st.warning(
                f"More than one upload matched **{DATASETS[key][0]}**. "
                f"Using the last matching file: `{uploaded_file.name}`."
            )
        uploaded_keys.add(key)
        frame = read_csv_source(uploaded_file, DATASETS[key][0])
        if frame is not None:
            loaded[key] = frame
            origins[key] = uploaded_file.name

    missing = [DATASETS[key][0] for key in DATASETS if key not in loaded]
    if missing:
        st.info(
            "Some datasets are unavailable. Upload the corresponding CSVs to "
            "include them: " + ", ".join(missing)
        )
    return loaded, origins


def add_gene_key(frame: pd.DataFrame, *column_names: str) -> pd.DataFrame:
    result = frame.copy()
    column = find_column(result, *column_names)
    result["_gene_key"] = (
        result[column].map(normalise_symbol) if column else ""
    )
    return result


def trait_text_for_row(row: pd.Series) -> str:
    columns = (
        "Obesity Phenotype",
        "GWAS Evidence",
        "Obesity Phenotypes",
        "trait",
    )
    values = []
    for key, value in row.items():
        if str(key).strip().lower() in {
            name.lower() for name in columns
        } and pd.notna(value):
            values.append(str(value))
    return " | ".join(values)


def trait_gene_set(top50: pd.DataFrame, selected_traits: list[str]) -> set[str]:
    gene_column = find_column(top50, "Gene Symbol", "gene_symbol")
    if not gene_column:
        return set()
    matched: set[str] = set()
    for _, row in top50.iterrows():
        text = trait_text_for_row(row)
        for trait in selected_traits:
            if any(re.search(pattern, text, re.IGNORECASE) for pattern in TRAIT_PATTERNS[trait]):
                matched.add(normalise_symbol(row[gene_column]))
                break
    return matched


def get_pathway_frame(data: dict[str, pd.DataFrame]) -> pd.DataFrame:
    sections: list[pd.DataFrame] = []
    for key in PATHWAY_KEYS:
        if key not in data:
            continue
        frame = data[key].copy()
        frame["Enrichment source"] = DATASETS[key][0]
        sections.append(frame)
    if not sections:
        return pd.DataFrame()

    result = pd.concat(sections, ignore_index=True, sort=False)
    aliases = {
        "rank": ("rank", "Rank"),
        "term": ("term", "Term", "pathway", "name"),
        "pvalue": ("pvalue", "p_value", "p-value"),
        "adjusted_pvalue": ("adj_pvalue", "adjusted_pvalue", "padj", "FDR"),
        "zscore": ("zscore", "z_score"),
        "combined_score": ("combined_score", "combined score"),
        "genes": ("genes", "Genes", "overlapping genes"),
        "gene_count": ("gene_count", "Count", "count"),
    }
    for canonical, candidates in aliases.items():
        column = find_column(result, *candidates)
        result[canonical] = result[column] if column else pd.NA

    result["rank"] = pd.to_numeric(result["rank"], errors="coerce")
    result["pvalue"] = pd.to_numeric(result["pvalue"], errors="coerce")
    result["adjusted_pvalue"] = pd.to_numeric(
        result["adjusted_pvalue"], errors="coerce"
    )
    result["combined_score"] = pd.to_numeric(
        result["combined_score"], errors="coerce"
    )
    result["GO ID"] = result["term"].astype(str).str.findall(GO_ID_PATTERN).str.join(", ")
    result["Gene symbols"] = result["genes"].fillna("").astype(str)
    return result


def pathway_link(term: str, source: str) -> str:
    go_match = GO_ID_PATTERN.search(term)
    if go_match:
        go_id = go_match.group(0).upper()
        return f"https://www.ebi.ac.uk/QuickGO/term/{go_id}"
    reactome_match = re.search(r"\bR-[A-Z]{3}-\d+(?:\.\d+)?\b", term)
    if reactome_match:
        return f"https://reactome.org/content/detail/{reactome_match.group(0)}"
    kegg_match = re.search(r"\b(?:hsa)?\d{5}\b", term, re.IGNORECASE)
    if "KEGG" in source and kegg_match:
        pathway_id = kegg_match.group(0).lower()
        if not pathway_id.startswith("hsa"):
            pathway_id = f"hsa{pathway_id}"
        return f"https://www.kegg.jp/pathway/{pathway_id}"
    return ""


def classification_counts(frame: pd.DataFrame) -> pd.DataFrame:
    classification_column = find_column(frame, "Classification", "classification")
    gene_column = find_column(frame, "Gene Symbol", "gene_symbol")
    if not classification_column or not gene_column:
        return pd.DataFrame(columns=["Classification", "Gene count", "Label"])
    counts = (
        frame.dropna(subset=[classification_column, gene_column])
        .assign(
            _classification=lambda current: current[classification_column].astype(str).str.strip()
        )
        .groupby("_classification")[gene_column]
        .nunique()
        .rename_axis("Classification")
        .reset_index(name="Gene count")
    )
    order = {"Established": 0, "Emerging": 1, "Novel": 2}
    counts["_slice_order"] = counts["Classification"].map(order).fillna(3)
    counts = counts.sort_values(["_slice_order", "Classification"])
    total_genes = counts["Gene count"].sum()
    counts["Label position"] = (
        counts["Gene count"].cumsum() - counts["Gene count"] / 2
    )
    counts["Slice label"] = counts.apply(
        lambda row: f"{row['Classification']} genes\n({row['Gene count']})",
        axis=1,
    )
    return counts.reset_index(drop=True)


def ontology_rows_for_gene(pathways: pd.DataFrame, gene: str) -> pd.DataFrame:
    if pathways.empty:
        return pd.DataFrame()
    gene_key = normalise_symbol(gene)
    go_pathways = pathways[
        pathways["Enrichment source"].astype(str).str.startswith("GO ")
    ].copy()
    if go_pathways.empty:
        return go_pathways

    def contains_gene(value: Any) -> bool:
        return gene_key in {
            normalise_symbol(symbol)
            for symbol in re.split(r"[|,;]", str(value))
            if symbol.strip()
        }

    matched = go_pathways[go_pathways["Gene symbols"].map(contains_gene)].copy()
    matched["Pathway link"] = [
        pathway_link(str(term), str(source))
        for term, source in zip(matched["term"], matched["Enrichment source"])
    ]
    columns = [
        "Enrichment source",
        "rank",
        "term",
        "GO ID",
        "pvalue",
        "adjusted_pvalue",
        "Gene symbols",
        "Pathway link",
    ]
    return matched[[column for column in columns if column in matched]]


def reference_rows(frame: pd.DataFrame, gene: str) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    gene_column = find_column(frame, "gene_symbol", "Gene Symbol")
    selected_frame = frame
    if gene_column:
        selected_frame = frame[
            frame[gene_column].map(normalise_symbol).eq(normalise_symbol(gene))
        ]
    for column in selected_frame.columns:
        if any(token in str(column).lower() for token in ("reference", "source", "url")):
            for value in selected_frame[column].dropna().astype(str):
                for url in URL_PATTERN.findall(value):
                    output.append(
                        {"Source": str(column), "Link": url.rstrip(".)")}
                    )
                for pmid in PMID_PATTERN.findall(value):
                    output.append(
                        {
                            "Source": f"PMID {pmid}",
                            "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                        }
                    )
    # PubMed search results often contain bare semicolon-separated PMID identifiers.
    if gene_column:
        for column in selected_frame.columns:
            if "pmid" in str(column).lower():
                for value in selected_frame[column].dropna().astype(str):
                    for pmid in re.findall(r"\b\d{5,10}\b", value):
                        output.append(
                            {
                                "Source": f"PMID {pmid}",
                                "Link": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                            }
                        )
    unique: dict[str, dict[str, str]] = {}
    for item in output:
        unique[item["Link"]] = item
    return list(unique.values())


def make_gene_link(gene_id: Any) -> str:
    identifier = "" if pd.isna(gene_id) else str(gene_id).strip()
    if re.fullmatch(r"ENSG\d{11}(?:\.\d+)?", identifier):
        return f"https://www.ensembl.org/id/{identifier}"
    return ""


def filter_gene_rows(
    frame: pd.DataFrame,
    gene_column_names: tuple[str, ...],
    active_genes: set[str] | None,
) -> pd.DataFrame:
    result = add_gene_key(frame, *gene_column_names)
    if active_genes is not None:
        result = result[result["_gene_key"].isin(active_genes)]
    return result.drop(columns=["_gene_key"], errors="ignore")


def render_gene_details(
    gene: str,
    top50: pd.DataFrame,
    data: dict[str, pd.DataFrame],
) -> None:
    gene_column = find_column(top50, "Gene Symbol", "gene_symbol")
    if not gene_column:
        st.info("The Top 50 file has no gene-symbol column.")
        return
    selected = top50[
        top50[gene_column].map(normalise_symbol).eq(normalise_symbol(gene))
    ]
    if selected.empty:
        st.info(f"{gene} is in the candidate data but not the Top 50 table.")
    else:
        linked_selected = selected.copy()
        linked_selected[gene_column] = linked_selected[gene_column].map(
            lambda symbol: (
                "https://www.ncbi.nlm.nih.gov/gene/?term="
                + quote(normalise_symbol(symbol), safe="")
            )
        )
        st.dataframe(
            linked_selected,
            width="stretch",
            hide_index=True,
            column_config={
                gene_column: st.column_config.LinkColumn(
                    gene_column,
                    display_text=r".*\?term=([^&]+)",
                )
            },
        )
        row = selected.iloc[0]
        ensembl_column = find_column(selected, "Ensembl Gene ID", "ensembl_gene_id")
        ensembl_url = make_gene_link(row[ensembl_column]) if ensembl_column else ""
        if ensembl_url:
            st.markdown(f"[Open {gene} in Ensembl]({ensembl_url})")

    refs: list[dict[str, str]] = []
    for key in ("top50", "hallmark", "pubmed", "evidence"):
        if key in data:
            refs.extend(reference_rows(data[key], gene))
    if refs:
        st.markdown("**Linked references**")
        st.dataframe(
            pd.DataFrame(refs),
            width="stretch",
            hide_index=True,
            column_config={
                "Link": st.column_config.LinkColumn(
                    "Link", display_text="Open source"
                )
            },
        )


def main() -> None:
    st.logo(str(LOGO_PATH), size="large", icon_image=str(LOGO_PATH))
    st.title("Obesity Target Prioritization Platform")
    st.caption(
        "Explore the 434-gene candidate pool, prioritized targets, pathway "
        "enrichment, scores, and supporting references."
    )

    with st.sidebar:
        st.header("Data")
        st.caption(
            "Bundled dashboard data loads automatically. Upload CSV or XLSX "
            "files to replace matching data."
        )
        uploads = st.file_uploader(
            "Upload CSV files or an Excel workbook",
            type=["csv", "xlsx"],
            accept_multiple_files=True,
            help="Recognized filenames include candidate_genes_union, top50, "
            "go_bp, go_cc, go_mf, kegg, reactome, ensembl_gene_records, "
            "gene_specific_pubmed, hallmark_mapping, and score_breakdown. "
            "Excel workbooks can provide TOP_50_GENES, HALLMARK_MAPPING, "
            "EVIDENCE_TABLE, and SCORE_BREAKDOWN sheets.",
        )
        data, data_origins = load_datasets(uploads)

        st.divider()
        st.header("Filters")
        candidate_symbols: set[str] = set()
        for key, names in (
            ("candidates", ("gene_symbol", "Gene Symbol")),
            ("top50", ("Gene Symbol", "gene_symbol")),
        ):
            if key in data:
                symbol_column = find_column(data[key], *names)
                if symbol_column:
                    candidate_symbols.update(
                        symbol
                        for symbol in data[key][symbol_column].map(normalise_symbol)
                        if symbol
                    )
        gene_options = sorted(candidate_symbols)
        selected_genes = st.multiselect(
            "Gene name",
            options=gene_options,
            placeholder="Search and select genes",
            help="Limits gene tables and pathway rows to these genes.",
        )
        selected_traits = st.multiselect(
            "Trait",
            options=list(TRAIT_PATTERNS),
            placeholder="Select phenotype terms",
            help="Matches the phenotype/evidence text in the prioritized Top 50. "
            "Candidate-only genes without phenotype annotations are not inferred.",
        )

        pathway_data = get_pathway_frame(data)
        all_go_ids = sorted(
            {
                go_id.upper()
                for values in pathway_data.get("GO ID", pd.Series(dtype=str)).dropna()
                for go_id in str(values).split(", ")
                if go_id.strip()
            }
        )
        selected_go_ids = st.multiselect(
            "GO ID",
            options=all_go_ids,
            placeholder="Search and select GO IDs",
            help="Limits enrichment results to the selected GO identifiers.",
        )

    top50 = data.get("top50", pd.DataFrame())
    trait_genes = trait_gene_set(top50, selected_traits) if selected_traits else None
    explicit_genes = set(selected_genes) if selected_genes else None
    active_genes: set[str] | None = None
    if explicit_genes is not None and trait_genes is not None:
        active_genes = explicit_genes & trait_genes
    elif explicit_genes is not None:
        active_genes = explicit_genes
    elif trait_genes is not None:
        active_genes = trait_genes

    candidate_count = (
        len(candidate_symbols)
        if candidate_symbols
        else 0
    )
    filtered_top50 = (
        filter_gene_rows(top50, ("Gene Symbol", "gene_symbol"), active_genes)
        if not top50.empty
        else pd.DataFrame()
    )
    top50_gene_column = find_column(top50, "Gene Symbol", "gene_symbol")

    metric_columns = st.columns(4)
    metric_columns[0].metric(
        "Candidate genes",
        f"{candidate_count:,}" if candidate_count else "—",
    )
    metric_columns[1].metric(
        "Prioritized genes",
        f"{top50[top50_gene_column].nunique():,}"
        if top50_gene_column and not top50.empty
        else "—",
    )
    pathway_count = len(pathway_data)
    metric_columns[2].metric(
        "Enrichment terms", f"{pathway_count:,}" if pathway_count else "—"
    )
    score_column = find_column(top50, "Overall Score", "overall_score")
    if score_column and not filtered_top50.empty:
        scores = pd.to_numeric(filtered_top50[score_column], errors="coerce").dropna()
        median_score = f"{scores.median():.0f}" if not scores.empty else "—"
    else:
        median_score = "—"
    metric_columns[3].metric("Median selected-set score", median_score)

    if not data:
        st.warning(
            "No datasets were loaded. Place the supplied CSVs in their configured "
            "paths or upload them using the sidebar."
        )
        return

    overview_tab, genes_tab, pathways_tab, sources_tab = st.tabs(
        ["Overview", "Gene explorer", "Pathway enrichment", "Sources & data"]
    )

    with overview_tab:
        st.subheader("Gene evidence classification")
        classification_column = find_column(
            filtered_top50, "Classification", "classification"
        )
        gene_column = find_column(filtered_top50, "Gene Symbol", "gene_symbol")
        classification_table = classification_counts(filtered_top50)
        selected_go_gene = None
        if classification_column and gene_column and not classification_table.empty:
            chart_column, genes_column = st.columns([1, 1])
            class_selection = alt.selection_point(
                name="classification_pick",
                fields=["Classification"],
                on="click",
                clear=False,
            )
            palette = alt.Scale(
                domain=["Established", "Emerging", "Novel"],
                range=["#178f79", "#df9a31", "#7a67b5"],
            )
            base_pie = alt.Chart(classification_table).encode(
                theta=alt.Theta("Gene count:Q", stack=True),
                order=alt.Order("_slice_order:Q", sort="ascending"),
            )
            pie = (
                base_pie.mark_arc(innerRadius=55, outerRadius=145, stroke="white")
                .encode(
                    color=alt.Color(
                        "Classification:N",
                        scale=palette,
                        legend=alt.Legend(title="Evidence classification"),
                    ),
                    opacity=alt.condition(
                        class_selection,
                        alt.value(1),
                        alt.value(0.55),
                    ),
                    tooltip=[
                        alt.Tooltip("Classification:N"),
                        alt.Tooltip("Gene count:Q", title="Genes"),
                    ],
                )
                .add_params(class_selection)
            )
            labels = base_pie.mark_text(
                radius=98,
                color="white",
                fontSize=13,
                fontWeight="bold",
                lineBreak="\n",
            ).encode(
                theta=alt.Theta(
                    "Label position:Q",
                    scale=alt.Scale(
                        domain=[
                            0,
                            float(classification_table["Gene count"].sum()),
                        ],
                        range=[0, 2 * math.pi],
                    ),
                ),
                text=alt.Text("Slice label:N"),
                tooltip=[
                    alt.Tooltip("Classification:N"),
                    alt.Tooltip("Gene count:Q", title="Genes"),
                ],
            )
            pie_chart = (pie + labels).properties(
                height=350,
                title="Top 50 genes by evidence classification",
            )
            with chart_column:
                pie_event = st.altair_chart(
                    pie_chart,
                    width="stretch",
                    on_select="rerun",
                    selection_mode="classification_pick",
                    key="classification_pie",
                )
            selected_classes = (
                pie_event.selection.get("classification_pick", [])
                if pie_event and pie_event.selection
                else []
            )
            chosen_class = (
                selected_classes[0].get("Classification")
                if selected_classes
                else None
            )
            genes_in_class = filtered_top50.copy()
            if chosen_class:
                genes_in_class = genes_in_class[
                    genes_in_class[classification_column]
                    .astype(str)
                    .eq(str(chosen_class))
                ]
            rank_column_for_gene = find_column(genes_in_class, "Rank", "rank")
            score_column_for_gene = find_column(
                genes_in_class, "Overall Score", "overall_score"
            )
            if rank_column_for_gene:
                genes_in_class[rank_column_for_gene] = pd.to_numeric(
                    genes_in_class[rank_column_for_gene], errors="coerce"
                )
                genes_in_class = genes_in_class.sort_values(rank_column_for_gene)
            gene_table_columns = [
                column
                for column in (
                    rank_column_for_gene,
                    gene_column,
                    classification_column,
                    score_column_for_gene,
                )
                if column is not None
            ]
            with genes_column:
                st.markdown(
                    f"**{chosen_class or 'All classifications'} genes** — "
                    "click a row to show that gene's GO ontology terms."
                )
                gene_selection = st.dataframe(
                    genes_in_class[gene_table_columns].reset_index(drop=True),
                    width="stretch",
                    hide_index=True,
                    on_select="rerun",
                    selection_mode="single-row",
                    key="classification_gene_table",
                )
            selected_rows = gene_selection.selection.rows
            if selected_rows:
                selected_row = genes_in_class.reset_index(drop=True).iloc[
                    selected_rows[0]
                ]
                selected_go_gene = normalise_symbol(selected_row[gene_column])

            if selected_go_gene:
                st.markdown(f"#### GO ontology results for {selected_go_gene}")
                gene_ontology = ontology_rows_for_gene(
                    pathway_data, selected_go_gene
                )
                if gene_ontology.empty:
                    st.info(
                        f"No GO Biological Process, Cellular Component, or "
                        f"Molecular Function rows in the supplied enrichment "
                        f"tables include {selected_go_gene}."
                    )
                else:
                    st.dataframe(
                        gene_ontology,
                        width="stretch",
                        hide_index=True,
                        column_config={
                            "Pathway link": st.column_config.LinkColumn(
                                "GO ontology",
                                display_text="Open in QuickGO",
                            ),
                            "pvalue": st.column_config.NumberColumn(
                                "P-value", format="%.3g"
                            ),
                            "adjusted_pvalue": st.column_config.NumberColumn(
                                "Adjusted p-value", format="%.3g"
                            ),
                        },
                    )
                st.caption(
                    "These rows show which GO terms in the supplied set-level "
                    "enrichment results contain the selected gene; they are not "
                    "independent gene-specific enrichment tests."
                )
            else:
                st.info(
                    "Click a classification slice, then click a gene row to "
                    "display its associated GO terms."
                )
        else:
            st.info(
                "Load a workbook Top 50 sheet with Gene Symbol and "
                "Classification columns to display the classification pie chart."
            )

        st.subheader("Prioritized gene scores")
        if score_column and not filtered_top50.empty:
            rank_column = find_column(filtered_top50, "Rank", "rank")
            gene_column = find_column(filtered_top50, "Gene Symbol", "gene_symbol")
            chart = filtered_top50.copy()
            chart[score_column] = pd.to_numeric(chart[score_column], errors="coerce")
            if rank_column and gene_column:
                chart[rank_column] = pd.to_numeric(chart[rank_column], errors="coerce")
                chart = chart.sort_values(rank_column)
                chart["Rank · Gene"] = (
                    chart[rank_column].astype("Int64").astype(str)
                    + " · "
                    + chart[gene_column].map(normalise_symbol)
                )
                chart = chart.set_index("Rank · Gene")[[score_column]].dropna()
            elif rank_column:
                chart[rank_column] = pd.to_numeric(chart[rank_column], errors="coerce")
                chart = chart.sort_values(rank_column)
                chart = chart.set_index(rank_column)[[score_column]].dropna()
                chart.index.name = "Rank"
            elif gene_column:
                chart = chart.set_index(gene_column)[[score_column]].dropna()
                chart.index.name = "Gene"
            else:
                chart = chart[[score_column]].dropna()
            st.line_chart(chart, y_label="Overall score", color="#7c3aed")
            st.caption(
                "Scores use the supplied Top 50 CSV. This is a ranked comparison, "
                "not a time series."
            )
        else:
            st.info("Load a Top 50 CSV with an Overall Score column to show this chart.")

        st.subheader("Highest-ranked genes in the active filter")
        if not filtered_top50.empty:
            rank_column = find_column(filtered_top50, "Rank", "rank")
            if rank_column:
                filtered_top50 = filtered_top50.sort_values(
                    rank_column, key=lambda values: pd.to_numeric(values, errors="coerce")
                )
            st.dataframe(filtered_top50.head(15), width="stretch", hide_index=True)
        else:
            st.info("No prioritized genes match the selected filters.")

    with genes_tab:
        st.subheader("Candidate pool")
        if "candidates" in data:
            candidate_table = filter_gene_rows(
                data["candidates"], ("gene_symbol", "Gene Symbol"), active_genes
            )
            st.caption(
                f"Showing {len(candidate_table):,} of "
                f"{len(data['candidates']):,} candidate rows."
            )
            st.dataframe(candidate_table, width="stretch", hide_index=True)
        else:
            st.info("Upload the candidate gene union CSV to browse the full pool.")

        st.subheader("Gene evidence and references")
        if candidate_symbols:
            detail_options = sorted(
                active_genes if active_genes is not None else candidate_symbols
            )
            if detail_options:
                selected_detail_gene = st.selectbox(
                    "Select a gene for details",
                    options=detail_options,
                )
                render_gene_details(selected_detail_gene, top50, data)
            else:
                st.info("No genes match the current sidebar filters.")

        for key, title, gene_cols in (
            ("scores", "Score breakdown", ("Gene Symbol", "gene_symbol")),
            ("hallmark", "Hallmark mapping", ("Gene Symbol", "gene_symbol")),
            ("evidence", "Evidence table", ("Gene Symbol", "gene_symbol")),
            ("ensembl", "Ensembl gene records", ("gene_symbol", "Gene Symbol")),
            ("pubmed", "PubMed search summary", ("gene_symbol", "Gene Symbol")),
        ):
            if key in data:
                with st.expander(title):
                    table = filter_gene_rows(data[key], gene_cols, active_genes)
                    url_column = find_column(table, "URL", "Source", "Supporting Reference")
                    column_config = (
                        {
                            url_column: st.column_config.LinkColumn(
                                url_column,
                                display_text="Open source",
                            )
                        }
                        if url_column
                        else None
                    )
                    st.dataframe(
                        table,
                        width="stretch",
                        hide_index=True,
                        column_config=column_config,
                    )

    with pathways_tab:
        st.subheader("Enrichment results")
        if pathway_data.empty:
            st.info("Upload at least one pathway enrichment CSV.")
        else:
            filtered_pathways = pathway_data.copy()
            if active_genes is not None:
                wanted = {normalise_symbol(gene) for gene in active_genes}

                def overlaps_genes(value: Any) -> bool:
                    genes = {
                        normalise_symbol(item)
                        for item in re.split(r"[|,;]", str(value))
                        if item.strip()
                    }
                    return bool(genes & wanted)

                filtered_pathways = filtered_pathways[
                    filtered_pathways["Gene symbols"].map(overlaps_genes)
                ]
            if selected_go_ids:
                selected_set = {value.upper() for value in selected_go_ids}
                filtered_pathways = filtered_pathways[
                    filtered_pathways["GO ID"].map(
                        lambda value: bool(
                            selected_set
                            & {item.upper() for item in str(value).split(", ")}
                        )
                    )
                ]

            source_filter = st.multiselect(
                "Pathway database",
                options=sorted(filtered_pathways["Enrichment source"].dropna().unique()),
                placeholder="All loaded databases",
            )
            if source_filter:
                filtered_pathways = filtered_pathways[
                    filtered_pathways["Enrichment source"].isin(source_filter)
                ]
            st.caption(
                f"{len(filtered_pathways):,} pathway rows match the active filters. "
                "Gene and trait filters are applied by overlapping genes from the "
                "supplied Top 50 evidence; GO ID filtering applies to GO terms."
            )

            adjusted = filtered_pathways.dropna(subset=["adjusted_pvalue"]).copy()
            if not adjusted.empty:
                adjusted = adjusted[adjusted["adjusted_pvalue"] > 0]
                adjusted["-log10 adjusted p"] = adjusted["adjusted_pvalue"].map(
                    lambda value: -math.log10(value)
                )
                adjusted = adjusted.sort_values(["rank", "-log10 adjusted p"])
                chart = adjusted.head(25).copy()
                chart["Chart term"] = (
                    chart["Enrichment source"].astype(str)
                    + " — "
                    + chart["term"].astype(str).str.slice(0, 70)
                )
                chart = chart.set_index("Chart term")[["-log10 adjusted p"]]
                st.line_chart(chart, y_label="-log10 adjusted p-value", color="#0891b2")
                st.caption(
                    "Top 25 terms by source rank. Lower adjusted p-values appear higher."
                )

            display_columns = [
                "Enrichment source",
                "rank",
                "term",
                "GO ID",
                "pvalue",
                "adjusted_pvalue",
                "zscore",
                "combined_score",
                "gene_count",
                "Gene symbols",
            ]
            display = filtered_pathways[
                [column for column in display_columns if column in filtered_pathways]
            ].copy()
            display["Pathway link"] = [
                pathway_link(str(term), str(source))
                for term, source in zip(
                    filtered_pathways["term"],
                    filtered_pathways["Enrichment source"],
                )
            ]
            st.dataframe(
                display,
                width="stretch",
                hide_index=True,
                column_config={
                    "Pathway link": st.column_config.LinkColumn(
                        "Pathway link", display_text="Open source"
                    ),
                    "adjusted_pvalue": st.column_config.NumberColumn(
                        "Adjusted p-value", format="%.3g"
                    ),
                    "pvalue": st.column_config.NumberColumn(
                        "P-value", format="%.3g"
                    ),
                },
            )

    with sources_tab:
        st.subheader("Loaded files")
        loaded_rows = [
            {
                "Dataset": DATASET_LABELS.get(key, key),
                "Rows": len(frame),
                "Columns": len(frame.columns),
                "Source": data_origins.get(key, "Uploaded file"),
            }
            for key, frame in data.items()
        ]
        st.dataframe(
            pd.DataFrame(loaded_rows),
            width="stretch",
            hide_index=True,
        )
        st.markdown(
            "- Gene identifiers link to [Ensembl](https://www.ensembl.org/).\n"
            "- PubMed references link to [PubMed](https://pubmed.ncbi.nlm.nih.gov/).\n"
            "- GO terms link to [QuickGO](https://www.ebi.ac.uk/QuickGO/).\n"
            "- Reactome and KEGG identifiers link to their respective pathway records."
        )
        st.caption(
            "Input CSV and Excel contents are parsed locally by this app. No "
            "external API calls are made; external links open only when clicked."
        )


main()
