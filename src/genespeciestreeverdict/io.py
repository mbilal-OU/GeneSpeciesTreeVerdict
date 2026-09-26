from __future__ import annotations

import csv
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


TREE_SUFFIXES = {".nwk", ".newick", ".tree", ".treefile", ".tre"}


@dataclass(frozen=True, slots=True)
class MappingRecord:
    gene_id: str
    species_id: str
    family: str | None = None


def read_tree(path: str | Path) -> Tree:
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Tree file not found: {path}")
    try:
        tree = Phylo.read(str(path), "newick")
    except Exception as exc:  # Biopython raises parser-specific exceptions
        raise ValueError(f"Could not parse Newick tree '{path}': {exc}") from exc
    names = [terminal.name for terminal in tree.get_terminals()]
    if any(name is None or str(name).strip() == "" for name in names):
        raise ValueError(f"Tree contains an unnamed terminal: {path}")
    if len(names) != len(set(names)):
        duplicates = sorted(name for name, n in Counter(names).items() if n > 1)
        raise ValueError(
            f"Tree contains duplicate terminal labels: {', '.join(duplicates)}. "
            "Gene copies must have unique sequence IDs."
        )
    return tree


def terminal_names(tree: Tree) -> list[str]:
    return [str(node.name) for node in tree.get_terminals()]


def read_mapping_table(path: str | Path) -> list[MappingRecord]:
    """Read two- or three-column TSV mapping data.

    Supported forms:
      gene_id <tab> species_id
      family  <tab> gene_id <tab> species_id

    A header row using these field names is optional.
    """

    path = Path(path)
    rows: list[list[str]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or not any(cell.strip() for cell in row):
                continue
            if row[0].lstrip().startswith("#"):
                continue
            rows.append([cell.strip() for cell in row])

    if not rows:
        raise ValueError(f"Mapping file is empty: {path}")

    header = [cell.lower() for cell in rows[0]]
    has_header = "gene_id" in header and "species_id" in header
    records: list[MappingRecord] = []

    if has_header:
        idx_gene = header.index("gene_id")
        idx_species = header.index("species_id")
        idx_family = header.index("family") if "family" in header else None
        for row in rows[1:]:
            required_idx = max(idx_gene, idx_species, idx_family or 0)
            if len(row) <= required_idx:
                raise ValueError(f"Malformed mapping row in {path}: {row}")
            records.append(
                MappingRecord(
                    gene_id=row[idx_gene],
                    species_id=row[idx_species],
                    family=row[idx_family] if idx_family is not None else None,
                )
            )
    else:
        width = len(rows[0])
        if width not in {2, 3}:
            raise ValueError(
                "Mapping table must contain 2 columns (gene_id, species_id) or "
                "3 columns (family, gene_id, species_id)."
            )
        for row in rows:
            if len(row) != width:
                raise ValueError(f"Inconsistent mapping row width in {path}: {row}")
            if width == 2:
                records.append(MappingRecord(gene_id=row[0], species_id=row[1]))
            else:
                records.append(MappingRecord(family=row[0], gene_id=row[1], species_id=row[2]))

    genes = [record.gene_id for record in records]
    duplicate_gene_ids = sorted(gene for gene, n in Counter(genes).items() if n > 1)
    if duplicate_gene_ids and all(record.family is None for record in records):
        raise ValueError(
            "A two-column mapping file contains duplicate gene IDs: "
            + ", ".join(duplicate_gene_ids[:10])
        )
    return records


def mapping_for_family(
    gene_tree: Tree,
    species_tree: Tree,
    records: list[MappingRecord] | None = None,
    family: str | None = None,
    delimiter: str | None = None,
    species_field: int = 0,
) -> dict[str, str]:
    gene_ids = terminal_names(gene_tree)
    species_ids = set(terminal_names(species_tree))

    if records is not None:
        relevant = [
            record
            for record in records
            if record.family is None or family is None or record.family == family
        ]
        mapping: dict[str, str] = {}
        for record in relevant:
            if record.gene_id in mapping and mapping[record.gene_id] != record.species_id:
                raise ValueError(
                    f"Gene '{record.gene_id}' maps to more than one species in the mapping table."
                )
            mapping[record.gene_id] = record.species_id
        missing = sorted(set(gene_ids) - set(mapping))
        if missing:
            raise ValueError(
                f"{len(missing)} gene-tree leaves are missing from the mapping table: "
                + ", ".join(missing[:10])
            )
    elif set(gene_ids).issubset(species_ids):
        mapping = {gene_id: gene_id for gene_id in gene_ids}
    elif delimiter is not None:
        mapping = {}
        for gene_id in gene_ids:
            parts = gene_id.split(delimiter)
            try:
                mapping[gene_id] = parts[species_field]
            except IndexError as exc:
                raise ValueError(
                    f"Cannot infer species from '{gene_id}' using delimiter={delimiter!r} "
                    f"and species_field={species_field}."
                ) from exc
    else:
        raise ValueError(
            "Gene labels do not match species-tree labels. Provide --mapping or an explicit "
            "--delimiter/--species-field rule."
        )

    unknown_species = sorted(set(mapping.values()) - species_ids)
    if unknown_species:
        raise ValueError(
            "Mapped species absent from the species tree: " + ", ".join(unknown_species[:10])
        )
    return {gene_id: mapping[gene_id] for gene_id in gene_ids}


def family_name_from_path(path: str | Path) -> str:
    name = Path(path).name
    for suffix in sorted(TREE_SUFFIXES, key=len, reverse=True):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(path).stem


def discover_gene_trees(directory: str | Path, recursive: bool = False) -> list[Path]:
    directory = Path(directory)
    if not directory.is_dir():
        raise NotADirectoryError(f"Gene-tree directory not found: {directory}")
    iterator = directory.rglob("*") if recursive else directory.glob("*")
    files = [path for path in iterator if path.is_file() and path.suffix.lower() in TREE_SUFFIXES]
    return sorted(files)
