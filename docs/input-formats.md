# Input formats

## Species tree
Newick with unique terminal labels:

```text
((Genome_1,Genome_2),(Genome_3,Genome_4));
```

## Gene tree
Newick with unique sequence/copy labels:

```text
((seqA,seqB),(seqC,seqD));
```

Internal numeric labels are read as branch/node support by Biopython when represented in standard Newick form.

## Mapping table

### Two-column form
Useful for one family:

```text
gene_id	species_id
seqA	Genome_1
seqB	Genome_2
```

### Three-column form
Useful for batch mode:

```text
family	gene_id	species_id
OG000001	seqA	Genome_1
OG000001	seqB	Genome_2
OG000002	seqC	Genome_1
```

Headers are optional, but recommended.

## Delimiter-based mapping
If labels follow a guaranteed convention such as `Genome_1|copyA`, you may use:

```bash
--delimiter '|' --species-field 0
```

Explicit mapping is preferred for publication workflows because it is auditable and avoids silent parsing assumptions.
