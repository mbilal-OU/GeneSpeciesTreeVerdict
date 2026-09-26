# FAQ

## My gene tree differs from my species tree. Which one is correct?
The difference alone cannot answer that. The species tree is a hypothesis and the gene tree is a locus-specific history estimate. Investigate the process and uncertainty that could produce the conflict.

## Does one sequence per species guarantee orthology?
No. Ancient duplication followed by differential loss or inconsistent copy selection can leave one sampled copy per species while mixing paralogous lineages.

## Should I remove every discordant gene from concatenation?
No automatic rule is scientifically justified from discordance alone. Consider copy number, orthology evidence, support, HGT/ILS/recombination, alignment quality and the downstream model.

## Can I include paralogs in a species-tree analysis?
Not as if they were ordinary single-copy orthologs. Some methods explicitly model or accept multi-copy gene-family trees, for example ASTRAL-Pro/ASTER and reconciliation-based frameworks.

## Is a core gene automatically a single-copy ortholog?
No. Core status describes prevalence; copy number and orthology are separate properties.

## Why does v0.1 use PASS/REVIEW/RESOLVE instead of KEEP/EXCLUDE?
Because a locus can be inappropriate for conventional concatenation yet still be biologically informative or usable with another method.

## Can I use the tool on bacteria?
Yes, but HGT and homologous recombination can be particularly important. Treat DL-only reconciliation as a diagnostic, not a complete microbial evolutionary model.
