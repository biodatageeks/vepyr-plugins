# PhenotypeOrthologous

Source: Ensembl release-116 FTP,
`variation/PhenotypeOrthologous/PhenotypesOrthologous_homo_sapiens_112_GRCh38.gff3.gz`
(BGZF + `.tbi`, 3.4 MB, 16,409 gene features on chr1–22, X, MT and four
unplaced contigs; no chrY). GRCh38 only, as the upstream plugin enforces.

The file is used as published — no preprocessing — so `verify_source="strict"`
applies. Ensembl's `CHECKSUMS` there is BSD `sum` output; the manifest `md5`
was computed on the downloaded file.

Match rule (from `PhenotypeOrthologous.pm`): tabix query by the variant's span,
first record whose `gene_id` equals the transcript's gene stable id. Fields
are emitted alphabetically (VEP sorts plugin header keys). RefSeq transcripts
in a merged cache never match (their gene id is not an ENSG id). Licence:
Ensembl data, free for any use (see the Ensembl FTP README).
