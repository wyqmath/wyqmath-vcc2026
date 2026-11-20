# Final count overrides

The final recipe replaced 15 of 900 count blocks. Indices are zero-based. Each block has 400 prediction cells on the full gene axis.

| Context / target | Index | Source and effect input | Evidence or scope |
| --- | ---: | --- | --- |
| D / ANP32B | 11 | [GSE242849](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE242849), MV4-11, relative FPKM | Deposited NC/S3 contrast; viral-control equivalence is unresolved |
| D / CARM1 | 35 | [GSE103528](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103528), MOLM13, RSEM expected reads | Day 3 contrast; both hairpins retained |
| D / PHF6 | 160 | [GSE277149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277149), THP-1, featureCounts | Clone/guide assignments for individual repeats are unresolved |
| D / RBM5 | 183 | [GSE225633](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225633), MOLM13, RSEM estimates | Three replicate-label holdouts |
| D / SMARCD1 | 209 | [GSE200201](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE200201), MOLM13, single-cell gene-expression counts | The final effect uses channels 1 to 10; channels 11 to 15 supplied comparison data |
| E / PHF23 | 459 | [GSE300095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE300095), basal HeLa, featureCounts | KO/WT contrast; RNA clone independence and WT editing equivalence are unresolved |
| E / PRMT6 | 470 | [GSE221971](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221971), U87, deposited gene counts | Three library holdouts; count columns used |
| E / VIM | 560 | [GSE286944](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE286944), U251-MG, featureCounts | Nine correlated holdout views |
| F / CUX1 | 649 | [GSE129487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129487), four RA synovial fibroblast donors, kallisto estimated counts | Refines supported coordinates of the earlier [GSE154148](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE154148) author-effect block |
| F / E2F1 | 665 | [GSE77260](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77260), Wharton's jelly mesenchymal cells, Partek estimates | Six whole-donor holdouts |
| F / GATAD2B | 680 | [GSE154148-associated workbook](https://doi.org/10.1038/s41556-021-00725-7), author-corrected differential effect | Consistency was checked against the supplied author effect |
| F / HMGB2 | 691 | [GSE269856](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269856), proliferating IMR90, gene counts | Library holdouts; WT transfection identity is unresolved |
| F / PHF14 | 758 | [GSE131576](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131576), proliferating WI38, HTSeq counts | Paired library holdouts |
| F / PRDM15 | 768 | [GSE154148-associated workbook](https://doi.org/10.1038/s41556-021-00725-7), author-corrected differential effect | Consistency was checked against the supplied author effect |
| F / PTBP2 | 773 | [GSE210131](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210131), AG22056 conversion day 3, featureCounts | AMnp/AMp conditional comparison; viral-burden balance is unresolved |
