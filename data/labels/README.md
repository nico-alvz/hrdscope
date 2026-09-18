# data/labels

`tcga_ov_hrd.tsv` is produced by `hrdscope labels` and is not committed (regenerate it, or download the released copy from Zenodo).

| column | meaning |
|---|---|
| patient | TCGA barcode (12 chars) |
| n_ascat_profiles | number of ASCAT3 tumour profiles averaged |
| hrd_loh, ntai, lst | scar components (Abkevich 2012, Birkbak 2012, Popova 2012) |
| hrd_sum | hrd_loh + ntai + lst |
| hrd_ge33, hrd_ge42, hrd_ge63 | binary labels at the thresholds used in the literature |
| knijnenburg2018_hrd | published HRD score when available (validation only) |
| n_dx_slides, n_ts_slides | diagnostic (FFPE) and tissue (frozen) slides available in GDC |
