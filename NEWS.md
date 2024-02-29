---
output:
  html_document: default
  pdf_document: default
---
# MSoutcomes 0.2.0

* New functions

  * Added the `CDW()` function.
  * Added the `CDI()` function.
  * Added the `CDEseq()` function.
  * Added the `PIRA()` function.
  
  
* Changes

  * Replaced SampleData with a new dataset.
  * Changes in `spmsDx()` and `spmsDx_no_fss()` functions
    * If a baseline EDSS score and corresponding FSS are not provided, these are now determined as the first EDSS score and corresponding FSS recorded in the dataset, outside 30 days (the default) of a relapse. Previously, the first EDSS score and corresponding FSS recorded in the dataset were used as the baseline scores, regardless of the recency of a relapse.
    
  
* Minor bug fixes and improvements

  * `spmsDx()` function: `ungroup` correctly placed (122).
  * `spmsDx()` function: sustained worsening in EDSS is now calculated including EDSS scores recorded beyond 30 days (the default) following a relapse (148).
  * `spmsDx_no_fss()` function: `ungroup` correctly placed (113).
  * `spmsDx_no_fss()` function: sustained worsening in EDSS is now calculated including EDSS scores recorded beyond 30 days (the default) following a relapse (131).
