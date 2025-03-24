# Tavarez_2025
Code For Figures in Tavarez 2025 Publication
Below are the packages used when to create figures, as time progressed I updated the packages and R studio itself. I recommend using the most recent version of all packages and Rstudio.
If a line of code has an error, you can look it up or email me at jrtavarez@albany.edu with the line that is not functioning and I will update the code.
# R version 4.4.1 (2024-06-14 ucrt) -- "Race for Your Life"
# Copyright (C) 2024 The R Foundation for Statistical Computing
# Platform: x86_64-w64-mingw32/x64
# 
# R is free software and comes with ABSOLUTELY NO WARRANTY.
# You are welcome to redistribute it under certain conditions.
# Type 'license()' or 'licence()' for distribution details.
# 
# Natural language support but running in an English locale
# 
# R is a collaborative project with many contributors.
# Type 'contributors()' for more information and
# 'citation()' on how to cite R or R packages in publications.
# 
# Type 'demo()' for some demos, 'help()' for on-line help, or
# 'help.start()' for an HTML browser interface to help.
# Type 'q()' to quit R.
# 
# 
# > sessionInfo()
# R version 4.4.1 (2024-06-14 ucrt)
# Platform: x86_64-w64-mingw32/x64
# Running under: Windows 10 x64 (build 19044)
# 
# Matrix products: default
# 
# 
# locale:
#   [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8    LC_MONETARY=English_United States.utf8
# [4] LC_NUMERIC=C                           LC_TIME=English_United States.utf8    
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      rstudioapi_0.16.0       jsonlite_1.8.8          magrittr_2.0.3          farver_2.1.2           
# [6] fs_1.6.4                zlibbioc_1.50.0         vctrs_0.6.5             memoise_2.0.1           ggtree_3.12.0          
# [11] htmltools_0.5.8.1       usethis_3.0.0           gridGraphics_0.5-1      htmlwidgets_1.6.4       plyr_1.8.9             
# [16] httr2_1.0.4             cachem_1.1.0            igraph_2.0.3            mime_0.12               lifecycle_1.0.4        
# [21] pkgconfig_2.0.3         gson_0.1.0              Matrix_1.7-0            R6_2.5.1                fastmap_1.2.0          
# [26] GenomeInfoDbData_1.2.12 shiny_1.9.1             digest_0.6.37           aplot_0.2.3             enrichplot_1.24.4      
# [31] colorspace_2.1-1        patchwork_1.3.0         AnnotationDbi_1.66.0    S4Vectors_0.42.1        pkgload_1.4.0          
# [36] RSQLite_2.3.7           fansi_1.0.6             httr_1.4.7              polyclip_1.10-7         compiler_4.4.1         
# [41] remotes_2.5.0           bit64_4.0.5             withr_3.0.1             BiocParallel_1.38.0     viridis_0.6.5          
# [46] DBI_1.2.3               pkgbuild_1.4.4          ggforce_0.4.2           R.utils_2.12.3          MASS_7.3-61            
# [51] rappdirs_0.3.3          sessioninfo_1.2.2       tools_4.4.1             scatterpie_0.2.4        ape_5.8                
# [56] httpuv_1.6.15           R.oo_1.26.0             glue_1.7.0              nlme_3.1-166            GOSemSim_2.30.2        
# [61] promises_1.3.0          shadowtext_0.1.4        grid_4.4.1              reshape2_1.4.4          fgsea_1.30.0           
# [66] generics_0.1.3          gtable_0.3.5            R.methodsS3_1.8.2       tidyr_1.3.1             data.table_1.16.0      
# [71] tidygraph_1.3.1         utf8_1.2.4              XVector_0.44.0          BiocGenerics_0.50.0     ggrepel_0.9.6          
# [76] pillar_1.9.0            stringr_1.5.1           yulab.utils_0.1.7       later_1.3.2             splines_4.4.1          
# [81] dplyr_1.1.4             tweenr_2.0.3            treeio_1.28.0           lattice_0.22-6          bit_4.0.5              
# [86] tidyselect_1.2.1        GO.db_3.19.1            Biostrings_2.72.1       miniUI_0.1.1.1          gridExtra_2.3          
# [91] IRanges_2.38.1          stats4_4.4.1            graphlayouts_1.1.1      Biobase_2.64.0          devtools_2.4.5         
# [96] stringi_1.8.4           UCSC.utils_1.0.0        lazyeval_0.2.2          ggfun_0.1.6             pacman_0.5.1           
# [101] codetools_0.2-20        ggraph_2.2.1            tibble_3.2.1            qvalue_2.36.0           BiocManager_1.30.25    
# [106] ggplotify_0.1.2         cli_3.6.3               xtable_1.8-4            munsell_0.5.1           Rcpp_1.0.13            
# [111] GenomeInfoDb_1.40.1     png_0.1-8               parallel_4.4.1          ellipsis_0.3.2          ggplot2_3.5.1          
# [116] blob_1.2.4              clusterProfiler_4.12.6  profvis_0.3.8           DOSE_3.30.5             urlchecker_1.0.1       
# [121] viridisLite_0.4.2       tidytree_0.4.6          scales_1.3.0            purrr_1.0.2             crayon_1.5.3           
# [126] rlang_1.1.4             cowplot_1.1.3           fastmatch_1.1-4         KEGGREST_1.44.1      
###########
