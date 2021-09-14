
DimRed <- setClass('DimRed', 
                   slots = c(
                     met_mat = 'matrix',
                     methylation_type = 'character',
                     vmr_count = 'numeric', 
                     max_ratio_of_na_cells = 'numeric',
                     pca_coords = 'matrix',
                     umap_coords = 'matrix'
                     
                   ) 
                   
                   )