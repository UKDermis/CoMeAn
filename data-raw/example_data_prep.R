
# healthy data
CTRL_example = as.matrix(read.csv("inst/extdata/CTRL_example.tsv", sep="\t"))
usethis::use_data(CTRL_example, overwrite = TRUE)

# atopic dermatitis data
AD_example = as.matrix(read.csv("inst/extdata/AD_example.tsv", sep="\t"))
usethis::use_data(AD_example, overwrite = TRUE)

# psoriasis data
PSO_example = as.matrix(read.csv("inst/extdata/PSO_example.tsv", sep="\t"))
usethis::use_data(PSO_example, overwrite = TRUE)

# gene annotation data
SkinSig_annotation = read.csv("inst/extdata/SkinSig_annotation.txt", sep="\t")
usethis::use_data(SkinSig_annotation, overwrite = TRUE)

#---

#roxygen2::roxygenise(clean = T)
