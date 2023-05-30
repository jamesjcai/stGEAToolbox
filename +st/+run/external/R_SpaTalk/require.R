if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools", repo="http://cran.rstudio.com/")
}
if (!requireNamespace("NNLM", quietly = TRUE)) {
    devtools::install_github("linxihui/NNLM")
}
if (!requireNamespace("SpaTalk", quietly = TRUE)) {
    devtools::install_github("ZJUFanLab/SpaTalk")
}