library(xtable)

# ---------------- Mom hiv -------------
momhiv_table = readr::read_csv("supptable_momhiv.csv")
print(xtable(momhiv_table), include.rownames = F)


# ---------------- Correlation analysis -----------------
corr_table = readr::read_csv("supptable_correlations.csv")
print(xtable(corr_table), include.rownames = F)

# ---------------- Cubic regression table ---------------
cubic_table = readr::read_csv("supptable_cubic_regression.csv")
print(xtable(cubic_table), include.rownames = F)
