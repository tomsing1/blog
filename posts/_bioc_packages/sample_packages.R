# retrieve a data.frame with all packages from a Bioconductor release
kVersion <- "3.18"
kDcfUrl <- "https://bioconductor.org/packages/%s/bioc/src/contrib/PACKAGES"
kPackageUrl <- "https://bioconductor.org/packages/%s/bioc/html/%s.html"

# read the full list of packages into a data.frame
df <- read.dcf(url(sprintf(kDcfUrlkVersion))) |> as.data.frame()

# shuffle the data.frame
set.seed(1234)
packages <- df[sample(seq.int(nrow(df)), nrow(df), replace = FALSE), "Package"]
head(packages)

# open the web site for the nth package
n = 1
if (isTRUE(interactive())) {
  browseURL(sprintf(kPackageUrl, kVersion, packages[n]))
}
