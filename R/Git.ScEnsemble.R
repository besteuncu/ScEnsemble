usethis::use_git()
usethis::use_github()
usethis::use_readme_rmd()
devtools::document()
usethis::use_gpl3_license()
usethis::use_roxygen_md()
usethis::use_data_raw("df.pollen")

usethis::use_testthat()
usethis::use_test("run_individual_algorithms")

BiocCheck::BiocCheck()
devtools::check()

# ScEnsemble Proje Kurulum Scripti


# Gerekli paketleri yükleyin
#install.packages(c("devtools", "usethis", "roxygen2", "testthat", "pkgdown"))

# Proje dizinine gidin
# setwd("path/to/ScEnsemble")

# 1. Temel klasör yapısını oluşturun
usethis::use_r("run_individual_algorithms")
usethis::use_r("calculate_all_validation_indices") 
usethis::use_r("generate_all_hypergraphs")
usethis::use_r("ensemble_clustering_algorithms")


# 2. Test yapısını kurun
usethis::use_testthat()
usethis::use_test("ensemble-methods")
usethis::use_test("clustering-algorithms")
usethis::use_test("validation-metrics")

# 3. Vignette oluşturun
usethis::use_vignette("ScEnsemble")
usethis::use_vignette("ScEnsemble-advanced")

# 4. GitHub yapılandırması
usethis::use_github_actions("R-CMD-check")
usethis::use_github_actions("pkgdown")
usethis::use_github_actions("test-coverage")

# 5. Paket dokümantasyonu
usethis::use_pkgdown()

# 6. News dosyası oluşturun
usethis::use_news_md()

# 7. .Rbuildignore'u güncelleyin
usethis::use_build_ignore(c("^.*\\.Rproj$", "^\\.Rproj\\.user$", "^data-raw$"))

# 8. Git yapılandırması
usethis::use_git_ignore(c("*.Rproj", ".Rproj.user", ".DS_Store", "*.RData", "*.Rhistory"))


print("✅ Proje yapısı başarıyla oluşturuldu!")
print("📝 Sonraki adım: S4 sınıflarını tanımlama")


pollen <- data.frame(df.pollen)
ann <- data.frame(numeric_labels_pann)

save(pollen, ann, file = "pollen.rda")


library(usethis)
use_git_ignore(c("*.Rproj", ".Rproj.user/", "README.Rmd"))
