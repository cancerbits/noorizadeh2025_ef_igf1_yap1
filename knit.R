# this top-level script may be used to run (knit) the individual
# Rmd files

# set up output path for the reports
config <- yaml::read_yaml("config.yaml")
rmd_dir <- file.path(config$project_root, "Rmd")
report_dir <- file.path(rmd_dir, "knit_html")

for(rmd in list.files(rmd_dir, pattern = ".Rmd$", full.names = TRUE)) {
 rmarkdown::render(input = rmd,
                   output_dir = report_dir,
                   knit_root_dir = config$project_root,
                   envir = new.env())
}