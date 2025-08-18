# Clear the current R environment and run garbage collection
rm(list = ls())
gc(full = TRUE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output/output_", timestamp)
outputdir <- paste0(getwd(),"/",output)
output_base_dir <- paste0("../../", output, "/")
dir.create(outputdir, showWarnings = FALSE)
print(getwd())
print(output)
# copy the config file into the output directory
config_path <- file.path(outputdir, "config.json")
file.copy(file.path(getwd(),"sc_pipeline/src/config.json"), config_path)

cat(paste0("Starting pipeline at ", timestamp, " with output directory ", output, "\n"))

rmarkdown::render(
    input = file.path(getwd(), "sc_pipeline", "src", "sc_pipeline.rmd"),
    output_format = "pdf_document",
    output_dir = output,
    intermediates_dir = output,
    params = list(
        output = output,
        config_path = config_path,
        output_base_dir = output_base_dir
    )
)
