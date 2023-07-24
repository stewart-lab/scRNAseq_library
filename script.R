timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output/output_", timestamp)
dir.create(output, showWarnings = FALSE)

# copy the config file into the output directory
config_path <- file.path(getwd(), output, "config.json")
file.copy(file.path(getwd(), "config.json"), config_path)

rmarkdown::render(
    input = file.path(getwd(), "src", "script.rmd"),
    output_format = "pdf_document",
    output_dir = output,
    intermediates_dir = output,
    params = list(
        output = output,
        config_path = config_path
    )
)
