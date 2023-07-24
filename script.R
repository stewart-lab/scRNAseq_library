timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output/output_", timestamp)
dir.create(output, showWarnings = FALSE)
file.copy("config.json", file.path(output, "config.json"))


rmarkdown::render(
    input = "src/script.rmd",
    output_format = "pdf_document",
    output_dir = output,
    intermediates_dir = output,
    params = list(
        output = output
    )
)
