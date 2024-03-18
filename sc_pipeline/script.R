timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output <- paste0("output/output_", timestamp)
dir.create(output, showWarnings = FALSE)

# copy the config file into the output directory
config_path <- file.path(getwd(), output, "config.json")
file.copy(file.path(getwd(), "/src/config.json"), config_path)
# print statement to confirm the start of the pipeline and the output directory
cat(paste0("Starting pipeline at ", timestamp, " with output directory ", output, "\n"))

rmarkdown::render(
    input = file.path(getwd(), "src", "sc_pipe_new.rmd"),
    output_format = "pdf_document",
    output_dir = output,
    intermediates_dir = output,
    params = list(
        output = output,
        config_path = config_path
    )
)
