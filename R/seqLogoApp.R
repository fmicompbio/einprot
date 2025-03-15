#' Generate a shiny app to explore sequence logos
#'
#' Launch a shiny app that displays the provided csv file as an interactive
#' table, and generates a sequence logo for the sequences in the
#' \code{seqWindow} column.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A shiny app object.
#'
#' @param seqTableCsv Character scalar, the path to a csv file. Typically
#'     this will be exported by the PTM einprot workflow. The file must have
#'     a column named \code{seqWindow}, containing sequences of a fixed
#'     width.
#' @param exportName Character scalar defining the default file name of the
#'     files exported from the shiny app.
#' @param ... Additional arguments passed to \code{utils::read.csv}.
#'
#' @examples
#' tbl <- data.frame(id = c("seq1", "seq2", "seq3"),
#'                   seqWindow = c("LITKDHE", "RQMKQPE", "LSNKVHG"))
#' tfile <- tempfile(fileext = "_seq.csv")
#' write.table(tbl, file = tfile, row.names = FALSE,
#'             col.names = TRUE, quote = TRUE, sep = ",")
#' if (interactive()) {
#'     seqLogoApp(tfile)
#' }
#'
#' @importFrom shiny fluidPage fluidPage sidebarLayout sidebarPanel
#'     mainPanel downloadButton plotOutput renderPlot reactive downloadHandler
#'     shinyApp textOutput renderText
#' @importFrom motifStack pcm2pfm colorset
#' @importFrom Biostrings consensusMatrix AAStringSet
#' @importFrom DT DTOutput renderDT
#' @importFrom writexl write_xlsx
#' @importFrom ggplot2 ggsave
#' @importFrom tools file_ext
#' @importFrom utils read.csv zip
#' @importFrom methods new
#'
seqLogoApp <- function(seqTableCsv,
                       exportName = sub("\\.csv$", "", basename(seqTableCsv)),
                       ...) {
    .assertScalar(x = seqTableCsv, type = "character")
    stopifnot(file.exists(seqTableCsv))
    stopifnot(tools::file_ext(seqTableCsv) == "csv")
    .assertScalar(x = exportName, type = "character")

    df <- utils::read.csv(seqTableCsv, ...)
    stopifnot("seqWindow" %in% colnames(df))
    if (length(unique(nchar(df$seqWindow))) != 1) {
        stop("All sequences in the seqWindow column must have the same width.")
    }

    ui <- shiny::fluidPage(
        shiny::titlePanel("Sequence logo generation"),

        shiny::sidebarLayout(
            shiny::sidebarPanel(width = 2,
                                shiny::textOutput("description"),
                                shiny::downloadButton("dl", "Download")),

            shiny::mainPanel(
                width = 10,
                DT::DTOutput(
                    outputId = "seqtable"
                ),
                shiny::plotOutput(
                    outputId = "seqlogo"
                )
            )
        )
    )

    #nocov start
    server <- function(input, output, session) {
        output$seqtable <- DT::renderDT(
            df, rownames = FALSE,
            extensions = "FixedColumns",
            filter = list(position = "top", clear = FALSE),
            options = list(scrollX = TRUE, pageLength = 5,
                           search = list(regex = TRUE, caseInsensitive = FALSE),
                           fixedColumns = list(leftColumns = 1))
        )

        output$description <- shiny::renderText(
            paste0("This is a simple interactive application to generate ",
                   "sequence logos. The table can be filtered using any of ",
                   "the columns, and the plot below it will show a ",
                   "sequence logo generated from the values in the ",
                   "'seqWindow' column, for all retained rows. Clicking on ",
                   "the 'Download' button will download a zip file with a ",
                   "pdf version of the sequence logo, as well as csv and ",
                   "xlsx files containing the retained rows of the table.\n\n"))

        seqlogo <- shiny::reactive({
            if (is.null(input$seqtable_rows_all)) {
                NULL
            } else {
                ## Replace ggseqlogo (removed from CRAN) with motifStack
                seqs <- df[input$seqtable_rows_all, "seqWindow"]
                pfm <- motifStack::pcm2pfm(Biostrings::consensusMatrix(
                    Biostrings::AAStringSet(seqs)))
                pfm <- new("pfm", mat = pfm, name = "",
                           color = motifStack::colorset(alphabet = "AA",
                                                        colorScheme = "chemistry"))
                pfm
                # ggseqlogo::ggseqlogo(df[input$seqtable_rows_all, "seqWindow"],
                #                      seq_type = "aa") +
                #     ggseqlogo::theme_logo(base_size = 15)
            }
        })
        output$seqlogo <- shiny::renderPlot({
            if (!is(seqlogo(), "pfm")) {
                NULL
            } else {
                motifStack::plotMotifLogo(seqlogo(), font = "sans", fontface = "plain")
                ## With ggseqlogo
                # seqlogo()
            }
        })

        output$dl <- shiny::downloadHandler(
            paste0(exportName, "_seqlogo_export.zip"),
            content = function(file) {
                dumptmp <- tempfile()
                dir.create(dumptmp)
                oldwd <- getwd()
                setwd(dumptmp)

                on.exit({
                    setwd(oldwd)
                    unlink(dumptmp, recursive = TRUE)
                })

                timestamp <- format(Sys.time(), "%Y-%m-%d-%H-%M")
                plotfile <- paste0(exportName, "-seqlogo-", timestamp, ".pdf")
                xlsxfile <- paste0(exportName, "-seqlogo-", timestamp, ".xlsx")
                csvfile <- paste0(exportName, "-seqlogo-", timestamp, ".csv")
                pdf(plotfile, width = 8, height = 5)
                motifStack::plotMotifLogo(seqlogo(), font = "sans", fontface = "plain")
                dev.off()
                # ggplot2::ggsave(seqlogo(), file = plotfile,
                #                 width = 8, height = 5)
                writexl::write_xlsx(df[input$seqtable_rows_all, ],
                                    path = xlsxfile)
                utils::write.table(df[input$seqtable_rows_all, ],
                                   file = csvfile,
                                   sep = ",", quote = TRUE,
                                   row.names = TRUE, col.names = TRUE)

                utils::zip(file, files = c(plotfile, xlsxfile, csvfile))
            }
        )
    }
    #nocov end

    shiny::shinyApp(ui, server)
}
