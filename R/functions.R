#' This function deals GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array) RAW.tar file using the robust multi-array average expression measure.
#'
#' @param file a file name specified by either a variable of mode character, or a double-quoted string, which is of 'GSE*_RAW.tar' which is of platform GPL570.
#' @param type a character string indicating which type of data frame is to be outputted. One of 'probeIDAndGeneSymbol'(default), 'probeID', or 'geneSymbol'.
#' @return a data frame which is a expression profile.
#' @import dplyr
#' @importFrom Biobase exprs
#' @importFrom GEOquery gunzip
#' @importFrom stats aggregate median
#' @importFrom utils untar
#' @export
#' @examples
#'
#You can put your own `GSE*_RAW.tar` under the working directory
#' file <- system.file("extdata", "GSE104683_RAW.tar", package = "DealGPL570")
#' file
#' #Next step would run for about 15s, so you can try it yourself or view the vignettes
#' \dontrun{
#' result <- DealGPL570(file = file)
#' }
DealGPL570 <- function(file, type = "probeIDAndGeneSymbol") {
    exdir <- file %>% stringr::str_extract("GSE\\d+")
    untar(file, exdir = exdir)
    gz <- list.files(path = exdir, pattern = "[gz]")
    sapply(paste(exdir, gz, sep = "/"), gunzip)
    cels <- affy::list.celfiles(exdir, full.names = T)
    l <- do.call(c, lapply(cels, FUN = function(x) {
        affy::ReadAffy(filenames = x) %>% BiocGenerics::annotation()
    }))
    cels <- cels[which(l == "hgu133plus2")]
    affy_data <- affy::ReadAffy(filenames = cels)
    unlink(exdir, recursive = T)
    affy_rma <- affy::rma(affy_data)
    exprSet <- exprs(affy_rma) %>% as.data.frame() %>% tibble::rownames_to_column(var = "probe_id")
    exprSet$probe_id <- as.character(exprSet$probe_id)
    if (type == "probeID") {
      return(exprSet)
    } else if (type == "geneSymbol") {
      exprSet_gene <- exprSet %>% dplyr::inner_join(probe2symbol_df, by = "probe_id") %>% dplyr::select(-"probe_id") %>%
        dplyr::select("symbol", dplyr::everything())
      names.col <- names(exprSet_gene)
      exprSet_gene <- aggregate(exprSet_gene[ , -1], by = list(exprSet_gene$symbol), FUN = median) %>%
        as.data.frame()
      names(exprSet_gene) <- names.col
      return(exprSet_gene)
    } else if (type == "probeIDAndGeneSymbol") {
      exprSet_both <- probe2symbol_df %>% dplyr::inner_join(exprSet)
      return(exprSet_both)
    }
}
