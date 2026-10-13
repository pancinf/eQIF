#This script extracts fine-mapped GWAS and control variants

###
### Libraries
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(reshape2)
library(optparse)
options(scipen=999)
set.seed(1234)

###
### Arguments
option_list = list(
  make_option("--finemap", type="character", default=NULL, help="path to input gwas file"),
  make_option("--chain", type="character", default=NULL, help="path to chain file"),
  make_option("--commonAf", type="character", default=NULL, help="path to file with common allele frequencies"),
  make_option("--caseIds", type="character", default=NULL, help="path to output case variant IDs"),
  make_option("--controlIds", type="character", default=NULL, help="path to output control variant IDs")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$finemap) | is.null(opt$chain) |
    is.null(opt$commonAf) | is.null(opt$caseIds) | is.null(opt$controlIds)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
### Main

##
## Read GWAS variants
gwas <- fread(
  opt$finemap,
  data.table = FALSE,
  header = TRUE,
  select = c(1,2,5,6,18)
)
gwas <- gwas[complete.cases(gwas),]

colnames(gwas)[colnames(gwas) == "ea"] <- "ref"
colnames(gwas)[colnames(gwas) == "nea"] <- "alt"

gwas <- gwas[
  nchar(gwas$ref) == 1 &
  nchar(gwas$alt) == 1,
]

gwas$chr <- paste0("chr", gwas$chr)

##
## Liftover
chain <- import.chain(opt$chain)

gwas <- makeGRangesFromDataFrame(
  gwas,
  seqnames.field = "chr",
  start.field = "bp",
  end.field = "bp",
  keep.extra.columns = TRUE,
  starts.in.df.are.0based = FALSE
)
gwas <- liftOver(gwas, chain)
uniqueLiftover <- lengths(gwas) == 1
gwas <- gwas[uniqueLiftover]
gwas <- unlist(gwas)

gwas <- data.frame(gwas)

gwas$width <- NULL
gwas$strand <- NULL
gwas$end <- NULL

colnames(gwas)[colnames(gwas) == "start"] <- "pos"
colnames(gwas)[1] <- "chr"

gwas$id <- paste0(
  gwas$chr,
  "_",
  gwas$pos,
  "_",
  gwas$ref,
  "_",
  gwas$alt
)

#Filter by susie
gwas <- gwas[order(gwas$id, -gwas$susie), ]
gwas <- gwas[!duplicated(gwas$id), ]
gwas <- gwas[gwas$susie >= 0.01, ]

##
## Assign AF
afs <- fread(
  opt$commonAf,
  header = FALSE,
  data.table = FALSE
)

colnames(afs) <- c("id","afEurNfe")

##
##Merge GWAS and afs considering flips

#Merge standard
gwas$mergeIndex <- 1:nrow(gwas)
gwasUnmerged <- gwas
gwas <- merge(gwas, afs)

#Merge allele flipped
gwasAlleleFlip <- gwasUnmerged[!(gwasUnmerged$mergeIndex %in% gwas$mergeIndex),]
gwasAlleleFlip[c("ref","alt")] <- gwasAlleleFlip[c("alt","ref")]
gwasAlleleFlip$id <- paste0(
  gwasAlleleFlip$chr, "_",
  gwasAlleleFlip$pos, "_",
  gwasAlleleFlip$ref, "_",
  gwasAlleleFlip$alt
)
gwasAlleleFlip <- merge(gwasAlleleFlip, afs)
gwas <- rbind(gwas, gwasAlleleFlip)

#Merge strand flipped
gwasStrandFlip <- gwasUnmerged[!(gwasUnmerged$mergeIndex %in% gwas$mergeIndex),]
gwasStrandFlip$ref <- chartr(
  "ACGT",
  "TGCA",
  gwasStrandFlip$ref
)
gwasStrandFlip$alt <- chartr(
  "ACGT",
  "TGCA",
  gwasStrandFlip$alt
)
gwasStrandFlip$id <- paste0(
  gwasStrandFlip$chr, "_",
  gwasStrandFlip$pos, "_",
  gwasStrandFlip$ref, "_",
  gwasStrandFlip$alt
)
gwasStrandFlip <- merge(gwasStrandFlip, afs)
gwas <- rbind(
  gwas,
  gwasStrandFlip
)

#Merge strand flipped + allele flipped
gwasStrandAlleleFlip <- gwasUnmerged[
  !(gwasUnmerged$mergeIndex %in% gwas$mergeIndex),
]
gwasStrandAlleleFlip$ref <- chartr(
  "ACGT",
  "TGCA",
  gwasStrandAlleleFlip$ref
)
gwasStrandAlleleFlip$alt <- chartr(
  "ACGT",
  "TGCA",
  gwasStrandAlleleFlip$alt
)
gwasStrandAlleleFlip[c("ref", "alt")] <-
  gwasStrandAlleleFlip[c("alt", "ref")]
gwasStrandAlleleFlip$id <- paste0(
  gwasStrandAlleleFlip$chr, "_",
  gwasStrandAlleleFlip$pos, "_",
  gwasStrandAlleleFlip$ref, "_",
  gwasStrandAlleleFlip$alt
)
gwasStrandAlleleFlip <- merge(
  gwasStrandAlleleFlip,
  afs,
  by = "id"
)
gwas <- rbind(
  gwas,
  gwasStrandAlleleFlip
)
gwas <- gwas[!duplicated(gwas$id), ]
enhancedVariantId <- unique(gwas$id)
gwas <- gwas[gwas$susie >= 0.95, ]

##
## Make controls
controls <- cbind(
  afs,
  colsplit(
    string = afs$id,
    pattern = "_",
    names = c("chr","pos","ref","alt")
  )
)

controls <- controls[
  !(controls$id %in% enhancedVariantId),
]

##
## Flatten and deduplicate
allControlIds <- unique(controls$id, use.names = FALSE)

##
## Write output
write.table(
  data.frame(caseVariantId = gwas$id),
  file = opt$caseIds,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  data.frame(controlVariantId = allControlIds),
  file = opt$controlIds,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
