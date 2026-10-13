#This script extracts fine-mapped eQTL and matching control variants

###
###Libraries
library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(reshape2)
library(optparse)
options(scipen = 999)
set.seed(1234)

###
###Arguments
option_list = list(

  make_option(
    "--inGtf",
    type = "character",
    default = NULL,
    help = "path to gtf file"
  ),

  make_option(
    "--caviar",
    type = "character",
    default = NULL,
    help = "path to input caviar file"
  ),

  make_option(
    "--commonAf",
    type = "character",
    default = NULL,
    help = "path to file with common allele frequencies"
  ),

  make_option(
    "--gwasCase",
    type = "character",
    default = NULL,
    help = "path to file with common allele frequencies"
  ),

  make_option(
    "--gwasControl",
    type = "character",
    default = NULL,
    help = "path to case gwas variants"
  ),

  make_option(
    "--caseIds",
    type = "character",
    default = NULL,
    help = "path to output case variant IDs"
  ),

  make_option(
    "--controlIds",
    type = "character",
    default = NULL,
    help = "path to output control variant IDs"
  )

)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if (is.null(opt$inGtf) |
    is.null(opt$caviar) |
    is.null(opt$commonAf) |
    is.null(opt$caseIds) |
    is.null(opt$controlIds)) {

  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###Functions
assignDist <- function(myDf) {
  variants <- makeGRangesFromDataFrame(
    df = myDf,
    seqnames = "chr",
    start.field = "pos",
    end.field = "pos",
    starts.in.df.are.0based = FALSE
  )

  gtf <- import(opt$inGtf)

  tx <- gtf[
    gtf$type == "transcript" &
    mcols(gtf)$transcript_type == "protein_coding" &
    grepl("basic", mcols(gtf)$tag)
  ]

  tss_pos <- ifelse(
    strand(tx) == "+",
    start(tx),
    end(tx)
  )

  tss_gr <- GRanges(
    seqnames = seqnames(tx),
    ranges = IRanges(tss_pos, tss_pos),
    strand = strand(tx),
    transcript_id = mcols(tx)$transcript_id,
    gene_id = mcols(tx)$gene_id
  )

  hits <- distanceToNearest(
    variants,
    tss_gr
  )

  res <- data.frame(
    chr = as.character(seqnames(variants)[queryHits(hits)]),
    pos = start(variants)[queryHits(hits)],
    nearest_tss_chr = as.character(seqnames(tss_gr)[subjectHits(hits)]),
    nearest_tss_pos = start(tss_gr)[subjectHits(hits)],
    strand = as.character(strand(tss_gr)[subjectHits(hits)]),
    distance = mcols(hits)$distance
  )

 res$relative_position <- ifelse(
  res$pos == res$nearest_tss_pos,
  "TSS",
  ifelse(
    res$strand == "+",
    ifelse(res$pos < res$nearest_tss_pos,
           "upstream_5p",
           "downstream_3p"),
    ifelse(res$pos > res$nearest_tss_pos,
           "upstream_5p",
           "downstream_3p")
		)
	)

  return(res)
}

###
###Main

##
##Read eQTL variants
eQTLs <- fread(
  opt$caviar,
  header = TRUE,
  data.table = FALSE
)
eQTLs <- eQTLs[complete.cases(eQTLs),]


eQTLs$GENE <- gsub(
  "\\..*",
  "",
  eQTLs$GENE
)

eQTLs$eQTL <- gsub(
  "_b38",
  "",
  eQTLs$eQTL
)

##
##Enhanced variants (PIP >= 0.01)
enhancedVariantId <- unique(
  eQTLs$eQTL[
    eQTLs$Probability >= 0.01
  ]
)

##
##Highly credible fine-mapped eQTLs (PIP >= 0.95)
eQTLs <- unique(
  eQTLs$eQTL[
    eQTLs$Probability >= 0.95
  ]
)

##
##optionally filter for GWAS variants
if(!is.null(opt$gwasCase)){
	gwasCaseDf <- fread(
	opt$gwasCase,
	header = FALSE,
	data.table = FALSE
	)
	eQTLs <- eQTLs[eQTLs %in% gwasCaseDf$V1]
}

eQTLs <- cbind(
  eQTLs,
  colsplit(
    string = eQTLs,
    pattern = "_",
    names = c(
      "chr",
      "pos",
      "ref",
      "alt"
    )
  )
)

colnames(eQTLs)[1] <- "id"


##
##SNPs only
eQTLs <- eQTLs[
  nchar(eQTLs$ref) == 1 &
  nchar(eQTLs$alt) == 1,
]

##
##Distance to nearest TSS
eQTLsDist <- assignDist(
  myDf = eQTLs
)

eQTLs <- merge(
  eQTLs,
  eQTLsDist
)

##
##Assign AF
afs <- fread(
  opt$commonAf,
  header = FALSE,
  data.table = FALSE
)

colnames(afs) <- c(
  "id",
  "afEurNfe"
)

eQTLs <- merge(
  eQTLs,
  afs
)

##
##Assign matching borders
eQTLs$lAf <- ifelse(
  eQTLs$afEurNfe <= 0.5,
  eQTLs$afEurNfe - 0.1 * eQTLs$afEurNfe,
  eQTLs$afEurNfe - 0.1 * (1 - eQTLs$afEurNfe)
)

eQTLs$uAf <- ifelse(
  eQTLs$afEurNfe <= 0.5,
  eQTLs$afEurNfe + 0.1 * eQTLs$afEurNfe,
  eQTLs$afEurNfe + 0.1 * (1 - eQTLs$afEurNfe)
)

eQTLs$ldistance <- eQTLs$distance - 250
eQTLs$udistance <- eQTLs$distance + 250

##
##Build control pool
controls <- cbind(
  afs,
  colsplit(
    string = afs$id,
    pattern = "_",
    names = c(
      "chr",
      "pos",
      "ref",
      "alt"
    )
  )
)

##
##Remove cases from controls
controls <- controls[
  !(controls$id %in% enhancedVariantId),
]

if(!is.null(opt$gwasControl)){
	gwasControlDf <- fread(
	opt$gwasControl,
	header = FALSE,
	data.table = FALSE
	)
	controls <- controls[controls$id %in% gwasControlDf$V1,]
}

##
##Assign distances to controls
controlsDist <- assignDist(
  myDf = controls
)

controls <- merge(
  controls,
  controlsDist
)

##
##Collect all matching controls
allControlIds <- vector(
  "list",
  nrow(eQTLs)
)

for (i in seq_len(nrow(eQTLs))) {
  controlSubset <- controls[
    controls$afEurNfe >= eQTLs$lAf[i] &
    controls$afEurNfe <= eQTLs$uAf[i] &
    controls$relative_position ==
      eQTLs$relative_position[i] &
    controls$distance >=
      eQTLs$ldistance[i] &
    controls$distance <=
      eQTLs$udistance[i],
  ]

  allControlIds[[i]] <- controlSubset$id
}

##
##Flatten and deduplicate
allControlIds <- unique(
  unlist(
    allControlIds,
    use.names = FALSE
  )
)

##
##Write case IDs
write.table(
  data.frame(
    caseVariantId = eQTLs$id
  ),
  file = opt$caseIds,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

##
##Write control IDs
write.table(
  data.frame(
    controlVariantId = allControlIds
  ),
  file = opt$controlIds,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
