# This script extracts all CDS regions for basic transcripts

###
###
###Libraries
library(data.table)
library(GenomicRanges)
library(optparse)
options(scipen=999)
set.seed(1234)

###
### Arguments
option_list = list(
  make_option("--inGtf", type="character", default=NULL, help="path to input gtf file"),
  make_option("--outBed", type="character", default=NULL, help="path to output bed file")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$inGtf) | is.null(opt$outBed)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###
###Main

##
##Read in GENCODE data
gencode <- fread(opt$inGtf,data.table = FALSE,header = FALSE)
colnames(gencode) <- c("Chr","source","type","start","end","dot1","strand","dot2","Anno")
gencode <- gencode[,c(1,3,4,5,9)]
gencode <- gencode[(grepl("tag \"basic\";",gencode$Anno)==TRUE),]

##
##Get exon positions for coding
gencodeCoding <- gencode[(grepl("gene_type \"protein_coding\";",gencode$Anno)==TRUE),]
gencodeExon <- gencodeCoding[gencodeCoding$type == "CDS",]
gencodeExon$Anno <- gsub("\\;.*","",gencodeExon$Anno)
gencodeExon$Anno <- gsub("\\..*","",gencodeExon$Anno)
gencodeExon$Anno <- gsub(".* \"", "", gencodeExon$Anno)
gencodeExon <- makeGRangesFromDataFrame(df = gencodeExon,seqnames.field = "Chr", start.field = "start", end.field = "end",starts.in.df.are.0based = FALSE, keep.extra.columns = TRUE)

gencodeExonTrans <- GenomicRanges::reduce(gencodeExon, ignore.strand=TRUE)

finalRemReg <- data.frame(gencodeExonTrans)
finalRemReg$width <- NULL
finalRemReg$strand <- NULL
finalRemReg$start <- finalRemReg$start -1

#
###Write output
write.table(finalRemReg,opt$outBed, append = FALSE, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
