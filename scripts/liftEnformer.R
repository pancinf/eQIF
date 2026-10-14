#This script lifts enformer to hg38

###
### Libraries
library(data.table)
library(rtracklayer)
library(GenomicRanges)
library(optparse)
options(scipen=999)
set.seed(1234)

###
### Arguments
option_list = list(
  make_option("--chain", type="character", default=NULL, help="path to chain file"),
  make_option("--enfDir", type="character", default=NULL, help="path to enformer directory"),
  make_option("--out", type="character", default=NULL, help="path to hg38 enformer output")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$chain) | is.null(opt$enfDir) |
    is.null(opt$out)) {
  print("You have to specify all inputs. Here is the help:")
  print_help(opt_parser)
  quit()
}

###
###Main
chain <- import.chain(opt$chain)
myPath <- opt$enfDir
myFiles <- list.files(path = myPath, pattern = ".tsv.bgz")
enf <- data.frame(chr = c(),pos = c(), ref = c(), alt = c(),PC0 = c(), PC1 = c(),PC2 = c(),PC3 = c(),PC4 = c(),PC5 = c(),PC6 = c(),PC7 = c(),PC8 = c(),PC9 = c(),PC10 = c(),PC11 = c(),PC12 = c(),PC13 = c(),PC14 = c(),PC15 = c(),PC16 = c(),
	PC17 = c(),PC18 = c(),PC19 = c())
for(i in 1:length(myFiles)){
	enfShort <- fread(paste0(myPath,myFiles[i]), select = c(1,2,4:25), data.table = FALSE, header = TRUE)
	enfShort <- enfShort[complete.cases(enfShort),]
	colnames(enfShort)[1:4] <- c("chr","pos","ref","alt")
	enfShort <- enfShort[nchar(enfShort$ref) == 1 & nchar(enfShort$alt) == 1,]
	enfShort <- makeGRangesFromDataFrame(enfShort, seqnames.field = "chr", start.field = "pos",
	                                 end.field = "pos",keep.extra.columns = TRUE, starts.in.df.are.0based = FALSE)

	enfShort <- liftOver(enfShort, chain)
	uniqueLiftover <- lengths(enfShort) == 1
	enfShort <- enfShort[uniqueLiftover]
	enfShort <- unlist(enfShort)
	enfShort <- data.frame(enfShort)
	enfShort$width <- NULL
	enfShort$strand <- NULL
	enfShort$end <- NULL
	colnames(enfShort)[1:2] <- c("chr","pos")
	enfShort$chr <- gsub("chr","",enfShort$chr)
	enf <- rbind(enf,enfShort)
}


##
##Write output
write.table(enf, file = opt$out,sep = "\t", col.names = TRUE, row.names = FALSE, append = FALSE, quote = FALSE)
