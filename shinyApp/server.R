#List of Cancer Gene Census (downloaded on 12.07.2016)
cgc <- c("ABI1","ABL1","ABL2","ACSL3","ACVR1","CASC5","MLLT11","AKAP9","AKT1","AKT2","ALDH2","ALK","RNF213","AMER1",
         "APC","ARHGEF12","RHOH","ARID1A","ARID1B","ARID2","ARNT","ASPSCR1","ASXL1","ATF1","ATIC","ATM","ATP1A1",
         "ATP2B3","ATR","ATRX","AXIN1","AXIN2","BAP1","BCL10","BCL11A","BCL11B","BCL2","BCL3","BCL5","BCL6","BCL7A",
         "BCL9","BCOR","BCR","BHD","BIRC3","BLM","BMPR1A","BRAF","BRCA1","BRCA2","BRD3","BRD4","BRIP1","BTG1","BUB1B",
         "HMGN2P46","NUTM1","RMI2","C2orf44","CACNA1D","CALR","CAMTA1","CANT1","CARD11","CARS","CASP8","RUNX1T1",
         "CBFA2T3","CBFB","CBL","CBLB","CBLC","CCDC6","CCNB1IP1","CCND1","CCND2","CCND3","CCNE1","PDCD1LG2","CD274",
         "CD74","CD79A","CD79B","CDC73","CDH1","CDH11","CDK12","CDK4","CDK6","CDKN1B","CDKN2A","CDKN2C","CDX2","CEBPA",
         "CNTRL","CEP89","CHCHD7","CHEK2","CHIC2","CHN1","CIC","CIITA","CLIP1","CLTC","CLTCL1","ACKR3","CNOT3","COL1A1",
         "COL2A1","KLF6","COX6C","CREB1","CREB3L1","CREB3L2","CREBBP","CRLF2","CRTC1","CRTC3","CSF3R","CTNNB1","CUX1",
         "CYLD","DAXX","DCTN1","DDB2","DDIT3","DDX10","DDX5","DDX6","DEK","DICER1","DNM2","DNMT3A","DUX4","EBF1",
         "ECT2L","EGFR","EIF3E","EIF4A2","ELF4","ELK4","ERC1","ELL","ELN","EML4","EP300","EPS15","ERBB2","ERC1","ERCC2",
         "ERCC3","ERCC4","ERCC5","ERG","ETNK1","ETV1","ETV4","ETV5","ETV6","MECOM","EWSR1","EXT1","EXT2","EZH2","EZR",
         "ACSL6","FAM131B","FAM46C","FANCA","FANCC","FANCD2","FANCE","FANCF","FANCG","FAS","FBXO11","FBXW7","FCGR2B",
         "FEV","FGFR1","FGFR1OP","FGFR2","FGFR3","FGFR4","FH","FHIT","FIP1L1","FLI1","FLT3","FNBP1","FOXA1","FOXL2",
         "FOXO1","FOXO3","FOXO4","FOXP1","FSTL3","FUBP1","FUS","KDSR","GAS7","GATA1","GATA2","GATA3","GMPS","GNA11",
         "GNAQ","GNAS","GOLGA5","GOPC","GPC3","GPHN","ARHGAP26","GRIN2A","H3F3A","H3F3B","SPECC1","CLP1","HERPUD1",
         "HEY1","HIP1","HIST1H3B","HIST1H4I","HLA-A","HLF","MNX1","HMGA1","HMGA2","HNRNPA2B1","HOOK3","HOXA11","HOXA13",
         "HOXA9","HOXC11","HOXC13","HOXD11","HOXD13","HRAS","HSP90AA1","HSP90AB1","IDH1","IDH2","IGH","IGK","IGL",
         "IKBKB","IKZF1","IL2","IL21R","IL6ST","IL7R","IRF4","FCRL4","ITK","JAK1","JAK2","JAK3","JAZF1","JUN","KCNJ5",
         "KDM5A","KDM5C","KDM6A","KDR","KIAA1549","SHTN1","KIF5B","KIT","KLF4","KLK2","KMT2D","KRAS","KTN1","AFF3",
         "LASP1","LCK","LCP1","TET1","LHFP","LIFR","LMNA","LMO1","LMO2","LPP","LRIG3","LSM14A","LYL1","MAF","MAFB",
         "MALAT1","MALT1","MAML2","MAP2K1","MAP2K2","MAP2K4","MAP3K1","MAP3K13","MAX","MDM2","MDM4","MECOM","MDS2",
         "CRTC1","MED12","MEN1","MET","MITF","MKL1","MLF1","MLH1","KMT2A","KMT2C","MLLT1","MLLT10","AFF1","MLLT3",
         "MLLT4","MLLT6","FOXO4","MN1","MPL","FN1","MSH2","MSH6","MSI2","MSN","MTCP1","MUC1","MUTYH","MYB","MYC",
         "MYCL","MYCN","MYD88","MYH11","MYH9","MYO5A","MYOD1","KAT6B","NAB2","NACA","NBN","NCOA1","NCOA2","NCOA4",
         "NCOR1","NDRG1","NF1","NF2","NFATC2","NFE2L2","NFIB","NFKB2","NFKBIE","NIN","NKX2-1","NONO","NOTCH1","NOTCH2",
         "NPM1","NR4A3","NRAS","NRG1","NSD1","NT5C2","NTRK1","NTRK3","NUMA1","NUP214","NUP98","NUTM2A","NUTM2B",
         "OLIG2","OMD","P2RY8","PAFAH1B2","PALB2","PAX3","PAX5","PAX7","PAX8","PBRM1","PBX1","PCM1","PCSK7","PDE4DIP",
         "PDGFB","PDGFRA","PDGFRB","PER1","PHF6","PHOX2B","PICALM","PIK3CA","PIK3R1","PIM1","PLAG1","PLCG1","PML",
         "PMS1","PMS2","PRRX1","SEPT5","POLE","POT1","POU2AF1","POU5F1","PPARG","PPFIBP1","PPP2R1A","PPP6C","PRCC",
         "PRDM1","PRDM16","PRF1","PRKAR1A","PSIP1","PTCH1","PTEN","PTPN11","PTPRB","PTPRC","PTPRK","PWWP2A","RABEP1",
         "RAC1","RAD21","RAD51B","RAF1","RALGDS","RANBP17","RANBP2","RAP1GDS1","RARA","RB1","RBM15","RECQL4","REL",
         "RET","RHOA","RNF43","ROS1","RPL10","RPL22","RPL5","RPN1","RSPO2","RSPO3","SNX29","RUNX1","KAT6A","SBDS",
         "SDC4","SDHAF2","SDHB","SDHC","SDHD","SET","SETBP1","SETD2","SF3B1","SFPQ","SRSF3","SH2B3","SH3GL1","PMEL",
         "SLC34A2","SLC45A3","SMAD4","SMARCA4","SMARCB1","SMARCD1","SMARCE1","SMO","SND1","SOCS1","SOX2","SPEN",
         "SPOP","SRGAP3","SRSF2","SS18","SS18L1","SSX1","SSX2","SSX4","STAG2","STAT3","STAT5B","STAT6","STK11",
         "RNF217-AS1","STRN","SUFU","SUZ12","SYK","TAF15","TAL1","TAL2","TBL1XR1","TBX3","TCEA1","HNF1A","TCF12",
         "TCF3","TCF7L2","TCL1A","TCL6","TERT","TET2","TFE3","TFEB","TFG","TFPT","TFRC","THRAP3","TRIM24","TLX1",
         "TLX3","TMPRSS2","TNFAIP3","TNFRSF14","TNFRSF17","TOP1","TP53","TPM3","TPM4","TPR","TRA","TRAF7","TRB",
         "TRD","TRIM27","TRIM33","TRIP11","TRRAP","TSC1","TSC2","TSHR","TTL","U2AF1","UBR5","USP6","VHL","VTI1A",
         "WAS","WHSC1","WHSC1L1","WIF1","WRN","WT1","WWTR1","XPA","XPC","XPO1","YWHAE","ZCCHC8","ZBTB16","ZMYM2",
         "PATZ1","ZNF331","ZNF384","ZNF521","CNBP","ZRSR2")

#Known amplified or deleted genes in cancer (Zack et al., Nature 2013)
known.ampl <- c('PAX8','NEDD9','E2F3','CDK6','AKT1','ZNF217','KRAS','BCL2L1','IGF1R','MCL1','TERT','MDM4','PDGFRA',
                'SOX2','CDK4','MDM2','MCL1','CCNE1','ERBB2','MYC','MYCL','NKX2-1','EGFR','CCND1','CCND2','MET','ERBB3')
known.del <- c('SMAD4','MYCN','CDKN1B','NOTCH1','BRCA1','ATM','IKZF2','KMT2C','APC','NF1','FAT1','RB1','PTEN',
               'ARID1A','PARK2','STK11','CDKN2A','TP53','CDKN2A')

#Actionable amplification/deletion in CIVIC (as of 11.07.2016)
civic <- read.table('CIVIC-ampl_del-20160711.txt', sep='\t', header = TRUE, text = '"', stringsAsFactors = FALSE)
civic$Type <- rep('Gain', length(civic$Status))
civic$Type[grep('DELETION', civic$Variant)] <- 'Loss'
civic.ampl <- unique(civic$Gene[civic$Type=='Gain'])
civic.del <- unique(civic$Gene[civic$Type=='Loss'])

#Function to cross a list of genes with the CGC genes, CIVIC and the ampl/del from Zack et al.
#Input is a list of items fomr a line of a text export from ChAS. Expects to find the gene list in positition 4 (1-based).
#Returns the positions 1 (Type), 2 (CN State), 3 (Size) and 5 (Cytoband Start). Appends three extra columns with the matching
#genes from CGC, Zack et al., and CIVIC.
crossLists <- function(tup){
  genes <- unlist(strsplit(tup[4], ', '))
  genes.cgc <- genes[genes %in% cgc]
  if(tup[1]=='Gain'){
    genes.known <- genes[genes %in% known.ampl]
    genes.civic <- genes[genes %in% civic.ampl]
  } else {
    genes.known <- genes[genes %in% known.del]
    genes.civic <- genes[genes %in% civic.del]
  }
  return(c(tup[1], tup[2], tup[3], tup[5], CGC=paste(genes.cgc, collapse=','), known=paste(genes.known, collapse=','), civic=paste(genes.civic, collapse=',')))
}

#Reads in an export file from ChAS. Expects the following columns: 
# 'Chromosome', 'Cytoband Start', 'Type', 'CN State', 'Size (kbp)' and 'Genes'.
#Filters out all lines that do not have a gene in either CGC, Zack et al., or CIVIC.
#Returns a data frame.
readResults <- function(fn){
  dat <- read.table(fn, sep='\t', header = TRUE, text = '"', stringsAsFactors = FALSE)
  dat.short <- dat[, c('Type','CN.State','Size..kbp.','Genes')]
  dat.short$Position <- paste0(dat$Chromosome, dat$Cytoband.Start)
  
  dat.xed <- data.frame(t(apply(dat.short, 1, crossLists)))
  colnames(dat.xed) <- c('Type', 'CN.state', 'Size','Cytoband.start', 'Cancer.genes', 'Known.gain_loss', 'CIVIC')
  
  dat.xed.short <- dat.xed[dat.xed$Cancer.genes!="" | dat.xed$Known.gain_loss!="" | dat.xed$CIVIC!="", ]
  return(dat.xed.short)
}


# By default, the file size limit is 5MB. It can be changed by
# setting this option. Here we'll raise limit to 9MB.
options(shiny.maxRequestSize = 9*1024^2)

#Basic shiny server function to render the page (calls readResults)
shinyServer(function(input, output) {
  output$contents <- renderTable({
    # input$file1 will be NULL initially. After the user selects
    # and uploads a file, it will be a data frame with 'name',
    # 'size', 'type', and 'datapath' columns. The 'datapath'
    # column will contain the local filenames where the data can
    # be found.
    
    fn <- input$file1
    
    if (is.null(fn))
      return(NULL)
    
    readResults(fn$datapath)
  })
})

