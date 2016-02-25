#!/usr/bin/env Rscript

file = commandArgs(TRUE)[1]
correctForTumorContent = (length(commandArgs(TRUE)) > 1)

genepos = c(1,2,3,4,8,9,13,14,15,16,17,18,19,23,24,25,26,27,31,32,33,34,38,39,40,41,45,46,50,51,52,53,54,55,59,60,61,65,66,67,68,69,70,74,75,76,77,81,82,83,84,88,89,90,91,92,96,100,101,105,106,107,108,109,110,111,115,119,120,124,125,126,130,131,135)

chrpos = c("1,6","2,11","3,21","4,29","5,36","7,43","8,48","9,57","10,63","11,72","12,79","13,86","14,94","15,98","16,103","17,113","18,117","19,122","20,128","22,133","X,137")

orderedGenes = c("MYCL","BCL9","MCL1","MDM4","MYCN","MSH2","VHL","PPARG","BAP1","PIK3CA","SOX2","ATP11B","DCUN1D1","FGFR3","PDGFRA","KIT","TET2","FBXW7","TERT","PIK3R1","APC","FGFR4","IL6","EGFR","CDK6","MET","FGFR1","MYC","CD274","PDCD1LG2","CDKN2A","PTCH1","TSC1","NOTCH1","GATA3","PTEN","FGFR2","WT1","CD44","CCND1","BIRC3","BIRC2","ATM","KRAS","ACVRL1","CDK4","MDM2","FLT3","BRCA2","RB1","GAS6","APEX1","PNP","NKX2-1","NKX2-8","AKT1","IGF1R","TSC2","CDH1","TP53","TIAF1","MYO18A","NF1","ERBB2","BRCA1","RPS6KB1","SMAD4","STK11","CCNE1","CSNK2A1","BCL2L1","ZNF217","SMARCB1","NF2","AR")

boxColors = c("red","red","red","red","red","blue","blue","red","blue","red","red","red","red","red","red","red","blue","blue","red","blue","blue","red","red","red","red","red","red","red","red","red","blue","blue","blue","blue","blue","blue","red","blue","red","red","red","red","blue","red","red","red","red","red","blue","blue","red","red","red","red","red","red","red","blue","blue","blue","red","red","blue","red","blue","red","blue","blue","red","red","red","red","blue","blue","red")


rows = readLines(file)

mapd = "N/A"
for (row in rows) {
	if (grepl("##mapd=", row)) {
		mapd = strsplit(row,"=")[[1]][2]
	}
}

cellularity = "N/A"
for (row in rows) {
	if (grepl("##Cellularity", row)) {
		cellularity = as.numeric(strsplit(row,"=")[[1]][2])
	}
}


vals = vector()
genes = vector()
for (row in rows[grepl("<CNV>", rows) & grepl("gene", rows)]) {
	cols = strsplit(row,";")[[1]]
	col = cols[grepl("CDF_MAPD", cols)]
	pairs = strsplit(col,"=")[[1]][2]
	percentiles = vector()
	for (pctval in strsplit(pairs, ",")[[1]]) {
		percentiles = c(percentiles, as.numeric(strsplit(pctval, ":")[[1]][2]))
	}

	gene = gsub("(.*gene':'([A-Z0-9-]+)'.*)", "\\2", row)

	for (x in percentiles[2:(length(percentiles) - 1)]) {
		genes = c(genes, gene)
		# correct for cellularity:
		if (correctForTumorContent)  {
			x = max(0, (x - 2*(1-cellularity)) / cellularity)
		}
		vals = c(vals, x)
	}
}

genes = factor(genes, levels=orderedGenes)
df = data.frame(genes=genes, vals=vals)

if (correctForTumorContent) {
	fn = paste(gsub("((.*).vcf)", "\\2", file), "_cellularity_corrected.png", sep="")
} else {
	fn = paste(gsub("((.*).vcf)", "\\2", file), ".png", sep="")
}

png(file=fn, width=2000, height=1000)
par(mai=c(1.5,1.5,2,0.42))
yrange = range(c(0,max(6, max(vals))), vals)
boxplot(vals~genes, data=df, at=genepos, las=2, cex.axis=1.2, cex.lab=2.0, boxlwd = 1, ylab='Copy Number', ylim=yrange, pin=c(4,8), border=boxColors, lwd=2)

abline(h = 2, col = "black", lty=2, lwd=2)
abline(h = 1, col = "blue", lwd=2)
abline(h = 4, col = "red", lwd=2)
#abline(h = 7, col = "red", lwd=2)

lastpos = 1
for (c in chrpos) {
	chr = strsplit(c,",")[[1]][1]
	pos = as.numeric(strsplit(c,",")[[1]][2])
	abline(v = pos, col = "black", lty=2)
	text(x = (pos+lastpos)/2.0, y = 0, chr, cex=2)
	lastpos = pos
}


shortname = gsub(".vcf", "", tail(strsplit(file, '/')[[1]], 1))
title(paste(shortname, "\nMAPD: ", mapd, "\nCellularity: ", cellularity), cex.main=2)
x = dev.off() # suppress output by assigning it to a varianble
