#!/usr/bin/env Rscript
## synplot.R by Brandon Seah (kbseah@mpi-bremen.de)
## Use with synplot.pl

# Rscript synplot.R --args tab0file tab1file tab2file tab3file outname cdstype colmax

## Get arguments from command line
args         <- commandArgs(trailingOnly=TRUE)
tab0file <- args[2]
tab1file <- args[3]
tab2file <- args[4]
tab3file <- args[5]
outname <- args[6]
cdstype <- args[7]
colmax <- args[8]

# Read input tables

tab0 <- read.table(file=tab0file,header=T)
tab1 <- read.table(file=tab1file,header=T)
tab2 <- read.table(file=tab2file,header=T)
tab3 <- read.table(file=tab3file,header=T)

pdf(file=outname,width=10,height=7)

par(mar=c(5,12,4,2)+0.1)
plot("0",
     type="n",
     xlim=c(0,max(tab1$stop)),
     ylim=c(0,max(tab1$y)),
     ylab="",
     xlab="Position (bp)",
     yaxt="n"
     )
axis(side=2,
     at=tab0$y,
     labels=tab0$label,
     las=2)

# Define color palette for ID values
colfunc <- colorRampPalette(c("white",colmax))

# Plot polygons corresponding to pairs of best Blast hits
for (i in 1:length(tab0$genome)) { 
    tab3.subset <- subset(tab3,
                          genome1 == as.character(tab0$genome[i])
                          )
    xses <- tab3.subset[,5:9]
    yses <- tab3.subset[,10:14]
    for (j in 1:dim(xses)[1]) {
        #polygon (xses[j,], yses[j,],col="pink",lty=0)
        polygon(xses[j,],
                yses[j,],
                col=colfunc(100)[round(tab3$pid[j])],
                lty=0)
    }
}

# Plot CDS regions
if (cdstype=="color") {
    # Plot CDS regions on top so that the colors will show
    for (i in 1:length(tab0$genome)) {      
        tab2.subset <- subset(tab2,
                              genome==as.character(tab0$genome[i])
                              )
        segments (tab2.subset$cumulstart,
                  rep(tab0$y[i],
                      dim(tab2.subset)[1]
                      ),
                  tab2.subset$cumulstop,
                  rep(tab0$y[i],
                      dim(tab2.subset)[1]
                      ),
                  lwd=5,
                  lend="butt", # Prevent round caps
                  col=as.character(tab2.subset$color)
                  )
    }
} else if (cdstype=="arrow") {
    # Plot CDS regions with arrows instead of line segments
    for (i in 1:length(tab0$genome)) {
        tab2.subset <- subset(tab2,
                              genome==as.character(tab0$genome[i])
                              )
        arrows (tab2.subset$cumulstart,
                rep(tab0$y[i],
                    dim(tab2.subset)[1]
                    ),
                tab2.subset$cumulstop,
                rep(tab0$y[i],
                    dim(tab2.subset)[1]
                    ),
                lwd=5,
                lend="butt",
                col=as.character(levels(tab2.subset$color)[1]),
                length=0.04
                )
    }
}

# Draw tick marks indicating contig boundaries
for (i in 1:length(tab1$contig)) { 
    segments(x0=tab1$start[i],
             x1=tab1$start[i],
             y0=tab1$y[i]-0.2,
             y1=tab1$y[i]+0.2,
             lend="butt"
             )
    segments(x0=tab1$stop[i],
             x1=tab1$stop[i],
             y0=tab1$y[i]-0.2,
             y1=tab1$y[i]+0.2,
             lend="butt"
             )
}

#text(x=rep(0,length(tab0$genome)),y=tab0$y,pos=2,labels=tab0$genome)
legend(x="topright",
       fill=colfunc(100)[seq(10,100,10)],
       legend=seq(10,100,10),
       title="Percent ID",
       y.intersp=0.65,
       cex=0.9,
       border="white",
       bty="n"
       )

dev.off()

