library(plot3D)


#CHeck means of quintiles within 2D scatter plot
means<-read.table("./output/Quantitative/MergedControls_Quintiles.txt",sep="\t",header=TRUE)
pdf(file="./output/Quantitative/TotalMeansQuintiles.pdf")
cols=rainbow(5)
plot(data$Broad,data$Small,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(0,30),ylim=c(0,30),col=rgb(0,0,0,0.1))
points(means$X2K,means$NDR,col=cols,bg=cols,pch=21)
legend(0,30,c("Quintile1 (Mean FPKM: 0.001)",
"Quintile2 (Mean FPKM: 0.173)",
"Quintile3 (Mean FPKM: 1.789)",
"Quintile4 (Mean FPKM: 7.637)",
"Quintile5 (Mean FPKM: 137.9)"),fill=cols)
dev.off()

pdf(file="./output/Quantitative/TotalMeansZoomQuintiles.pdf")
cols=rainbow(5)
plot(data$Broad,data$Small,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(10,20),ylim=c(5,15),col=rgb(0,0,0,0.1))
points(means$X2K,means$NDR,col=cols,bg=cols,pch=21)
legend(10,15,c("Quintile1 (Mean FPKM: 0.001)",
"Quintile2 (Mean FPKM: 0.173)",
"Quintile3 (Mean FPKM: 1.789)",
"Quintile4 (Mean FPKM: 7.637)",
"Quintile5 (Mean FPKM: 137.9)"),fill=cols)
dev.off()



#Check the same for deciles
means<-read.table("./output/Quantitative/MergedControls_Deciles.txt",sep="\t",header=TRUE)
pdf(file="./output/Quantitative/TotalMeansDeciles.pdf")
cols=rainbow(10)
plot(data$Broad,data$Small,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(0,30),ylim=c(0,30),col=rgb(0,0,0,0.1))
points(means$X2K,means$NDR,col=cols,bg=cols,pch=21)
legend(0,30,c("Deciles1 (Mean FPKM: 0.001)",
"Deciles2  (Mean FPKM: 0.001)",
"Deciles3  (Mean FPKM: 0.044)",
"Deciles4  (Mean FPKM: 0.303)",
"Deciles5  (Mean FPKM: 1.037)",
"Deciles6  (Mean FPKM: 2.546)",
"Deciles7  (Mean FPKM: 5.121)",
"Deciles8  (Mean FPKM: 10.17)",
"Deciles9  (Mean FPKM: 22.15)",
"Deciles10 (Mean FPKM: 253.9)"),fill=cols)
dev.off()

pdf(file="./output/Quantitative/TotalMeansZoomDeciles.pdf")
cols=rainbow(10)
plot(data$Broad,data$Small,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(10,20),ylim=c(5,15),col=rgb(0,0,0,0.1))
points(means$X2K,means$NDR,col=cols,bg=cols,pch=21)
legend(10,15,c("Decile1 (Mean FPKM: 0.001)",
"Decile2  (Mean FPKM: 0.001)",
"Decile3  (Mean FPKM: 0.044)",
"Decile4  (Mean FPKM: 0.303)",
"Decile5  (Mean FPKM: 1.037)",
"Decile6  (Mean FPKM: 2.546)",
"Decile7  (Mean FPKM: 5.121)",
"Decile8  (Mean FPKM: 10.17)",
"Decile9  (Mean FPKM: 22.15)",
"Decile10 (Mean FPKM: 253.9)"),fill=cols)
dev.off()

#add central points for percentiles
pdf(file="./output/Quantitative/MergedControls_PercentileMeans.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
data_2K_TSS<-rep(0,100)
data_NDR<-rep(0,100)
for( i in c(1:100)){
    data_2K_TSS[i]<-mean(data$X2K_Cov[which(data$FPKM_percentile < i & data$FPKM_percentile > i-1)])
    data_NDR[i]<-mean(data$NDR_Cov[which(data$FPKM_percentile < i & data$FPKM_percentile > i-1)])
}
plot(data$X2K_Cov,data$NDR_Cov,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(0,30),ylim=c(0,30),col=rgb(0,0,0,0.1))
colfunc <- colorRampPalette(c("green", "red"))
points(data_2K_TSS,data_NDR,pch=21,col="black",bg=colfunc(100))
legend(0,30,c("Percentile1",
"Percentile10",
"Percentile20",
"Percentile30",
"Percentile40",
"Percentile50",
"Percentile60",
"Percentile70",
"Percentile80",
"Percentile90",
"Percentile100"),fill=c(colfunc(100)[1],colfunc(100)[10],colfunc(100)[20],colfunc(100)[30],colfunc(100)[40],colfunc(100)[50],colfunc(100)[60],colfunc(100)[70],colfunc(100)[80],colfunc(100)[90],colfunc(100)[100]))
dev.off()

#Check variability of each quintile
pdf(file="./output/Quantitative/MergedControls_QuntileDistribution.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
data_2K_TSS<-rep(0,5)
data_NDR<-rep(0,5)
colfunc <- rainbow(5)
col_distribution<-rep(rgb(0,0,0,0),length(data$X2K_Cov))
for( i in c(1:5)){
    data_2K_TSS[i]<-mean(data$X2K_Cov[which(data$FPKM_percentile < i*20 & data$FPKM_percentile > (i-1)*20)])
    data_NDR[i]<-mean(data$NDR_Cov[which(data$FPKM_percentile < i*20 & data$FPKM_percentile > (i-1)*20)])
    col_distribution[which(data$FPKM_percentile < i*20 & data$FPKM_percentile > (i-1)*20)]<-colfunc[i]
}
plot(data$X2K_Cov,data$NDR_Cov,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(0,30),ylim=c(0,30),pch=".",col=col_distribution)
points(data_2K_TSS,data_NDR,pch=21,col="black",bg=colfunc)
legend(0,30,c("Quintile1 (Mean FPKM: 0.001)",
"Quintile2 (Mean FPKM: 0.173)",
"Quintile3 (Mean FPKM: 1.789)",
"Quintile4 (Mean FPKM: 7.637)",
"Quintile5 (Mean FPKM: 137.9)"),fill=colfunc)
dev.off()


#Check variability of each decile
pdf(file="./output/Quantitative/MergedControls_DecileDistribution.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
data_2K_TSS<-rep(0,10)
data_NDR<-rep(0,10)
colfunc <- rainbow(10)
col_distribution<-rep(rgb(0,0,0,0),length(data$X2K_Cov))
for( i in c(1:10)){
    data_2K_TSS[i]<-mean(data$X2K_Cov[which(data$FPKM_percentile < i*10 & data$FPKM_percentile > (i-1)*10)])
    data_NDR[i]<-mean(data$NDR_Cov[which(data$FPKM_percentile < i*10 & data$FPKM_percentile > (i-1)*10)])
    col_distribution[which(data$FPKM_percentile < i*10 & data$FPKM_percentile > (i-1)*10)]<-colfunc[i]
}
plot(data$X2K_Cov,data$NDR_Cov,xlab="2K-TSS Coverage",ylab="NDR Coverage",xlim=c(0,30),ylim=c(0,30),pch=".",col=col_distribution)
points(data_2K_TSS,data_NDR,pch=21,col="black",bg=colfunc)
legend(0,30,c("Decile1 (Mean FPKM: 0.001)",
"Decile2  (Mean FPKM: 0.001)",
"Decile3  (Mean FPKM: 0.044)",
"Decile4  (Mean FPKM: 0.303)",
"Decile5  (Mean FPKM: 1.037)",
"Decile6  (Mean FPKM: 2.546)",
"Decile7  (Mean FPKM: 5.121)",
"Decile8  (Mean FPKM: 10.17)",
"Decile9  (Mean FPKM: 22.15)",
"Decile10 (Mean FPKM: 253.9)"),fill=colfunc)
dev.off()

#Plot distribution of mean FPKM values of four non-pregnant women
fpkms<-read.table("./ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",sep="\t",header=TRUE)
pdf(file="./output/Quantitative/fpkm_distribution.pdf")
hist(fpkms$Mean.FPKM, breaks=10000,xlim=c(0,100))
dev.off()

pdf(file="./output/Quantitative/fpkm_distribution_zoom.pdf")
hist(fpkms$Mean.FPKM, breaks=100000,xlim=c(0,10))
dev.off()
pdf(file="./output/Quantitative/fpkm_distribution_zoom_omitzero.pdf")
hist(fpkms$Mean.FPKM[which(fpkms$Mean.FPKM>0)], breaks=100000,xlim=c(0,10))
dev.off()

pdf(file="./output/Quantitative/MergedControls_FPKM_3D_Scatter.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
scatter3D(data$X2K_Cov,data$NDR_Cov,data$FPKM,xlab="2K-TSS Coverage",ylab="NDR Coverage", zlab="FPKM",xlim=c(0,30),ylim=c(0,30),zlim=c(0,100),clim=c(0,150))
dev.off()


pdf(file="./output/Quantitative/MergedControls_FPKM_3D_Scatter_Ranks.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
scatter3D(data$X2K_Cov,data$NDR_Cov,data$FPKM_percentile,xlab="2K-TSS Coverage",ylab="NDR Coverage", zlab="FPKM_percentile",xlim=c(0,30),ylim=c(0,30),zlim=c(0,100),clim=c(0,100))
dev.off()

#3D Histogram of binned mean fpkms
pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM_binned.txt",sep="\t",header=TRUE)
fpkm_masked<-data$FPKM
fpkm_masked[is.na(fpkm_masked)]<-0
fpkm_values<-matrix(fpkm_masked,nrow=30,byrow=TRUE)
hist3D(z=fpkm_values,border="black",bty="g",space=0.3,shade=0.2,ticktype = "detailed",xlab="2K-TSS Coverage",ylab="NDR Coverage",zlab="FPKM")
dev.off()


#3D Histogram of binned mean fpkms. FPKM values > 100 are capped to 100
pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot_100FPKMcap.pdf")
data<-read.table("MergedControls_FPKM_binned.txt",sep="\t",header=TRUE)
fpkm_masked<-data$FPKM
fpkm_masked[is.na(fpkm_masked)]<-0
fpkm_masked[which(fpkm_masked > 100)]<-100
fpkm_values<-matrix(fpkm_masked,nrow=30,byrow=TRUE)
hist3D(z=fpkm_values,border="black",bty="g",space=0.3,shade=0.2,ticktype = "detailed",zlim=c(0,100),xlab="2K-TSS Coverage",ylab="NDR Coverage",zlab="FPKM")
dev.off()

#3D Histogram of binned mean fpkms. FPKM values > 100 are capped to 100; filter out bins with < 10 samples
pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot_100FPKMcap_masked.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM_binned.txt",sep="\t",header=TRUE)
fpkm_masked<-data$FPKM
fpkm_masked[is.na(fpkm_masked)]<-0
fpkm_masked[which(fpkm_masked > 100)]<-100
fpkm_masked[which(data$Count < 11)]<-0
fpkm_values<-matrix(fpkm_masked,nrow=30,byrow=TRUE)
hist3D(z=fpkm_values,border="black",bty="g",space=0.3,shade=0.2,ticktype = "detailed",zlim=c(0,100),xlab="2K-TSS Coverage",ylab="NDR Coverage",zlab="FPKM")
dev.off()


#Make 3D Histograms of percentiles
pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot_Percentiles.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM_binned.txt",sep="\t",header=TRUE)
perc_masked<-data$FPKM_percentile
perc_masked[is.na(perc_masked)]<-0
perc_values<-matrix(perc_masked,nrow=30,byrow=TRUE)
hist3D(x=c(0:29),y=c(0:29),z=perc_values,border="black",xlab="2K-TSS Coverage",ylab="NDR Coverage",zlab="FPKM percentiles",bty="g",space=0.3,shade=0.2)
dev.off()

#Make 3D Histograms of percentiles
pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot_Percentiles_masked.pdf")
data<-read.table("./output/Quantitative/MergedControls_FPKM_binned.txt",sep="\t",header=TRUE)
perc_masked<-data$FPKM_percentile
perc_masked[which(data$Count < 11)]<-0
perc_masked[is.na(perc_masked)]<-0
perc_values<-matrix(perc_masked,nrow=30,byrow=TRUE)
hist3D(x=c(0:29),y=c(0:29),z=perc_values,border="black",opaque.top=TRUE,alpha=0.5,facets=TRUE,xlab="2K-TSS Coverage",ylab="NDR Coverage",zlab="FPKM percentiles",bty="g",space=0.3,shade=FALSE,ticktype = "detailed")
dev.off()

pdf(file="./output/Quantitative/MergedControls_FPKM_3D_BinnedBarplot_Percentiles_masked_2DImage.pdf")
perc_masked<-data$FPKM_percentile
perc_masked[which(data$Count < 11)]<-NA
perc_values<-matrix(perc_masked,nrow=30,byrow=TRUE)
image2D(x=c(0:29),y=c(0:29),z=perc_values, border="black",xlab="2K-TSS Coverage",ylab="NDR Coverage")
dev.off()

#Correlation tests
data<-read.table("./output/Quantitative/MergedControls_FPKM.txt",sep="\t",header=TRUE)
cor.test(data$X2K_Cov,data$FPKM,method="pearson")
cor.test(data$X2K_Cov,data$FPKM,method="spearman")
cor.test(data$X2K_Cov,data$FPKM_percentile,method="pearson")
cor.test(data$X2K_Cov,data$FPKM_percentile,method="spearman")

cor.test(data$NDR_Cov,data$FPKM,method="pearson")
cor.test(data$NDR_Cov,data$FPKM,method="spearman")
cor.test(data$NDR_Cov,data$FPKM_percentile,method="pearson")
cor.test(data$NDR_Cov,data$FPKM_percentile,method="spearman")

library(Cairo)
CairoPDF(file="./output/Quantitative/Correlation_2KTSS_vs_FPKM.pdf")
plot(data$X2K_Cov,data$FPKM,xlab="2K-TSS Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100),pch=20,col=rgb(0,0,0,0.1))
dev.off()
CairoPDF(file="./output/Quantitative/Correlation_2KTSS_vs_FPKMpercentile.pdf")
plot(data$X2K_Cov,data$FPKM_percentile,xlab="2K-TSS Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100),pch=20,col=rgb(0,0,0,0.1))
dev.off()
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKM.pdf")
plot(data$NDR_Cov,data$FPKM,xlab="NDR Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100),pch=20,col=rgb(0,0,0,0.1))
dev.off()
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKMpercentile.pdf")
plot(data$NDR_Cov,data$FPKM_percentile,xlab="NDR Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100),pch=20,col=rgb(0,0,0,0.1))
dev.off()

library(MASS)
kd_estimation<-kde2d(data$X2K_Cov,data$FPKM,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_2KTSS_vs_FPKM_heatmap.pdf")
image(kd_estimation,xlab="2K-TSS Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$X2K_Cov,data$FPKM_percentile,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_2KTSS_vs_FPKMpercentile_heatmap.pdf")
image(kd_estimation,xlab="2K-TSS Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$NDR_Cov,data$FPKM,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKM_heatmap.pdf")
image(kd_estimation,xlab="NDR Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$NDR_Cov,data$FPKM_percentile,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKMpercentile_heatmap.pdf")
image(kd_estimation,xlab="NDR Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100))
dev.off()


kd_estimation<-kde2d(data$X2K_Cov,data$FPKM,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_2KTSS_vs_FPKM_perspective.pdf")
persp(kd_estimation,xlab="2K-TSS Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$X2K_Cov,data$FPKM_percentile,n=1000)
CairoPDF(file="Correlation_2KTSS_vs_FPKMpercentile_perspective.pdf")
persp(kd_estimation,xlab="2K-TSS Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$NDR_Cov,data$FPKM,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKM_perspective.pdf")
persp(kd_estimation,xlab="NDR Coverage",ylab="FPKM",xlim=c(0,30),ylim=c(0,100))
dev.off()
kd_estimation<-kde2d(data$NDR_Cov,data$FPKM_percentile,n=1000)
CairoPDF(file="./output/Quantitative/Correlation_NDR_vs_FPKMpercentile_perspective.pdf")
persp(kd_estimation,xlab="NDR Coverage",ylab="FPKM percentile",xlim=c(0,30),ylim=c(0,100))
dev.off()

