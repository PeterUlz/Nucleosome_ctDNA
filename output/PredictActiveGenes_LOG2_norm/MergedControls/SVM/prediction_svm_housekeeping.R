library(e1071)
data_broad<-read.csv("../MergedControls_TSSCoverage_Broad_PlasmaRNASeq_NMonly_NormLog2.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data_small<-read.csv("../MergedControls_TSSCoverage_Small_PlasmaRNASeq_NMonly_NormLog2.txt",header=TRUE,sep="\t",stringsAsFactors=FALSE)
data<-data.frame(Gene=data_broad$Gene,TSS=data_broad$TSS,Broad=data_broad$Coverage,Small=data_small$Coverage)
#############################################################################################################################################################################
analyze_and_plot_clustering<-function(data,output_dir) {
    cols<-rep("black",length(data$Type))
    cols[which(data$Clustering == "Expressed")]<-"red"
    cols[which(data$Clustering == "Unexpressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_clustering.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",xlim=c(0,30),ylim=c(0,30))
    dev.off()

    for (xx in which(data$Clustering == "Unexpressed")) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in which(data$Clustering == "Expressed")) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

}
#############################################################################################################################################################################
sensitivity_top5000<-function(data,output_dir){
    #Check in more detail whether Top5000 genes are in expressde cluster and bottom 1000 in unexpressed cluster:
    top5000<-read.csv("../../../../ref/Plasma-RNASeq/Top5000_NMonly.txt",header=FALSE)
    bottom5000<-read.csv("../../../../ref/Plasma-RNASeq/Bottom5000_NMonly.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top5000$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top5000"
    }
    for (bottom in bottom5000$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom5000"
    }
    true_expressed<-which(expressed == "Top5000" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom5000" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom5000" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top5000" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top5000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top5000" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom5000" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top5000_.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top5000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top5000_Highlight.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top5000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom5000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top5000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom5000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top5000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom5000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top5000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom5000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top5000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
sensitivity_top1000<-function(data,output_dir){
    #Check in more detail whether Top1000 genes are in expressde cluster and bottom 1000 in unexpressed cluster:
    top1000<-read.csv("../../../../ref/Plasma-RNASeq/Top1000_NMonly.txt",header=FALSE)
    bottom1000<-read.csv("../../../../ref/Plasma-RNASeq/Bottom1000_NMonly.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top1000$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top1000"
    }
    for (bottom in bottom1000$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom1000"
    }
    true_expressed<-which(expressed == "Top1000" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom1000" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom1000" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top1000" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top1000" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom1000" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top1000.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top1000_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom1000_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top1000.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
checkSpecialGenes<-function(data,output_dir){
    #Check in more detail whether Top1000 genes are in expressed cluster and bottom 1000 in unexpressed cluster:
    special<-c("CCND1","ERBB2","MYC","FGFR1","JAK1")
    special_genes<-rep("unknown",length(data$Gene))
    for (top in special) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        special_genes[index] <- "Special"
    }
    special_index_expressed<-which(special_genes == "Special" &  data$Clustering == "Expressed")
    special_index_unexpressed<-which(special_genes == "Special" & data$Clustering == "Unexpressed")

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Special_genes_Highlight.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[special_index_expressed],data$Small[special_index_expressed],col="red")
    text(data$Broad[special_index_expressed]+1,data$Small[special_index_expressed],labels=data$Gene[special_index_expressed],col="red")
    points(data$Broad[special_index_unexpressed],data$Small[special_index_unexpressed],col="green")
    text(data$Broad[special_index_unexpressed]+1,data$Small[special_index_unexpressed],labels=data$Gene[special_index_unexpressed],col="green")
    dev.off()
}
#############################################################################################################################################################################
cancer_driver_genes<-function(data,output_dir){
    #Check in more detail whether cancer driver genes are called expressed or not:
    vogelstein_driver<-read.csv("../../../../ref/vogelstein_driver_genes.txt",header=FALSE)

    driver_in_expressed_cluster<-0
    driver_in_unexpressed_cluster<-0
    driver<-rep("unknown",length(data$Gene))
    for (top in vogelstein_driver$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        driver[index] <- "Driver"
    }
    driver_expressed<-which(driver == "Driver" & data$Clustering == "Expressed")
    driver_unexpressed<-which(driver == "Driver" & data$Clustering == "Unexpressed")

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[driver_expressed]<-"red"
    cols[driver_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Driver_genes.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[driver_expressed],data$Small[driver_expressed],col=cols[driver_expressed])
    points(data$Broad[driver_unexpressed],data$Small[driver_unexpressed],col=cols[driver_unexpressed])
    dev.off()

    for (xx in driver_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_driver_assigned_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in driver_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_driver_assigned_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Driver in expressed cluster",length(driver_expressed))), file=paste(output_dir,"/Statistics_Vogelstein_Driver.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Driver in unexpressed cluster",length(driver_unexpressed))), file=paste(output_dir,"/Statistics_Vogelstein_Driver.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_housekeeping<-function(data,output_dir){
    #Check in more detail whether Housekeeping genes are assigned into expressed cluster:
    housekeeping<-read.csv("../../../../ref/Housekeeping/HK_gene_names.txt",header=FALSE)

    housekeeping_in_expressed_cluster<-0
    housekeeping_in_unexpressed_cluster<-0
    data$Housekeeping<-rep("unknown",length(data$Gene))
    for (top in housekeeping$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        data$Housekeeping[index] <- "Housekeeping"
    }
    true_expressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering == "Expressed")
    false_unexpressed<-which(data$Housekeeping == "Housekeeping" & data$Clustering == "Unexpressed")

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Housekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_Housekeeping_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_Housekeeping_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Housekeeping (Eisenberg) in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Housekeeping (Eisenberg) in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Houskeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}

#############################################################################################################################################################################
clustering_vs_rma<-function(data,output_dir){
    #check RMA values for genes and plot "Heatmap"
    data$rma<-rep(0,length(data$Gene))
    rma_values<-read.csv("../../../../ref/Plasma-RNASeq/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",header=TRUE,sep="\t")
    for (gene in rma_values$Gene) {
        index<-which(data$Gene == gene)
        rma_index<-which(rma_values$Gene == gene)
        if (length(index) == 0) {
            next
        }
        if (length(rma_index) == 1) {
            rma_value <- rma_values$Mean[rma_index]
        }
        else  {
            rma_value <- mean(rma_values$Mean[rma_index])
        }
        data$rma[index] <- rep(rma_value,length(index))
    }
    colors<-colorRampPalette(c("green","red"))(3)
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    level<-which(data$rma > 7)
    cols[level]<-colors[3]
    level<-which(data$rma < 4)
    cols[level]<-colors[1]
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Heatmap.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()

    #compare distance to expressed cluster center to RMA
    colors<-colorRampPalette(c("green","red"))(max(floor(data$rma)))
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    for (i in min(floor(data$rma)):max(floor(data$rma))) {
       level<-which(floor(data$rma) == i)
       cols[level]<-colors[i]
    }
}

#############################################################################################################################################################################
clustering_vs_fpkm<-function(data,output_dir) {
    #check FPKM values for genes and plot "Heatmap"
    data$fpkm<-rep(0,length(data$Gene))
    fpkm_values<-read.csv("../../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
    for (gene in fpkm_values$tracking_id) {
        index<-which(data$Gene == gene)
        fpkm_index<-which(fpkm_values$tracking_id == gene)
        if (length(index) == 0) {
            next
        }
        if (length(fpkm_index) == 1) {
            fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
        }
        else  {
            fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
        }
        data$fpkm[index] <- rep(fpkm_value,length(index))
    }

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Boxplot.png",sep=""))
    boxplot(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
    dev.off()
    true_expressed<-which(data$fpkm > 8 & data$Clustering == "Expressed")
    true_unexpressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering == "Unexpressed")
    false_expressed<-which(data$fpkm > 0 & data$fpkm < 0.1 & data$Clustering ==  "Expressed")
    false_unexpressed<-which(data$fpkm > 8 & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-rgb(0.4,0,0.77)
    cols[true_unexpressed]<-rgb(0.4,0.6,1)
    cols[false_unexpressed]<-rgb(1,0,0)
    cols[false_expressed]<-rgb(1,0,0)

    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Heatmap.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Highlight.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],col=cols[false_expressed])
    points(data$Broad[true_expressed],data$Small[true_expressed],col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],col=cols[true_unexpressed])
    dev.off()


    write.table(t(c("Top1000 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom1000 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top1000 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom1000 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_FPKM.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_top100<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../../ref/Plasma-RNASeq/Top100_NMonly.txt",header=FALSE)
    bottom100<-read.csv("../../../../ref/Plasma-RNASeq/Bottom100_NMonly.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    for (bottom in bottom100$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom100" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top100" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom100" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom100 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom100 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top100.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)

}

#############################################################################################################################################################################
sensitivity_top100_amplified_B7<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Plasma/B7_Top100_1q_PlasmaLower1.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100_B7.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100_B7.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100_B7.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_B7_Top100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_B7_Top100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_top100_amplified_B13<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Plasma/B13_Top100_8pq_PlasmaLower1.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100_B13.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100_B13.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100_B13.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    dev.off()

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_B13_Top100_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_B13_Top100_in_unexpressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_top100_amplified_noHousekeeping_B13<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Housekeeping/B13_Top100_8pq_noHousekeeping.txt",header=FALSE)
    bottom100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/B13_Bottom100.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    for (bottom in bottom100$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom100" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100_8pq_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top100" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom100" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top100_8pq_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100_8pq_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100_8pq_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

	data$fpkm<-rep(0,length(data$Gene))
    fpkm_values<-read.csv("../../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
	for (gene in fpkm_values$tracking_id) {
		index<-which(data$Gene == gene)
		fpkm_index<-which(fpkm_values$tracking_id == gene)
		if (length(index) == 0) {
		    next
		}
		if (length(fpkm_index) == 1) {
		    fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
		}
		else  {
		    fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
		}
		    data$fpkm[index] <- rep(fpkm_value,length(index))
	}

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_expressed_cluster_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_unexpressed_cluster_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_unexpressed_cluster_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_expressed_cluster_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    fpkm_true_expressed<-mean(data$fpkm[true_expressed])
    fpkm_false_unexpressed<-mean(data$fpkm[false_unexpressed])
    fpkm_true_expressed_min<-min(data$fpkm[true_expressed])
    fpkm_false_unexpressed_min<-min(data$fpkm[false_unexpressed])
    fpkm_true_expressed_max<-max(data$fpkm[true_expressed])
    fpkm_false_unexpressed_max<-max(data$fpkm[false_unexpressed])

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom100 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom100 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Mean FPKM Top100 expressed cluster",fpkm_true_expressed)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Mean FPKM Top100 unexpressed cluster",fpkm_false_unexpressed)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Min FPKM Top100 expressed cluster",fpkm_true_expressed_min)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Min FPKM Top100 unexpressed cluster",fpkm_false_unexpressed_min)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Max FPKM Top100 expressed cluster",fpkm_true_expressed_max)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Max FPKM Top100 unexpressed cluster",fpkm_false_unexpressed_max)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top100_8pq_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)

}
#############################################################################################################################################################################
sensitivity_top100_amplified_noHousekeeping_B7<-function(data,output_dir){
    #Check if snesitivity is higher for Top100 vs Bottom 100:
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Housekeeping/B7_Top100_1q_noHousekeeping.txt",header=FALSE)
    bottom100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/B7_Bottom100.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    for (bottom in bottom100$V1) {
        index<-which(data$Gene == bottom)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Bottom100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    true_unexpressed<-which(expressed == "Bottom100" & data$Clustering == "Unexpressed")
    false_expressed<-which(expressed == "Bottom100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))
    accuracy<-(length(true_expressed)+length(true_unexpressed) )/ (length(true_expressed)+length(false_unexpressed)+length(true_unexpressed)+length(false_expressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_True_Top100_1pq_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[which(expressed == "Top100" & data$Clustering == "Unexpressed")]<-"red"
    cols[which(expressed == "Bottom100" & data$Clustering == "Expressed")]<-"green"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_Wrong_Top100_1q_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Top100_1q_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[true_unexpressed]<-"green"
    cols[false_expressed]<-"blue"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/MergedControls_Broad_vs_Small_TSS_coverage_All_Highlight_Top100_1q_noHousekeeping.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,10),ylim=c(0,10))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[true_unexpressed],data$Small[true_unexpressed],,col=cols[true_unexpressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    points(data$Broad[false_expressed],data$Small[false_expressed],,col=cols[false_expressed])
    dev.off()

	data$fpkm<-rep(0,length(data$Gene))
    fpkm_values<-read.csv("../../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
	for (gene in fpkm_values$tracking_id) {
		index<-which(data$Gene == gene)
		fpkm_index<-which(fpkm_values$tracking_id == gene)
		if (length(index) == 0) {
		    next
		}
		if (length(fpkm_index) == 1) {
		    fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
		}
		else  {
		    fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
		}
		    data$fpkm[index] <- rep(fpkm_value,length(index))
	}

    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_expressed_cluster_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in true_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_unexpressed_cluster_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Top100_in_unexpressed_cluster_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/MergedControls_assigned_in_Bottom100_in_expressed_cluster_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }


    fpkm_true_expressed<-mean(data$fpkm[true_expressed])
    fpkm_false_unexpressed<-mean(data$fpkm[false_unexpressed])
    fpkm_true_expressed_min<-min(data$fpkm[true_expressed])
    fpkm_false_unexpressed_min<-min(data$fpkm[false_unexpressed])
    fpkm_true_expressed_max<-max(data$fpkm[true_expressed])
    fpkm_false_unexpressed_max<-max(data$fpkm[false_unexpressed])

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Bottom100 in unexpressed cluster",length(true_unexpressed))), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Bottom100 in expressed cluster",length(false_expressed))), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Mean FPKM Top100 expressed cluster",fpkm_true_expressed)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Mean FPKM Top100 unexpressed cluster",fpkm_false_unexpressed)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Min FPKM Top100 expressed cluster",fpkm_true_expressed_min)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Min FPKM Top100 unexpressed cluster",fpkm_false_unexpressed_min)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Max FPKM Top100 expressed cluster",fpkm_true_expressed_max)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Max FPKM Top100 unexpressed cluster",fpkm_false_unexpressed_max)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Accuracy",accuracy)), file=paste(output_dir,"/Statistics_Top100_1q_noHousekeeping.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)

}
#############################################################################################################################################################################
sensitivity_log2over1_notPlasma_B7<-function(data,output_dir){
    #Check if sensitivity is high for genes not expressed in plasma RNAseq but in Tumor RNASeq of B7
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Plasma/B7_Log2_over1_FPKMover5_PlasmaLower01.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B7_Broad_vs_Small_TSS_coverage_True_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B7_Broad_vs_Small_TSS_coverage_All_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B7_Broad_vs_Small_TSS_coverage_All_Highlight_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    dev.off()


    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/B7_assigned_in_B7_log2over1_notPlasma_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/B7_assigned_in_B7_log2over1_notPlasma_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B7.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
#############################################################################################################################################################################
sensitivity_log2over1_notPlasma_B13<-function(data,output_dir){
    #Check if sensitivity is high for genes not expressed in plasma RNAseq but in Tumor RNASeq of B7
    top100<-read.csv("../../../../ref/Tumor_RNASeq/10mio_reads/Tumor_minus_Plasma/B13_Log2_over1_FPKMover5_PlasmaLower01.txt",header=FALSE)

    top_in_expressed_cluster<-0
    top_in_unexpressed_cluster<-0
    bottom_in_expressed_cluster<-0
    bottom_in_unexpressed_cluster<-0
    expressed<-rep("unknown",length(data$Gene))
    for (top in top100$V1) {
        index<-which(data$Gene == top)
        if (length(index) == 0) {
            next
        }
        expressed[index] <- "Top100"
    }
    true_expressed<-which(expressed == "Top100" & data$Clustering == "Expressed")
    false_unexpressed<-which(expressed == "Top100" & data$Clustering == "Unexpressed")

    sensitivity<-length(true_expressed) / (length(true_expressed)+length(false_unexpressed))

    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B13_Broad_vs_Small_TSS_coverage_True_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",,col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B13_Broad_vs_Small_TSS_coverage_All_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    dev.off()
    cols<-rep(rgb(0,0,0,0.2),length(data$Gene))
    cols[true_expressed]<-"red"
    cols[false_unexpressed]<-"blue"
    png(paste(output_dir,"/B13_Broad_vs_Small_TSS_coverage_All_Highlight_log2over1_notPlasma.png",sep=""))
    plot(data$Broad,data$Small,xlab="Broad TSS Coverage",ylab="Small TSS Coverage",col=cols,xlim=c(0,30),ylim=c(0,30))
    points(data$Broad[true_expressed],data$Small[true_expressed],,col=cols[true_expressed])
    points(data$Broad[false_unexpressed],data$Small[false_unexpressed],col=cols[false_unexpressed])
    dev.off()


    for (xx in true_expressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/B13_assigned_in_B13_log2over1_notPlasma_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    for (xx in false_unexpressed) {
        write.table(t(c(as.character(data$Gene[xx]),as.character(data$TSS[xx]))), file=paste(output_dir,"/B13_assigned_in_B13_log2over1_notPlasma_in_expressed_cluster.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    }

    write.table(t(c("Top100 in expressed cluster",length(true_expressed))), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")
    write.table(t(c("Top100 in unexpressed cluster",length(false_unexpressed))), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
    write.table(t(c("Sensitivity",sensitivity)), file=paste(output_dir,"/Statistics_log2over1_notPlasma_B13.txt",sep=""), quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t",append=TRUE)
}
##############################################################################################
#use housekeeping and unexpressed genes for evaluation
top1000<-read.csv("../../../../ref/Housekeeping/HK_gene_names.txt",header=FALSE)
bottom1000<-read.csv("../../../../ref/FANTOM5/Fantom5_all_lower0.1.txt",header=FALSE)

top_in_expressed_cluster<-0
top_in_unexpressed_cluster<-0
bottom_in_expressed_cluster<-0
bottom_in_unexpressed_cluster<-0
data$Expressed<-rep("unknown",length(data$Gene))
for (top in top1000$V1) {
    index<-which(data$Gene == top)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Top1000"
}
for (bottom in bottom1000$V1) {
    index<-which(data$Gene == bottom)
    if (length(index) == 0) {
        next
    }
    data$Expressed[index] <- "Bottom1000"
}

data$Prediction_Expr<-rep(0,length(data$Gene))
data$Prediction_Unexpr<-rep(0,length(data$Gene))
for (i in 1:1000) {
    set.seed(i)
    data$Training<-rep(FALSE,length(data$Gene))
    training_top_index<-sample(which(data$Expressed == "Top1000"),300)
    training_bottom_index<-sample(which(data$Expressed == "Bottom1000"),300)
    data$Training[training_top_index]<-TRUE
    data$Training[training_bottom_index]<-TRUE
    training_data<-data.frame(Broad=c(data$Broad[training_top_index],data$Broad[training_bottom_index]),Small=c(data$Small[training_top_index],data$Small[training_bottom_index]))
    type<-rep("Unexpressed",length(training_data$Broad))
    type[1:300]<-"Expressed"
    training_data$Type<-as.factor(type)
    svm_model<-svm(Type ~ .,data=training_data,type="C-classification")

    test_data_top_index<-which(data$Training == FALSE)
    test_data_bottom_index<-which(data$Training == FALSE)
    test_data<-data.frame(Broad=c(data$Broad[test_data_top_index],data$Broad[test_data_bottom_index]),Small=c(data$Small[test_data_top_index],data$Small[test_data_bottom_index]),Index=c(test_data_top_index,test_data_bottom_index))
    type<-rep("Unexpressed",length(test_data$Broad))
    type[1:length(test_data_top_index)]<-"Expressed"
    test_data$Type<-as.factor(type)
    prediction<-predict(svm_model,test_data)
    data$Prediction_Expr[test_data$Index[which(prediction == "Expressed")]]<-data$Prediction_Expr[test_data$Index[which(prediction == "Expressed")]]+1
    data$Prediction_Unexpr[test_data$Index[which(prediction == "Unexpressed")]]<-data$Prediction_Unexpr[test_data$Index[which(prediction == "Unexpressed")]]+1
}

true_expressed<-which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Top1000")
false_unexpressed<-which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Top1000")
true_unexpressed<-which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Bottom1000")
false_expressed<-which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Bottom1000")


col=rep(rgb(0,0,0,0.2),length(data$Gene))
col[true_expressed]<-"red"
col[true_unexpressed]<-"green"
col[false_unexpressed]<-"blue"
col[false_expressed]<-"blue"
png("/Housekeeping/MergedControls_Prediction_Top1000.png")
plot(data$Broad,data$Small,col=col,xlim=c(0,30),ylim=c(0,30),xlab="Broad Coverage",ylab="Small Coverage")
dev.off()

length(which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Top1000"))
length(which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Top1000"))
length(which(data$Prediction_Expr < data$Prediction_Unexpr & data$Expressed == "Bottom1000"))
length(which(data$Prediction_Expr > data$Prediction_Unexpr & data$Expressed == "Bottom1000"))

data$Clustering<-rep("Unexpressed",length(data$Gene))
expressed_index<-which(data$Prediction_Expr > data$Prediction_Unexpr)
data$Clustering[expressed_index]<-"Expressed"

#analyze_and_plot_clustering(data,".")
#sensitivity_top5000(data,"./Housekeeping/Majority")
#sensitivity_top1000(data,"./Housekeeping/Majority")
#sensitivity_top100(data,"./Housekeeping/Majority")
#sensitivity_top100_amplified_B7(data,"./Housekeeping/Majority")
#sensitivity_top100_amplified_B13(data,"./Housekeeping/Majority")
#checkSpecialGenes(data,"./Housekeeping/Majority")
#cancer_driver_genes(data,"./Housekeeping/Majority")
#sensitivity_housekeeping(data,"./Housekeeping/Majority")
#clustering_vs_rma(data,"./Housekeeping/Majority")
#clustering_vs_fpkm(data,"./Housekeeping/Majority")

#Check up whether predicted expressed genes have higher RMA
#data$rma<-rep(0,length(data$Gene))
#rma_values<-read.csv("../../../../ref/Plasma-RNASeq/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",header=TRUE,sep="\t")
#for (gene in rma_values$Gene) {
#    index<-which(data$Gene == gene)
#    rma_index<-which(rma_values$Gene == gene)
#    if (length(index) == 0) {
#        next
#    }
#    if (length(rma_index) == 1) {
#        rma_value <- rma_values$Mean[rma_index]
#    }
#    else  {
#        rma_value <- mean(rma_values$Mean[rma_index])
#    }
#    data$rma[index] <- rep(rma_value,length(index))
#}
#w_maj_rma<-wilcox.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
#t_maj_rma<-t.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
#med_expressed_maj_rma<-median(data$rma[which(data$Clustering=="Expressed")])
#med_unexpressed_maj_rma<-median(data$rma[which(data$Clustering=="Unexpressed")])

#png("./Housekeeping/Majority/MergedControls_Broad_vs_Small_TSS_coverage_RMA_Boxplot.png")
#boxplot(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")],names=c("Expressed","Unexpressed"),col=c("red","green"))
#dev.off()


#Check up whether expressed genes have higher FPKM
#data$fpkm<-rep(0,length(data$Gene))
#fpkm_values<-read.csv("../../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
#for (gene in fpkm_values$tracking_id) {
#    index<-which(data$Gene == gene)
#    fpkm_index<-which(fpkm_values$tracking_id == gene)
#    if (length(index) == 0) {
#        next
#    }
#    if (length(fpkm_index) == 1) {
#        fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
#    }
#    else  {
#        fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
#    }
#        data$fpkm[index] <- rep(fpkm_value,length(index))
#}
#w_maj_fpkm<-wilcox.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
#t_maj_fpkm<-t.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
#med_expressed_maj_fpkm<-median(data$fpkm[which(data$Clustering=="Expressed")])
#med_unexpressed_maj_fpkm<-median(data$fpkm[which(data$Clustering=="Unexpressed")])##

#png("./Housekeeping/Majority/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Boxplot.png")
#boxplot(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")],ylim=c(0,100),names=c("Expressed","Unexpressed"),col=c("red","green"))
#dev.off()


#############################################################################################
# Only analyze genes which are consistently called in >75% of the SVM runs
expr_consistent<-which(data$Prediction_Expr > 0.75 * (data$Prediction_Expr+data$Prediction_Unexpr)) 
unexpr_consistent<-which(data$Prediction_Unexpr > 0.75 * (data$Prediction_Expr+data$Prediction_Unexpr)) 
data$Clustering<-rep("unknown",length(data$Gene))
data$Clustering[expr_consistent]<-"Expressed"
data$Clustering[unexpr_consistent]<-"Unexpressed"

col=rep(rgb(0,0,0,0.2),length(data$Gene))
col[expr_consistent]<-"red"
col[unexpr_consistent]<-"green"
png("./Housekeeping/Consistency/Consistent_calls.png")
plot(data$Broad,data$Small,col=col,xlim=c(0,30),ylim=c(0,30))
dev.off()

col=rep("blue",length(data$Gene))
col[expr_consistent]<-rgb(0,0,0,0.2)
col[unexpr_consistent]<-rgb(0,0,0,0.2)
png("./Housekeeping/Consistency/Inconsistent_calls.png")
plot(data$Broad,data$Small,col=col,xlim=c(0,30),ylim=c(0,30))
dev.off()

analyze_and_plot_clustering(data,"./Housekeeping/Consistency")
sensitivity_top5000(data,"./Housekeeping/Consistency")
sensitivity_top1000(data,"./Housekeeping/Consistency")
sensitivity_top100(data,"./Housekeeping/Consistency")
sensitivity_top100_amplified_B7(data,"./Housekeeping/Consistency")
sensitivity_top100_amplified_B13(data,"./Housekeeping/Consistency")
sensitivity_top100_amplified_noHousekeeping_B7(data,"./Housekeeping/Consistency")
sensitivity_top100_amplified_noHousekeeping_B13(data,"./Housekeeping/Consistency")
sensitivity_log2over1_notPlasma_B7(data,"./Housekeeping/Consistency")
sensitivity_log2over1_notPlasma_B13(data,"./Housekeeping/Consistency")
checkSpecialGenes(data,"./Housekeeping/Consistency")
cancer_driver_genes(data,"./Housekeeping/Consistency")
sensitivity_housekeeping(data,"./Housekeeping/Consistency")
clustering_vs_rma(data,"./Housekeeping/Consistency")
clustering_vs_fpkm(data,"./Housekeeping/Consistency")

# Check distribution of Top1000 and Bottom1000 Genes
top_bottom<-which((data$Expressed == "Top1000" | data$Expressed == "Bottom1000") & (data$Broad < 25) & (data$Small < 20))
top<-which(data$Expressed[top_bottom] == "Top1000")
bottom<-which(data$Expressed[top_bottom] == "Bottom1000")
col=rep("black",length(data$Broad[top_bottom]))
col[top]<-"red"
col[bottom]<-"green"
png("./Housekeeping/KDE/distribution_top_and_bottom1000.png")
plot(data$Broad[top_bottom],data$Small[top_bottom],col=col,xlim=c(0,25),ylim=c(0,20))
dev.off()
#Kernel density estimation
library(MASS)
kd_estimation<-kde2d(data$Broad[top_bottom],data$Small[top_bottom],n=50)

png("./Housekeeping/KDE/KDE_Contour.png") 
contour(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,25),ylim=c(0,20))
dev.off()

png("./Housekeeping/KDE/KDE_Contour_data.png") 
plot(data_numeric)
contour(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",add=TRUE,col="green",xlim=c(0,25),ylim=c(0,20))
dev.off()
png("./Housekeeping/KDE/Heatmap.png") 
image(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,25),ylim=c(0,20))
dev.off()
png("./Housekeeping/KDE/Perspective.png") 
persp(kd_estimation, phi = 45, theta = 30,xlab="Broad Normalized Coverage",ylab="Small Coverage",zlab="Density",xlim=c(0,25),ylim=c(0,20))
dev.off()

png("./Housekeeping/KDE/Heatmap_contour.png")
# fancy contour with image
image(kd_estimation,xlab="Broad Normalized Coverage",ylab="Small Coverage",xlim=c(0,25),ylim=c(0,20)); contour(kd_estimation, add = T)
dev.off()

# fancy perspective
png("./Housekeeping/KDE/Perspective_angle.png")
persp(kd_estimation, phi = 30, theta = 35, shade = .1, border = NA,xlab="Broad Normalized Coverage",ylab="Small Coverage",zlab="Density")
dev.off()



#Check up whether expressed genes have higher FPKM
data$fpkm<-rep(0,length(data$Gene))
fpkm_values<-read.csv("../../../../ref/Plasma-RNASeq/FPKM/filtered.merged.genes.fpkm_tracking",header=TRUE,sep="\t")
for (gene in fpkm_values$tracking_id) {
    index<-which(data$Gene == gene)
    fpkm_index<-which(fpkm_values$tracking_id == gene)
    if (length(index) == 0) {
        next
    }
    if (length(fpkm_index) == 1) {
        fpkm_value <- fpkm_values$Mean.FPKM[fpkm_index]
    }
    else  {
        fpkm_value <- mean(fpkm_values$Mean.FPKM[fpkm_index])
    }
        data$fpkm[index] <- rep(fpkm_value,length(index))
}
w_cons<-wilcox.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
t_cons<-t.test(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")])
med_expressed_cons<-median(data$fpkm[which(data$Clustering=="Expressed")])
med_unexpressed_cons<-median(data$fpkm[which(data$Clustering=="Unxpressed")])
sd_expressed_cons<-sd(data$fpkm[which(data$Clustering=="Expressed")])
sd_unexpressed_cons<-sd(data$fpkm[which(data$Clustering=="Unxpressed")])

png("./Housekeeping/Consistency/MergedControls_Broad_vs_Small_TSS_coverage_FPKM_Boxplot.png")
boxplot(data$fpkm[which(data$Clustering=="Expressed")],data$fpkm[which(data$Clustering=="Unexpressed")],ylim=c(0,100),names=c("Expressed","Unexpressed"),col=c("red","green"))
dev.off()

#Check up whether predicted expressed genes have higher RMA
data$rma<-rep(0,length(data$Gene))
rma_values<-read.csv("../../../../ref/Plasma-RNASeq/NonPregnant_annotated_noChrM_dedup_onlyNM.txt.csv",header=TRUE,sep="\t")
for (gene in rma_values$Gene) {
    index<-which(data$Gene == gene)
    rma_index<-which(rma_values$Gene == gene)
    if (length(index) == 0) {
        next
    }
    if (length(rma_index) == 1) {
        rma_value <- rma_values$Mean[rma_index]
    }
    else  {
        rma_value <- mean(rma_values$Mean[rma_index])
    }
    data$rma[index] <- rep(rma_value,length(index))
}
w_cons_rma<-wilcox.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
t_cons_rma<-t.test(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")])
med_expressed_cons_rma<-median(data$rma[which(data$Clustering=="Expressed")])
med_unexpressed_cons_rma<-median(data$rma[which(data$Clustering=="Unexpressed")])

png("./Housekeeping/Consistency/MergedControls_Broad_vs_Small_TSS_coverage_RMA_Boxplot.png")
boxplot(data$rma[which(data$Clustering=="Expressed")],data$rma[which(data$Clustering=="Unexpressed")],names=c("Expressed","Unexpressed"),col=c("red","green"))
dev.off()

