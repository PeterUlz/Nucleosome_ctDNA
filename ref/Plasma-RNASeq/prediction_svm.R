library(e1071)
data<-read.csv("MergedControls_allGenes_prediction.csv",header=TRUE,sep="\t")

top1000<-read.csv("../../../ref/Plasma-RNASeq/Top1000_NMonly.txt",header=FALSE)
bottom1000<-read.csv("../../../ref/Plasma-RNASeq/Bottom1000_NMonly.txt",header=FALSE)

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

top1000<-which(data$Expressed == "Top1000")
bottom1000<-which(data$Expressed == "Bottom1000")
data$Training<-rep(FALSE,length(data$Gene))
set.seed(1234)
training_top1000_index<-sample(top1000,300)
set.seed(1234)
training_bottom1000_index<-sample(bottom1000,300)
training<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[training_top1000_index],data$Broad.TSS.Coverage.Norm[training_bottom1000_index]),Small=c(data$Small.TSS.coverage[training_top1000_index],data$Small.TSS.coverage[training_bottom1000_index]))
training$Type<-rep(0,length(training$Broad))
training$Type[1:300]="Expressed"
training$Type[301:600]="Unexpressed"
training$Type<-as.factor(training$Type)
data$Training[training_top1000_index]<-TRUE
data$Training[training_bottom1000_index]<-TRUE

test_data_top_index<-which(data$Expressed == "Top1000" & data$Training == FALSE)
test_data_bottom_index<-which(data$Expressed == "Bottom1000" & data$Training == FALSE)
test_data<-data.frame(Broad=c(data$Broad.TSS.Coverage.Norm[test_data_top_index],data$Broad.TSS.Coverage.Norm[test_data_bottom_index]),Small=c(data$Small.TSS.coverage[test_data_top_index],data$Small.TSS.coverage[test_data_bottom_index]))
test_data$Type<-rep("Unexpressed",length(test_data$Broad))
test_data$Type[1:length(test_data_top_index)]<-"Expressed"
svm.model<-svm(Type ~ .,data=training,type="C-classification")
svm.pred<-predict(svm.model, test_data)
table(pred = svm.pred, true = test_data$Type)


