## version 0.2
# source("https://bioconductor.org/biocLite.R") ͨ��bioconductor����������Ӧ��R��
# biocLite("edgeR")
args = commandArgs(TRUE)
a=length(args)
if(a<4){
	cat("Usage: Rscript edgeR.r <count.xls> <col_num1,col_num2> <outpfx> <dispersion>\n")
	q()
}
## input prepare
infile = args[1]
cols=as.numeric((unlist(strsplit(args[2],','))))
group_colnames=(unlist(strsplit(args[3],',')))
#cols
allcols=sum(cols)
outpfx = args[4]
dispersion = args[5]  # û���ظ�ʱʹ�ã���Ϊ�趨����ѧ����ϵ����edgeR���飺����ʵ��0.4���Ŵ��������Ƶ�ģʽ����0.1�������ظ�0.01 
pq = as.numeric(args[6])  ## 1=> p , 2=>q
cutp = as.numeric(args[7])
cutfc = as.numeric(args[8])
cutfc = log2(cutfc)
norm = as.character(args[9])
fpkm = as.character(args[10])

#cutfc = 1 ## log2(fc) filter cutoff
#cutp = 0.05 ## pvalue filter cutoff

pq1=c("PValue","FDR")
pq=pq1[pq]
#pq
##1. �����ļ���������У������
library(edgeR,quietly=TRUE)   
data = read.table(args[1], header=T, row.names=1, sep='\t',check.names = F,quote = "",comment.char = "")
x = data[, 1:ncol(data)] #���������У�
all_list = DGEList(counts=x, group=colnames(x)) # ����ÿ���������������
all_list <- calcNormFactors(all_list,method="TMM")    #����������Ԥ��У������
#all_list  # �鿴��y��������ݣ�ע��У�����ӵ���ֵ�� 
#plotMDS(all_list)  # ʹ��mutlidimensional scaling�ķ������о��࣬MDS���㷨����PCA������Ч��

##2. ׼�����ݣ������в������
C_number=cols[1];   # ������������
T_number=cols[2];   # ������������
col_ordering = 1:allcols # �趨 ���ڲ�������������ı�ţ�����ѡ��A���C��
#head(data)
DiffMatrix = data[,col_ordering] # ���ܱ���ɸѡ���ڲ������������
## filter row and num
filter = 0.001
f1=(rowSums(DiffMatrix)>filter)
f2=(rowSums(DiffMatrix)<=filter)
f1n=nrow(DiffMatrix[rowSums(DiffMatrix)>filter,])
f2n=nrow(DiffMatrix[rowSums(DiffMatrix)<=filter,])
#f2n
f2name = rownames(DiffMatrix[rowSums(DiffMatrix)<=filter,])

DiffMatrix = DiffMatrix[rowSums(DiffMatrix)>filter,] # ֻ��ѡreadcount��ֵ����1�Ļ��򣬼�ȥ�������Ļ���;
conditions = factor(c(rep("control", C_number), rep("treat", T_number))) # �趨������Ϣ���ƣ�ע��������������ֵ����ǰ�趨�õ�
Diff_list = DGEList(counts=DiffMatrix, group=conditions) # ���������Ƚ���Ĳ��������
#Diff_list$samples[3]=rep(1,ncol(DiffMatrix))  # У�����Ӷ�ʹ��Ĭ�ϵ�1
# # ����2�е�ѭ��������ǽ�У�������滻ΪedgeR�������У�����ӣ� 
Diff_list = calcNormFactors(Diff_list,method="TMM")
if(norm == "no"){
	Diff_list$samples[3]=rep(1,ncol(DiffMatrix))
}
Diff_list$samples[3]
# name= rownames(Diff_list$samples)
# for (i in 1:length(name)){Diff_list$samples[name[i],3]= all_list$samples[name[i],3]}

#if (dispersion == 'no' & (C_number > 1 || T_number > 1)){             #��Ĺ�ϵ��ֻҪ��һ�����ظ���
if(C_number > 1 && T_number > 1){
    Diff_list = estimateCommonDisp(Diff_list)  # ���������ɢ��
    Diff_list = estimateTagwiseDisp(Diff_list) # ����ÿ���������ɢ��
#    plotBCV(Diff_list)
    etest = exactTest(Diff_list)                  # ��ȷ����
	cat("group\n")
}else{
    etest = exactTest(Diff_list, pair = c("control", "treat"), dispersion=as.numeric(dispersion)^2)  # ѡ�������ƣ���Ϊ�趨bcv
#	cat("dispersion\n")
}
tTags = topTags(etest,n=NULL) # �������,������FDR


#head(data)

## calculate cpm
data2=data

if(fpkm == "none"){
	for(i in 1:allcols)
	{
		data2[,i] = round(data[,i]*1000000/sum(data[,i]),2)
		colnames(data2)[i] = paste(colnames(data)[i],"_CPM",sep="")
	}
}else{
	data2 = read.table(fpkm, header=T, row.names=1, sep='\t',check.names = F,quote = "",comment.char = "")
}
#head(data2)

## mean 
data3=data2[,1:2]
data3[,1] = apply(data.frame(data2[,1:C_number]), 1, mean )
data3[,2] = apply(data.frame(data2[,(C_number+1):allcols]), 1, mean)
data3[which(data3[,1]==0),1]=0.001
data3[which(data3[,2]==0),2]=0.001
log2fc=log2(data3[,2]/data3[,1])
colnames(data3)=paste(group_colnames,"_mean",sep="")
#head(data3)

## ouput file, id,count,log2(fc),pvalue,qvalue
data4=data.frame(rownames(data),data,data2,data3,check.names=FALSE)
colnames(data4)[1] = "id"
data4$'log2(fc)' = log2fc
datat = tTags@.Data[[1]][,3:4]
datat$id = rownames(datat)
if(f2n != 0){
	n = nrow(datat)
	datat[(n+1):(n+f2n),1:3]=rep(1,f2n)
	datat[(n+1):(n+f2n),2]=rep(1,f2n)
	datat[(n+1):(n+f2n),3]=f2name
}

data5 = merge(data4,datat,by="id")
#head(data5)

#write.table(data2, file=paste(outpfx,".edgeR.cpm.xls",sep=""), sep='\t', quote=F, row.names=T)  # 写出数据
#write.table(data3, file=paste(outpfx,".edgeR.mean.xls",sep=""), sep='\t', quote=F, row.names=T)  # 写出数据
#write.table(tTags, file=paste(outpfx,".edgeR.DE_results.xls",sep=""), sep='\t', quote=F, row.names=T)  # 写出数据
write.table(data5, file=paste(outpfx,".edgeR.all.xls",sep=""), sep='\t', quote=F, row.names=F)  # д������
data6 = data5[which(abs(data5[,'log2(fc)'])>cutfc & data5[,pq]<cutp),]
#head(data6)
write.table(data6, file=paste(outpfx,".edgeR.filter.xls",sep=""), sep='\t', quote=F, row.names=F)

1
## 3. ��ͼ 
#summary(de <- decideTestsDGE(etest,p.value=0.05,lfc=1)) # �趨������� FDR �� log2FC����ֵ��lfc������ֻ��R 3.2.3�İ汾��Ч��
#detags <- rownames(Diff_list)[as.logical(de)]  # �趨һ���Ƿ�����������߼��ж�����
#plotSmear(etest, de.tags=detags) # ���ƻ�ɽͼ
#abline(h=c(-1, 1), col="blue")   # ������
