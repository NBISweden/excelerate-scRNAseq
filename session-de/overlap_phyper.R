overlap_phyper<-function(L,bg=length(unique(unlist(L))),with.unique=TRUE,plot=FALSE){
   # takes a list with indices and creates a matrix with overlap between elements
   # can also plot as a heatmap with coloring according to significance in overlap.
   # phyper test uses all entries in L as background if bg is not specified.

   nL<-length(L)
   M<-mat.or.vec(nL,nL)
   P<-mat.or.vec(nL,nL)
   P[,]<-1
   nU<-mat.or.vec(nL,1)
   for (i in 1:nL){
      nU[i]<-length(setdiff(L[[i]],unlist(L[-i])))
      for (j in  i:nL){
        M[i,j]<-length(intersect(L[[i]],L[[j]]))
	P[i,j]<-1-phyper(M[i,j],length(L[[i]]),bg-length(L[[i]]),length(L[[j]]))
	if (i==j) {P[i,j]<-NA}
      }
   }
   if (with.unique){
      M<-cbind(M,nU)
      colnames(M)<-c(names(L),"unique")
      P<-cbind(P,rep(1,length(nU)))
      colnames(P)<-c(names(L),"unique")
   }else {
      colnames(M)<-names(L)
   }

   rownames(M)<-names(L)
   rownames(P)<-names(L)
   if (plot){
      library(gplots)
      lab<-matrix(as.character(M),nrow=nL)
      lab[is.na(lab)]<-''
      par(oma=c(4,1,2,4),xpd=T,cex=0.5,mfrow=c(1,1))
      cexC = 0.2 + 0.2/log10(nL)
      notecex = 0.2 + 0.2/log10(nL)
      h<-heatmap.2(-log10(P+1e-16),cellnote=lab,scale="none",trace="none",density.info="none",notecex=notecex,notecol="black",dendrogram="none",Colv=F,Rowv=F,key=T,cexRow=cexC,cexCol=cexC,key.xlab="-log10(p.value)",key.title='')
   }
   return(list(overlap=M,pval=P))
}
