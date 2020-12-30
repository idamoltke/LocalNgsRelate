plot.localngsrelate.posterior <- function(x,col=1:3,lwd=2,ylab="Posterier probability",xlab="Position (Mb)",chr=NULL,legend=TRUE,...){
    pos<-x$Pos/1e6
    post0<-x$Post0
    post1<-x$Post1
    post2<-x$Post2
    post = cbind(post0,post1,post2)

    par(mar=c(5.1, 4.1, 2.1, 8.1), xpd=TRUE)
    if(!is.null(chr)){
        pos<-pos[x$Chr==chr]
        post<-post[x$Chr==chr,]
        plot(pos,post[,1],type="l",ylim=c(0,1),col=col[1],lwd=lwd,ylab=ylab,xlab=xlab,las=1,main="Posterior decoding",...)
        lines(pos,post[,2],col=col[2],lwd=lwd)
        lines(pos,post[,3],col=col[3],lwd=lwd)
    }
    else{
        C<-names(table(x$Chr))
        m<-c(0,cumsum(tapply(pos,x$Chr,max)))
        pos2<-rep(NA,length(pos))
        for(tal in 1:length(C))
            pos2[x$Chr==C[tal]]<-pos[x$Chr==C[tal]]+m[tal]
        plot(pos2,post[,1],type="l",ylim=c(-.05,1),col=col[1],lwd=lwd,ylab=ylab,xlab=xlab,las=1,main="Posterior decoding",...)
        lines(pos2,post[,2],col=col[2],lwd=lwd)
        lines(pos2,post[,3],col=col[3],lwd=lwd)
        par(xpd=FALSE)
        abline(v=m)
       
        for(tal in 1:length(C))
            text(m[tal]+diff(m)[tal]/2,-0.05,C[tal],col="gray45")
    }

    if(legend)
        legend("topright",paste("IBD=",2:0),col=rev(1:3),lty=1,inset=c(-0.11,0),xpd=TRUE)
}

plot.localngsrelate.viterbi <- function(x,col=1:3,cex=0.2,ylab="Viterbi inferred IBD state",xlab="Position (Mb)",chr=NULL,legend=TRUE,...){
    pos<-x$Pos/1e6
    vit<-x$Viterbi
    par(mar=c(5.1, 4.1, 2.1, 8.1), xpd=TRUE)
    if(!is.null(chr)){
        pos<-pos[x$Chr==chr]
        vit<-vit[x$Chr==chr]
        plot(pos,y=vit,type="p",ylim=c(0,2),col=(vit+1),cex=cex,ylab=ylab,xlab=xlab,main="Viterbi",las=1,yaxt="n",...)
        axis(2, seq(0,2,1),las=1)
    }
    else{
        C<-names(table(x$Chr))
        m<-c(0,cumsum(tapply(pos,x$Chr,max)))
        pos2<-rep(NA,length(pos))
        for(tal in 1:length(C))
            pos2[x$Chr==C[tal]]<-pos[x$Chr==C[tal]]+m[tal]
        plot(x=pos2,y=vit,col=(vit+1),type="p",ylim=c(-0.1,2),cex=cex,ylab=ylab,xlab=xlab,main="Viterbi",las=1,yaxt="n",...)#ylab=ylab,xlab=xlab,yaxt="n",main="Viterbi",...)
        axis(2,seq(0,2,1),las=1)
        par(xpd=FALSE)
        abline(v=m)
       
        for(tal in 1:length(C))
            text(m[tal]+diff(m)[tal]/2,-0.1,C[tal],col="gray45")
    }

    if(legend)
        legend("topright",paste("IBD=",2:0),col=rev(1:3),lty=1,inset=c(-0.11,0),xpd=TRUE)
}

## Example of use
if(FALSE){
    
    ## Read in results
    nam = "testoutput_a2_b3" # Insert prefix of the name of results file here
    nam = "testoutput_a3_b4" # Insert prefix of the name of results file here
    nam = "testoutput_a0_b1" # Insert prefix of the name of results file here
    
    res = read.table(paste(nam,".IBDtractinference.gz",sep=""),header=T)
    
    ## Plot posterior
    pdf(paste(nam,"_posterior_allchrs.pdf",sep=""),h=4,w=16)
    plot.localngsrelate.posterior(res)
    graphics.off()
    
    pdf(paste(nam,"_posterior_chr1.pdf",sep=""),h=4,w=16)
    plot.localngsrelate.posterior(res,chr=1)
    graphics.off()
    
    ## Plot viterbi
    pdf(paste(nam,"_viterbi_allchrs.pdf",sep=""),h=4,w=16)
    plot.localngsrelate.viterbi(res)
    graphics.off()
    
    pdf(paste(nam,"_viterbi_chr1.pdf",sep=""),h=4,w=16)
    plot.localngsrelate.viterbi(res,chr=1)
    graphics.off()
    
}

