#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly=TRUE)
out <- normalizePath(args[[1]])
d <- fread(file.path(out,'selected_gene_features.tsv'), na.strings=c('','NA'))
nums <- c('transcript_length','CDS_length','UTR5_length','UTR3_length','GC_content','A_content',
          'half_life_min','published_expression','sum_editing','CDS_A_count')
for (nm in nums) set(d,j=nm,value=as.numeric(d[[nm]]))
d[, edited := as.integer(editing_binary_status=='Edited')]
d[, decay_class := factor(decay_class, levels=c('II','III','IV','V','I'))]
d[, group := ifelse(edited==1,'Edited','Unedited')]
d[, log2_half_life := log2(half_life_min)]
z <- function(x) as.numeric((x-mean(x,na.rm=TRUE))/sd(x,na.rm=TRUE))
d[, `:=`(z_log2_transcript_length=z(log2(transcript_length)),
    z_log2_CDS_length=z(log2(CDS_length)), z_log2_UTR5_length=z(log2(UTR5_length+1)),
    z_log2_UTR3_length=z(log2(UTR3_length+1)),
    z_log2_published_expression=z(log2(published_expression+1)),
    z_GC_content=z(GC_content), z_A_content=z(A_content))]
decay <- d[decay_class %in% c('II','III','IV','V') & is.finite(log2_half_life)]
tab <- table(factor(d$edited,levels=c(1,0)),factor(d$decay_class!='I',levels=c(TRUE,FALSE)))
ft <- fisher.test(tab,alternative='greater')

violin <- function(values,lo,hi) {
    plotted <- log10(values[values>=lo & values<=hi])
    den <- density(plotted,n=512)
    q <- quantile(plotted,c(.25,.5,.75))
    iqr <- q[3]-q[1]
    whisk <- range(plotted[plotted>=q[1]-1.5*iqr & plotted<=q[3]+1.5*iqr])
    list(n=length(values), mean=mean(values), sd=sd(values),
         density_y=den$x,density=den$y,q1=unname(q[1]),median=unname(q[2]),q3=unname(q[3]),
         whisk_lo=whisk[1],whisk_hi=whisk[2])
}
a <- lapply(c('Edited','Unedited'),function(g) c(list(group=g),violin(d[group==g,transcript_length],200,100000)))
cc <- lapply(c('Edited','Unedited'),function(g) c(list(group=g),violin(decay[group==g,half_life_min],5,50000)))
forms <- list(H3=log2_half_life~edited+z_log2_transcript_length+z_log2_published_expression+decay_class,
 H4=log2_half_life~edited+z_log2_CDS_length+z_log2_UTR5_length+z_log2_UTR3_length+z_log2_published_expression+decay_class,
 H5=log2_half_life~edited+z_log2_transcript_length+z_log2_published_expression+z_GC_content+z_A_content+decay_class)
models <- lapply(forms,function(f)lm(f,data=decay))
reg <- rbindlist(lapply(names(models),function(nm){
 m<-models[[nm]];s<-summary(m)$coefficients
 data.table(model=nm,variable=rownames(s),estimate=s[,1],SE=s[,2],P=s[,4],
            lower=s[,1]-qnorm(.975)*s[,2],upper=s[,1]+qnorm(.975)*s[,2],N=nobs(m),
            N_edited=sum(model.frame(m)$edited==1),N_unedited=sum(model.frame(m)$edited==0))
}))
e <- decay[edited==1 & sum_editing>0 & CDS_A_count>0]
e[, editing_level:=sum_editing/CDS_A_count]
e[, x:=log2(editing_level)]
e[, y:=log2_half_life]
cor <- cor.test(e$editing_level,e$half_life_min,method='spearman',exact=FALSE)
lmfit <- lm(y~x,data=e)
x <- seq(min(e$x),max(e$x),length.out=101)
pred <- predict(lmfit,newdata=data.frame(x=x),interval='confidence')
eshape <- list(n=nrow(e),rho=unname(cor$estimate),P=cor$p.value,
 x=e$x,y=e$y,fit_x=x,fit=pred[,'fit'],lower=pred[,'lwr'],upper=pred[,'upr'])
result <- list(A=a,B=list(P=ft$p.value,edited_decay=unname(tab[1,1]),edited_stable=unname(tab[1,2]),
 control_decay=unname(tab[2,1]),control_stable=unname(tab[2,2])),C=cc,
 D=as.data.frame(reg),E=eshape)
jsonlite::write_json(result,file.path(out,'plot_data.json'),auto_unbox=TRUE,digits=15,pretty=FALSE)
fwrite(d,file.path(out,'all_gene_source_data.tsv'),sep='\t')
fwrite(reg,file.path(out,'regression_results.tsv'),sep='\t')
fwrite(e,file.path(out,'panel_E_source_data.tsv'),sep='\t')
fwrite(data.table(panel='E',n=nrow(e),rho=unname(cor$estimate),P=cor$p.value),file.path(out,'panel_E_statistics.tsv'),sep='\t')
cat('A:',nrow(d),'edited',sum(d$edited),'control',sum(d$edited==0),'\n')
cat('C:',nrow(decay),'edited',sum(decay$edited),'control',sum(decay$edited==0),'\n')
cat('E:',nrow(e),'rho',unname(cor$estimate),'P',cor$p.value,'\n')
print(tab)
