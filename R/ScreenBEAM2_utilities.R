###############################################################################
# Author: Jiyang Yu, Xinge Wang
# 2019.8.2
# Utility functions
###############################################################################


check.columns.lib<-function(lib){
  if(ncol(lib)<3){
    stop("Need more than 3 columns.")
  } else{
    lib[,2]<-as.character(lib[,2])
    if(is.na(str_match(lib[,2][1],"A")) & is.na(str_match(lib[,2][1],"T")) & is.na(str_match(lib[,2][1],"G")) & is.na(str_match(lib[,2][1],"C"))){
      stop("Second column need to be sequence.")
    }
  }
}

path2name<-function(path, file.extent){
  # Strip off sample name from file path
  path<-as.character(path)
  path.strip<-unlist(strsplit(path, split = "/"))
  name.extent<-path.strip[length(path.strip)]
  name<-gsub(file.extent, "", name.extent)
  return(name)
}

get_read_count<-function(path, every.n.line, file.extent){
  reads.path<-as.character(path)
  num.reads<-strsplit(system(paste0("wc -l <", reads.path), intern = T), split = " ")[[1]][1]
  num.reads<-as.numeric(num.reads)/every.n.line
  if(num.reads%%1!=0){
    print("ERROR: please check whether every.n.line is correct.")
    break
  }
  print(paste0(path2name(path,file.extent = ".fastq"), " has ", num.reads, " raw reads."))
  return(num.reads)
}

get_reads_info<-function(reads.path, file.extent, every.n.line){
  if(is.na(file.extent)) stop("Please give file extension name. Example: .fastq")
  if(is.na(every.n.line)) stop("Please give how many lines each read has. Example: 4")
  # Create a list object to store reads info
  reads.stat<-list()
  sample.names<-c()
  sample.reads.count<-c()
  for(path in reads.path){
    sample.names<-c(sample.names,path2name(path, file.extent=file.extent))
    print(paste0("Processing ", path2name(path, file.extent=file.extent)))
    sample.reads.count<-c(sample.reads.count, get_read_count(path, every.n.line=every.n.line, file.extent=file.extent))
  }
  reads.stat$sample.names<-sample.names
  reads.stat$sample.reads.count<-sample.reads.count
  reads.stat$fastq.path.list<-reads.path
  return(reads.stat)
}

create_unique_readfile<-function(fs, out.path){
  awk_com<-" | awk \'(NR%4==2){read=$1;total++;count[read]++} END {for(read in count) {print \">\" count[read] \"_\" num++ ; print read}}'"
  for(i in 1:length(fs)){
    sample.name<-path2name(fs[i], file.extent = ".fastq")
    message<-paste0("Converting ", sample.name , " to unique FASTA file.")
    print(message)
    comand.convert<-paste0("cat ", as.character(fs[i]), awk_com, " > ", out.path, sample.name,".fa")
    system(comand.convert)
    print(paste0("Done converting ",sample.name))
  }
}

run.blat<-function(analysis.par, query){
  blat.stat<-list()
  # Set parameter for blat (DON'T CHANGE)
  tileSize<-11
  stepSize<-floor(((analysis.par$RNA.len+1) - tileSize)/2)-1
  # Run blat
  blat.com.handler<-paste0("blat -tileSize=", tileSize, " -stepSize=", stepSize, " -oneOff=5 -minScore=0 -minIdentity=50 -maxGap=0 -repMatch=1000000 -dots=10000 ")
  for(path in analysis.par$fasta.path.list){
    sample.name<-path2name(path, file.extent = ".fa")
    print(paste0("Running blat for ", sample.name))
    blat.command<-paste0(blat.com.handler, path, paste0(" ", query)," -out=blast8 ", analysis.par$out.dir.output.mapping.step2, sample.name, ".txt")
    system(blat.command)
    print(paste0("Done blat mapping for ", sample.name))
  }
  blat.path.list<-list.files(analysis.par$out.dir.output.mapping.step2, '*.txt', recursive=T, full.names=T)
  blat.stat$blat.path.list <-blat.path.list
  blat.stat$blat.tileSize<-tileSize
  blat.stat$blat.stepSize<-stepSize
  return(blat.stat)
}

get_raw_count<-function(analysis.par, save.data.every.run = T){
  require(data.table)
  # Create a list to store count tables
  count.table.list<-list()
  if(save.data.every.run) save(count.table.list, file=paste0(analysis.par$out.dir.output.mapping.step3, "count_table_list.RData"))

  for(path in analysis.par$blat.path.list){
    # Get sample name
    sample.name<-path2name(path, file.extent = ".txt")
    print(paste0("Collecting raw count table for sample ", sample.name))

    # Read blat table from /output/mapping/Step2 folder
    print(paste0("Read sample's blat table ",sample.name))
    blat.df<-fread(path, col.names=c("sgRNA","unique_read_ID","identity","n.match","n.mid.mismatch","gap","start","end","ref.start","ref.end","e.value","bit.score"))

    # Remove duplicated/substring RNA
    if(!is.null(analysis.par$removed.RNA)){
      keep<-which(!blat.df$sgRNA %in% analysis.par$removed.RNA[,1])
      blat.df<-blat.df[keep,]
    }
    print(paste0("Done read sample's blat table ",sample.name))

    # Remove Gap mapping
    print(paste0("Remove gapped mapping for ",sample.name))
    blat.df<-blat.df[blat.df$gap==0, ]

    # Remove wrong direction mapping
    print(paste0("Get the right direction mapping for ",sample.name))
    blat.plus<-blat.df[blat.df$ref.end>blat.df$ref.start,]

    # Remove blat.df object to save memory
    rm(blat.df)

    # Order blat.plus by unique_read_ID, this is essential for pick reads
    setorder(blat.plus, unique_read_ID)

    # Get the most frequently trim position
    # Attention: there is a special case, when multiple trim position account for similar percentage, we still use only 1 position.
    # This may effect sh/sgRNA letter count number variates, but because similar sh/sgRNA targeting the same gene, this will not effect gene level
    ref.count<-setDT(blat.plus)[, .N, by=.(ref.start)]
    setorder(ref.count, -N)
    REF<-as.numeric(ref.count[1,1])
    print(paste0("The most common adapter starts from ", REF))

    # Get the most common sh/sgRNA length
    LEN<-as.numeric(tail(names(sort(table(blat.plus$ref.end - blat.plus$ref.start +1))), 1))
    print(paste0("The most common RNA length is ", LEN))

    # Add a column to the blat.plus table to mark how far the mapping is from the most common trim position
    blat.plus$ref.dis<--abs(blat.plus$ref.start - REF)

    # Calculate frequency of unique_read_ID to seperate 1-to-1 mapping and 1-to-multiple mapping
    tab<-setDT(blat.plus)[, .N, by=.(unique_read_ID)]

    # Store 1-to-1 mapping to blat.plus.one
    one.uniqueID<-tab$unique_read_ID[tab$N==1]
    print("Get unique mapping reads.")
    blat.plus.one<-blat.plus[which(blat.plus$unique_read_ID %in% one.uniqueID),] #Store unique match
    print("Deal with multiple mapping reads.")

    # Store 1-to-multiple mapping to blat.plus.more
    blat.plus.more<-blat.plus[which(!(blat.plus$unique_read_ID %in% one.uniqueID)),] #store multiple match, need to pick one

    # Sorting blat.plus.more to only maintain the first row, to get unique mapping
    setorder(blat.plus.more, -unique_read_ID, -n.match, -ref.dis, -bit.score)
    blat.plus.more.unique<-blat.plus.more %>% distinct(unique_read_ID, .keep_all = T )

    # Remove intermediate variable for memory use
    rm(blat.plus)
    rm(blat.plus.more)

    # Combine blat.plus.one and blat.plus.more.unique
    l<-list(blat.plus.more.unique, blat.plus.one)
    blat.clean<-rbindlist(l)

    # Remove intermediate variable for memory use
    rm(blat.plus.more.unique)
    rm(blat.plus.one)

    # Add n.total.mismatch = n.mid.mismatch + n.end.mismatch
    blat.clean$n.end.mismatch <- LEN - blat.clean$n.match
    blat.clean$n.total.mismatch <- blat.clean$n.mid.mismatch + blat.clean$n.end.mismatch

    # Split the unique_read_ID string to get the count number for each unique_read
    blat.clean$count<-as.numeric(str_split_fixed(blat.clean$unique_read_ID, "_", 2)[,1])

    ######### Collect raw count table
    # Count number of TOTAL mapped reads without n.mismatch
    print(paste0("Collecting total mapped reads for ", sample.name))
    count.table.list$total.mapping[[sample.name]]<-sum(blat.clean$count)

    # Count number of RAW COUNT with Maximum n.mismatch
    print(paste0("Collecting count table for ", sample.name))
    RNAcount = aggregate(count ~ sgRNA + count, data = blat.clean, sum)
    names(RNAcount)<-c("sgRNA","count")
    zero.RNA<-setdiff(analysis.par$RNA.list, RNAcount$sgRNA)
    if(length(zero.RNA)!=0){
      zero.RNA.count<-rep(0, length(zero.RNA))
      zero.df<-data.frame(sgRNA = zero.RNA, count = zero.RNA.count)
      RNAcount<-rbind(RNAcount, zero.df)
    }
    names(RNAcount)<-c("sgRNA",sample.name)
    count.table.list$RNA.count.df[[sample.name]]<-RNAcount

    # Count n.mismatch table
    print(paste0("Collecting n.mismatch table for ", sample.name))
    count.dist.df<-ddply(blat.clean, "sgRNA", function(df){
      n.vec<-rep(0,LEN+1)
      mis.vec<-as.numeric(unique(df$n.total.mismatch))
      for(n in mis.vec){
        # start from 0
        n.vec[n+1]<-sum(df[df$n.total.mismatch==n,"count"])
      }
      return(c(sample.name, n.vec))
    })
    colnames(count.dist.df)<-c("shRNAId","sample",paste0("X",seq(0,LEN)))
    zero.RNA<-setdiff(analysis.par$RNA.list, count.dist.df[,1])
    if(length(zero.RNA)!=0){
      zero.mat<-matrix(0L, nrow = length(zero.RNA), ncol = LEN+1)
      zero.df<-data.frame(zero.RNA, sample.name, zero.mat)
      colnames(zero.df)<-c("shRNAId","sample",paste0("X",seq(0,LEN)))
      count.dist.df<-rbind(count.dist.df, zero.df)
    }
    count.table.list$count.dist.df[[sample.name]]<-count.dist.df
    print(paste0("Done collecting raw counts for ", sample.name))
    if(save.data.every.run) save(count.table.list, file=paste0(analysis.par$out.dir.output.mapping.step3, "count_table_list.RData"))
  }
  return(count.table.list)
}

write_count_table<-function(raw.count.list, analysis.par){
  collect.all.list<-list()
  # get the count table
  count.df<-data.frame()
  count.df<-raw.count.list$RNA.count.df[[1]]
  for(i in 2:length(raw.count.list$RNA.count.df)){
    count.df<-merge(count.df, raw.count.list$RNA.count.df[[i]], by ="sgRNA")
  }
  rownames(count.df)<-count.df[,1]
  count.df<-count.df[,-1]
  collect.all.list[["raw.count.table"]]<-count.df
  write.csv(count.df, paste0(analysis.par$out.dir.output.mapping.step3, "count.csv"), quote = F)
  print(paste0("count.csv is ready in folder ", analysis.par$out.dir.output.mapping.step3))

  # get the count.dist.table
  count.dist.table<-rbindlist(raw.count.list$count.dist.df)
  collect.all.list[["raw.count.dist"]]<-count.dist.table
  write.csv(count.dist.table, paste0(analysis.par$out.dir.output.mapping.step3, "count.dist.csv"), quote = F)
  print(paste0("count.dist.csv is ready in folder ", analysis.par$out.dir.output.mapping.step3))

  # get the summary table
  total.df<-data.frame(sample = analysis.par$sample.names, total = analysis.par$sample.reads.count)
  mapped.df<-data.frame(sample = names(raw.count.list$total.mapping), n.matched = unlist(raw.count.list$total.mapping))
  rownames(mapped.df)<-c()
  summary.df<-merge(mapped.df, total.df, by = "sample")
  collect.all.list[["raw.summary"]]<-summary.df
  write.csv(summary.df, paste0(analysis.par$out.dir.output.mapping.step3, "summary.csv"), quote = F)
  print(paste0("summary.csv is ready in folder ", analysis.par$out.dir.output.mapping.step3))
  return(collect.all.list)
}

normalizeCountDist<-function(count.dist,total=1e6){

  options(digits=2+round(log10(total)))
  col.mm<-grep('^X[0-9]+$',names(count.dist),value=T)
  mismatch.table<-apply(count.dist[,col.mm],2,as.numeric)
  count.dist$total<-apply(mismatch.table,1,sum)
  names(count.dist)[2]<-'sampleName'

  count.dist<-merge(count.dist,plyr::ddply(count.dist,'sampleName',calT<-function(x){
    c(sample.total=sum(x$total))
  }),by='sampleName')

  fac<-total/count.dist$sample.total

  if('Undetermined'%in%count.dist$sampleName){
    count.dist$sampleName[grepl('Undetermined',count.dist$sampleName)]<-'Undetermined'
    ns<-nrow(subset(count.dist,!grepl('Undetermined',sampleName) & !duplicated(sampleName)))

    fac[count.dist$sampleName=='Undetermined']<-fac[count.dist$sampleName=='Undetermined']*ns*subset(count.dist,sampleName=='Undetermined')$sample.total[1]/sum(unique(subset(count.dist,sampleName!='Undetermined')$sample.total))
  }

  count.dist[,col.mm]<-round(mismatch.table*fac)

  count.dist[,-c(ncol(count.dist):(ncol(count.dist)-1))]
}

generateEset<-function(m,n.mismatch=NULL,normalize=TRUE,normalize.total=1e6){

  if(is.null(n.mismatch)){n.mismatch<-as.integer(gsub('X','',tail(names(m$count.dist),1)))}

  profiles<-getCountByMismatch(m$count.dist,n.mismatch=n.mismatch,normalize=normalize,normalize.total = normalize.total,annotation = NULL)

  lib<-m$features
  row.names(lib)<-lib[,1]
  group<-subset(m$samples,!grepl('Undetermined',sampleName))
  row.names(group)<-group$sampleName

  if(!all(row.names(group)%in%names(profiles))){
    stop('sampleNames are different between sampleAnnotation and profiles!')

  }
  if(!all(row.names(lib)%in%row.names(profiles)))
    stop('feature names are different between the lib and profile!\n')


  profiles<-profiles[row.names(lib),-1]
  profiles<-profiles[,row.names(group)]


  eset<-new("ExpressionSet",phenoData = new("AnnotatedDataFrame",group),
            featureData=new("AnnotatedDataFrame",lib),annotation='',
            exprs=as.matrix(profiles))

  eset
}

normalizeMiseqProfile<-function(d,total=NULL,pseudoCount=1,undertermined.col='Undetermined'){
  #special case of Undetermined, total number of Undetermined proporitial its percentage
  if(!is.null(undertermined.col)){
    if(is.numeric(undertermined.col)){
      if(!undertermined.col%in%1:ncol(d))
        stop('undertermined.col \'',undertermined.col, '\' doesn\'t exist!\n')
      dn<-data.frame(d[,!1:ncol(d)%in%undertermined.col],row.names=row.names(d))
      if(ncol(dn)==1) names(dn)<-setdiff(names(d),names(d)[undertermined.col])
    }else{
      if(!undertermined.col%in%names(d)){
        undertermined.col<-grep(undertermined.col,names(d),value=T)
        #stop('undertermined.col \'',undertermined.col, '\' doesn\'t exist!\n')
      }
      dn<-data.frame(d[,!names(d)%in%undertermined.col],row.names=row.names(d))
      if(ncol(dn)==1) names(dn)<-setdiff(names(d),undertermined.col)
    }


    dn<-normalize.scale(dn,total=total,pseudoCount = pseudoCount)

    du<-d[,undertermined.col]

    u.total<-total*ncol(data.frame(dn))*sum(d[,undertermined.col])/(sum(d)-sum(d[,undertermined.col]))

    du<-normalize.scale(du,total=u.total,pseudoCount = 0)

    dnew<-cbind(dn,du)
    names(dnew)[ncol(dnew)]<-setdiff(names(d),names(dn))

  }else{
    dnew<-normalize.scale(d,total=total,pseudoCount = pseudoCount)
  }

  dnew

}

normalize.scale<-function(d,total=NULL,pseudoCount=1){

  if(!is.data.frame(d)) {d<-data.frame(d)}
  if(!all(d>0)){ d<-d+pseudoCount}

  s<-apply(d,2,sum)
  m<-ifelse(is.null(total),as.integer(mean(s)),as.integer(total))
  options(digits=2+nchar(m))
  fac<-m/s
  for(i in 1:length(s)){
    d[,i]<-round(d[,i]*fac[i],0)
  }

  if(!all(d>0)) {d<-d+pseudoCount}
  return(d)
}

getCountByMismatch<-function(count.dist,n.mismatch,normalize=TRUE,normalize.total=1e6,annotation=NULL,undertermined.col=NULL,...){

  col.sel<-paste('X',0:n.mismatch,sep='')
  sel.df<-apply(count.dist[,col.sel],2, as.numeric)
  count.dist$total<-apply(sel.df,1,sum)
  dnew<-tidyr::spread(count.dist[,c(1:2,ncol(count.dist))],key='sampleName',value='total')
  dnew<-data.frame(dnew[,-1],row.names=dnew[,1])
  if(normalize)
    dnew<-normalizeMiseqProfile(dnew,total=normalize.total,undertermined.col=undertermined.col
                                ,...
    )
  dnew<-data.frame(shId=row.names(dnew),dnew)

  if(!is.null(annotation)){
    if(nrow(dnew)!=nrow(annotation))
      stop('diff shRNA number in annotation!\n')

    dnew<-merge(dnew,annotation,by.x=names(dnew)[1],by.y=names(annotation)[1],all.x=TRUE,sort=FALSE)
  }

  dnew
}

saveCountEset<-function(m,save.path,total=1e6, n.mismatch=2){

  cat('count.maxmm.raw.eset is saved...\n')
  count.maxmm.raw.eset<-generateEset(m,n.mismatch=NULL,normalize=FALSE,normalize.total=total)
  save(count.maxmm.raw.eset,file=file.path(save.path,'count.maxmm.raw.eset'))

  cat('count.maxmm.normalized.eset is saved...\n')
  count.maxmm.normalized.eset<-generateEset(m,n.mismatch=NULL,normalize=TRUE,normalize.total=total)
  save(count.maxmm.normalized.eset,file=file.path(save.path,paste('count.maxmm.normalized.',total/1e6,'M.eset',sep='')))

  prefix<-paste0("count.",n.mismatch,"mm")
  raw.name<-paste0(prefix,".raw.eset")
  cat(paste(raw.name, "is saved...\n"))
  count.Nmm.raw.eset<-generateEset(m,n.mismatch=n.mismatch,normalize=FALSE,normalize.total=total)
  save(count.Nmm.raw.eset,file=file.path(save.path,raw.name))

  norm.name<-paste0(prefix,".normalized.eset")
  cat(paste(norm.name," is saved...\n"))
  count.Nmm.normalized.eset<-generateEset(m,n.mismatch=n.mismatch,normalize=TRUE,normalize.total=total)
  save(count.Nmm.normalized.eset,file=file.path(save.path,paste(norm.name,".",total/1e6,'M.eset',sep='')))

}

generateEset.ScreenBEAM<-function(input.file,control.samples,case.samples,control.groupname='control',case.groupname='treatment',gene.columnId=2){
  # Create eset by only picking control and case columns
  require(Biobase)
  # Read the tsv table
  d<-read.table(input.file,sep='\t',header=T, check.names = F)
  # Check if sample names are duplicated
  if(sum(duplicated(names(d)))>0){
    print("ATTENTION! Sample names should not be duplicated!!")
    print("Duplicated sample columns will create wrong result!")
  }

  # Check the gene information column is the second one
  gene.column<-as.integer(gene.columnId)
  if(!gene.column%in%1:2) stop ('gene.columnId must be 1 or 2!\n')

  # sh/sgRNA column is the first 1, gene column is the second one, or vise versa
  id.column<-setdiff(1:2,gene.column)

  # row names of fd are sh/sgRNAs and get the first 2 columns
  fd<-data.frame(d[,c(id.column,gene.column)],row.names=d[,id.column])
  names(fd)<-c('id','gene')

  # Check if control.samples and case.samples are inside the count table
  if(is.character(c(control.samples,case.samples))){
    s.ctrl<-setdiff(control.samples,names(d))
    if(length(s.ctrl)>0) stop(paste('NO control samples: ',s.ctrl,'\n',collapse=', '))
    s.case<-setdiff(case.samples,names(d))
    if(length(s.case)>0) stop(paste('NO case samples: ',s.case,'\n',collapse=', '))
  }

  # To subset the raw count table by control and case
  d1<-data.frame(d[,c(control.samples,case.samples)],row.names=d[,id.column])

  # Add condition column for the Pdata(eset)
  pd<-data.frame(sampleName=names(d1),group=c(rep(control.groupname,length(control.samples)),rep(case.groupname,length(case.samples))),condition=c(rep('control',length(control.samples)),rep('treatment',length(case.samples))),row.names=names(d1))

  eset<-new("ExpressionSet",phenoData = new("AnnotatedDataFrame",pd),featureData=new("AnnotatedDataFrame",fd),exprs=as.matrix(d1))

  eset
}

DRAgeneLevel2<-function(eset,data.type=c('microarray','NGS'),do.normalization=FALSE, total=1e6, filterLowCount=TRUE,
                        filterBy='control',count.cutoff=4,nitt=15000,burnin=5000, thin=10, rna.size=6, sample.rna.time=100, method="Bayesian", pooling="partial",...){
  ### Check data type
  if(missing(data.type)){
    data.type<-ifelse(all(exprs(eset)>0) & is.integer(exprs(eset)),'NGS','microarray')
  }
  data.type<-match.arg(data.type,c('microarray','NGS'))
  cat(paste('data.type:',data.type,'data\n'))

  ## do.log2=TRUE, if data is NGS
  do.log2<-ifelse(data.type=='NGS',TRUE,FALSE)

  ###normalization
  if(do.normalization){
    if(data.type=='NGS'){
      cat('normalization: scale normlaization to NGS count data!\n')
      pseudoCount<-ifelse(all(exprs(eset)>0),0,1)
      if(pseudoCount>0){print("Some expression value are 0, pseudoCount=1 will be added.")}

      if(all(apply(exprs(eset),2,sum)<total)) {
        total<-total
      }else{
        print("Warning!!!Not all ColSum of expression data are larger than total!!! Please reset total. Now only return not normalized data.")
        total<-NULL
      }

      exprs(eset)<-as.matrix(normalize.scale(exprs(eset),total=total,pseudoCount = pseudoCount))
    }

    if(data.type=='microarray'){
      #quantile normalization
      cat('normalization: quantile normlaization to microarray log2(intensity) data!\n')
      exprs(eset)<-normalize.quantile(exprs(eset), ties=TRUE)
    }
  }

  ## Require fData column name to be gene or geneSymbol
  if(sum(names(fData(eset))=='gene')==0){
    if(sum(names(fData(eset))=='geneSymbol')==1){
      names(fData(eset))[names(fData(eset))=='geneSymbol']<-'gene'
    }else{
      stop('gene column doesn\'t exist!\n')
    }
  }

  # DR dataframe for each gene
  DR<-data.frame(gene=unique(fData(eset)$gene))

  condition<-'treatment'
  ctrl<-'control'

  if(length(condition)!=length(ctrl))
    stop('condion and ctrl have diff length!\n')

  if(!all(c(condition,ctrl)%in%pData(eset)$condition))
    stop('conditon, ctrl are not all in pData(eset)$condition!\n')

  #add n.shRNA column
  n.shRNAs<-table(fData(eset)$gene)
  n.shRNAs<-data.frame(gene=names(n.shRNAs),n.sh_sgRNAs.raw=as.integer(n.shRNAs))
  DR<-merge(n.shRNAs,DR,by='gene')

  eset.sel<-eset

  ## Filter out low count
  if(data.type=='NGS'){
    if(filterLowCount){
      sel<-grep(filterBy,pData(eset)$condition)
      if(length(sel)>=1){
        eset.sel<-eset[apply(exprs(eset[,sel]),1,median)>=count.cutoff,]
      }else{
        # For the case only 1 sample without replicates
        eset.sel<- eset[exprs(eset[,sel])>=count.cutoff,]
      }
    }
  }

  for(i in 1:length(condition)){

    cat(condition[i],'vs.', ctrl[i], 'in processs...\n')

    d.eset<-eset.sel[,pData(eset.sel)$condition %in% c(ctrl[i],condition[i])]

    condition.tag<-pData(d.eset)$group[pData(d.eset)$condition%in%condition[i]][1]
    ctrl.tag<-pData(d.eset)$group[pData(d.eset)$condition%in%ctrl[i]][1]

    if(do.log2){
      print("Perform log2 transformation for NGS data!")
      pseudoCount<-ifelse(all(exprs(d.eset)>0),0,1)
      if(pseudoCount>0){print("Some expression value are 0, pseudoCount=1 will be added.")}
      #take log2 for count data
      exprs(d.eset)<-log2(exprs(d.eset))
    }

    comp<-factor(gsub(condition[i],1,gsub(ctrl[i],0,pData(d.eset)$condition)))

    d<-data.frame(gene=as.character(fData(d.eset)$gene),exprs(d.eset),stringsAsFactors=FALSE)

    dr<-plyr::ddply(d,'gene', function(df){
      print(df$gene[1])
      # If number of RNA targeting the same gene is too large, we sample a subset of them and take median
      if(nrow(df)>rna.size){
        print(paste0(df$gene[1]," has more rna targeting 1 gene, now sample."))
        shuffle.list<-list()# to store each sample bid result
        for(iter in 1:sample.rna.time){
          #print(iter)
          sample_index<-sample(nrow(df),rna.size, replace = F)
          subdf<-df[sample_index,2:ncol(df)]
          df.bid<-combRowEvid.2grps(comp = comp, d = as.data.frame(subdf), family = gaussian, method = method, nitt = nitt,
                                    burnin = burnin, thin = thin, pooling = pooling, logTransformed = T, restand = F)
          shuffle.list[[iter]]<-df.bid
        }
        shuffle.df<-do.call(rbind, shuffle.list)
        shuffle.df<-shuffle.df[order(shuffle.df[,1], decreasing = F),]
        return(shuffle.df[round(nrow(shuffle.df)/2),])
      }else {
        subdf<-df[,2:ncol(df)]
        df.bid<-combRowEvid.2grps(comp = comp, d = as.data.frame(subdf), family = gaussian, method = method, nitt = nitt,
                                  burnin = burnin, thin = thin, pooling = pooling, logTransformed = T, restand = F)
        return(df.bid)
      }
    })

    if(do.log2){
      dr$AveSignal<-round(2^dr$AveSignal,0)
    }else{
      dr$AveSignal<-round(dr$AveSignal,0)
    }

    FDR.BH.partial<-p.adjust(dr$pval.partial,'BH')
    partial.id<-grep('pval.partial',names(dr))
    dr<-data.frame(dr[,1:partial.id],FDR.BH.partial=FDR.BH.partial,dr[,(partial.id+1):ncol(dr)])

    names(dr)<-gsub('.partial','',names(dr))

    m<-min(dr$pval[dr$pval>0])[1]
    dr$pval[dr$pval==0]<-ifelse(m<1e-7,m,1e-7)
    dr$z<-sign(dr$t)*qnorm(dr$pval/2,lower.tail = FALSE)

    #transform FC to log2FC
    dr$FC<-sign(dr$FC)*log2(abs(dr$FC))
    names(dr)<-gsub('^FC','log2FC',names(dr))
    names(dr)[-1]<-paste(names(dr)[-1],'.',sum(comp==1),'VS',sum(comp==0),'.',condition.tag,'_VS_',ctrl.tag,sep='')
    DR<-merge(DR,dr,by='gene')
  }


  names(DR)<-gsub('^FDR.BH','FDR',gsub('^sd','B.sd',gsub('^coef','B',gsub('n.levels','n.sh_sgRNAs.passFilter',names(DR)))))

  col.sel<-c('n.sh_sgRNAs.passFilter','log2FC','B','z','pval','FDR','B.sd')

  condition<-gsub('^z.','',grep('^z',names(DR),value=T));condition

  col.pre<-c(
    'gene'
    ,
    'n.sh_sgRNAs.raw'
    ,
    paste(rep(col.sel,each=length(condition)),rep(condition,length(col.sel)),sep='.')
  )


  DR<-DR[,col.pre]

  DR<-DR[order(DR[,grep('^B\\.',names(DR))[1]]),]

  return(DR)

}

FC <- function(x,cl,logTransformed=TRUE,log.base=2,average.method=c('geometric','arithmetic'),pseudoCount=0){
  x.class0 <- x[(cl == 0)]+pseudoCount
  x.class1 <- x[(cl == 1)]+pseudoCount
  if(missing(average.method))
    average.method<-'geometric'
  if(logTransformed){
    if(is.na(log.base)|log.base<0)
      stop('You must specify log.bsae !\n')
    logFC<-mean(x.class1)-mean(x.class0)
    FC.val<-sign(logFC)*log.base^abs(logFC)
  }else{
    logFC<-ifelse(average.method=='arithmetic',log(mean(x.class1))-log(mean(x.class0)),mean(log(x.class1)-mean(log(x.class0))))
    FC.val<-sign(logFC)*exp(abs(logFC))
  }
  FC.val[FC.val==0 | is.na(FC.val)]<-1
  FC.val
}

normalize.scale<-function(d,total=NULL,pseudoCount=1){

  if(!is.data.frame(d)) d<-data.frame(d)

  if(!all(d>0)) d<-d+pseudoCount
  s<-apply(d,2,sum)
  m<-ifelse(is.null(total),as.integer(mean(s)),as.integer(total))
  options(digits=2+nchar(m))
  fac<-m/s
  for(i in 1:length(s)){
    d[,i]<-round(d[,i]*fac[i],0)
  }
  if(!all(d>0)) d<-d+pseudoCount
  d

}


normalize.quantile<-function (M, ties = TRUE)
{
  n <- dim(M)
  if (is.null(n))
    return(M)
  if (n[2] == 1)
    return(M)
  O <- S <- array(, n)
  nobs <- rep(n[1], n[2])
  i <- (0:(n[1] - 1))/(n[1] - 1)
  for (j in 1:n[2]) {
    Si <- sort(M[, j], method = "quick", index.return = TRUE)
    nobsj <- length(Si$x)
    if (nobsj < n[1]) {
      nobs[j] <- nobsj
      isna <- is.na(M[, j])
      S[, j] <- approx((0:(nobsj - 1))/(nobsj - 1), Si$x,
                       i, ties = "ordered")$y
      O[!isna, j] <- ((1:n[1])[!isna])[Si$ix]
    }
    else {
      S[, j] <- Si$x
      O[, j] <- Si$ix
    }
  }
  m <- rowMeans(S)
  for (j in 1:n[2]) {
    if (ties)
      r <- rank(M[, j])
    if (nobs[j] < n[1]) {
      isna <- is.na(M[, j])
      if (ties)
        M[!isna, j] <- approx(i, m, (r[!isna] - 1)/(nobs[j] -
                                                      1), ties = "ordered")$y
      else M[O[!isna, j], j] <- approx(i, m, (0:(nobs[j] -
                                                   1))/(nobs[j] - 1), ties = "ordered")$y
    }
    else {
      if (ties)
        M[, j] <- approx(i, m, (r - 1)/(n[1] - 1), ties = "ordered")$y
      else M[O[, j], j] <- m
    }
  }
  M
}


combRowEvid.2grps<-function(d,comp,
                            family=gaussian
                            ,method=c('MLE','Bayesian'),pooling=c('full','no','partial'),

                            n.iter=1000
                            ,
                            prior.V.scale=0.02
                            ,
                            prior.R.nu=1
                            ,
                            prior.G.nu=2
                            ,
                            nitt = 13000
                            ,
                            burnin =3000
                            ,
                            thin=10
                            ,
                            restand=FALSE
                            ,
                            logTransformed=TRUE
                            ,
                            log.base=2
                            ,
                            average.method=c('geometric')
                            ,
                            pseudoCount=0
){

  if(!all(comp %in% c(1,0))){
    stop('comp only takes 1 or 0 !!! \n')
  }

  if(missing(pooling))
    pooling<-c('full','no','partial')
  else
    pooling<-tolower(pooling)

  pooling<-match.arg(pooling,several.ok=T)


  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (!family$family %in% c('gaussian','binomial', 'poisson')) {
    print(family)
    stop("Only Gaussian Poisson, and Binomial are supported!!! \n")
  }

  if(missing(method))
    method<-'Bayesian'
  if(grepl('Bayes',method,ignore.case = T)) method<-'Bayesian'
  if(grepl('MLE|MaxLikelihood',method,ignore.case = T)) method<-'MLE'
  method<-match.arg(method)

  d<-data.frame(t(d))


  dat<-data.frame(
    response=
      c(unlist(d[comp==levels(comp)[1],]),unlist(d[comp==levels(comp)[2],]))
    ,
    treatment=
      factor(c(rep(levels(comp)[1],sum(comp==levels(comp)[1])*ncol(d)),rep(levels(comp)[2],sum(comp==levels(comp)[2])*ncol(d))))
    ,
    probe=
      factor(c(rep(colnames(d),each=sum(comp==levels(comp)[1])),rep(colnames(d),each=sum(comp==levels(comp)[2]))))
  )


  #calculate FC
  FC.val<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='geometric',pseudoCount=pseudoCount)

  AveSignal<-mean(dat$response)
  n.levels<-nlevels(dat$probe)


  rs<-c(FC=FC.val,AveSignal=AveSignal,n.levels=as.integer(n.levels))

  #cat(rs,'\n')

  if('arithmetic' %in% average.method){
    FC.ari<-FC(dat$response,dat$treatment,logTransformed=logTransformed,log.base=log.base,average.method='arithmetic',pseudoCount=0)
    rs<-c(FC.ari_raw=FC.ari,rs)
  }


  #re-standarize the input data
  if(restand & sd(dat$response)>0)
    dat$response<-0.5*(dat$response-mean(dat$response))/sd(dat$response)

  #MLE approach to estimate parameters
  if(method=='MLE'){
    if(family$family=='gaussian'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat)
        else
          M.no<-glm(response ~ treatment, data=dat)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }


      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          #consider sd==0 case
          if(sd(dat$response)==0)
          {
            dat$response<-rnorm(nrow(dat),mean(dat$response),sd(dat$response)+0.001)
          }
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(	sum.tmp$devcomp$dims['n']-	sum.tmp$devcomp$dims['p']-	sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          #########################################
          #Another way to calculate pvalue
          #########################################
          #M.null<-glmer(response ~ (treatment + 1 | probe), data=dat)
          #anova(M.partial, M.null)
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='binomial'){
      if('full'%in%pooling){
        M.full<-glm(treatment ~ response, data=dat, family=family)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(treatment ~ response + probe, data=dat, family=family)
        else
          M.no<-glm(treatment ~ response, data=dat, family=family)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(treatment ~ response + (response + 1 | probe), data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            Dev.partial=as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(treatment ~ response, data=dat, family=family)
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }


    }else if(family$family=='poisson'){
      #comlete pooling
      if('full'%in%pooling){
        M.full<-glm(response ~ treatment, data=dat,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual,
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      #no pooling
      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-glm(response ~ treatment + probe, data=dat,family='poisson')
        else
          M.no<-glm(response ~ treatment, data=dat,family='poisson')

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual,
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      #partial pooling with multilevel model of varing slopes and intercepts
      if('partial'%in%pooling){
        if(n.levels>1){
          M.partial<-glmer(response ~ treatment + (treatment + 1 | probe), data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-as.numeric(sum.tmp$devcomp$dims['n']-sum.tmp$devcomp$dims['p']-sum.tmp$devcomp$dims['q'])
          pval.partial<-2*pt(abs(sum.tmp$coef[2,3]),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=as.numeric(sum.tmp$AICtab[1]),
            BIC.partial=as.numeric(sum.tmp$AICtab[2]),
            #REMLDev.partial=-as.numeric(sum.tmp$AICtab[5]),
            logLik=-as.numeric(sum.tmp$AICtab[3]),
            Dev.partial=-as.numeric(sum.tmp$AICtab[4])
          )
        }else{
          M.partial<-glm(response ~ treatment, data=dat,family='poisson')
          sum.tmp<-summary(M.partial)
          t.partial<-sum.tmp$coef[2,3]
          df.partial<-sum.tmp$df.residual
          pval.partial<-2*pt(abs(t.partial),lower.tail=FALSE,df=df.partial)
          z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
          rs.partial<-c(
            coef.partial=sum.tmp$coef[2,1],
            se.partial=sum.tmp$coef[2,2],
            t.partial=t.partial,
            pval.partial=pval.partial,
            z.partial=z.partial,
            df.partial=df.partial,
            AIC.partial=sum.tmp$aic,
            BIC.partial=AIC(M.partial,k=log(nrow(dat))),
            #REMLDev.partial=sum.tmp$deviance,
            logLik=logLik(M.partial),
            Dev.partial=sum.tmp$deviance
          )
        }
        rs<-c(rs,rs.partial)
      }

    }
    else{
      stop('Only liner model with gaussian or poisson distrn and binomial family model are supported !!! \n')
    }
  }

  #Bayesian approach
  else if(method=='Bayesian'){
    #if('full'%in%pooling | 'no'%in%pooling){
    #  require(arm)
    #}else{
    #  require(MCMCglmm)
    #}
    if(family$family=='gaussian'){
      if('full'%in%pooling){

        #comlete pooling
        M.full<-bayesglm(response ~ treatment, data=dat, maxit = n.iter)
        sum.tmp<-summary(M.full)

        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          M.no<-bayesglm(response ~ treatment + probe, data=dat, maxit=n.iter)
        else
          M.no<-bayesglm(response ~ treatment, data=dat, maxit =n.iter)
        sum.tmp<-summary(M.no)

        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)

      }

      if('partial'%in%pooling){

        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
          #if(df.partial<2)
          #	df.partial<-2
        }else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))

          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin)

          df.partial<-nrow(dat)-2
          #if(df.partial<2)
          #	df.partial<-2
          #cat(n.levels,df.partial,'\n')
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)

        if(is.na(pval.partial))
          pval.partial<-sum.tmp$sol[2,5]
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1]
          ,
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }

    }else if(family$family=='binomial'){

      if(family$link=='logit')
        prior.scale<-2.5
      else if(family$link=='probit')
        prior.scale<-2.5*1.6

      if('full'%in%pooling){
        M.full<-bayesglm(treatment ~ response, data=dat, family=family, n.iter = n.iter, prior.scale = prior.scale)
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          z.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #a bug in bayeglm for summary, fix: total obs - rank of the model
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        if(n.levels>1)
          M.no<-bayesglm(treatment ~ response + probe, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)
        else
          M.no<-bayesglm(treatment ~ response, data=dat, family=family, n.iter=n.iter,prior.scale = prior.scale)

        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          pval.no=sum.tmp$coef[2,4],
          z.no=sum.tmp$coef[2,3],
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }
      if('partial'%in%pooling){
        #categorial=logit, ordinal=probit
        if(family$link=='logit'){
          glmm.family<-'categorial'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels))
          }
        }
        else if(family$link=='probit'){
          glmm.family<-'ordinal'
          if(n.levels>1){
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = n.levels+1)))
          }else{
            prior<-list(R = list(V = prior.V.scale, nu=n.levels+1))
          }
        }
        else
          stop('For multilevel model with Binomial family and Bayeisan method, only logit and probit model are supported!!! \n')

        #			coef.scale<-ifelse(family$link=='logit', 2.5^2, (2.5*1.6)^2)
        #			intercept.scale<-ifelse(family$link=='logit', 10^2, (10*1.6)^2)
        #			scale<-var(dat$treatment)/var(dat$response)
        ################################# prior to be modified #################################
        #alpha.mu=rep(0,2),alpha.V=diag(2)*25^2)))
        #########################################################################################
        sd.threshold<-prior.scale*4
        #coef.sign<- -sign(logFC)
        #				while(sd.threshold>prior.scale*2){
        if(n.levels>1){
          M.partial<-MCMCglmm(treatment~response, random=~idh(1+response):probe,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-(n.levels+1)*2
        }else{
          M.partial<-MCMCglmm(treatment~response,family=glmm.family, prior=prior, data=dat, verbose=FALSE, nitt = nitt, burnin = burnin,thin=thin)
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #						z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          #effective sample size
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        #coef.sign<-sign(rs.partial['coef.partial'])
        sd.threshold<-as.double(rs.partial['sd.partial'])
        #				}
        rs<-c(rs,rs.partial)
      }
    }else if(family$family=='poisson'){
      if('full'%in%pooling){
        #comlete pooling
        M.full<-bayesglm(response ~ treatment, data=dat,n.iter = n.iter,family='poisson')
        sum.tmp<-summary(M.full)
        rs.full<-c(
          coef.full=sum.tmp$coef[2,1],
          se.full=sum.tmp$coef[2,2],
          t.full=sum.tmp$coef[2,3],
          pval.full=sum.tmp$coef[2,4],
          #z.full=sign(sum.tmp$coef[2,1])*abs(qnorm(sum.tmp$coef[2,4]/2)),
          df.full=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.full=sum.tmp$aic,
          BIC.full=AIC(M.full,k=log(nrow(dat))),
          Dev.full=sum.tmp$deviance
        )
        rs<-c(rs,rs.full)
      }

      if('no'%in%pooling){
        #no pooling
        if(n.levels>1)
          M.no<-bayesglm(response ~ treatment + probe, data=dat, n.iter=n.iter,family='poisson')
        else
          M.no<-bayesglm(response ~ treatment, data=dat, n.iter=n.iter,family='poisson')
        sum.tmp<-summary(M.no)
        rs.no<-c(
          coef.no=sum.tmp$coef[2,1],
          se.no=sum.tmp$coef[2,2],
          t.no=sum.tmp$coef[2,3],
          pval.no=sum.tmp$coef[2,4],
          #z.no=sign(coef.no)*abs(qnorm(pval.no/2)),
          df.no=sum.tmp$df.residual-sum.tmp$df[1],
          AIC.no=sum.tmp$aic,
          BIC.no=AIC(M.no,k=log(nrow(dat))),
          Dev.no=sum.tmp$deviance
        )
        rs<-c(rs,rs.no)
      }

      if('partial'%in%pooling){
        #partial pooling with multilevel model of varing slopes and intercepts
        prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu), G = list(G1 = list(V = diag(2)*prior.V.scale, nu = prior.G.nu )))

        if(n.levels>1){
          M.partial<-MCMCglmm(response ~ treatment, random=~idh(treatment+1):probe, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-(n.levels+1)*2
        }
        else{
          prior<-list(R = list(V = prior.V.scale, nu=prior.R.nu))
          M.partial<-MCMCglmm(response ~ treatment, data=dat, prior=prior,verbose=FALSE, nitt=nitt, burnin = burnin,thin=thin,family='poisson')
          df.partial<-nrow(dat)-2
        }
        sum.tmp<-summary(M.partial)
        t.partial<-sum.tmp$sol[2,1]/sd(M.partial$Sol[,2])
        pval.partial<-2*pt(-abs(t.partial),df=df.partial)
        z.partial<-sign(t.partial)*abs(qnorm(pval.partial/2))
        rs.partial<-c(
          coef.partial=sum.tmp$sol[2,1],
          'l-95% CI.partial'=sum.tmp$sol[2,2],
          'u-95% CI.partial'=sum.tmp$sol[2,3],
          t.partial=t.partial,
          pval.partial=pval.partial,
          z.partial=z.partial,
          df.partial=df.partial,
          #pMCMC
          pvalMCMC.partial=sum.tmp$sol[2,5],
          #					z.partial=sign(sum.tmp$sol[2,1])*abs(qnorm(sum.tmp$sol[2,5]/2)),
          n.effet.partial=sum.tmp$sol[2,4],
          n.MCMC.partial=sum.tmp$cstats[3],
          #Standard Deviation
          sd.partial=sd(M.partial$Sol[,2]),
          #only DIC saved, set AIC, BIC as DIC
          DIC.partial=sum.tmp$DIC,
          Dev.partial=mean(M.partial$Deviance)
        )
        rs<-c(rs,rs.partial)
      }
    }else{
      stop('Only liner model and binomial family model are supported !!! \n')
    }
  }

  rs

}
