ProbeMake <- function(fafile,LN=90,ln=60,TM=80,tm=60,CG=70,cg=30,gap=0,method="S2L",direction="3to5",
                      prohibitseq=NULL,nn_table="DNA_NN4",tmm_table ="DNA_TMM1",imm_table ="DNA_IMM1",
                      de_table="DNA_DE1",dnac1=25,dnac2=25,Na=50,K=0,Tris=0,Mg=0,dNTPs=0,saltcorr=5){
  if(method!="S2L" && method!="L2S"){
    stop("method must be one of 'S2L' and 'L2S'")
  }
  if(direction!='3to5' && direction!='5to3'){
    stop("direction must be one of '3to5' and '5to3'")
  }
  seqNameSet <- names(fafile)
  seqWid <- width(fafile)
  
  if(any(width(fafile)<ln)){
    warning(paste0("Sequence length in line ",paste(which(width(fafile)<ln),collapse=" ")," is less than the specified minimum length of probe"))
  }
  seqNum <- which(!1:length(fafile) %in% which(width(fafile)<ln))
  
  #Probe table list
  TotalProbe <- NULL
  ln <- ln - 1
  LN <- LN - 1
  #i=1
  for(i in seqNum){
    ProbeSet <- matrix(NA,1,7)
    colnames(ProbeSet) <- c("TargetID","Chr","Start","End","Sequence","Tm","GC")
    seqCont <- as.character(fafile[[i]])
    if(direction=='3to5'){
      seqCont <- seqCont
    }else{
      seqCont <- TmCalculator::complement(seqCont,reverse=TRUE)
    }
    SeqInfor <- unlist(strsplit(seqNameSet[i], "[ \t]+"))
    SeqID <- SeqInfor[1]
    SeqRange <- unlist(strsplit(SeqInfor[2], "[=:-]"))
    ChrID <- SeqRange[2]
    seqS <- as.numeric(SeqRange[3])
    seqE <- as.numeric(SeqRange[4])
    pstart <- 1
    
    if(method=="L2S"){
      n <- LN
      unitLen <- -1
      while(TRUE){
        # break if probe length is more than end of target
        if(pstart+n > seqWid[i]){
          n <- seqWid[i] - pstart
          if(n<ln){
            break
          }
        }
        SubSeqObj <- substr(seqCont, pstart, pstart+n)
        # remove prohibits sequence
        if(!is.null(prohibitseq)){
          AllMaxProhPos <- 0
          for(ps in prohibitseq){
            if(grepl(ps,SubSeqObj)){
              prohPosSet <- unlist(gregexpr(pattern=ps,SubSeqObj))
              prohPos <- max(prohPosSet)
              AllMaxProhPos <- append(AllMaxProhPos,prohPos)
            }
            maxpos <- max(AllMaxProhPos)
            pstart <- maxpos+pstart
          }
          if(maxpos != 0){
            next
          }
        }
        CGcont <- GC(SubSeqObj)
        Tm <- Tm_NN(SubSeqObj, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = nn_table,
                    tmm_table = tmm_table, imm_table = imm_table,de_table = de_table, dnac1 = dnac1,
                    dnac2 = dnac2, selfcomp = FALSE, Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
                    saltcorr = saltcorr)
        Index <- c((CGcont <= CG & CGcont >= cg) & (Tm >=tm & Tm <= TM))
        
        if(Index==FALSE){
          # extend probe start position if Tm > TM
          if(Tm < tm){
            pstart <- pstart + 1
            n <- LN
          }else{# improve Tm by extend probe length #Tm < tm
            n <- n + unitLen
          }
          # set probe sequence start position
          if(n < ln){
            pstart <- pstart + 1
            n <- LN
          }
        }else{
          if(direction=='3to5'){
            Start1 <- seqS+pstart-1
            End1 <- seqS+pstart-1+n
          }else{
            Start2 <- seqS+pstart-1
            End2 <- seqS+pstart-1+n
            Start1 <- seqS+seqE-End2
            End1 <- seqS+seqE-Start2
            SubSeqObj <- complement(SubSeqObj,reverse=TRUE)
          }
          ProbeSet <- rbind(ProbeSet,data.frame(TargetID=SeqID,Chr=ChrID,Start=Start1,End=End1,
                                                Sequence=SubSeqObj,Tm=Tm,GC=CGcont))
          pstart <- pstart + n + gap
          n <- LN
        }
      }
    }else if(method=="S2L"){
      n <- ln
      unitLen <- 1
      while(TRUE){
        # break if probe length is more than end of target
        if(pstart+n > seqWid[i]){
          break
        }
        SubSeqObj <- substr(seqCont, pstart, pstart+n)
        # remove prohibits sequence
        if(!is.null(prohibitseq)){
          AllMaxProhPos <- 0
          for(ps in prohibitseq){
            if(grepl(ps,SubSeqObj)){
              prohPosSet <- unlist(gregexpr(pattern=ps,SubSeqObj))
              prohPos <- max(prohPosSet)
              AllMaxProhPos <- append(AllMaxProhPos,prohPos)
            }
            maxpos <- max(AllMaxProhPos)
            pstart <- maxpos+pstart
          }
          if(maxpos != 0){
            next
          }
        }
        CGcont <- GC(SubSeqObj)
        Tm <- Tm_NN(SubSeqObj, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = nn_table,
                    tmm_table = tmm_table, imm_table = imm_table,de_table = de_table, dnac1 = dnac1,
                    dnac2 = dnac2, selfcomp = FALSE, Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
                    saltcorr = saltcorr)
        Index <- c((CGcont <= CG & CGcont >= cg) & (Tm >=tm & Tm <= TM))
        
        if(Index==FALSE){
          # extend probe start position if Tm > TM
          if(Tm > TM){
            pstart <- pstart + 1
            n <- ln
          }else{# improve Tm by extend probe length #Tm < tm
            n <- n + unitLen
          }
          # set probe sequence start position
          if(n > LN){
            pstart <- pstart + 1
            n <- ln
          }
        }else{
          if(direction=='3to5'){
            Start1 <- seqS+pstart-1
            End1 <- seqS+pstart-1+n
          }else{
            Start2 <- seqS+pstart-1
            End2 <- seqS+pstart-1+n
            Start1 <- seqS+seqE-End2
            End1 <- seqS+seqE-Start2
            SubSeqObj <- TmCalculator::complement(SubSeqObj,reverse=TRUE)
          }
          ProbeSet <- rbind(ProbeSet,data.frame(TargetID=SeqID,Chr=ChrID,Start=Start1,End=End1,
                                                Sequence=SubSeqObj,Tm=Tm,GC=CGcont))
          pstart <- pstart + n + gap
          n <- ln
        }
      }
    }
    TotalProbe <- rbind(TotalProbe,ProbeSet[-1,])
    message(paste("Target region ",i," (",SeqID,") is completed with ",nrow(ProbeSet)-1," probes at ",as.character(Sys.time()),"\n",sep=""))
  }
  rownames(TotalProbe) <- c(1:nrow(TotalProbe))
  return(TotalProbe)
}