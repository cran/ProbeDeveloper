#' Make probes
#'
#' Probes are made with a FASTA-formatted input file containing the target sequence. User can specify the allowable ranges of probe length, percent GC content, and adjust melting temperature calculated using nearest neighbor thermodynamics or empirical formulas based on GC content. Candidate probe sequences passing all checks output in BED format.
#'
#' @param fafile Input file with a FASTA format read by function readDNAStringSet in R package 'Biostrings'
#'
#' @param LN The maximum allowed probe length, default is 90
#'
#' @param ln The minimum allowed probe length, default is 60
#'
#' @param TM The maximum allowed melting temperature, default is 80
#'
#' @param tm The minimum allowed melting temperature, default is 60
#'
#' @param CG The maximum allowed percent GC content, default is 70
#'
#' @param cg The minimum allowed percent GC content, default is 30
#'
#' @param gap The minimum gap between adjacent probes, default is 0
#'
#' @param method 'S2L' is used to design probe extending from minimal length to the maximum until passing all checks, conversely 'L2S' make probe from maximal length to the minimum. Default is 'S2L'
#'
#' @param direction Design probes from 3 to 5 end or from 5 to 3 end of target sequence, default is '3to5'
#'
#' @param prohibitseq Prohibited sequence list, e.g prohibitseq=c("GGGGG","CCCCC"), default is NULL
#'
#' @param TmMethod The method used to calculate Tm, 'Tm_NN' and 'Tm_GC' can be seleted
#'
#' @param variant Empirical constants coefficient with 8 variant for 'Tm_GC' method: Chester1993, QuikChange, Schildkraut1965, Wetmur1991_MELTING, Wetmur1991_RNA, Wetmur1991_RNA/DNA, Primer3Plus and vonAhsen2001
#'
#' @param nn_table Thermodynamic NN values, eight tables are implemented.
#'
#' For DNA/DNA hybridizations:
#'   DNA_NN1,DNA_NN2,DNA_NN3,DNA_NN4
#'
#' For RNA/RNA hybridizations:
#'   RNA_NN1,RNA_NN2,RNA_NN3
#'
#' For RNA/DNA hybridizations:
#'   R_DNA_NN1
#'
#' @param tmm_table Thermodynamic values for terminal mismatches. Default: DNA_TMM1
#'
#' @param imm_table Thermodynamic values for internal mismatches, may include insosine mismatches. Default: DNA_IMM1
#'
#' @param de_table Thermodynamic values for dangling ends: DNA_DE1(default),RNA_DE1
#'
#' @param dnac1 Concentration of the higher concentrated strand [nM]. Typically this will be the primer (for PCR) or the probe. Default: 25
#'
#' @param dnac2 Concentration of the lower concentrated strand [nM]. Default: 25
#'
#' @param Na Millimolar concentration of Na, default is 0
#'
#' @param K Millimolar concentration of K, default is 0
#'
#' @param Tris Millimolar concentration of Tris, default is 0
#'
#' @param Mg Millimolar concentration of Mg, default is 0
#'
#' @param dNTPs Millimolar concentration of dNTPs, default is 0
#'
#' @param saltcorr Salt correction method. Options are "Schildkraut2010", "Wetmur1991","SantaLucia1996","SantaLucia1998-1","Owczarzy2004","Owczarzy2008". Note that "SantaLucia1998-2" is not available for this function.
#'
#' @param DMSO Percent of DMSO
#'
#' @param fmd Formamide concentration in percentage (fmdmethod="concentration") or molar (fmdmethod="molar")
#'
#' @param DMSOfactor Coeffecient of Tm decreases per percent DMSO. Default=0.75 von Ahsen N (2001) <PMID:11673362>. Other published values are 0.5, 0.6 and 0.675.
#'
#' @param fmdfactor Coeffecient of Tm decrease per percent formamide. Default=0.65. Several papers report factors between 0.6 and 0.72.
#'
#' @param fmdmethod "concentration" method for formamide concentration in percentage and "molar" for formamide concentration in molar
#'
#' @returns Returns a bed file in the format TargetID <tab> Chr <tab> Start <tab> End <tab> Sequence <tab> Tm <tab> GC
#'
#' @references
#'
#' Beliveau B J, Kishi J Y, Nir G, et al. (2017). OligoMiner: A rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. bioRxiv.
#'
#' Breslauer K J , Frank R , Blocker H , et al. Predicting DNA duplex stability from the base sequence.[J]. Proceedings of the National Academy of Sciences, 1986, 83(11):3746-3750.
#'
#' Sugimoto N , Nakano S , Yoneyama M , et al. Improved Thermodynamic Parameters and Helix Initiation Factor to Predict Stability of DNA Duplexes[J]. Nucleic Acids Research, 1996, 24(22):4501-5.
#'
#' Allawi, H. Thermodynamics of internal C.T mismatches in DNA[J]. Nucleic Acids Research, 1998, 26(11):2694-2701.
#'
#' Hicks L D , Santalucia J . The thermodynamics of DNA structural motifs.[J]. Annual Review of Biophysics & Biomolecular Structure, 2004, 33(1):415-440.
#'
#' Freier S M , Kierzek R , Jaeger J A , et al. Improved free-energy parameters for predictions of RNA duplex stability.[J]. Proceedings of the National Academy of Sciences, 1986, 83(24):9373-9377.
#'
#' Xia T , Santalucia , J , Burkard M E , et al. Thermodynamic Parameters for an Expanded Nearest-Neighbor Model for Formation of RNA Duplexes with Watson-Crick Base Pairs,[J]. Biochemistry, 1998, 37(42):14719-14735.
#'
#' Chen J L , Dishler A L , Kennedy S D , et al. Testing the Nearest Neighbor Model for Canonical RNA Base Pairs: Revision of GU Parameters[J]. Biochemistry, 2012, 51(16):3508-3522.
#'
#' Bommarito S, Peyret N, Jr S L. Thermodynamic parameters for DNA sequences with dangling ends[J]. Nucleic Acids Research, 2000, 28(9):1929-1934.
#'
#' Turner D H , Mathews D H . NNDB: the nearest neighbor parameter database for predicting stability of nucleic acid secondary structure[J]. Nucleic Acids Research, 2010, 38(Database issue):D280-D282.
#'
#' Sugimoto N , Nakano S I , Katoh M , et al. Thermodynamic Parameters To Predict Stability of RNA/DNA Hybrid Duplexes[J]. Biochemistry, 1995, 34(35):11211-11216.
#'
#' Allawi H, SantaLucia J: Thermodynamics and NMR of internal G-T mismatches in DNA. Biochemistry 1997, 36:10581-10594.
#'
#' Santalucia N E W J . Nearest-neighbor thermodynamics of deoxyinosine pairs in DNA duplexes[J]. Nucleic Acids Research, 2005, 33(19):6258-67.
#'
#' Peyret N , Seneviratne P A , Allawi H T , et al. Nearest-Neighbor Thermodynamics and NMR of DNA Sequences with Internal A-A, C-C, G-G, and T-T Mismatches, [J]. Biochemistry, 1999, 38(12):3468-3477.
#'
#' Marmur J , Doty P . Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature.[J]. Journal of Molecular Biology, 1962, 5(1):109-118.
#'
#' Schildkraut C . Dependence of the melting temperature of DNA on salt concentration[J]. Biopolymers, 2010, 3(2):195-208.
#'
#' Wetmur J G . DNA Probes: Applications of the Principles of Nucleic Acid Hybridization[J]. CRC Critical Reviews in Biochemistry, 1991, 26(3-4):33.
#'
#' Untergasser A , Cutcutache I , Koressaar T , et al. Primer3--new capabilities and interfaces[J]. Nucleic Acids Research, 2012, 40(15):e115-e115.
#'
#' von Ahsen N, Wittwer CT, Schutz E , et al. Oligonucleotide melting temperatures under PCR conditions: deoxynucleotide Triphosphate and Dimethyl sulfoxide concentrations with comparison to alternative empirical formulas. Clin Chem 2001, 47:1956-1961.
#'
#' @author Junhui Li
#'
#' @examples
#' data(samplefa)
#' ProbeMake(samplefa,LN=90,ln=60,TM=80,tm=70,CG=80,cg=20,TmMethod="Tm_NN",Na=50)
#'
#' @importFrom  TmCalculator complement GC Tm_NN Tm_GC
#' @importFrom  Biostrings width
#'
#' @export

ProbeMake <- function(fafile,
                      LN=90,
                      ln=60,
                      TM=80,
                      tm=60,
                      CG=70,
                      cg=30,
                      gap=0,
                      method=c("S2L","L2S"),
                      direction=c("3to5","5to3"),
                      prohibitseq=NULL,
                      TmMethod=c("Tm_GC","Tm_NN"),
                      variant=c("Primer3Plus",
                                "Chester1993",
                                "QuikChange",
                                "Schildkraut1965",
                                "Wetmur1991_MELTING",
                                "Wetmur1991_RNA",
                                "Wetmur1991_RNA/DNA",
                                "vonAhsen2001"),
                      nn_table=c("DNA_NN4",
                                 "DNA_NN1",
                                 "DNA_NN2",
                                 "DNA_NN3",
                                 "RNA_NN1",
                                 "RNA_NN2",
                                 "RNA_NN3",
                                 "R_DNA_NN1"),
                      tmm_table="DNA_TMM1",
                      imm_table="DNA_IMM1",
                      de_table=c("DNA_DE1",
                             "RNA_DE1"),
                      dnac1=25,
                      dnac2=25,
                      Na=0,
                      K=0,
                      Tris=0,
                      Mg=0,
                      dNTPs=0,
                      saltcorr=c("Schildkraut2010",
                                 "Wetmur1991",
                                 "SantaLucia1996",
                                 "SantaLucia1998-1",
                                 "SantaLucia1998-2",
                                 "Owczarzy2004",
                                 "Owczarzy2008"),
                      DMSO=0,
                      fmd=0,
                      DMSOfactor=0.75,
                      fmdfactor=0.65,
                      fmdmethod=c("concentration","molar")){
  method <- match.arg(method)
  direction <- match.arg(direction)
  nn_table <- match.arg(nn_table)
  tmm_table <- match.arg(tmm_table)
  imm_table <- match.arg(imm_table)
  de_table <- match.arg(de_table)
  saltcorr <- match.arg(saltcorr)
  fmdmethod <- match.arg(fmdmethod)

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
        if(TmMethod=="Tm_NN"){
          Tm <- Tm_NN(SubSeqObj, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = nn_table,
                      tmm_table = tmm_table, imm_table = imm_table,de_table = de_table, dnac1 = dnac1,
                      dnac2 = dnac2, selfcomp = FALSE, Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
                      saltcorr = saltcorr,DMSO=DMSO,fmd=fmd,DMSOfactor=DMSOfactor,
                      fmdfactor=fmdfactor,fmdmethod=fmdmethod)
        }else{
          Tm <- Tm_GC(SubSeqObj,
                       ambiguous=FALSE,
                       userset=NULL,
                       variant=variant,
                       Na=Na,
                       K=K,
                       Tris=Tris,
                       Mg=Mg,
                       dNTPs=dNTPs,
                       saltcorr=saltcorr,
                       mismatch=TRUE,
                       DMSO=DMSO,
                       fmd=fmd,
                       DMSOfactor=DMSOfactor,
                       fmdfactor=fmdfactor,
                       fmdmethod=fmdmethod)
        }
        Tm <- Tm$Tm
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
        if(TmMethod=="Tm_NN"){
          Tm <- Tm_NN(SubSeqObj, ambiguous = FALSE, comSeq = NULL, shift = 0, nn_table = nn_table,
                      tmm_table = tmm_table, imm_table = imm_table,de_table = de_table, dnac1 = dnac1,
                      dnac2 = dnac2, selfcomp = FALSE, Na = Na, K = K, Tris = Tris, Mg = Mg, dNTPs = dNTPs,
                      saltcorr = saltcorr,DMSO=DMSO,fmd=fmd,DMSOfactor=DMSOfactor,
                      fmdfactor=fmdfactor,fmdmethod=fmdmethod)
        }else{
          Tm <- Tm_GC(SubSeqObj,
                      ambiguous=FALSE,
                      userset=NULL,
                      variant=variant,
                      Na=Na,
                      K=K,
                      Tris=Tris,
                      Mg=Mg,
                      dNTPs=dNTPs,
                      saltcorr=saltcorr,
                      mismatch=TRUE,
                      DMSO=DMSO,
                      fmd=fmd,
                      DMSOfactor=DMSOfactor,
                      fmdfactor=fmdfactor,
                      fmdmethod=fmdmethod)
        }
        Tm <- Tm$Tm
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
