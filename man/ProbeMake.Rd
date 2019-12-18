\name{ProbeMake}
\alias{ProbeMake}
\title{
Make probes
}
\description{
ProbeMake searches for probes with a FASTA-formatted input file containing the target sequence. And it allows users to specify allowable ranges of probe length, percent GC content, and adjust melting temperature calculated using nearest neighbor thermodynamics. Candidate probe sequences passing all checks are outputted in BED format.
}
\usage{
ProbeMake(fafile, LN = 90, ln = 60, TM = 80, tm = 60, CG = 70, cg = 30, gap = 0,
method = "S2L",direction = "3to5",prohibitseq=NULL,nn_table = "DNA_NN4", 
tmm_table = "DNA_TMM1", imm_table = "DNA_IMM1",de_table = "DNA_DE1",
dnac1 = 25, dnac2 = 25, Na = 50,K = 0, Tris = 0, Mg = 0, dNTPs = 0, saltcorr = 5)
}
\arguments{
  \item{fafile}{
Input file with a FASTA format read by function readDNAStringSet in R package 'Biostrings'
}
  \item{LN}{
The maximum allowed probe length, default is 90
}
  \item{ln}{
The minimum allowed probe length, default is 60
}
  \item{TM}{
The maximum allowed melting temperature, default is 80
}
  \item{tm}{
The minimum allowed melting temperature, default is 60
}
  \item{CG}{
The maximum allowed percent GC content, default is 70
}
  \item{cg}{
The minimum allowed percent GC content, default is 30
}
  \item{gap}{
The minimum gap between adjacent probes, default is 0
}
  \item{method}{
'S2L' is used to design probe extending from minimal probe length to the maximum until passing all checks, conversely 'L2S' make probe from maximal probe length to the minimum. Default is 'S2L'
}
  \item{direction}{
Design probes from 3 to 5 end of target sequence or from 5 to 3 end, default is '3to5'
}
  \item{prohibitseq}{
Prohibited sequence list, e.g prohibitseq=c("GGGGG","CCCCC"), default is NULL
}
  \item{nn_table}{
Thermodynamic NN values, eight tables are implemented. For DNA/DNA hybridizations: DNA_NN1,DNA_NN2,DNA_NN3,DNA_NN4

For RNA/RNA hybridizations: RNA_NN1,RNA_NN2,RNA_NN3

For RNA/DNA hybridizations: R_DNA_NN1

Default: DNA_NN4
}
  \item{tmm_table}{
Thermodynamic values for terminal mismatches. Default: DNA_TMM1
}
  \item{imm_table}{
Thermodynamic values for internal mismatches, may include insosine mismatches. Default: DNA_IMM1
}
  \item{de_table}{
Thermodynamic values for dangling ends: DNA_DE1(default),RNA_DE1
}
  \item{dnac1}{
Concentration of the higher concentrated strand [nM]. Typically this will be the primer (for PCR) or the probe. Default: 25
}
  \item{dnac2}{
Concentration of the lower concentrated strand [nM]. Default: 25
}
  \item{Na}{
Millimolar concentration of Na. Default: 50
}
  \item{K}{
Millimolar concentration of K. Default: 0
}
  \item{Tris}{
Millimolar concentration of Tris. Default: 0
}
  \item{Mg}{
Millimolar concentration of Mg
}
  \item{dNTPs}{
Millimolar concentration of dNTPs. Default: 50
}
  \item{saltcorr}{
Type of salt correction. Default: 5.
}
}
\value{
Returns a bed file in the format TargetID <tab> Chr <tab> Start <tab> End <tab> Sequence <tab> Tm <tab> GC
}
\references{
Beliveau B J, Kishi J Y, Nir G, et al. (2017). OligoMiner: A rapid, flexible environment for the design of genome-scale oligonucleotide in situ hybridization probes. bioRxiv.
}
\author{
Junhui Li
}
\examples{
data(samplefa)
ProbeMake(samplefa,LN=90,ln=60,TM=80,tm=70,CG=80,cg=20,gap=0,method="S2L",direction='3to5')
}
\keyword{ probe }