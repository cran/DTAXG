\name{DTAXG2}
\alias{DTAXG2}
\title{Two Diagnostic Tests Assessment in the Absence of Gold Standard}
\usage{
DTAXG2(group1,group2,
    prior.se.group1,prior.sp.group1,
    prior.se.group2,prior.sp.group2,
    prior.pi,n.sample,n.burnin,SUM)
}
\description{
To calculate the sensitivity and specificity for two diagnostic tests in the absence of gold standard using the Bayesian method. The prior information of sensitivity and specificity must be provided.
}
\arguments{
  \item{group1}{vector of 0 and 1, and 1 indicates the positive}
  \item{group2}{vector of 0 and 1, and 1 indicates the positive}
  \item{prior.se.group1}{the prior range of sensitity of group1}
  \item{prior.sp.group1}{the prior range of specificity of group1}
  \item{prior.se.group2}{the prior range of sensitity of group2}
  \item{prior.sp.group2}{the prior range of specificity of group2}
  \item{prior.pi}{the prior range of detection rate}
  \item{n.sample}{the number of the Gibbs sampling, defalt is 12000}
  \item{n.burnin}{the number of the burn-in in Gibbs sampling, defalt is 2000}
  \item{SUM}{wheter to return the summary results, defalt is 'TRUE'}
}
\value{
  \item{table}{if 'SUM' is TRUE, the Q50, Q2.5, and Q97.5 of sensitivity and specificity will be shown}
  \item{PI}{detection rate}
  \item{S1}{sensitivity of group1}
  \item{S2}{sensitivity of group2}
  \item{C1}{specificity of group1}
  \item{C2}{specificity of group2}
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\references{
Haiyan Gu, Qiguang Chen. Diagnostic Test Assessment in the Absence of Gold Standard. Zhong Guo Wei Sheng Tong Ji (in Chinese). 1999. 16(4): 203-205.
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:

more functions will be included in 'DATXG' package!
}
\examples{
group1=c(rep(1,86),rep(0,64))
group2=c(rep(1,27),rep(0,86-27),rep(1,5),rep(0,64-5))
prior.se.group1=c(0.5,0.85)
prior.sp.group1=c(0.8,1)
prior.se.group2=c(0.1,0.5)
prior.sp.group2=c(0.1,0.5)
prior.pi=c(0.66,0.67)

rst=DTAXG2(
  group1,
  group2,
  prior.se.group1,
  prior.sp.group1,
  prior.se.group2,
  prior.sp.group2,
  prior.pi,
  n.sample=120,
  n.burnin=20)

print(rst)

}
