require(paleoPhylo)
tab=read.csv("paleoPhylo.tab")
data(berggren95)
p93 <- with(tab, as.paleoPhylo(Code, Parent, Start, End, label=Name))
drawPhylo(p93,l2r=T,cexLab=1)
