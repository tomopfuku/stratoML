require(ggplot2)
ymin=seq(from=0,to=(nrow(tab)-1)*3,by=3)
ymax = seq(from=1,to=(nrow(tab))*3,by=3)
tab$ymin=ymin
tab$ymax=ymax
p=ggplot(tab)
tab$text_pos_x=tab$End#((tab$Start+tab$End)/2)-.5
tab$text_pos_y=(tab$ymax+tab$ymin)/2
for (row in 1:nrow(tab)) {
	currow = tab[row,]
	if (is.na(currow$Parent)){
		next
	}
	parrow = tab[tab$Code ==currow$Parent,]
	curind = which(tab$Code==currow$Code)
	parind = which(tab$Code==parrow$Code)
	if (curind > parind) {
		ystart = currow$ymin
		yend = parrow$ymax
	} else {
		ystart = currow$ymax
		yend = parrow$ymin
	}
	if (parrow$End > currow$Start) {
		xend = parrow$End
	} else{
		xend = parrow$Start-((abs(parrow$Start-currow$Start)/5)*4)
	}
	p = p + annotate(geom = "segment", linetype=3,x = currow$Start,y = ystart,xend=xend,yend=yend)
}

p=p+geom_rect(aes(xmin=Start,xmax=End,ymin=ymin,ymax=ymax),fill="black")+theme_bw()+geom_text(aes(x=tab$text_pos_x,y=tab$text_pos_y),hjust=-0.01,label=tab$Name,size=5.2)
p=p+geom_rect(aes(xmin=FAD,xmax=LAD,ymin=ymin,ymax=ymax),fill="grey")+theme_bw()+geom_text(aes(x=tab$text_pos_x,y=tab$text_pos_y),hjust=-0.01,label=tab$Name,size=5.2)
p=p+labs(x="ma",y="")
p=p+theme(axis.text.y=element_blank(),axis.ticks.y=element_blank(),panel.grid.major.y = element_blank())
p=p+scale_x_reverse()
p=p+xlim(max(tab$Start),min(tab$text_pos_x)-3.01)
#traits = read.table("ectocyon.p3.tab",sep=" ")
#tab$ord=row.names(tab)
#alldat=merge(tab,traits,by=1,all.x=T)
