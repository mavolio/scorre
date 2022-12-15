library(tidyverse)

#run pie2_function first

trt<-c("CO2", "Drt.","Irg.", "Temp.", "N", "P", "Mult. Nuts.","Inter.")
value<-c(1,1,1,1,1,1,1,1)
pie<-data.frame(trt, value)

##doing this in baseR
collst = RColorBrewer::brewer.pal(n = 8, name = "Dark2")

#make Legend
pie(pie$value, labels=pie$trt, col=collst, border="white")

#for poaceae overall
labellistpoa<-c("", "","","", "","","+","+")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for poaceae sub1
labellistpoa<-c("", "","","", "","+","","")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for poaceae sub2
labellistpoa<-c("", "-","","", "","","","")
colactive<-collst
colactive[labellistpoa==""]="darkgray"
pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")



#for cyperaceae
labellistcyp<-c("", "","","", "","","+","+")
colactive<-collst
colactive[labellistcyp==""]="darkgray"
pie2(pie$value, labels=labellistcyp, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for orchidaceae
labellistorc<-c("", "","","", "","","-","")
colactive<-collst
colactive[labellistorc==""]="darkgray"
pie2(pie$value, labels=labellistorc, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for asteraceae sub1
labellistast<-c("", "","","", "-","","","")
colactive<-collst
colactive[labellistast==""]="darkgray"
pie2(pie$value, labels=labellistast, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for asteraceae sub2
labellistast<-c("", "","+","", "","","","")
colactive<-collst
colactive[labellistast==""]="darkgray"
pie2(pie$value, labels=labellistast, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Lamiacea and Orobanaceae and Polomonaceae
labellistLO<-c("", "","","", "-","","","")
colactive<-collst
colactive[labellistLO==""]="darkgray"
pie2(pie$value, labels=labellistLO, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Plantaginaceae
labellistpla<-c("", "","","", "","","","-")
colactive<-collst
colactive[labellistpla==""]="darkgray"
pie2(pie$value, labels=labellistpla, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Gentinaceae
labellistgen<-c("", "","","", "-","","-","-")
colactive<-collst
colactive[labellistgen==""]="darkgray"
pie2(pie$value, labels=labellistgen, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Solaneaceae
labellistsol<-c("", "","","", "+","","+","")
colactive<-collst
colactive[labellistsol==""]="darkgray"
pie2(pie$value, labels=labellistsol, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Amaranth
labellistama<-c("", "+","-","", "+","-","","")
colactive<-collst
colactive[labellistama==""]="darkgray"
pie2(pie$value, labels=labellistama, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Fabaceae
labellistfab<-c("", "-","+","", "-","","-","")
colactive<-collst
colactive[labellistfab==""]="darkgray"
pie2(pie$value, labels=labellistfab, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")

#for Euphorb
labellisteup<-c("", "","","+", "","","","")
colactive<-collst
colactive[labellisteup==""]="darkgray"
pie2(pie$value, labels=labellisteup, col=colactive, border="white", line_length = 1, text_center = 0.6, textcol = "white")


#for Brassicaceae
labellistbra<-c("", "","","", "","","","+")
colactive<-collst
colactive[labellistbra==""]="darkgray"
pie2(pie$value, labels=labellistbra, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")