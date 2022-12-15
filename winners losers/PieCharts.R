library(tidyverse)

trt<-c("CO2", "Drt.","Irg.", "Temp.", "N", "P", "Mult. Nuts.","Inter.")
value<-c(1,1,1,1,1,1,1,1)
fill<-c("Black", "Purple", "Darkgreen", "Green", "Blue", "Goldenro3", "Darkgray", "Orange")
pie<-data.frame(trt, value, fill)
#using ggplot - not sure i want to do this
ggplot(data=pie, aes(x=="",y=value, fill=trt))+
  geom_bar(stat = "identity")+
  coord_polar("x",start=0)+
  scale_fill_manual(values=fill)

##doing this in baseR
collst = RColorBrewer::brewer.pal(n = 8, name = "Dark2")

#make base pie chart 
pie(pie$value, labels=pie$trt, col=collst, border="white")

#for poaceae
labellistpoa<-c("+", "+","+","", "","-","","")
colactive<-collst
colactive[labellist==""]="darkgray"
pie2(pie$value, labels=labellistpoa, col=colactive, border="white",line_length = 1, text_center = 0.6, textcol = "white")


