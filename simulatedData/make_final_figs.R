load("res_bin.Rdata")
load("res_lin.Rdata")
load("res_pois.Rdata")


dfpallete<-c("black","red","cadetblue","blue","blueviolet","deeppink")
dfpallete2<-c("red","cadetblue","blue","blueviolet","deeppink")



df<-rbind(df_bin,df_lin,df_pois)

df2<-rbind(df_bin2,df_lin2,df_pois2)

df2$model<-factor(df2$model,levels=c("LMM","BMM","PMM"))
df$model<-factor(df$model,levels=c("LMM","BMM","PMM"))

df$model2<-ifelse(df$model=="LMM","Linear~Mixed~Model",ifelse(df$model=="BMM","Binomial~Mixed~Model", "Poisson~Mixed~Model"))
df2$model2<-ifelse(df2$model=="LMM","Linear~Mixed~Model",ifelse(df2$model=="BMM","Binomial~Mixed~Model", "Poisson~Mixed~Model"))

df2$model2<-factor(df2$model2,levels=c("Linear~Mixed~Model","Binomial~Mixed~Model","Poisson~Mixed~Model"))
df$model2<-factor(df$model2,levels=c("Linear~Mixed~Model","Binomial~Mixed~Model","Poisson~Mixed~Model"))



p0<-ggplot(df[df$rho==0.5,])+
  geom_point(aes(x=Nn2,y=bound_t*100,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=bound_t*100,colour =  method),size=.7)+
  facet_grid(model2~vars,labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=label_parsed(levels(factor(df$Nn))))+
  ylab("boundary estimate (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0) +
  scale_colour_manual(values=dfpallete)+
  theme(legend.title=element_blank())


pp<-ggplot(df2[df2$rho==0.5&!(df2$method=="BM"&df2$model=="LMM"),])+
  geom_point(aes(x=Nn2,y=prial,colour =  method),size=2)+
  geom_line(aes(x=Nn2,y=prial,colour =  method),size=0.7)+
  facet_grid(model2~vars,scales="free_y",labeller = label_parsed)+
  theme_light()+
  scale_x_continuous(name="", breaks=c(1,2,3,4), labels=label_parsed(levels(factor(df2$Nn))))+
  ylab("PRIAL (%)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  geom_hline(yintercept=0)+
  scale_colour_manual(values=dfpallete2)+
  theme(legend.title=element_blank())

tiff("Fig_bound_all.tif",height=6*3/2,width=12,units="in",res=1200,compression = "lzw")
p0
dev.off()


pdf("Fig_bound_all.pdf",height=6*3/2,width=12 )
p0
dev.off()


tiff("Fig_prial_all.tif",height=6*3/2,width=12,units="in",res=1200,compression = "lzw")
pp
dev.off()


pdf("Fig_prial_all.pdf",height=6*3/2,width=12 )
pp
dev.off()



