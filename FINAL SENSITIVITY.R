##########################################################################################
## R-CODE
## Dr. I Gede Nyoman Mindra Jaya, M.Si
##  
##########################################################################################

#(1) call libraries 
library(spatialreg)
library(maptools)
library(spdep)
library(INLA)
library(sp)
library(RColorBrewer)
library(lattice)
library(gstat)
library(raster)
require(splancs)
library(ggplot2)
library(dplyr)
library(tidyr)
library(brinla)
library(rgdal)
library(maps)
library(ggplot2)
library(grid)
library(ggsn)
library(sf)
library(rgeos)


 
setwd("/Users/mindra/Jabar") 
Jabar<-readShapePoly("Jabar.shp")
proj4string(Jabar) <- CRS("+proj=longlat +datum=WGS84")
Jabar$ID<-c(1:27) 
#project your data (I'm using California Albers here) and apply a zero-width buffer
Jabar <- spTransform(Jabar, CRS("+proj=longlat +datum=WGS84"))
Jabar <- gBuffer(Jabar, byid = T, width = 0)

#Attempting a UnionSpatialPolygons based on the COUNTY field
Jabar.df <- as(Jabar, "data.frame")
countycol <- Jabar.df$ID
Jabar.diss <- unionSpatialPolygons(Jabar, countycol)
Jabar1<-Jabar.diss
 

setwd("/Users/mindra/@Sensitivity")

DataDengue<-read.csv("DBD.csv", sep=";")			#All data with NA

DataDengue$id<-DataDengue$ID

JabarData <- fortify(Jabar1, region="ID")
DataDengue.shp<-merge(JabarData, DataDengue, by="id", all.x=TRUE)
DataDengueMap<-DataDengue.shp[order(DataDengue.shp$order), ] 
#DataDengueMap<-na.omit(DataDengueMap)


 
#####Map discrete



pretty_breaks <- c(0.5, 1, 1.50, 2) # find the extremes

minVal <- 0
maxVal <- 100

# compute labels
labels <- c()
brks <- c(-1, pretty_breaks, 180) # round the labels (actually, only the extremes)
for(idx in 1:length(brks)){
  labels <- c(labels,round(brks[idx + 1], 2))
}

labels <- labels[1:length(labels)-1]  # define a new variable on the data set just as above
DataDengueMap$brks <- cut(DataDengueMap$SMR, 
                          breaks = brks, 
                          include.lowest = TRUE, 
                          labels = labels)

brks_scale <- levels(DataDengueMap$brks)
labels_scale <- rev(brks_scale)


#DataDengueMap1<-subset(DataDengueMap, IT<4)

png("Fig2.png", units="in", width=14, height=8, res=600) 
 
ggplot() + 
  geom_polygon(data = DataDengueMap, aes(fill = brks, x = long,  y = lat, group = group),color="black",size=0.1) +
  facet_wrap(~Year, ncol=3,scale="free")+   
  theme_bw()+ ylab(" ")+xlab(" ")+
  theme(legend.position = "bottom")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(labels = c("(0.0-0.5]", "(0.5-1.0]","(1.0-1.5]", "(1.5-2.0]", ">2"), values =  c("cyan1","green", "yellow","orange", "red"))+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(fill = "SIR") +
  theme(legend.position="bottom")+theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        )

dev.off()


 

png("Fig1.png", units="in", width=8, height=8, res=600) 
 
ggplot() + 
  geom_line(data = DataDengue, aes(x = Year,  y = Cases, color = District),size=0.3) + 
  theme_bw()+ ylab(" ")+xlab(" ")+
  theme(legend.position = "bottom")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+  
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(fill = "Number of cases")  



dev.off()



  
#Load Dengue data
setwd("/Users/mindra/Jabar") 

Jabar<-readShapePoly("Jabar.shp")

proj4string(Jabar) <- CRS("+proj=longlat +datum=WGS84")

#Compute adjacency matrix, as nb object 'adj' and sparse matrix 'W'

adj <- poly2nb(Jabar)
W <- as(nb2mat(adj, style = "B"), "Matrix")

setwd("/Users/mindra/@Sensitivity")
 
# (4) Lakukan pemdoelan dengan INLA
###  SETUP INLA OUTPUT ###
control <- list(
  predictor = list(compute = TRUE,link=1),
results = list(return.marginals.random = TRUE,  return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals=TRUE, return.marginals.predictor=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
  po = TRUE, waic=TRUE, graph=TRUE, openmp.strategy="huge"))
  
  
  
  
prior.c1=c(0.01,0.01)
igprior1=list(theta=list(prior="loggamma",param=prior.c1))
 
prior.c2=c(0.1,0.1)
igprior2=list(theta=list(prior="loggamma",param=prior.c2))


prior.c3=c(1,1)
igprior3=list(theta=list(prior="loggamma",param=prior.c3))


prior.c4=c(1,0.1)
igprior4=list(theta=list(prior="loggamma",param=prior.c4))


prior.c5=c(1,0.01)
igprior5=list(theta=list(prior="loggamma",param=prior.c5))

prior.c6=c(1,0.001)
igprior6=list(theta=list(prior="loggamma",param=prior.c6))

 
prior.c7=c(1,0.0001)
igprior7=list(theta=list(prior="loggamma",param=prior.c7))

prior.c8=c(1,0.00001)
igprior8=list(theta=list(prior="loggamma",param=prior.c8))
 
 
#### HC Prior
halfcauchy = "expression:
lambda = 10;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
 
hcprior = list(prec = list(prior = halfcauchy))
hcprior1<-hcprior 


halfcauchy = "expression:
lambda = 15;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
 
hcprior = list(prec = list(prior = halfcauchy))
hcprior2<-hcprior 



halfcauchy = "expression:
lambda = 20;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
 
hcprior = list(prec = list(prior = halfcauchy))
hcprior3<-hcprior 

halfcauchy = "expression:
lambda = 25;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
 
hcprior = list(prec = list(prior = halfcauchy))
hcprior4<-hcprior 


halfcauchy = "expression:
lambda = 30;
precision = exp(log_precision);
logdens = -1.5*log_precision-log(pi*lambda)-log(1+1/(precision*lambda^2));
log_jacobian = log_precision;
return(logdens+log_jacobian);"
 
hcprior = list(prec = list(prior = halfcauchy))
hcprior5<-hcprior 

#### Uniform Prior

sdunif="expression:
logdens=-log_precision/2;
return(logdens)"

uprior=list(prec=list(prior=sdunif))  

###PC PRIOR




sdres<-sd(DataDengue$Cases)
pcprior1 <- list(prec = list(prior="pc.prec", param = c(sdres,0.01)))

pcprior2 <- list(prec = list(prior="pc.prec", param = c(sdres,0.1)))

pcprior3 <- list(prec = list(prior="pc.prec", param = c(sdres,0.5)))


DataDengue$IT1<-DataDengue$IT
DataDengue$ID1<-DataDengue$ID



MF_IG1<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+  
       f(ID1,model="besagproper2",hyper = igprior1,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior1)) 
	   	

RF_IG1<- inla(MF_IG1, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
MF_IG2<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior2,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior2)) 
	   	

RF_IG2<- inla(MF_IG2, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
MF_IG3<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior3,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior3)) 
	   	

RF_IG3<- inla(MF_IG3, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG4<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior4,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior4)) 
	   	

RF_IG4<- inla(MF_IG4, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG5<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior5,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior5)) 
	   	

RF_IG5<- inla(MF_IG5, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG6<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior6,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior6)) 
	   	

RF_IG6<- inla(MF_IG6, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG7<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior7,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior7)) 
	   	

RF_IG7<- inla(MF_IG7, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG8<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior8,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior8)) 
	   	

RF_IG8<- inla(MF_IG8, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			 


######			
			
MF_HC1<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior1,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC1<- inla(MF_HC1, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))



		
MF_HC2<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior2,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC2<- inla(MF_HC2, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
				
			
MF_HC3<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior3,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC3<- inla(MF_HC3, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
				
			
MF_HC4<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior4,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC4<- inla(MF_HC4, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
			
				
			
MF_HC5<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior5,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC5<- inla(MF_HC5, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
			
				 

DataDengue$IT1<-DataDengue$IT
DataDengue$ID1<-DataDengue$ID

MF_U<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = uprior,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_U<- inla(MF_U, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			 
###### PC




			 



MF_PC1<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = pcprior1,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC1<- inla(MF_PC1, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_PC2<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = pcprior2,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC2<- inla(MF_PC2, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			



MF_PC3<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+  
       f(ID1,model="besagproper2",hyper = pcprior3,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC3<- inla(MF_PC3, family="poisson",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


###NB
MF_IG1NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+  
       f(ID1,model="besagproper2",hyper = igprior1,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior1)) 
	   	

RF_IG1NB<- inla(MF_IG1, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
MF_IG2NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior2,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior2)) 
	   	

RF_IG2NB<- inla(MF_IG2, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
MF_IG3NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior3,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior3)) 
	   	

RF_IG3NB<- inla(MF_IG3, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG4NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior4,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior4)) 
	   	

RF_IG4NB<- inla(MF_IG4, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG5NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior5,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior5)) 
	   	

RF_IG5NB<- inla(MF_IG5, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG6NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior6,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior6)) 
	   	

RF_IG6NB<- inla(MF_IG6, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG7NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior7,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior7)) 
	   	

RF_IG7NB<- inla(MF_IG7, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_IG8NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = igprior8,  graph=W, group=IT1, control.group=list(model="ar1", hyper = igprior8)) 
	   	

RF_IG8NB<- inla(MF_IG8, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			 


######			
			
MF_HC1NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior1,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC1NB<- inla(MF_HC1, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))



		
MF_HC2NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior2,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC2NB<- inla(MF_HC2, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
				
			
MF_HC3NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior3,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC3NB<- inla(MF_HC3, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
				
			
MF_HC4NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior4,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC4NB<- inla(MF_HC4, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
			
				
			
MF_HC5NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = hcprior5,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_HC5NB<- inla(MF_HC5, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			
			
			
				 

MF_UNB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = uprior,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_UNB<- inla(MF_U, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			 
###### PC



MF_PC1NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = pcprior1,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC1NB<- inla(MF_PC1, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			


MF_PC2NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+ 
       f(ID1,model="besagproper2",hyper = pcprior2,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC2NB<- inla(MF_PC2, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			



MF_PC3NB<-Cases~ HealtyBehaviour+ f(ID,model="besagproper2",graph=W,hyper = uprior)+f(IT,model="ar1",hyper = uprior)+  
       f(ID1,model="besagproper2",hyper = pcprior3,  graph=W, group=IT1, control.group=list(model="ar1")) 
	   	

RF_PC3NB<- inla(MF_PC3, family="nbinomial",E=Expected, data=DataDengue, 
            control.compute = control$compute,
            control.predictor = control$predictor,
            control.inla = list(int.strategy = "eb", strategy = "simplified.laplace"),
			control.fixed=list(prec=.1,prec.intercept=.01))
			  
  
  
  
  
  
  
  
###########

  
  
DICRF_IG1<-RF_IG1$dic$dic
DICRF_IG2<-RF_IG2$dic$dic
DICRF_IG3<-RF_IG3$dic$dic
DICRF_IG4<-RF_IG4$dic$dic
DICRF_IG5<-RF_IG5$dic$dic
DICRF_IG6<-RF_IG6$dic$dic
DICRF_IG7<-RF_IG7$dic$dic
DICRF_IG8<-RF_IG8$dic$dic
DICRF_HC1<-RF_HC1$dic$dic
DICRF_HC2<-RF_HC2$dic$dic
DICRF_HC3<-RF_HC3$dic$dic
DICRF_HC4<-RF_HC4$dic$dic
DICRF_HC5<-RF_HC5$dic$dic
DICRF_U<-RF_U$dic$dic
DICRF_PC1<-RF_PC1$dic$dic
DICRF_PC2<-RF_PC2$dic$dic
DICRF_PC3<-RF_PC3$dic$dic
DICRF_IG1NB<-RF_IG1NB$dic$dic
DICRF_IG2NB<-RF_IG2NB$dic$dic
DICRF_IG3NB<-RF_IG3NB$dic$dic
DICRF_IG4NB<-RF_IG4NB$dic$dic
DICRF_IG5NB<-RF_IG5NB$dic$dic
DICRF_IG6NB<-RF_IG6NB$dic$dic
DICRF_IG7NB<-RF_IG7NB$dic$dic
DICRF_IG8NB<-RF_IG8NB$dic$dic
DICRF_HC1NB<-RF_HC1NB$dic$dic
DICRF_HC2NB<-RF_HC2NB$dic$dic
DICRF_HC3NB<-RF_HC3NB$dic$dic
DICRF_HC4NB<-RF_HC4NB$dic$dic
DICRF_HC5NB<-RF_HC5NB$dic$dic
DICRF_UNB<-RF_UNB$dic$dic
DICRF_PC1NB<-RF_PC1NB$dic$dic
DICRF_PC2NB<-RF_PC2NB$dic$dic
DICRF_PC3NB<-RF_PC3NB$dic$dic




DIC<-c(DICRF_IG1,DICRF_IG2,
DICRF_IG3,
DICRF_IG4,
DICRF_IG5,
DICRF_IG6,
DICRF_IG7,
DICRF_IG8,
DICRF_HC1,
DICRF_HC2,
DICRF_HC3,
DICRF_HC4,
DICRF_HC5,
DICRF_U,
DICRF_PC1,
DICRF_PC2,
DICRF_PC3,
DICRF_IG1NB,
DICRF_IG2NB,
DICRF_IG3NB,
DICRF_IG4NB,
DICRF_IG5NB,
DICRF_IG6NB,
DICRF_IG7NB,
DICRF_IG8NB,
DICRF_HC1NB,
DICRF_HC2NB,
DICRF_HC3NB,
DICRF_HC4NB,
DICRF_HC5NB,
DICRF_UNB,
DICRF_PC1NB,
DICRF_PC2NB,
DICRF_PC3NB)




waicRF_IG1<-RF_IG1$waic$waic
waicRF_IG2<-RF_IG2$waic$waic
waicRF_IG3<-RF_IG3$waic$waic
waicRF_IG4<-RF_IG4$waic$waic
waicRF_IG5<-RF_IG5$waic$waic
waicRF_IG6<-RF_IG6$waic$waic
waicRF_IG7<-RF_IG7$waic$waic
waicRF_IG8<-RF_IG8$waic$waic
waicRF_HC1<-RF_HC1$waic$waic
waicRF_HC2<-RF_HC2$waic$waic
waicRF_HC3<-RF_HC3$waic$waic
waicRF_HC4<-RF_HC4$waic$waic
waicRF_HC5<-RF_HC5$waic$waic
waicRF_U<-RF_U$waic$waic
waicRF_PC1<-RF_PC1$waic$waic
waicRF_PC2<-RF_PC2$waic$waic
waicRF_PC3<-RF_PC3$waic$waic
waicRF_IG1NB<-RF_IG1NB$waic$waic
waicRF_IG2NB<-RF_IG2NB$waic$waic
waicRF_IG3NB<-RF_IG3NB$waic$waic
waicRF_IG4NB<-RF_IG4NB$waic$waic
waicRF_IG5NB<-RF_IG5NB$waic$waic
waicRF_IG6NB<-RF_IG6NB$waic$waic
waicRF_IG7NB<-RF_IG7NB$waic$waic
waicRF_IG8NB<-RF_IG8NB$waic$waic
waicRF_HC1NB<-RF_HC1NB$waic$waic
waicRF_HC2NB<-RF_HC2NB$waic$waic
waicRF_HC3NB<-RF_HC3NB$waic$waic
waicRF_HC4NB<-RF_HC4NB$waic$waic
waicRF_HC5NB<-RF_HC5NB$waic$waic
waicRF_UNB<-RF_UNB$waic$waic
waicRF_PC1NB<-RF_PC1NB$waic$waic
waicRF_PC2NB<-RF_PC2NB$waic$waic
waicRF_PC3NB<-RF_PC3NB$waic$waic




waic<-c(waicRF_IG1,waicRF_IG2,
waicRF_IG3,
waicRF_IG4,
waicRF_IG5,
waicRF_IG6,
waicRF_IG7,
waicRF_IG8,
waicRF_HC1,
waicRF_HC2,
waicRF_HC3,
waicRF_HC4,
waicRF_HC5,
waicRF_U,
waicRF_PC1,
waicRF_PC2,
waicRF_PC3,
waicRF_IG1NB,
waicRF_IG2NB,
waicRF_IG3NB,
waicRF_IG4NB,
waicRF_IG5NB,
waicRF_IG6NB,
waicRF_IG7NB,
waicRF_IG8NB,
waicRF_HC1NB,
waicRF_HC2NB,
waicRF_HC3NB,
waicRF_HC4NB,
waicRF_HC5NB,
waicRF_UNB,
waicRF_PC1NB,
waicRF_PC2NB,
waicRF_PC3NB)


mlikRF_IG1<-RF_IG1$mlik[2]
mlikRF_IG2<-RF_IG2$mlik[2]
mlikRF_IG3<-RF_IG3$mlik[2]
mlikRF_IG4<-RF_IG4$mlik[2]
mlikRF_IG5<-RF_IG5$mlik[2]
mlikRF_IG6<-RF_IG6$mlik[2]
mlikRF_IG7<-RF_IG7$mlik[2]
mlikRF_IG8<-RF_IG8$mlik[2]
mlikRF_HC1<-RF_HC1$mlik[2]
mlikRF_HC2<-RF_HC2$mlik[2]
mlikRF_HC3<-RF_HC3$mlik[2]
mlikRF_HC4<-RF_HC4$mlik[2]
mlikRF_HC5<-RF_HC5$mlik[2]
mlikRF_U<-RF_U$mlik[2]
mlikRF_PC1<-RF_PC1$mlik[2]
mlikRF_PC2<-RF_PC2$mlik[2]
mlikRF_PC3<-RF_PC3$mlik[2]
mlikRF_IG1NB<-RF_IG1NB$mlik[2]
mlikRF_IG2NB<-RF_IG2NB$mlik[2]
mlikRF_IG3NB<-RF_IG3NB$mlik[2]
mlikRF_IG4NB<-RF_IG4NB$mlik[2]
mlikRF_IG5NB<-RF_IG5NB$mlik[2]
mlikRF_IG6NB<-RF_IG6NB$mlik[2]
mlikRF_IG7NB<-RF_IG7NB$mlik[2]
mlikRF_IG8NB<-RF_IG8NB$mlik[2]
mlikRF_HC1NB<-RF_HC1NB$mlik[2]
mlikRF_HC2NB<-RF_HC2NB$mlik[2]
mlikRF_HC3NB<-RF_HC3NB$mlik[2]
mlikRF_HC4NB<-RF_HC4NB$mlik[2]
mlikRF_HC5NB<-RF_HC5NB$mlik[2]
mlikRF_UNB<-RF_UNB$mlik[2]
mlikRF_PC1NB<-RF_PC1NB$mlik[2]
mlikRF_PC2NB<-RF_PC2NB$mlik[2]
mlikRF_PC3NB<-RF_PC3NB$mlik[2]




mlik<-c(mlikRF_IG1,mlikRF_IG2,
mlikRF_IG3,
mlikRF_IG4,
mlikRF_IG5,
mlikRF_IG6,
mlikRF_IG7,
mlikRF_IG8,
mlikRF_HC1,
mlikRF_HC2,
mlikRF_HC3,
mlikRF_HC4,
mlikRF_HC5,
mlikRF_U,
mlikRF_PC1,
mlikRF_PC2,
mlikRF_PC3,
mlikRF_IG1NB,
mlikRF_IG2NB,
mlikRF_IG3NB,
mlikRF_IG4NB,
mlikRF_IG5NB,
mlikRF_IG6NB,
mlikRF_IG7NB,
mlikRF_IG8NB,
mlikRF_HC1NB,
mlikRF_HC2NB,
mlikRF_HC3NB,
mlikRF_HC4NB,
mlikRF_HC5NB,
mlikRF_UNB,
mlikRF_PC1NB,
mlikRF_PC2NB,
mlikRF_PC3NB)

MAERF_IG1<-mean(na.omit(abs(RF_IG1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG2<-mean(na.omit(abs(RF_IG2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG3<-mean(na.omit(abs(RF_IG3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG4<-mean(na.omit(abs(RF_IG4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG5<-mean(na.omit(abs(RF_IG5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG6<-mean(na.omit(abs(RF_IG6$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG7<-mean(na.omit(abs(RF_IG7$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG8<-mean(na.omit(abs(RF_IG8$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC1<-mean(na.omit(abs(RF_HC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC2<-mean(na.omit(abs(RF_HC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC3<-mean(na.omit(abs(RF_HC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC4<-mean(na.omit(abs(RF_HC4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC5<-mean(na.omit(abs(RF_HC5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_U<-mean(na.omit(abs(RF_U$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC1<-mean(na.omit(abs(RF_PC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC2<-mean(na.omit(abs(RF_PC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC3<-mean(na.omit(abs(RF_PC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG1NB<-mean(na.omit(abs(RF_IG1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG2NB<-mean(na.omit(abs(RF_IG2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG3NB<-mean(na.omit(abs(RF_IG3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG4NB<-mean(na.omit(abs(RF_IG4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG5NB<-mean(na.omit(abs(RF_IG5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG6NB<-mean(na.omit(abs(RF_IG6NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG7NB<-mean(na.omit(abs(RF_IG7NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_IG8NB<-mean(na.omit(abs(RF_IG8NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC1NB<-mean(na.omit(abs(RF_HC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC2NB<-mean(na.omit(abs(RF_HC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC3NB<-mean(na.omit(abs(RF_HC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC4NB<-mean(na.omit(abs(RF_HC4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_HC5NB<-mean(na.omit(abs(RF_HC5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_UNB<-mean(na.omit(abs(RF_UNB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC1NB<-mean(na.omit(abs(RF_PC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC2NB<-mean(na.omit(abs(RF_PC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))
MAERF_PC3NB<-mean(na.omit(abs(RF_PC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)))


MAE<-c(MAERF_IG1,MAERF_IG2,
MAERF_IG3,
MAERF_IG4,
MAERF_IG5,
MAERF_IG6,
MAERF_IG7,
MAERF_IG8,
MAERF_HC1,
MAERF_HC2,
MAERF_HC3,
MAERF_HC4,
MAERF_HC5,
MAERF_U,
MAERF_PC1,
MAERF_PC2,
MAERF_PC3,
MAERF_IG1NB,
MAERF_IG2NB,
MAERF_IG3NB,
MAERF_IG4NB,
MAERF_IG5NB,
MAERF_IG6NB,
MAERF_IG7NB,
MAERF_IG8NB,
MAERF_HC1NB,
MAERF_HC2NB,
MAERF_HC3NB,
MAERF_HC4NB,
MAERF_HC5NB,
MAERF_UNB,
MAERF_PC1NB,
MAERF_PC2NB,
MAERF_PC3NB)


MSERF_IG1<-mean((na.omit(RF_IG1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG2<-mean((na.omit(RF_IG2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG3<-mean((na.omit(RF_IG3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG4<-mean((na.omit(RF_IG4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG5<-mean((na.omit(RF_IG5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG6<-mean((na.omit(RF_IG6$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG7<-mean((na.omit(RF_IG7$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG8<-mean((na.omit(RF_IG8$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC1<-mean((na.omit(RF_HC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC2<-mean((na.omit(RF_HC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC3<-mean((na.omit(RF_HC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC4<-mean((na.omit(RF_HC4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC5<-mean((na.omit(RF_HC5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_U<-mean((na.omit(RF_U$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC1<-mean((na.omit(RF_PC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC2<-mean((na.omit(RF_PC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC3<-mean((na.omit(RF_PC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG1NB<-mean((na.omit(RF_IG1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG2NB<-mean((na.omit(RF_IG2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG3NB<-mean((na.omit(RF_IG3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG4NB<-mean((na.omit(RF_IG4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG5NB<-mean((na.omit(RF_IG5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG6NB<-mean((na.omit(RF_IG6NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG7NB<-mean((na.omit(RF_IG7NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_IG8NB<-mean((na.omit(RF_IG8NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC1NB<-mean((na.omit(RF_HC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC2NB<-mean((na.omit(RF_HC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC3NB<-mean((na.omit(RF_HC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC4NB<-mean((na.omit(RF_HC4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_HC5NB<-mean((na.omit(RF_HC5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_UNB<-mean((na.omit(RF_UNB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC1NB<-mean((na.omit(RF_PC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC2NB<-mean((na.omit(RF_PC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)
MSERF_PC3NB<-mean((na.omit(RF_PC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases))^2)


MSE<-c(MSERF_IG1,MSERF_IG2,
MSERF_IG3,
MSERF_IG4,
MSERF_IG5,
MSERF_IG6,
MSERF_IG7,
MSERF_IG8,
MSERF_HC1,
MSERF_HC2,
MSERF_HC3,
MSERF_HC4,
MSERF_HC5,
MSERF_U,
MSERF_PC1,
MSERF_PC2,
MSERF_PC3,
MSERF_IG1NB,
MSERF_IG2NB,
MSERF_IG3NB,
MSERF_IG4NB,
MSERF_IG5NB,
MSERF_IG6NB,
MSERF_IG7NB,
MSERF_IG8NB,
MSERF_HC1NB,
MSERF_HC2NB,
MSERF_HC3NB,
MSERF_HC4NB,
MSERF_HC5NB,
MSERF_UNB,
MSERF_PC1NB,
MSERF_PC2NB,
MSERF_PC3NB)



CORRF_IG1<-cor(na.omit(data.frame(yhat=RF_IG1$summary.fitted.values[,1]*DataDengue$Expected,y=DataDengue$Cases)))[1,2]
CORRF_IG2<-cor(na.omit(data.frame(yhat=RF_IG2$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG3<-cor(na.omit(data.frame(yhat=RF_IG3$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG4<-cor(na.omit(data.frame(yhat=RF_IG4$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG5<-cor(na.omit(data.frame(yhat=RF_IG5$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG6<-cor(na.omit(data.frame(yhat=RF_IG6$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG7<-cor(na.omit(data.frame(yhat=RF_IG7$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG8<-cor(na.omit(data.frame(yhat=RF_IG8$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC1<-cor(na.omit(data.frame(yhat=RF_HC1$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC2<-cor(na.omit(data.frame(yhat=RF_HC2$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC3<-cor(na.omit(data.frame(yhat=RF_HC3$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC4<-cor(na.omit(data.frame(yhat=RF_HC4$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC5<-cor(na.omit(data.frame(yhat=RF_HC5$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_U<-cor(na.omit(data.frame(yhat=RF_U$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC1<-cor(na.omit(data.frame(yhat=RF_PC1$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC2<-cor(na.omit(data.frame(yhat=RF_PC2$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC3<-cor(na.omit(data.frame(yhat=RF_PC3$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG1NB<-cor(na.omit(data.frame(yhat=RF_IG1NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG2NB<-cor(na.omit(data.frame(yhat=RF_IG2NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG3NB<-cor(na.omit(data.frame(yhat=RF_IG3NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG4NB<-cor(na.omit(data.frame(yhat=RF_IG4NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG5NB<-cor(na.omit(data.frame(yhat=RF_IG5NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG6NB<-cor(na.omit(data.frame(yhat=RF_IG6NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG7NB<-cor(na.omit(data.frame(yhat=RF_IG7NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_IG8NB<-cor(na.omit(data.frame(yhat=RF_IG8NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC1NB<-cor(na.omit(data.frame(yhat=RF_HC1NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC2NB<-cor(na.omit(data.frame(yhat=RF_HC2NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC3NB<-cor(na.omit(data.frame(yhat=RF_HC3NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC4NB<-cor(na.omit(data.frame(yhat=RF_HC4NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_HC5NB<-cor(na.omit(data.frame(yhat=RF_HC5NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_UNB<-cor(na.omit(data.frame(yhat=RF_UNB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC1NB<-cor(na.omit(data.frame(yhat=RF_PC1NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC2NB<-cor(na.omit(data.frame(yhat=RF_PC2NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]
CORRF_PC3NB<-cor(na.omit(data.frame(yhat=RF_PC3NB$summary.fitted.values[,1]*DataDengue$Expected,DataDengue$Cases)))[1,2]


COR<-c(CORRF_IG1,CORRF_IG2,
CORRF_IG3,
CORRF_IG4,
CORRF_IG5,
CORRF_IG6,
CORRF_IG7,
CORRF_IG8,
CORRF_HC1,
CORRF_HC2,
CORRF_HC3,
CORRF_HC4,
CORRF_HC5,
CORRF_U,
CORRF_PC1,
CORRF_PC2,
CORRF_PC3,
CORRF_IG1NB,
CORRF_IG2NB,
CORRF_IG3NB,
CORRF_IG4NB,
CORRF_IG5NB,
CORRF_IG6NB,
CORRF_IG7NB,
CORRF_IG8NB,
CORRF_HC1NB,
CORRF_HC2NB,
CORRF_HC3NB,
CORRF_HC4NB,
CORRF_HC5NB,
CORRF_UNB,
CORRF_PC1NB,
CORRF_PC2NB,
CORRF_PC3NB)





MAPERF_IG1<-mean(na.omit(abs(RF_IG1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG2<-mean(na.omit(abs(RF_IG2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG3<-mean(na.omit(abs(RF_IG3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG4<-mean(na.omit(abs(RF_IG4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG5<-mean(na.omit(abs(RF_IG5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG6<-mean(na.omit(abs(RF_IG6$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG7<-mean(na.omit(abs(RF_IG7$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG8<-mean(na.omit(abs(RF_IG8$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC1<-mean(na.omit(abs(RF_HC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC2<-mean(na.omit(abs(RF_HC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC3<-mean(na.omit(abs(RF_HC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC4<-mean(na.omit(abs(RF_HC4$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC5<-mean(na.omit(abs(RF_HC5$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_U<-mean(na.omit(abs(RF_U$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC1<-mean(na.omit(abs(RF_PC1$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC2<-mean(na.omit(abs(RF_PC2$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC3<-mean(na.omit(abs(RF_PC3$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG1NB<-mean(na.omit(abs(RF_IG1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG2NB<-mean(na.omit(abs(RF_IG2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG3NB<-mean(na.omit(abs(RF_IG3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG4NB<-mean(na.omit(abs(RF_IG4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG5NB<-mean(na.omit(abs(RF_IG5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG6NB<-mean(na.omit(abs(RF_IG6NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG7NB<-mean(na.omit(abs(RF_IG7NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_IG8NB<-mean(na.omit(abs(RF_IG8NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC1NB<-mean(na.omit(abs(RF_HC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC2NB<-mean(na.omit(abs(RF_HC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC3NB<-mean(na.omit(abs(RF_HC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC4NB<-mean(na.omit(abs(RF_HC4NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_HC5NB<-mean(na.omit(abs(RF_HC5NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_UNB<-mean(na.omit(abs(RF_UNB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC1NB<-mean(na.omit(abs(RF_PC1NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC2NB<-mean(na.omit(abs(RF_PC2NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100
MAPERF_PC3NB<-mean(na.omit(abs(RF_PC3NB$summary.fitted.values[,1]*DataDengue$Expected-DataDengue$Cases)/DataDengue$Cases))*100


MAPE<-c(MAPERF_IG1,MAPERF_IG2,
MAPERF_IG3,
MAPERF_IG4,
MAPERF_IG5,
MAPERF_IG6,
MAPERF_IG7,
MAPERF_IG8,
MAPERF_HC1,
MAPERF_HC2,
MAPERF_HC3,
MAPERF_HC4,
MAPERF_HC5,
MAPERF_U,
MAPERF_PC1,
MAPERF_PC2,
MAPERF_PC3,
MAPERF_IG1NB,
MAPERF_IG2NB,
MAPERF_IG3NB,
MAPERF_IG4NB,
MAPERF_IG5NB,
MAPERF_IG6NB,
MAPERF_IG7NB,
MAPERF_IG8NB,
MAPERF_HC1NB,
MAPERF_HC2NB,
MAPERF_HC3NB,
MAPERF_HC4NB,
MAPERF_HC5NB,
MAPERF_UNB,
MAPERF_PC1NB,
MAPERF_PC2NB,
MAPERF_PC3NB)



ResultAp<-data.frame(DIC=DIC, WAIC=waic,MLIK=mlik,MAE=MAE,MSE=MSE,MAPE=MAPE,COR=COR)

write.csv(ResultAp,"ResultAp.csv")

bri.hyperpar.plot(RF_IG2)


  

##### Final SIR

DataDengue$SIR<-RF_U$summary.fitted.values[,1]

Risk<-1

DataDengue$Prob<-unlist(lapply(RF_U$marginals.fitted.values, function(X){
  1-inla.pmarginal(Risk, X)
}))
 



JabarData <- fortify(Jabar1, region="ID")
DataDengue.shp<-merge(JabarData, DataDengue, by="id", all.x=TRUE)
DataDengueMap<-DataDengue.shp[order(DataDengue.shp$order), ] 
#DataDengueMap<-na.omit(DataDengueMap)


pretty_breaks <- c(0.5, 1, 1.50, 2) # find the extremes

minVal <- 0
maxVal <- 100

# compute labels
labels <- c()
brks <- c(-1, pretty_breaks, 180) # round the labels (actually, only the extremes)
for(idx in 1:length(brks)){
  labels <- c(labels,round(brks[idx + 1], 2))
}

labels <- labels[1:length(labels)-1]  # define a new variable on the data set just as above
DataDengueMap$brks <- cut(DataDengueMap$SIR, 
                          breaks = brks, 
                          include.lowest = TRUE, 
                          labels = labels)

brks_scale <- levels(DataDengueMap$brks)
labels_scale <- rev(brks_scale)


#DataDengueMap1<-subset(DataDengueMap, IT<4)

png("Fig4.png", units="in", width=14, height=8, res=600) 
 
ggplot() + 
  geom_polygon(data = DataDengueMap, aes(fill = brks, x = long,  y = lat, group = group),color="black",size=0.1) +
  facet_wrap(~Year, ncol=3,scale="free")+   
  theme_bw()+ ylab(" ")+xlab(" ")+
  theme(legend.position = "bottom")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(labels = c("(0.0-0.5]", "(0.5-1.0]","(1.0-1.5]", "(1.5-2.0]", ">2"), values =  c("cyan1","green", "orange","red1", "red3"))+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(fill = "Estimated SIR") +
  theme(legend.position="bottom")+theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        )

dev.off()





pretty_breaks <- c(0.20, 0.5, 0.80, 0.9) # find the extremes

minVal <- 0
maxVal <- 100

# compute labels
labels <- c()
brks <- c(-1, pretty_breaks, 180) # round the labels (actually, only the extremes)
for(idx in 1:length(brks)){
  labels <- c(labels,round(brks[idx + 1], 2))
}

labels <- labels[1:length(labels)-1]  # define a new variable on the data set just as above
DataDengueMap$brks <- cut(DataDengueMap$Prob, 
                          breaks = brks, 
                          include.lowest = TRUE, 
                          labels = labels)

brks_scale <- levels(DataDengueMap$brks)
labels_scale <- rev(brks_scale)


#DataDengueMap1<-subset(DataDengueMap, IT<4)

png("Fig5.png", units="in", width=14, height=8, res=600) 
 
ggplot() + 
  geom_polygon(data = DataDengueMap, aes(fill = brks, x = long,  y = lat, group = group),color="black",size=0.1) +
  facet_wrap(~Year, ncol=3,scale="free")+   
  theme_bw()+ ylab(" ")+xlab(" ")+
  theme(legend.position = "bottom")  +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(labels = c("(0.0-0.2]", "(0.2-0.5]","(0.5-0.8]", "(0.8-0.9]", ">0.9"), values =  c("cyan1","green", "orange","red1", "red3"))+  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(fill = "Exceedance probability") +
  theme(legend.position="bottom")+theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + theme(axis.text.x=element_blank(), #remove x axis labels
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.text.y=element_blank(),  #remove y axis labels
        axis.ticks.y=element_blank()  #remove y axis ticks
        )

dev.off()


ModelFinal<-RF_U

bri.hyperpar.summary(ModelFinal)
summary(ModelFinal)





bri.hyperpar.plot(RF_IG1)

bri.hyperpar.plot(RF_IG2)

bri.hyperpar.plot(RF_IG3)

bri.hyperpar.plot(RF_IG4)

bri.hyperpar.plot(RF_IG5)

bri.hyperpar.plot(RF_IG6)

bri.hyperpar.plot(RF_IG7)

bri.hyperpar.plot(RF_IG8)

bri.hyperpar.plot(RF_HC1)

bri.hyperpar.plot(RF_HC2)

bri.hyperpar.plot(RF_HC3)

bri.hyperpar.plot(RF_HC4)

bri.hyperpar.plot(RF_HC5)

bri.hyperpar.plot(RF_U)

bri.hyperpar.plot(RF_PC1)

bri.hyperpar.plot(RF_PC2)

bri.hyperpar.plot(RF_PC3) 



Temporal<-data.frame(Mean=exp(ModelFinal$summary.random$IT[,2]),
					 Min=exp(ModelFinal$summary.random$IT[,4]),
					 Max=exp(ModelFinal$summary.random$IT[,6]), Year=c(2016,2017,2018,2019,2020,2021))


png("Temporal.png", units="in", width=6, height=5, res=100)

ggplot(Temporal)+
  geom_ribbon(aes(ymin = Min, ymax = Max,x=Year, fill="95% Credible Interval"),alpha = 0.25) +
  geom_line(aes(x = Year, y = Mean), colour = "black",  size=1.5 ) +
  xlab("Year")+ylab("(b) Temporal effect")+
  scale_colour_manual(" ",values=c("blue"))+
  scale_fill_manual("",values="blue")+theme_bw()+
  theme(legend.position="bottom")+theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

 dev.off()




Districts_names <- c("1"= "Kab Bogor",       
  "2"= "Kab Sukabumi",     
  "3"= "Kab Cianjur",      
  "4"= "Kab Bandung",      
  "5"= "Kab Garut",        
  "6"= "Kab Tasikmalaya",  
  "7"= "Kab Ciamis",       
  "8"= "Kab Kuningan",     
  "9"= "Kab Cirebon",      
 "10"= "Kab Majalengka",   
 "11"= "Kab Sumedang",     
 "12"= "Kan Indramayu",    
 "13"= "Kab Subang",       
 "14"= "Kab Purwakarta",   
 "15"= "Kab Karawang",     
 "16"= "Kab Bekasi",       
 "17"= "Kab Bandung Barat",
 "18"= "Kab Ciamis",       
 "19"= "Kota Bogor",       
 "20"= "Kota Sukabumi",    
 "21"= "Kota Bandung",     
 "22"= "Kota Cirebon",     
 "23"= "Kab Bekasi",       
 "24"= "Kota depok",       
 "25"= "Kota Cimahi",      
 "26"= "Kota Tasikmalaya", 
 "27"= "Kota Banjar"  )
 


Interaction <-data.frame(Mean=exp(ModelFinal$summary.random$ID1[,2]), 
                        Min=exp(ModelFinal$summary.random$ID1[,4]), 
                        Max=exp(ModelFinal$summary.random$ID1[,6]),
                        Districts=ModelFinal$summary.random$ID1[,1], 
                        Year=rep(c(2016,2017,2018,2019,2020,2021),each=27))        



png("Interaction.png", units="in", width=8, height=6, res=100)

ggplot(Interaction)+
  geom_ribbon(aes(ymin = Min, ymax = Max,x=Year, fill="95% Credible Interval"),alpha = 0.25) +
  geom_line(aes(x = Year, y = Mean), colour = "black",  size=1.5 ) +
  xlab("Year")+ylab("(c) Interaction effect")+
  scale_colour_manual(" ",values=c("blue"))+
  scale_fill_manual("",values="blue")+theme_bw()+
  facet_wrap(~Districts,labeller = as_labeller(Districts_names), ncol=5)+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position="bottom")+theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

 dev.off()

 
 
 
 
 
setwd("/Users/mindra/@Sensitivity")

DataDengueEst<-data.frame(RR=exp(ModelFinal$summary.random$ID[,2]),id=ModelFinal$summary.random$ID[,1])


JabarData <- fortify(Jabar1, region="ID")
DataDengueEst.shp<-merge(JabarData, DataDengueEst, by="id", all.x=TRUE)
DataDengueEstMap<-DataDengueEst.shp[order(DataDengueEst.shp$order), ] 
#DataDengueEstMap<-na.omit(DataDengueEstMap)


 
#####Map discrete



pretty_breaks <- c(0.5, 1, 1.50, 2) # find the extremes

minVal <- 0
maxVal <- 100

# compute labels
labels <- c()
brks <- c(-1, pretty_breaks, 180) # round the labels (actually, only the extremes)
for(idx in 1:length(brks)){
  labels <- c(labels,round(brks[idx + 1], 2))
}

labels <- labels[1:length(labels)-1]  # define a new variable on the data set just as above
DataDengueEstMap$brks <- cut(DataDengueEstMap$RR, 
                          breaks = brks, 
                          include.lowest = TRUE, 
                          labels = labels)

brks_scale <- levels(DataDengueEstMap$brks)
labels_scale <- rev(brks_scale)

 
png("Spatial.png", units="in", width=6, height=5, res=100) 
 
ggplot() + 
  geom_polygon(data = DataDengueEstMap, aes(fill = brks, x = long,  y = lat, group = group),color="black",size=0.1) + 
  theme_bw()+ ylab(" ")+xlab(" ")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  scale_fill_manual(labels = c("(0.0-0.5]", "(0.5-1.0]","(1.0-1.5]", "(1.5-2.0]", ">2"), values =  c("cyan1","green", "orange","red", "red3"))+  
  theme(text = element_text(size=12), axis.text.x = element_text(size=8),axis.text.y = element_text(size=8)) + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
        labs(fill = "(a) Spatial effect")+theme(legend.position="bottom")
dev.off()
 
 
 



