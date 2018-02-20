###################################################################
#Soil gas budget, gas diffusion, soil respiration
#
#by Raphael Menke & Laurin Freiberg
###################################################################


# Respiration ####
hi
#create function to calculate soil respiration for different depths
novak <- function(l_s=30, 
                  depth, 
                  tree = "beech", 
                  temp=10){
  
  #parameters for spruce 1.Q10 2.t_ref 3.res_t_ref
  spruce_para <- c( 3.3, 9.3, 1.69)
  #parameters for beech 1.Q10 2.t_ref 3.res_t_ref
  beech_para <- c(3.8, 9.5,1.27 )
  
  #choose parameter set for beech or spruce
  if (tree == "beech") {params <- beech_para 
  }else if (tree=="spruce"){ params <-  spruce_para
  #stop function if tree is not beech or spruce
  }else{
    stop("tree has to be 'beech' or 'spruce'")}
  
  # formula for temperature dependent respiration (vesterdal et al. 2011)
  resp <- 10^((log10(params[1])*(temp-params[2]))/10 + log10(params[3]))#mymol/m2/s
  
  #chage unit from mymol/m2/s to cm3/cm2/s with ideal gas law
  #universal gas constant R
  R<-8.314#J/mol/K=Nm/mol/K
  #amount of molecules n
  n<-resp/1e6/1e4#mol/cm2/s
  #pressure p
  p<-101.3*1000#Pa#N/m2
  #calculate Volume by converting p*V=n*R*t 
  V<-n*R*(temp+273.15)/p #m3/cm2/s
  #respiration rate vco2
  vco2<-V*1e6#cm3/cm2/s
  
  #formula for depths dependent respiration rate in soil (novak 2007)
  s_10 <-  vco2/l_s
  s_c<-matrix(nrow=length(depth),ncol=length(s_10))
  for (i in 1:length(depth)){
    s_c[i,] <- s_10 *exp(-depth[i]/l_s)*params[1]^((temp - 10)/10)
  }
  return(s_c)
}


##############################################################
#Diffusion function
##############################################################

#sinoidal oscillating temperature with an average of 10?C and an amplitude of 10
t<-sinpi(seq(-0.5,1.5,len=365))*10+10

co2_soil_depth <- function(epsilon=.2, 
                           timestep=30, 
                           max_depth=250, 
                           z=10, 
                           ambient=0.0004, 
                           temp=t, 
                           tree_c="beech", 
                           outputresolution=86400, #in sek=one day
                           total_t= 365, #time of interest in days
                           ls=30){
  
  ps <- novak(tree=tree_c, depth = seq(z,max_depth,z)-z/2, temp=temp,l_s=ls)
  
  J<-total_t*86400/timestep
  n <- max_depth/z
  output<-matrix(nrow=J*timestep/outputresolution,ncol=n) #result matrix
  result <- rep(ambient, n)
  resultcounter<-1
  Ds20<- 0.496*eps^1.661 #transfer-function according to Schack-Kirchner et al. 
  Ds_t<-Ds20*((temp+273)/293)^1.72 #temperature correction Currie ..
  
  for (j in seq(1,J)){
    time<-j*timestep
    day<-rep(1:365,ceiling(total_t/365))
    p<-ps[,day[ceiling(time/86400)]]
    Ds<-Ds_t[day[ceiling(time/86400)]]
    c1 <- result[1]  + timestep * Ds/epsilon*
      (result[2]- 3*result[1]+2*ambient)/(.75*z^2) + p[1]*z *timestep/epsilon
    c2 <- result[2:(n-1)] + timestep * Ds/epsilon*
      (result[3:n]- 2*result[2:(n-1)]+result[1:(n-2)])/(.75*z^2) + 
      p[2:(n-1)] *timestep*z/epsilon
    c3 <- result[n] + timestep * Ds/epsilon*
      (0- result[n] + result [n-1])/(z^2) + p[n] *timestep*z/epsilon
    result <- c(c1,c2,c3)
    
    if(time%%outputresolution==0){
      output[resultcounter,]<-result
      resultcounter<-resultcounter+1}
    if ((j/J*100)%%1==0){
      print(paste((j/J)*100,"% complete"))}
  }
  return(output)}

##############################################################
#run modell

out_30<-co2_soil_depth(total_t = 365*2,ls=30)
out_20<-co2_soil_depth(total_t = 365*2,ls=20)
out_50<-co2_soil_depth(total_t = 365*2,ls=50)

###############################################################
#plot output
###############################################################

matplot(t(out[seq(1,length(out[1,]),len=15),]),type="l")
matplot(rep(1:365,2),out,type="l")

image(out_30)
image(out_20)
image(out_50)

library(ggplot2)
data_30<-data.frame(co2=as.vector(out_30[366:730,]),day=rep(1:365,25),depth=-rep(seq(10,250,10)-10/2,each=365))
ggplot(data_30,aes(day,depth,fill=co2))+geom_tile()+geom_contour(aes(z=co2),col="white")+scale_fill_gradient2(low="blue",mid = "yellow",high = "red",midpoint = mean(out_30)+4e-2)+theme_minimal()

data_20<-data.frame(co2=as.vector(out_20[366:730,]),day=rep(1:365,25),depth=-rep(seq(10,250,10)-10/2,each=365))
ggplot(data_20,aes(day,depth,fill=co2))+geom_tile()+geom_contour(aes(z=co2),col="white")+scale_fill_gradient2(low="blue",mid = "yellow",high = "red",midpoint = mean(out_30)+4e-2)+theme_minimal()

data_50<-data.frame(co2=as.vector(out_50[366:730,]),day=rep(1:365,25),depth=-rep(seq(10,250,10)-10/2,each=365))
ggplot(data_50,aes(day,depth,fill=co2))+geom_tile()+geom_contour(aes(z=co2),col="white")+scale_fill_gradient2(low="blue",mid = "yellow",high = "red",midpoint = mean(out_30)+4e-2)+theme_minimal()
