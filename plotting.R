#####################################
## This file has all the plotting functions
#####################################

# libraries
require(graphics)
require(ggplot2)
library(colorspace)
#require(xlsx)
library(reshape2)
#alpha<-1860


# This function plots the TB mortality data from the Dawson 1950 paper
plotDawMortality<-function(thiswd){

  # plot daw data
  daw<-getDAWtable()
  daw_total<-daw$t1Total
  daw_years<-seq(from=1850, to=1930, by=10)
  daw_total<-cbind(daw_years, daw_total)
  daw_total_df<-melt(daw_total, id.vars = "daw_years", measure.vars =  c("age0", "age5", "age15", "age25", "age35", "age45", "age55", "age65", "age75" ) )
  colnames(daw_total_df)<-c("daw_years", "Age", "Mortality") 
  
  #cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  daw_plot<-ggplot(daw_total_df, aes(x=Age, y=Mortality, group=daw_years) ) +
    geom_line(aes(color=daw_years))+
    geom_point(aes(color=daw_years))+
    scale_colour_gradientn(colours=rainbow(4))+
    labs(title="TB Death Rates in England and Wales",x="Lower bound of dge groups (years)", y = "Death rate ( per million)")
  
  
  #save the plot as a pdf
  pdf(paste0(thiswd,"/Plots/Daw1950.pdf") )
  print(daw_plot)
  dev.off()
  
}


# plot the decaying exponential
# ### initial starting values
# beta_o<-0.1
# tau<-0.001
# gamma<-0.00001
# alpha<-1860
#flag indicates if we save a new plot as a pdf (=1) or overlay on existing plot(=0)
plotBeta_decayExp<-function(beta_o, tau,gamma, flag){
  t_values<-seq(from=1850, to=1940, by=10)
  t_values_timeshift_daily<-(t_values-1850)*365
  beta<-(beta_o * exp( - tau * ( t_values_timeshift_daily))) + gamma
  beta_func<-approxfun(t_values_timeshift_daily, beta, method="linear")
  
  if(flag==1){
    #save the plot as a pdf
    pdf(paste0(thiswd,"/Plots/betaPlot_beta_o",beta_o,"_tau",tau,".pdf") )
    print(plot(t_values_timeshift_daily, beta_func(t_values_timeshift_daily), type='b', lwd=2, pch=16, lty=1, main="time varying beta", ylim=c(0,beta_o)))
    dev.off()
    
  }
  if(flag==0){
    lines(t_values_timeshift_daily, beta_func(t_values_timeshift_daily),type='b', lwd=2, lty=2)
  }
  
}


### plot the data in table per calendar year
plotTablePerCY<-function(table, type){
  startYears<-seq(from=1851, to=1931, by=10)
  endYears<-seq(from=1860, to=1940, by=10)
  
  startAge<-c(0, 5, 15, 25,35,45,55,65,75)
  quartz()
  par(mfrow=c(1,1), oma=c(0,0,2,0))
  for(i in 1:dim(table)[1]){
    if(i==1)
      plot(startAge, table[i,], type='b', pch=16, col=i, ylab="Deaths per Million", xlab="Age (mid point)", ylim=range(table) )
    if(i>1)
      lines(startAge, table[i,], type='b', pch=16, col=i)
  }
  legend(x=55,y=max(table)-1000, legend=as.character(startYears), col=1:dim(table)[1], lwd=2 )
  title(paste("DAW TABLE of", type, "by calendar year", sep=" "), outer=TRUE)
  
}

# plot the data in the table by mean birth year
plotTableMeanBirthYear<-function(table, type){
  meanBirthYears<-seq(from=1855, to=1935, by=10)
  tableByBirthYear<-matrix(0, nrow=9, ncol=9)
  for(i in 1:dim(table)[1]){
    tmp<-table[i, 1]
    if(i<9){
      for(t in 1:(9-i)){
        tmp<-c(tmp, table[i+t,t+1 ])
      }#end of t loop
    }#end of if
    tableByBirthYear[i,1:length(tmp)]<-tmp
  }#end of i loop
  
  startAge<-c(0, 5, 15, 25,35,45,55,65,75)
  quartz()
  par(mfrow=c(1,1), oma=c(0,0,2,0))
  for(i in 1:dim(tableByBirthYear)[1]){
    if(i==1)
      plot(startAge, tableByBirthYear[i,], type='b', pch=16, col=i, ylab="Deaths per Million", xlab="Age (mid point)", ylim=range(tableByBirthYear) )
    if(i>1)
      lines(startAge, tableByBirthYear[i,], type='b', pch=16, col=i)
  }
  legend(x=55,y=max(tableByBirthYear)-1000, legend=as.character(meanBirthYears), col=1:dim(tableByBirthYear)[1], lwd=2 )
  title(paste("DAW TABLE of", type, "by mean birth year", sep=" "), outer=TRUE)
  
}

plotOutputbyDiseaseClass<-function(output_pre){
  quartz()
  par(mfrow=c(3,2), oma=c(0,0,2,0))
  
  for(c in 2:dim(output_pre)[2]){
    
    if(c %in% 2:10){
      yrange<-range(output_pre[,2:10])
      if(c==2){
        plot(output_pre$time, output_pre[,c], type='l', lty=c, col=c, lwd=2, ylab=" Numbers", ylim=yrange, main="Susceptibles")
        legend(x=output_pre$time[2]-100, y=yrange[2]-10, legend=c("0-5", "5-15", "15-25", "25-35", "35-45", "45-55", "55-65","65-75", ">75"), lty=2:10, col=2:10)
      }
      if(c>2)  
        lines(output_pre$time, output_pre[,c], lty=c, lwd=1, col=c)
    }
    
    
    if(c %in% 11:19){
      yrange<-range(output_pre[,11:19])
      if(c==11){
        plot(output_pre$time, output_pre[,c], type='l', lty=c, col=c, lwd=2, ylab=" Numbers", ylim=yrange, main="Latent_a")
        legend(x=output_pre$time[2]-100, y=yrange[2], legend=c("0-5", "5-15", "15-25", "25-35", "35-45", "45-55", "55-65","65-75", ">75"), lty=11:19, col=11:19)
      }
      if(c>11)  
        lines(output_pre$time, output_pre[,c], lty=c, lwd=1, col=c)
    }
    
    if(c %in% 20:28){
      yrange<-range(output_pre[,20:28])
      if(c==20){
        plot(output_pre$time, output_pre[,c], type='l', lty=c, col=c, lwd=2, ylab=" Numbers", ylim=yrange, main="Latent_b")
        legend(x=output_pre$time[2]-100, y=yrange[2], legend=c("0-5", "5-15", "15-25", "25-35", "35-45", "45-55", "55-65","65-75", ">75"), lty=20:28, col=20:28)
      }
      if(c>20)  
        lines(output_pre$time, output_pre[,c], lty=c, lwd=1, col=c)
    }
    
    if(c %in% 29:37){
      yrange<-range(output_pre[,29:37])
      if(c==29){
        plot(output_pre$time, output_pre[,c], type='l', lty=c, col=c, lwd=2, ylab=" Numbers", ylim=yrange, main="Infected")
        legend(x=output_pre$time[2]-100, y=yrange[2], legend=c("0-5", "5-15", "15-25", "25-35", "35-45", "45-55", "55-65","65-75", ">75"), lty=29:37, col=29:37)
      }
      if(c>29)  
        lines(output_pre$time, output_pre[,c], lty=c, lwd=1, col=c)
    }
    
    if(c %in% 38:46){
      yrange<-range(output_pre[,38:46])
      if(c==38){
        plot(output_pre$time, output_pre[,c], type='l', lty=c, col=c, lwd=2, ylab=" Numbers", ylim=yrange, main="Recovered")
        legend(x=output_pre$time[2]-100, y=yrange[2], legend=c("0-5", "5-15", "15-25", "25-35", "35-45", "45-55", "55-65","65-75", ">75"), lty=38:46, col=38:46)
      }
      if(c>38)  
        lines(output_pre$time, output_pre[,c], lty=c, lwd=1, col=c)
    }
    
    
  }# end of 'for' loop
  
  title("Output by Disease Class", outer=TRUE) 
}

plotOutputByAge<-function(output_pre){
  quartz()
  par(mfrow=c(3,3), oma=c(0,0,2,0))
  maintext<-""
  for(c in 2:10){
    if(c==2) maintext<-"Age 0-5"
    if(c==3) maintext<-"Age 5-15"
    if(c==4) maintext<-"Age 15-25"
    if(c==5) maintext<-"Age 25-35"
    if(c==6) maintext<-"Age 35-45"
    if(c==7) maintext<-"Age 45-55"
    if(c==8) maintext<-"Age 55-65"
    if(c==9) maintext<-"Age 65-75"
    if(c==10) maintext<-"Age >75"
    yrange<-range(output_pre[,c(c,c+9,c+18,c+27,c+36)])
    plot(output_pre$time, output_pre[,c],type='l', lty=c, lwd=1, col=1, ylim=yrange, ylab="Numbers", main=maintext)
    legend(x=output_pre$time[2]-100, y=yrange[2], legend=c("S", "L_a", "L_b", "I", "S_r"), lty=c, col=1:5 )
    for(i in 1:4){
      lines(output_pre$time, output_pre[,c+(i*9)], lty=c, lwd=2, col=i+1)
    }
    
  }# end of 'for' loop
  title("Output by Age Class", outer=TRUE)
}

plotAllAcceptedBeta<-function(M_accepted){
  pdf(paste0(thiswd,"/Plots/betaAllAccepted.pdf") )
  print(
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10) )
      )
  for(m in 1:dim(M_accepted)[1]){
    print(plotBeta_decayExp(M_accepted$beta_o[m], M_accepted$tau[m],gamma, 0) )
  }
  dev.off()
  
  # # plot histogram of accepted parameter values
  # pdf(paste0(thiswd,"/Plots/histAllAccepted.pdf") )
  # print(
  # par(mfrow=c(2,2), oma=c(0,0,2,0))
  # hist(M_accepted$mult)
  # hist(M_accepted$beta_o)
  # hist(M_accepted$tau)
  # )
  # dev.off()
  # 
  
}

### Function to plot the total population over time. 
plotTotalPopln<-function(output){
  total_popln<-cbind(output[,1],rowSums(output[,2:46]) )
  quartz()
  plot(total_popln[,1], total_popln[,2], type='l', lwd=2, ylim=c(max(total_popln[,2])-1000,max(total_popln[,2])+1000) )
}

# yearsFrom1940 is the number of years backwards from 1940 we want to view the plot
# yearsFrom1940=10,20, etc
plotCompartmentalDisbnByAgeGroup<-function(output_all, t_year){
  library(reshape2) # to put the dataframe in the right form for ggplot
  
  
  output_last<-output_all[ (dim(output_all)[1]), ]
  AgeClass<-1:9
  
  output_mx<-matrix(as.numeric(output_last[-1]), nrow=5, ncol=9, byrow=TRUE)
  output_total<-colSums(output_mx)
  output_proportion<-matrix(0, nrow=5, ncol=9)
  
  for(diseaseClassIndex in 1:5){
    output_proportion[diseaseClassIndex,]<-output_mx[diseaseClassIndex,]/output_total
  }
  DF<-as.data.frame( cbind(1:9, t(output_proportion) ) ) 
  colnames(DF)<-c("Ageclass", "S", "L_a", "L_b", "I", "R")
  DF1<-melt(DF, id.var="Ageclass")
  plot1<-ggplot(DF1, aes(x = Ageclass, y = value, fill = variable)) + 
    geom_bar(stat = "identity") +
    scale_x_continuous(breaks=seq(1,9,1))+
    ggtitle(paste0("year ", t_year))
  
  # quartz()
  # par(mfrow=c(1,1), oma=c(0,0,2,0))
  # plot1
  # 
  #save the plot
  pdf(paste0(thiswd,"/Plots/DiseaseDisbnByAgeGroup.pdf"))
  print(plot1)
  dev.off()

}

#function to plot each disease class across age at a point in time
#disease class Index varies from 1 to 5 ( S, L_a, L_b, I, R)
plotDiseaseClassbyAge<-function(output_all){
  
  output_last<-output_all[dim(output_all)[1], ]
  AgeClass<-1:9
  
  output_mx<-matrix(as.numeric(output_last[-1]), nrow=5, ncol=9, byrow=TRUE)
  output_total<-colSums(output_mx)
  
  
  quartz()
  par(mfrow=c(3,2), oma=c(0,0,2,0))
  for(diseaseClassIndex in 1:5){
    
    titleString<-character()
    if(diseaseClassIndex==1) titleString<-"S"
    if(diseaseClassIndex==2) titleString<-"L_a"
    if(diseaseClassIndex==3) titleString<-"L_b"
    if(diseaseClassIndex==4) titleString<-"I"
    if(diseaseClassIndex==5) titleString<-"R"
    
  #  plot(AgeClass, output_mx[diseaseClassIndex,]/output_total, type='b', lty=1, col="blue", pch=16, ylab="Numbers", xlab="Age Classes", main=paste0("Disease Class ",titleString) )
      plot(AgeClass, output_mx[diseaseClassIndex,], type='b', lty=1, col="blue", pch=16, ylab="Numbers", xlab="Age Classes", main=paste0("Disease Class ",titleString) )
    
  }
  titleString<-" ALL "
  plot(AgeClass, output_total, type='b', lty=1, col="blue", pch=16, ylab="Numbers", xlab="Age Classes", main=paste0("Disease Class ",titleString) )
  
  
  
}

