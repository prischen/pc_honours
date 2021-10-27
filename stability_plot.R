library(stringr)
library(dplyr)
library(plotly)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(scales)

for (c in c(1:2,3)){
  mod_num = c
  #packages = c("glint","hiern","lmaic","lmsub","lmbic")
  packages = c("glint","hiern","lmbic") # for erful
  mod_name = packages[mod_num]
  #package_name = c("glinternet","hierNet","lmIntSubsets w/ AIC","lmSubsets","lmIntSubsets w/ BIC")
  package_name = c("glinternet","hierNet","lmIntSubsets w/ BIC") # for erful
  
  #glint : 1
  #hiern : 2
  #lmaic : 3
  #lmsub : 4
  #lmbic: 5 / 3 for erful
  
  mod = c()
  for (i in 1:10){
    mod = cbind(mod, get(paste0(mod_name,"_iteration-", i)))
  }
  
  err = c()
  for (i in 1:100){
    #err = cbind(err, get(paste0("toter_iteration-",i))[mod_num,])
    err = cbind(err, get(paste0("erful_iteration-",i))[mod_num,])
  }
  
  pack_res = data.frame(Model = as.vector(mod), RSS = as.numeric(err))
  pack = pack_res %>% group_by(Model) %>%
    summarise(AvgRSS = mean(RSS))
  
  df = data.frame(table(mod))
  colnames(df) = c("Model", "Frequency")
  df$Size = as.integer(str_count(df$Model, ",") + 1)
  
  if(c == 3){
    df$Size = df$Size - 1
  }
  
  count = df %>% group_by(Size) %>%
    summarise(TotalOccurences = sum(Frequency))
  
  df = merge(df, count, by = "Size")
  df$SelProb = df$Frequency/sum(df$TotalOccurences)
  df = df %>% select(-c('Frequency','TotalOccurences'))
  
  df = merge(df, pack, by = "Model")
  
  df = df[c('Model','Size','AvgRSS','SelProb')]
  
  titl = paste0("Model Stability Plot - Package: ", package_name[c])
  
  fig2 = plot_ly(df, y =~ AvgRSS, x =~ as.factor(Size), type = "scatter", 
                text = ~paste("Parameters: ", Model,"\n", "Prob: ", SelProb), 
                hoverinfo = 'text',
                marker = list(size = ~SelProb*200, 
                              color = ~SelProb,
                              reversescale = TRUE,
                              opacity = 0.5,
                              sizemin = 5, 
                              line = list(width = 0),
                              sizemax = 10,
                              colorscale = "Viridis")) %>% 
    layout(title = titl,
           xaxis = list(title = "Number of Parameters"),
           yaxis = list(title = "Information Criterion")) 

fig = ggplot(df, aes(x = Size, y = AvgRSS, col = SelProb)) + 
  geom_point(aes(size = SelProb, color = SelProb), alpha = 0.7) + 
  scale_colour_gradient(low = "orange", high = "cornflowerblue", guide = "legend", labels = scales::percent) +
  scale_size_continuous(range = c(6,20), labels = scales::percent) + 
  scale_x_continuous(breaks = seq(min(df$Size), max(df$Size),1)) +
  theme_minimal() +
  ggtitle(titl) +
  xlab("Number of Parameters") +
  ylab("Residual Sum of Squares") +
  labs(color = "Selection\nProbability", size = "Selection\nProbability") +
  theme(text = element_text(size = 15)) + 
  theme(legend.text=element_text(size=10), legend.title = element_text(size = 11))

  assign(paste0("stab_plot-",mod_num), fig)
  }
