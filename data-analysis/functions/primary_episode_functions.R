# Author: Bryan Mayer
# Created: 10.28.2014

## This script contains functions to:
#1.) extract summary information from the shedding episodes
#2.) 


# -------- Primary episode processing ----------------

#episodeData makes a data frame listing the episodes and basic attributes
#otherfunction turns episodeData into time series for further analysis
#new three phases (growth, transition, clearance)
primary_episode_fun = function(dataIn){
  dataIn = arrange(dataIn, days2)
  
  peak = max(dataIn$count)
  
  peak_range = which((peak - dataIn$count) <= 1) #find counts within window of measured peak
  peak_days = dataIn$days2[peak_range]
  
  outdata = with(dataIn, data.frame(
    PatientID2 = PatientID2[1],
    Virus = Virus[1],
    startTime = days2[1],
    endTime = tail(days2, 1),
    firstpos = count[1],
    lastpos = tail(count, 1),
    totalSwabs = sum(!is.na(count2)),
    peak = max(count2, na.rm = T),
    peakphase_start = min(peak_days),
    whenmax = days2[which(count2 == max(count2, na.rm = T))],
    peakphase_end = days2[which(count2 == max(count2, na.rm = T))],
    peak_shedding_days = max(peak_days) - min(peak_days),
    peak_shedding_days_end = max(peak_days),
    duration = tail(days2, 1) - days2[1],
    relative_DOB = days2[1] - as.numeric(days_dob[1])
  ))
  
  outdata
  
}

#this takes all of the raw time series data and creates episode data by patient ID
create_primary_episodes_fun = function(inData){
  
  episodeOut = plyr::ldply(unique(inData$PatientID2), function(id){
    
    tempdata = subset(inData, PatientID2 == as.character(id))
    temp_episode = primary_episode_fun(tempdata)

    temp_episode$momhiv = tempdata$momhiv[1]
    
    return(temp_episode)
  })
  
  episodeOut
  
}

change_firstpos_fun = function(episodeDataTS){

  ldply(unique(episodeDataTS$episodeGroup), function(episodeN){
    tempdata = subset(episodeDataTS, episodeGroup == episodeN)
    
    #only work with uncensored data 
    if(tempdata$censor[1] != "uncensored" | dim(tempdata)[1] == 1) return(NULL)
    
    timeinterval_fromfirst = tail(tempdata$days, -1) - min(tempdata$days)
    deltacount = tail(tempdata$count, -1) -  head(tempdata$count, -1)
    
    data.frame(
      PatientID2 = tempdata$PatientID2[1],
      idpar = tempdata$idpar2[1],
      EpisodeID = episodeN,
      starting_count = head(tempdata$count, -1),
      time_from_first = timeinterval_fromfirst,
      deltacount = deltacount)
    
  })
  
  change_firstpos = change_firstpos %>% group_by(EpisodeID) %>% dplyr::mutate(lastpos = (time_from_first == max(time_from_first)))
  change_firstpos$lastpos2 = ifelse(change_firstpos$lastpos, 1, NA)

return(change_firstpos)
}

change_lastpos_fun = function(episodeDataTS){
  
  ldply(unique(episodeDataTS$episodeGroup), function(episodeN){
    tempdata = subset(episodeDataTS, episodeGroup == episodeN)
    
    #only work with uncensored data 
    if(tempdata$censor[1] != "uncensored" | dim(tempdata)[1] == 1) return(NULL)
    
    timeinterval_fromlast = -rev(head(tempdata$days, -1) - max(tempdata$days))
    deltacount = tail(rev(tempdata$count), -1) -  head(rev(tempdata$count), -1)
    
    data.frame(
      PatientID2 = tempdata$PatientID2[1],
      EpisodeID = episodeN,
      idpar = tempdata$idpar2[1],
      starting_count = rev(tail(tempdata$count, -1)),
      time_from_last = timeinterval_fromlast,
      deltacount = deltacount) 
  })
  
  change_lastpos = change_lastpos %>% group_by(EpisodeID) %>% dplyr::mutate(firstpos = (time_from_last == max(time_from_last)))
  change_lastpos$firstpos2 = ifelse(change_lastpos$firstpos, 1, NA)

  return(change_lastpos)
}

plot_first_peak_last_blip = function(inputEpisodeData, hivsplit = F, virus = "CMV", gtitle = NULL, ylimits = c(1.8, 5.75)){
  
  first_last_posdata = melt(inputEpisodeData,
                            measure.vars = c("firstposR", "peak", "lastposR", "blipValue"), 
                            value.name = "count", variable.name = "category")
  
  if(hivsplit) {
    first_last_posdataTXT = ddply(first_last_posdata, .(category, momhiv), summarize, total = sum(!is.na(count)))
    plot = ggplot(data = first_last_posdata, aes(x = category, y = count, fill = momhiv)) + 
      geom_boxplot(alpha = 0.5) +
      scale_fill_manual("Mom HIV", values = c("khaki1", "grey")) +
      scale_y_continuous(paste(virus, "DNA concentration"), breaks = 2:9, limits = ylimits) +
      scale_x_discrete("", breaks = c("firstposR", "peak", "lastposR", "blipValue"),
                       labels = c("First +", "Peak", "Last +", "Blip")) +
      geom_vline(x = 3.5) +
      geom_hline(yintercept = 2.176, linetype = "dashed", size = 1) + 
      geom_text(data = first_last_posdataTXT, aes(x = category, y = 2, label = total), 
                family = "Times", position = position_dodge(width = 0.75)) +
      ggtitle(gtitle)
      theme(text = element_text(family = "Times", size = 20))
  } else {
    first_last_posdataTXT = ddply(first_last_posdata, .(category), summarize, total = sum(!is.na(count)))
    plot = ggplot(data = first_last_posdata, aes(x = category, y = count)) + geom_boxplot(fill = "grey", alpha = 0.5) +
      scale_y_continuous(paste(virus, "DNA concentration"), breaks = 2:9, limits = ylimits) +
      scale_x_discrete("", breaks = c("firstposR", "peak", "lastposR", "blipValue"),
                       labels = c("First +", "Peak", "Last +", "Blip")) +
      geom_vline(x = 3.5) +
      geom_hline(yintercept = 2.176, linetype = "dashed", size = 1) + 
      geom_text(data = first_last_posdataTXT, aes(x = category, y = 2, label = total), family = "Times") +
      theme(text = element_text(family = "Times", size = 20))
  }
  
  return(plot)
  
}

biphasicFit_function = function(dataIn){
  #sort by days so things are in chronological order
  dataIn = arrange(dataIn, days2)
  
  peak = max(dataIn$count)
  whenPeakInd = with(dataIn, which(count ==  peak))
  peak_range = which((peak - dataIn$count) <= 1) #find counts within one log of measured peak
  
  dataIn$count[dataIn$count == 0] = NA
  
  #Alt models
  model1A = lm(count ~ days2, data = dataIn[1:whenPeakInd, ])
  model2A = lm(count ~ days2, data = dataIn[max(peak_range):dim(dataIn)[1], ])
  
  model1 = lm(count ~ days2, data = dataIn[1:min(peak_range), ])
  model2 = lm(count ~ days2, data = dataIn[whenPeakInd:dim(dataIn)[1], ])
  
  with(dataIn, data.frame(
    PatientID2 = PatientID2[1],
    PatientID22 = PatientID22[1],
    rsq_topeak = round(summary(model1)$r.sq, 2),
    rsq_growth = round(summary(model1A)$r.sq, 2),
    rsq_GrowthToPeak = round(summary(model1)$r.sq/summary(model1A)$r.sq, 2),
    rsq_frompeak = round(summary(model2)$r.sq, 2),
    rsq_clr = round(summary(model2A)$r.sq, 2),
    rsq_ClrToPeak = round(summary(model2)$r.sq/summary(model2A)$r.sq, 2)
  ))
}
  
triphasicFit_function = function(dataIn, old = F){
  
  #sort by days so things are in chronological order
  dataIn = arrange(dataIn, days2)
  
  if(old){temp_episode = primary_episode_fun_old(dataIn)
  }else temp_episode = primary_episode_fun(dataIn)  
  
  totaltimes = dim(dataIn)[1]
  peak = max(dataIn$count)
  first_peak = with(temp_episode, which(dataIn$days2 == peakphase_start))
  last_peak =  with(temp_episode, which(dataIn$days2 == peakphase_end))
  
  #remove 0
  dataIn$count[dataIn$count == 0] = NA
    
  model1 = lm(count ~ days2, data = dataIn[1:first_peak, ])
  model2 = lm(count ~ days2, data = dataIn[first_peak:last_peak, ])
  model3 = lm(count ~ days2, data = dataIn[last_peak:totaltimes, ])
  
  output = with(dataIn, data.frame(
    PatientID2 = PatientID2[1],
    PatientID2 = PatientID2[1],
    momhiv = momhiv[1],
    age = days_dob[1],
    growth_i = first_peak,
    middle_i = last_peak - first_peak,
    decay_i = totaltimes - last_peak,
    peak = peak,
    first_peak_day = days2[first_peak],
    peak_day = temp_episode$whenmax,
    last_peak_day = days2[last_peak],
    high_shed_length =  days2[last_peak] - days2[first_peak],
    decay_length = days2[totaltimes] - days2[last_peak],
    total_days = days2[totaltimes],
    growth = unname(model1$coef[2]),
    middle = unname(model2$coef[2]),
    decay = unname(model3$coef[2]),
    rsq_1 = summary(model1)$r.sq, 
    rsq_2 = summary(model2)$r.sq,
    rsq_3 = summary(model3)$r.sq,
    double = log10(2)/unname(model1$coef[2]),
    half = log10(0.5)/unname(model3$coef[2])
    #rsq_avg = (min(peak_range)*summary(model1)$r.sq + (max(peak_range) - min(peak_range) + 1)*summary(model2)$r.sq + (totaltimes - max(peak_range) + 1)*summary(model3)$r.sq) / (totaltimes + 2)
  ))

  output
}

triphasicCI_function = function(dataIn, old = F){
  
  #sort by days so things are in chronological order
  dataIn = arrange(dataIn, days2)
  
  if(old){temp_episode = primary_episode_fun_old(dataIn)
  }else temp_episode = primary_episode_fun(dataIn)  
  
  totaltimes = dim(dataIn)[1]
  peak = max(dataIn$count)
  first_peak = with(temp_episode, which(dataIn$days2 == peakphase_start))
  last_peak =  with(temp_episode, which(dataIn$days2 == peakphase_end))
  

  #remove 0
  dataIn$count[dataIn$count == 0] = NA
  dataIn$weeks = dataIn$days2/7
  
  model1 = lm(count ~ weeks, data = dataIn[1:first_peak, ])
  model2 = lm(count ~ weeks, data = dataIn[first_peak:last_peak, ])
  model3 = lm(count ~ weeks, data = dataIn[last_peak:totaltimes, ])
  
  growth = round(unname(model1$coef[2]), 2)
  growthL = round(confint(model1)[2,1], 2)
  growthH = round(confint(model1)[2,2], 2)
  
  high = round(unname(model2$coef[2]), 2)
  highL = round(confint(model2)[2,1], 2)
  highH = round(confint(model2)[2,2], 2)
  
  decay = round(unname(model3$coef[2]), 2)
  decayL = round(confint(model3)[2,1], 2)
  decayH = round(confint(model3)[2,2], 2)
  
  growthF = paste(growth, " (", growthL, ", ", growthH, ")", sep = "")
  highF = paste(high, " (", highL, ", ", highH, ")", sep = "")
  decayF = paste(decay, " (", decayL, ", ", decayH, ")", sep = "")
  
  output = with(dataIn, data.frame(
    PatientID2 = PatientID2[1],
    ng = as.character(first_peak),
    growth = growthF,
    rsq_1 = round(summary(model1)$r.sq, 2),
    nh = as.character(last_peak - first_peak + 1),
    high = highF,
    rsq_2 = round(summary(model2)$r.sq, 2),
    nd = as.character(totaltimes - last_peak + 1),
    decay = decayF,
    rsq_3 = round(summary(model3)$r.sq, 2)
    #rsq_avg = (min(peak_range)*summary(model1)$r.sq + (max(peak_range) - min(peak_range) + 1)*summary(model2)$r.sq + (totaltimes - max(peak_range) + 1)*summary(model3)$r.sq) / (totaltimes + 2)
  ))
  
  output
}

swab2swab_function = function(dataIn, old = F){
  #sort by days so things are in chronological order
  dataIn = arrange(dataIn, days2)
  
  if(old){temp_episode = primary_episode_fun_old(dataIn)
  }else temp_episode = primary_episode_fun(dataIn)  
  
  totaltimes = dim(dataIn)[1]
  peak = max(dataIn$count)
  
  #remove 0
  dataIn$count[dataIn$count == 0] = NA
  
  
  output = with(dataIn, data.frame(
    PatientID = PatientID[1],
    momhiv = momhiv[1],
    age = days_dob[1],
    days = head(days2, -1),
    days_change = tail(days2, -1) - head(days2, -1),
    change = tail(count, -1) - head(count, -1),
    change_std = 7 * (tail(count, -1) - head(count, -1))/(tail(days2, -1) - head(days2, -1))
    ))
  
  output$phase = sapply(output$days, function(time){
    if(time < temp_episode$peakphase_start) return("growth")
    else if (time < temp_episode$peakphase_end) return("middle")
    else return("decay")
  })
  
  output$phase = factor(output$phase, levels = c("growth", "middle", "decay"), ordered = T)
    
  output
}





tri_plots_fun = function(triphasicFit) {
  llply(unique(triphasicFit$PatientID), function(id){
    temp = subset(CMVPrimaryEpisodes, PatientID == id)
    #recode zeros as missing
    temp$count[temp$count == 0] = NA
    
    epfit = subset(triphasicFit, PatientID == id)
    
    #using max avg r2
    model1 = lm(count ~ days2, data = temp[1:epfit$growth_i, ])
    model2 = lm(count ~ days2, data = temp[epfit$growth_i:epfit$middle_i, ])
    model3 = lm(count ~ days2, data = temp[epfit$middle_i:dim(temp)[1], ])
    
    testdata1 = data.frame(
      days2 = temp$days2[1:epfit$growth_i])
    testdata2 = data.frame(
      days2 = temp$days2[epfit$growth_i:epfit$middle_i])
    testdata3 = data.frame(
      days2 = temp$days2[epfit$middle_i:dim(temp)[1]])
    
    testdata1$count = predict(model1, testdata1)
    testdata2$count = predict(model2, testdata2)
    testdata3$count = predict(model3, testdata3)
    
    max_rsq = with(epfit, round(c(rsq_1, rsq_2, rsq_3, rsq_avg), 2))
    
    ggplot(data = temp, aes(x = days2, y = count)) + geom_point() + geom_line() +
      geom_smooth(method = "loess", se = F, colour = "black") +
      geom_line(data = testdata1, aes(x = days2, y = count), colour = "blue", alpha = 0.8, size = 2) +
      geom_line(data = testdata2, aes(x = days2, y = count), colour = "blue", alpha = 0.8, size = 2) +
      geom_line(data = testdata3, aes(x = days2, y = count), colour = "blue", alpha = 0.8, size = 2) +
      annotate("text", x = -Inf, y = -Inf, label = paste("best fit", max_rsq[1], max_rsq[2], max_rsq[3], max_rsq[4], sep = "\n"), 
               hjust = -0.1, vjust = -0.2, colour = "blue", size = 5) +
      ggtitle(id)
  })
  #ggsave(paste(figureDir, "trifit.pdf", sep = ""), do.call(marrangeGrob, c(tri_plots, list(nrow = 2, ncol = 2))))
}


quadratic_decay_fit = function(pid){
  
  
  if(with(subset(episodeSummary, PatientID == pid), duration == peakphase_end)) return(NULL)
  
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  decay_data$days_sq = decay_data$days2 * decay_data$days2
  model1 = lm(count ~ days2, data = decay_data)
  model2 = lm(count ~ days2 + days_sq, data = decay_data)

  last_day =  max(decay_data$days2)
  minPoint = -unname(model2$coef[2])/(2 * unname(model2$coef[3]))
  last_model_count = predict(model2, data.frame(days2 = last_day, days_sq = last_day^2))
  
  data.frame(
    PatientID = pid,
    simple_decay = unname(model1$coef[2]),
    simple_rsq = summary(model1)$r.sq,
    quad_decay1 = unname(model2$coef[2]),
    quad_decay2 = unname(model2$coef[3]),
    quad_pvalue = as.data.frame(summary(model2)$coef)$Pr[3],
    quad_rsq = summary(model2)$r.sq,
    peak_day = with(subset(episodeSummary, PatientID == pid), whenmax),
    peak_end = with(subset(episodeSummary, PatientID == pid), peakphase_end),
    minPoint = minPoint,
    last_day = last_day,
    min_count = min(decay_data$count),
    last_count = tail(decay_data$count, 1),
    min_model_count = ifelse(minPoint < last_day, 
                             predict(model2, data.frame(days2 = minPoint, days_sq = minPoint^2)), 
                             last_model_count),
    last_model_count = last_model_count
  )

}

quadratic_decay_fit_table =  function(pid){
  
  
  if(with(subset(episodeSummary, PatientID == pid), duration == peakphase_end)) return(NULL)
  
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  
  decay_data$weeks = decay_data$days2/7
  decay_data$weeks_sq = decay_data$weeks * decay_data$weeks
  model1 = lm(count ~ weeks, data = decay_data)
  model2 = lm(count ~ weeks + weeks_sq, data = decay_data)
  
  last_day =  max(decay_data$days2)
  minPoint = -unname(model2$coef[2])/(2 * unname(model2$coef[3])) * 7
  #last_model_count = predict(model2, data.frame(days2 = last_day, days_sq = last_day^2))
  quad_est = round(1000*unname(model2$coef[3]), 2)
  quad_estL = round(1000*confint(model2)[3,1], 2)
  quad_estH = round(1000*confint(model2)[3,2], 2)
  
  quadF = paste(quad_est, " (", quad_estL, ", ", quad_estH, ")", sep = "")
    
    
  data.frame(
    PatientID2 = decay_data$PatientID2[1],
    quad_est = quadF,
    #quad_lower = confint(model2)[3,1],
    #quad_upper = confint(model2)[3,2],
    quad_rsq = round(summary(model2)$r.sq, 2),
    minPoint = as.character(round(minPoint, 0)),
    last_day = as.character(last_day),
    diff = if(unname(model2$coef[3]) > 0) round(last_day - minPoint, 0) else NA,
    rebound = minPoint < last_day
  )
  
}

plot_quad_fit = function(pid){
  
  if(with(subset(episodeSummary, PatientID == pid), duration == peakphase_end)) return(NULL)
  
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  decay_data$days_sq = decay_data$days2 * decay_data$days2
  model1 = lm(count ~ days2, data = decay_data)
  model2 = lm(count ~ days2 + days_sq, data = decay_data)
  
  days = seq(min(decay_data$days2), max(decay_data$days2), length = 100)
  temp = data.frame(
    days2 = days,
    days_sq = days * days)
  
  predData1 = data.frame(
    days2 = days,
    count = predict(model1, temp),
    model = "simple")
  predData2 = data.frame(
    days2 = days,
    count = predict(model2, temp),
    model = "quadratic")
  
  predData = rbind(predData1, predData2)
  
  ggplot(data = decay_data, aes(x = days2, y = count)) + geom_point() +
    geom_line(data = predData, aes(colour = model)) +
    geom_vline(xintercept = with(decay_data, days2[which(count == min(count))]),
              linetype = "dotted") +
    geom_vline(xintercept = with(subset(predData, model == "quadratic"), days2[which(count == min(count))])) +
    scale_colour_discrete(guide = F) + ggtitle(pid)
  
}

#cubic fit 
cubic_decay_fit = function(pid){
  
  if(with(subset(episodeSummary, PatientID == pid), duration == peakphase_end)) return(NULL)
  
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  decay_data$days_sq = decay_data$days2 * decay_data$days2
  decay_data$days_cb = decay_data$days2 * decay_data$days2 * decay_data$days2
  
  model2 = lm(count ~ days2 + days_sq + days_cb, data = decay_data)
    
  last_day =  max(decay_data$days2)
  
  quad_test = function(a, b, c) (b^2 - 4*a*c) > 0
  
  quad_form = function(a,b,c) (-b + c(-1, 1) * sqrt(b^2 - 4*a*c)) / (2 * a)
  
  a = 3 * unname(model2$coef[4])
  b = 2 * unname(model2$coef[3])
  c = unname(model2$coef[2])
  
  if(quad_test(a,b,c)) critPoints = quad_form(a,b,c) else critPoints = NA
  
  inf_point = -b / (2 * a)
  
  last_model_count = predict(model2, data.frame(days2 = last_day, days_sq = last_day^2, days_cb = last_day^3))
  
    data.frame(
    PatientID = pid,
    PatientID2 =  PatientID_key$newID[ which(PatientID_key$originalID == as.character(pid))],
    cube_decay1 = unname(model2$coef[2]),
    cube_decay2 = unname(model2$coef[3]),
    cube_decay3 = unname(model2$coef[4]),
    cube_rsq = summary(model2)$r.sq,
    peak_day = with(subset(episodeSummary, PatientID == pid), whenmax),
    peak_end = with(subset(episodeSummary, PatientID == pid), peakphase_end),
    critPoint1 = critPoints[1],
    minPoint1 = (2*a*critPoints[1] + b) > 0,
    critPoint2 = critPoints[2],
    minPoint2 = (2*a*critPoints[2] + b) > 0,
    inf_point = inf_point,
    concave_after = (2 * a * (inf_point + 1) +  b) < 0,
    last_day = last_day,
    min_count = min(decay_data$count),
    last_count = tail(decay_data$count, 1),
    last_model_count = last_model_count
  )
  
}

cubic_decay_fit_table =  function(pid){

  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  
  when_peak = subset(episodeSummary, PatientID == pid)$whenmax
  
  decay_data$days_sq = decay_data$days2 * decay_data$days2
  decay_data$days_cb = decay_data$days2 * decay_data$days2 * decay_data$days2
  model2 = lm(count ~ days2 + days_sq + days_cb, data = decay_data)
  
  quad_test = function(a, b, c) (b^2 - 4*a*c) > 0
  quad_form = function(a,b,c) (-b + c(-1, 1) * sqrt(b^2 - 4*a*c)) / (2 * a)
  second_der = function(x, a, b) 2 * a * x + b
  
  #note the cubic polynomial coef are already included
  a = 3 * unname(model2$coef[4])
  b = 2 * unname(model2$coef[3])
  c = unname(model2$coef[2])
  inf_point = -b / (2 * a)
  if(quad_test(a,b,c)) {
    critPoints = quad_form(a,b,c) 
    which_max = which((2*a*critPoints + b) > 0)
    }else {
      critPoints = NA
      which_max = NA
    }
  
  rebound = if(second_der(inf_point + 10, a, b) > 0) {
    if(max(decay_data$days2) - critPoints[which_max] > 0) "**" else "*"
  }else " "
  
  data.frame(
    PatientID2 = decay_data$PatientID2[1],
    decel =  "*",
    rebound = rebound, 
    inflection_point = as.character(round(inf_point, 0)),
    switch_pos = second_der(inf_point + 10, a, b) > 0,
    min_point = round(critPoints[which_max], 0),
    min_point_peak = round(critPoints[which_max], 0) - round(when_peak, 0),
    final_day = (round(max(decay_data$days2), 0))
    )
  
}


plot_cube_fit = function(pid){
  
  if(with(subset(episodeSummary, PatientID == pid), duration == peakphase_end)) return(NULL)
  
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)
  decay_data$days_sq = decay_data$days2 * decay_data$days2
  decay_data$days_cb = decay_data$days2 * decay_data$days2 * decay_data$days2
  
  model2 = lm(count ~ days2 + days_sq + days_cb, data = decay_data)
  
  a = 3 * unname(model2$coef[4])
  b = 2 * unname(model2$coef[3])
  
  days = seq(min(decay_data$days2), max(decay_data$days2), length = 100)
  temp = data.frame(
    days2 = days,
    days_sq = days * days,
    days_cb = days * days * days)
  
  predData = data.frame(
    days2 = days,
    count = predict(model2, temp)
    )
    
  ggplot(data = decay_data, aes(x = days2, y = count)) + geom_point() +
    geom_line(data = predData) +    
    geom_vline(xintercept = -b / (2 * a), linetype = "dotted") +
    geom_vline(xintercept = with(predData, days2[which(count == min(count))])) +
    scale_colour_discrete(guide = F) + ggtitle(pid)
  
}

plot_loess_fit = function(pid){
  decay_data = subset(CMVPrimaryEpisodes, PatientID == pid & count > 0 &
                        days2 >= subset(episodeSummary, PatientID == pid)$peakphase_end)

  lo_spline = with(decay_data, loess(count ~ days2, span = 1, degree = 2))
  
  predData = data.frame(
    days2 = with(decay_data, seq(min(days2), max(days2), 1)),
    count = predict(lo_spline, with(decay_data, seq(min(days2), max(days2), 1)))
  )
  
  ggplot(data = decay_data, aes(x = days2, y = count)) + geom_point() +
    geom_line(data = predData) +
    geom_vline(xintercept = with(decay_data, days2[which(count == min(count))]),
               linetype = "dotted") +
    geom_vline(xintercept = with(predData, days2[which(count == min(count))])) +
    scale_colour_discrete(guide = F) + ggtitle(pid)
}


## correlation boostrap ##

boot_fun = function(r, dataIn, conf = 0.95){
  
  pair_cor = function(d){
    cor(d[,1], d[,2], method = "spearman")
  }
  
  size = dim(dataIn)[1]
  output = rep(NA, length = r)
  for(i in 1:r){
    sample_set = dataIn[sample(size, replace = T) ,]     
    output[i] = pair_cor(sample_set)
  }
  lowerCI = quantile(output, probs = (1 - conf)/2)
  upperCI = quantile(output, probs = (1 - (1 - conf)/2))
  sig = if((lowerCI < 0 & upperCI < 0) | (lowerCI > 0 & upperCI > 0)) "*" else " "
  
  group1 = as.numeric(unlist(dataIn[,1]))
  group2 =  as.numeric(unlist(dataIn[,2]))
  base_test= cor.test(group1, group2, method = "spearman")$p.value
  
  data.frame(
    pair1 = names(dataIn)[1],
    pair2 = names(dataIn)[2],
    base_cor = pair_cor(dataIn),
    base_pvalue = as.numeric(base_test),
    base_test = if(base_test < 0.05) "*" else " ",
    boot_mean = mean(output),
    boot_median = median(output),
    confidence = conf,
    boot_lowerCI = lowerCI,
    boot_upperCI = upperCI,
    sig = sig
  )
}
