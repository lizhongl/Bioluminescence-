library(ggplot2)
library(scales)
library(dplyr)
library(tidyr)
library(tibble)
# library(cowplot)
# library(MetaCycle)
library(lomb) #lsp



files <- list.files(pattern='_Raw.csv')  # data file to read
fnames <- gsub('_Raw.csv', '', files)


# read csv data exported from Lumicycle Analysis & normalize --------
dats <- list()  #normalized data
for (i in 1:length(files)){
  df<-read.csv(files[i], skip=3, col.names = c("Date","Time",	"Days",	"CountsPerSec",	"Baseline")) %>%
    transmute(Time=as.POSIXct(paste(Date, Time), format="%m/%d/%Y %H:%M"), #combine date time
              Hours=(Days-Days[1])*24,  #Hours
              Lum=CountsPerSec) %>% 
    filter(Hours>1.5 & Hours<36)  #subset data, remove the first 1.5 hour@@@@@@@@@@@@@@@@@@@@@@

  #detrend lum by subtracting the loess line
  loess_fit <- loess(Lum ~ Hours, data.frame(df), span=3*6/diff(range(df$Hours)))
  dats[[fnames[i]]] <- df %>% 
    mutate(Lum_detrd=Lum-predict(loess_fit)) #+mean(predict(loess_fit))) #detrend
}


#compute period by Lomb-Scargle Periodogram ----------
pers <- {}  #periods
pval <- {}
for (i in 1:length(dats)){
  d=dats[[i]]%>% select(Hours, Lum_detrd) %>% filter(Hours<24)  #subset data @@@@@@@@@@@@@@@@@@@@@@
  ls <- lsp(d, from=2, to=10, type='period', plot=F) #Lomb-Scargle Periodogram
  pval[i] <- ls$p.value[1]  #pvalue
  pers[i] <- ifelse(pval[i]<0.01, ls[['peak.at']][1], NA)  #period
  #plot Lomb-Scargle Periodogram -----------
  p <- plot(ls, main=paste0(fnames[i]), xlabel="Period (h)", ylabel="Power") +
    annotate('text', x=2, y=ls$peak, label=paste0('period=',signif(pers[i],3)), hjust=0, size=6)
  ggsave(paste0('LS_periodogram_', fnames[i], '.png'), plot=p, width=5, height=5, units='in', dpi=600)
}
d.pers <- data.frame(Sample=fnames, Period=signif(pers,4), pvalue=signif(pval, 4))
write.csv(d.pers, 'Periods.csv')

#put data in one df
dat <- bind_rows(dats, .id='name') %>% mutate(sample=gsub('_(.*)', '', name), rep=gsub('(.*)_rep', '', name))


#plot detrended data ----------
scaleFUN <- function(x) sprintf("%.f", x)  # digits

for (i in unique(dat$sample)){
  d <- dat%>%filter(sample==i)
  ymax <- max(abs(d$Lum_detrd))
  pl <- ggplot(data=d, aes(x=Hours, y=Lum_detrd, color=rep)) +  
    geom_line(linetype=1, size=0.2)+
    ylab("Luminescence") +
    xlab("Hours") + 
    ggtitle(i)+
    theme_classic() +       # Removes gridlines & background
    theme(plot.title = element_text(hjust = 0.5, size=7, margin=margin(0,0,2,0)),
          axis.title = element_text(face = "plain", color = "black", size = 6), 
          # axis.title.x=element_blank(),
          axis.text = element_text(face = "plain", color = "black", size = 5), 
          # axis.text.x = element_blank(),
          axis.line = element_line(colour = "black", size = 0.1),
          axis.ticks = element_line(colour = "black", size = 0.1),
          axis.ticks.length = unit(1, "pt"),
          # panel.border = element_rect(colour = "black", fill=NA, size=0.1),
          plot.margin = margin(1, 1, 1, 1, "pt"),
          legend.position = c(0.1,0.95),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size=5),
          legend.key.height = unit(0.01,"line"),
          legend.key.width = unit(0.2, "line"),
          legend.spacing.x = unit(0.2, "line")
          # legend.spacing.y = unit(0.005, "line")
    )+
    guides(color = ifelse(length(unique(d$rep))==1, 'none', 'legend')) +
    # scale_color_manual(values=c)+
    scale_x_continuous(expand=c(0.1,0.1))+
    scale_y_continuous(expand=c(0.1,0.1), limits=c(-ymax, ymax), breaks=pretty_breaks(n=4), labels=scaleFUN) # set y-axis labels
  # output pdf
  # ggsave(filename=paste0('Plot-Line_', i, "_detrd.pdf"), plot=pl, useDingbats=F, width=2, height=1) 
  ggsave(filename=paste0('Plot-Line_', i, "_detrd.png"), plot=pl, width=2, height=1, dpi=600) 
}

#plot raw data ----------
for (i in unique(dat$sample)){
  d <- dat%>%filter(sample==i)
  ymax <- max(abs(d$Lum))
  pl <- ggplot(data=d, aes(x=Hours, y=Lum, color=rep)) +  
    geom_line(linetype=1, size=0.2)+
    ylab("Luminescence") +
    xlab("Hours") + 
    ggtitle(i)+
    theme_classic() +       # Removes gridlines & background
    theme(plot.title = element_text(hjust = 0.5, size=7, margin=margin(0,0,2,0)),
          axis.title = element_text(face = "plain", color = "black", size = 6), 
          # axis.title.x=element_blank(),
          axis.text = element_text(face = "plain", color = "black", size = 5), 
          # axis.text.x = element_blank(),
          axis.line = element_line(colour = "black", size = 0.1),
          axis.ticks = element_line(colour = "black", size = 0.1),
          axis.ticks.length = unit(1, "pt"),
          # panel.border = element_rect(colour = "black", fill=NA, size=0.1),
          plot.margin = margin(1, 1, 1, 1, "pt"),
          legend.position = c(0.1,0.95),
          legend.title = element_blank(),
          legend.background = element_blank(),
          legend.text = element_text(size=5),
          legend.key.height = unit(0.01,"line"),
          legend.key.width = unit(0.2, "line"),
          legend.spacing.x = unit(0.2, "line")
          # legend.spacing.y = unit(0.005, "line")
    )+
    guides(color = ifelse(length(unique(d$rep))==1, 'none', 'legend')) +
    # scale_color_manual(values=c)+
    scale_x_continuous(expand=c(0.1,0.1))+
    scale_y_continuous(expand=c(0.1,0.1), limits=c(0, ymax), breaks=pretty_breaks(n=4), labels=scaleFUN) # set y-axis labels
  # output pdf
  # ggsave(filename=paste0('Plot-Line_', i, "_detrd.pdf"), plot=pl, useDingbats=F, width=2, height=1) 
  ggsave(filename=paste0('Plot-Line_', i, "_raw.png"), plot=pl, width=2, height=1, dpi=600) 
}


# #Wavelet analysis ====================
# library(WaveletComp)
# 
# for (i in 1:length(dats)){
#   d <- dats[[i]] %>% column_to_rownames(var='Hours')
#   
#   my.w = analyze.wavelet(d, "Lum",
#                          loess.span = 0,
#                          dt = 1/6,        # Time unit = Hour
#                          dj = 1/250,      # y-axis resolution
#                          lowerPeriod = 1,
#                          upperPeriod = 24,
#                          make.pval = T, 
#                          n.sim = 100)
#   
#   # plot wavelet as pdf
#   pdf(paste0("Plot-Wvlt_", fnames[i], ".pdf"), width = 3.2, height = 2)
#   par(mar=c(2.5,2.5,0.5,1), mgp=c(2,.8,0), cex.axis=.8)
#   pwt<- wt.image(my.w, n.levels = 250, # colors
#                  color.key = 'i',
#                  maximum.level = 8,
#                  plot.coi = F, 
#                  plot.contour = F,
#                  plot.ridge = F,
#                  siglvl = 0.01,
#                  timelab = "Time (Hours)",
#                  periodlab = "Period (hour)",
#                  timetcl = -0.2,
#                  periodtcl = -0.2,
#                  spec.time.axis = list(at = seq(0,24,4)*6, labels = seq(0,24,4), las = 1, hadj = NA, padj = 0),
#                  spec.period.axis = list(at = c(1,6,24), labels = T, las = 1, hadj = 1, padj = NA),
#                  # lwd = 2, 
#                  lwd.axis = .1, 
#                  legend.params = list(width=0.5, shrink=0.7, mar=2, n.ticks=2, label.digits=0, lab.line=0.1))
#   dev.off()
# }

# # plot wvlt & line (failed )
# library(ggplotify)
# pcow<-plot_grid(pl, as.ggplot(pwt), nrow = 2, rel_heights = c(1,1.2), align = 'v')
# save_plot(filename=paste0(plotname, 'WVLT-line.pdf'), pcow, base_width = 1.5, base_height = 1.2)

