# p004
# reads in data.csv, runs ANOVA for FFA/PPA, runs post-hocs, makes plots
# packages (run separately) ----------------------------------------------------------------
library('reshape2'); library('plotrix'); library('ggplot2'); library('ez'); library('stringr'); library('grid')
library('plyr'); library('lsr'); library('nlme'); library('multcomp'); library('broom')
library(extrafont)
library(extrafontdb)
library(Rttf2pt1)
loadfonts()

# subject ages, useful variables and functions ------------------------------------------------------------------

projDir = '/Users/david/Work/research/p004'
data <- read.csv(file.path(projDir,'data/data.csv'), header = T)
  
# subject ages are stored here
ages = c(25,34,23,20,23,30,25,16,21,20,22,28,21,21,56,44,28) # excluded subjects are males aged 66, 33 and 24
mean_age = mean(ages)
sd_age = sd(ages)

# set up condition levels
imtypes = c('intact', 'local', 'global')
data$imtype <- factor(data$imtype, levels = imtypes)
imtypes_names = c('intact', 'locally\nscrambled', 'globally\nscrambled')
imtypes_legend = c("intact", "locally scrambled", "globally scrambled")
categories = c('face', 'house')
sequences = c('same', 'diff')
pref_cats = c('preferred', 'nonpreferred')

sig_code <- function(x){
  y <- c()
  y[x > .05] = ''
  y[x <= .05] = '*'
  y[x <= .01] = '**'
  y[x <= .001] = '***'
  return(y)
}

# FFA/PPA
regions = c('ffa', 'ppa') 
resdir <- file.path(projDir,'results/FFA_PPA_100')
subjects <- levels(data$subject)
facet_labels = c(ffa = 'FFA', ppa = 'PPA')
ylims = c(-.2,2.3)
temp_data <- dplyr::filter(data, grepl('_100', region)) # reduce to ROIs of 100 voxels
temp_data <- dplyr::filter(temp_data, !grepl('bilat|ofa|opa|psts', region)) # remove bilateral and other regions
DF = 16

  
# initial 5-way ANOVA -----------------------------------------------------
temp_data$hemi[grep("right", temp_data$region)] <- 'right' # set up hemisphere variable
temp_data$hemi[grep("left", temp_data$region)] <- 'left'
temp_data$hemi <- factor(temp_data$hemi)
temp_data$region <- factor(substr(temp_data$region, 0, 3)) # remove hemi and voxel info from region variable
temp_data <- droplevels(temp_data[temp_data$subject %in% subjects,]) # restrict to relevant subjects
DF <- length(levels(temp_data$subject))-1

# recode category based on the regions preference for that category
temp_data$pref_cat='nonpreferred'
for (cond in 1:length(temp_data$subject)){
  if (temp_data$region[cond] %in% c('ffa','ofa') && temp_data$category[cond] == 'face'){ temp_data$pref_cat[cond] = 'preferred'}
  if (temp_data$region[cond] %in% c('ppa','opa') && temp_data$category[cond] == 'house'){ temp_data$pref_cat[cond] = 'preferred'}
}
temp_data$pref_cat <- factor(temp_data$pref_cat, levels = pref_cats)

# add the labels you wish to use in the plots
temp_data$pref_cat_names <- revalue(temp_data$pref_cat, c('preferred' = 'P','nonpreferred'='NP'))
temp_data$imtype_names <- revalue(temp_data$imtype, c('global' = 'globally\nscrambled','local' = 'locally\nscrambled'))

# run 5-way anova
model_all <- ezANOVA(data = temp_data, dv = signal_change, wid = subject, within = .(region, hemi, imtype, pref_cat, sequence), detailed = TRUE)
model_all$ANOVA$pes = model_all$ANOVA$SSn/(model_all$ANOVA$SSn + model_all$ANOVA$SSd) # add effect sizes

# print results to file
dir.create(resdir, showWarnings = F)
anovaFile <- file.path(resdir, 'ANOVA.txt')
sink(anovaFile, append = F)
cat('# ANOVA #\n')
cat(capture.output(model_all), sep = '\n')
sink() # unsink text file

# examine interaction (image type * preferred category) -----------

# collapse across hemisphere, region and sequence 
temp_data_ip = aggregate(signal_change ~ subject + imtype + pref_cat + imtype_names + pref_cat_names, data = temp_data, FUN = mean)
temp_data_ip$condition <- factor(paste(temp_data_ip$imtype, temp_data_ip$pref_cat, sep='_'))

## pairwise
pwp_unc <- pairwise.t.test(temp_data_ip$signal_change, temp_data_ip$condition, paired = T, p.adjust.method = "none")
pwp_cor <- pairwise.t.test(temp_data_ip$signal_change, temp_data_ip$condition, paired = T, p.adjust.method = "fdr")
df <- data.frame(contrast = imtypes, meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
for (i in 1:length(imtypes)){
  pref <- temp_data_ip$signal_change[temp_data_ip$imtype == imtypes[i] & temp_data_ip$pref_cat == 'preferred']
  nonpref <- temp_data_ip$signal_change[temp_data_ip$imtype == imtypes[i] & temp_data_ip$pref_cat == 'nonpreferred']
  df$meandiff[i] <- mean(pref-nonpref)
  test <- t.test(pref,nonpref,paired=T)
  df$CI_lower[i] <- test$conf.int[1]
  df$CI_upper[i] <- test$conf.int[2]
  df$t[i] <- test$statistic[[1]]
  df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == sprintf('%s_preferred', imtypes[i]), colnames(pwp_unc$p.value) == sprintf('%s_nonpreferred', imtypes[i])]
  df$p.fdr[i] <- pwp_cor$p.value[rownames(pwp_unc$p.value) == sprintf('%s_preferred', imtypes[i]), colnames(pwp_unc$p.value) == sprintf('%s_nonpreferred', imtypes[i])]
  df$d[i] <- cohensD(pref,nonpref,method='paired')
}
df$p.sig = sig_code(df$p)
df$p.fdr.sig = sig_code(df$p.fdr)

pwFile <- file.path(resdir, 'pairwise.txt')
sink(pwFile, append = F)
cat(sprintf('# PW contrasts - preferred > non-preferred within image type #\n'))
cat(capture.output(df), sep = '\n')
sink() # unsink text file

# barplot
temp_plot = aggregate(signal_change ~ imtype + pref_cat, data = temp_data_ip, FUN = mean) # get means for each condition
temp2 = aggregate(signal_change ~ imtype + pref_cat, data = temp_data_ip, FUN = std.error) # get standard errors for each condition
temp_plot$se = temp2$signal_change # add standard errors to plot data
colnames(temp_plot)[3] = 'mean' # change the column name to mean
temp_plot$imtype = factor(temp_plot$imtype, levels = imtypes)
temp_plot$pref_cat = factor(temp_plot$pref_cat, levels = pref_cats)

pdf(file=file.path(resdir,'imtypeXprefcat.pdf'),
    bg='transparent',
    width=3.5, 
    height=3.5)
print(ggplot(temp_plot, aes(x=imtype, y=mean, fill=pref_cat)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.text.x  = element_text(colour = 'black', size=11), 
              axis.text.y  = element_text(colour = 'black', size=11), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size=11),
              legend.title = element_blank(),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        scale_fill_manual(values=c('white', 'gray')) +
        scale_x_discrete(labels = imtypes_names) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.5), size = .5) +
        geom_bar(aes(fill = pref_cat), stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black') +
        labs(y = 'signal change (%)') +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=seq(-.5,2,.5), expand = c(0,0)))
dev.off()

# smarties plot
mean_data = aggregate(signal_change ~ pref_cat_names + imtype_names, data = temp_data_ip, FUN = mean)
pdf(file=file.path(resdir,'imtypeXprefcat_smarties.pdf'),
    bg='transparent',
    width=4, 
    height=3.5)
print(ggplot(temp_data_ip, aes(x=pref_cat_names, y=signal_change)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.ticks.length=unit(-0.15, "cm"),
              axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
              axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
              strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = 14),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        facet_grid(.~imtype_names) +
        geom_violin(colour = NA, fill = 'black', alpha = .25) +
        labs(x = c('P','NP')) +
        geom_point(aes(x=pref_cat_names, y=signal_change, colour=subject)) +
        geom_line(aes(x = pref_cat_names, y=signal_change, colour = subject, group = subject)) +
        geom_point(data = mean_data, aes(x = pref_cat_names, y=signal_change, group = imtype_names), colour = 'black', shape = 15, size = 3) +
        geom_line(data = mean_data, aes(x = pref_cat_names, y=signal_change, group = imtype_names), colour = 'black') +
        ylab("signal change (%)") +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=seq(-.5,2,.5), expand = c(0,0)))
dev.off()
embed_fonts(file.path(resdir,'imtypeXprefcat_smarties.pdf'))

# examine interaction (region * image type * preferred category) for supplementary figure ----

temp_data_rip <- aggregate(signal_change ~ subject + region + imtype + pref_cat + imtype_names + pref_cat_names, data = temp_data, FUN = mean) # collapse x sequence
temp_data_rip$condition <- factor(paste(temp_data_rip$region, temp_data_rip$imtype, temp_data_rip$pref_cat, sep='_'))
temp_data_rip$contrast <- factor(paste(temp_data_rip$region, temp_data_rip$imtype, sep='_'))
cons <- levels(temp_data_rip$contrast)[c(2,3,1,5,6,4)]

## pairwise
pwp_unc <- pairwise.t.test(temp_data_rip$signal_change, temp_data_rip$condition, paired = T, p.adjust.method = "none")
pwp_cor <- pairwise.t.test(temp_data_rip$signal_change, temp_data_rip$condition, paired = T, p.adjust.method = "fdr")
df <- data.frame(contrast = cons, meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
for (i in 1:length(cons)){
  pref <- temp_data_rip$signal_change[temp_data_rip$contrast == cons[i] & temp_data_rip$pref_cat == 'preferred']
  nonpref <- temp_data_rip$signal_change[temp_data_rip$contrast == cons[i] & temp_data_rip$pref_cat == 'nonpreferred']
  df$meandiff[i] <- mean(pref-nonpref)
  test <- t.test(pref,nonpref,paired=T)
  df$CI_lower[i] <- test$conf.int[1]
  df$CI_upper[i] <- test$conf.int[2]
  df$t[i] <- test$statistic[[1]]
  df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == sprintf('%s_preferred', cons[i]), colnames(pwp_unc$p.value) == sprintf('%s_nonpreferred', cons[i])]
  df$p.fdr[i] <- pwp_cor$p.value[rownames(pwp_unc$p.value) == sprintf('%s_preferred', cons[i]), colnames(pwp_unc$p.value) == sprintf('%s_nonpreferred', cons[i])]
  df$d[i] <- cohensD(pref,nonpref,method='paired')
}
df$p.sig = sig_code(df$p)
df$p.fdr.sig = sig_code(df$p.fdr)

sink(pwFile, append = T)
cat(sprintf('\n# PW contrasts - preferred > non-preferred within region and image type #\n'))
cat(capture.output(df), sep = '\n')
sink() # unsink text file

# compare magnitude of preferred v non-preferred between image types

temp_data_ri.p <- aggregate(data=temp_data_rip, signal_change ~ subject + region + imtype, FUN = diff) # get preferred > nonpreferred
temp_data_ri.p <- aggregate(data=temp_data_ri.p, signal_change ~ subject + imtype, FUN = mean) # get mean across region
consX <- imtypes[c(1,1,2)]
consY <- imtypes[c(2,3,3)]

## pairwise
pwp_unc <- pairwise.t.test(temp_data_ri.p$signal_change, temp_data_ri.p$imtype, paired = T, p.adjust.method = "none")
pwp_unc_melt <- melt(pwp_unc$p.value)
pwp_cor <- pairwise.t.test(temp_data_ri.p$signal_change, temp_data_ri.p$imtype, paired = T, p.adjust.method = "fdr")
pwp_cor_melt <- melt(pwp_cor$p.value)
df <- data.frame(contrast = 1:length(consX), meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
for (i in 1:length(consX)){
  df$contrast[i] <- paste(consX[i],consY[i],sep=' v ')
  x <- temp_data_ri.p$signal_change[temp_data_ri.p$imtype == consX[i]]
  y <- temp_data_ri.p$signal_change[temp_data_ri.p$imtype == consY[i]]
  df$meandiff[i] <- mean(x-y)
  test <- t.test(x,y,paired=T)
  df$CI_lower[i] <- test$conf.int[1]
  df$CI_upper[i] <- test$conf.int[2]
  df$t[i] <- test$statistic[[1]]
  df$p[i] <- max(pwp_unc_melt$value[pwp_unc_melt$Var1 == consX[i] & pwp_unc_melt$Var2 == consY[i]], pwp_unc_melt$value[pwp_unc_melt$Var1 == consY[i] & pwp_unc_melt$Var2 == consX[i]], na.rm=T)
  df$p.fdr[i] <- max(pwp_cor_melt$value[pwp_cor_melt$Var1 == consX[i] & pwp_cor_melt$Var2 == consY[i]], pwp_cor_melt$value[pwp_unc_melt$Var1 == consY[i] & pwp_cor_melt$Var2 == consX[i]], na.rm=T)
  df$d[i] <- cohensD(x,y,method='paired')
}
df$p.sig = sig_code(df$p)
df$p.fdr.sig = sig_code(df$p.fdr)

sink(pwFile, append = T)
cat(sprintf('\n# PW contrasts - preferred > non-preferred between image type #\n'))
cat(capture.output(df), sep = '\n')
sink() # unsink text file

# plot
data_rip_plot = aggregate(signal_change ~ region + imtype + pref_cat, data = temp_data_rip, FUN = mean) # get means for each condition
temp = aggregate(signal_change ~ region + imtype + pref_cat, data = temp_data_rip, FUN = std.error) # get standard errors for each condition
data_rip_plot$se = temp$signal_change # add standard errors to plot data
colnames(data_rip_plot)[4] = 'mean' # change the column name to mean
data_rip_plot$imagetype = factor(data_rip_plot$imtype, levels = imtypes) # reorder factors
data_rip_plot$pref_cat = factor(data_rip_plot$pref_cat, levels = pref_cats)
pdf(file=file.path(resdir, 'regionXimtypeXprefcat.pdf'), 
    bg='transparent',
    width=6, 
    height=3.8)
print(ggplot(data_rip_plot, aes(x=imagetype, y=mean, fill=pref_cat)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.text.x  = element_text(colour = 'black', size=11), 
              axis.text.y  = element_text(colour = 'black', size=11), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size=11),
              legend.title = element_blank(),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2),
              strip.text.x = element_text(colour = 'black', size=13),
              strip.background = element_blank()) +
        facet_grid(. ~ region, labeller=labeller(region = facet_labels)) +
        scale_fill_manual(values=c('white', 'gray')) +
        scale_x_discrete(labels = imtypes_names) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.5), size = .5) +
        geom_bar(aes(fill = pref_cat), stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black') +
        labs(y = 'signal change (%)') +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=c(0.0,.5,1,1.5,2), expand = c(0,0)))
dev.off()


# smarties plot
mean_data = aggregate(signal_change ~ pref_cat_names + imtype_names + region, data = temp_data_rip, FUN = mean)
pdf(file=file.path(resdir,'regionXimtypeXprefcat_smarties.pdf'),
    bg='transparent',
    width=4, 
    height=3.5)
print(ggplot(temp_data_rip, aes(x=pref_cat_names, y=signal_change)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.ticks.length=unit(-0.15, "cm"),
              axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
              axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
              strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = 14),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        facet_grid(region ~ imtype_names, labeller=labeller(region = facet_labels)) +
        geom_hline(yintercept = 0) +
        geom_violin(colour = NA, fill = 'black', alpha = .25) +
        geom_point(aes(x=pref_cat_names, y=signal_change, colour=subject)) +
        geom_line(aes(x = pref_cat_names, y=signal_change, colour = subject, group = subject)) +
        geom_point(data = mean_data, aes(x = pref_cat_names, y=signal_change), colour = 'black', shape = 15, size = 3) +
        geom_path(data = mean_data, aes(x = rep(1:2, 6), y = signal_change)) +
        ylab("signal change (%)") +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=seq(0,2,1), expand = c(0,.2)))
dev.off()

# examine interaction (image type * preferred category * sequence) --------

temp_data_ips <- aggregate(signal_change ~ subject + imtype + pref_cat + sequence + pref_cat_names + imtype_names, data = temp_data, FUN = mean) # collapse x region
temp_data_ips$condition <- factor(paste(temp_data_ips$imtype, temp_data_ips$pref_cat, temp_data_ips$sequence, sep = '_'))
temp_data_ips$contrast <- factor(paste(temp_data_ips$imtype, temp_data_ips$pref_cat, sep='_'))
cons <- levels(temp_data_ips$contrast)[c(4,3,6,5,2,1)]

## pairwise
pwp_unc <- pairwise.t.test(temp_data_ips$signal_change, temp_data_ips$condition, paired = T, p.adjust.method = "none")
pwp_cor <- pairwise.t.test(temp_data_ips$signal_change, temp_data_ips$condition, paired = T, p.adjust.method = "fdr")
df <- data.frame(contrast = cons, meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
for (i in 1:length(cons)){
  different <- temp_data_ips$signal_change[temp_data_ips$contrast == cons[i] & temp_data_ips$sequence == 'diff']
  same <- temp_data_ips$signal_change[temp_data_ips$contrast == cons[i] & temp_data_ips$sequence == 'same']
  df$meandiff[i] <- mean(different-same)
  test <- t.test(different,same,paired=T)
  df$CI_lower[i] <- test$conf.int[1]
  df$CI_upper[i] <- test$conf.int[2]
  df$t[i] <- test$statistic[[1]]
  df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == sprintf('%s_same', cons[i]), colnames(pwp_unc$p.value) == sprintf('%s_diff', cons[i])]
  df$p.fdr[i] <- pwp_cor$p.value[rownames(pwp_unc$p.value) == sprintf('%s_same', cons[i]), colnames(pwp_unc$p.value) == sprintf('%s_diff', cons[i])]
  df$d[i] <- cohensD(different,same,method='paired')
}
df$p.sig = sig_code(df$p)
df$p.fdr.sig = sig_code(df$p.fdr)

sink(pwFile, append = T)
cat(sprintf('\n# PW contrasts - different > same within imtype and within category #\n'))
cat(capture.output(df), sep = '\n')
sink() # unsink text file

# plot
temp_data_plot <- aggregate(signal_change ~ subject + imtype + pref_cat + sequence, data = temp_data_ips, FUN = mean) # collapse across region
temp_data_plot <- aggregate(signal_change ~ imtype + pref_cat + sequence, data = temp_data_plot, FUN = mean) # collapse across subject
temp = aggregate(signal_change ~ imtype + pref_cat + sequence, data = temp_data, FUN = std.error) # get standard errors across subjects
temp_data_plot$se = temp$signal_change # add standard errors to plot data
colnames(temp_data_plot)[4] = 'mean' # change the column name to mean
temp_data_plot$imtype = factor(temp_data_plot$imtype, levels = imtypes)
temp_data_plot$pref_cat = factor(temp_data_plot$pref_cat, levels = pref_cats)

pdf(file=file.path(resdir,'imtypeXprefcatXsequence.pdf'), 
    bg='transparent',
    width=6, 
    height=3.8)
print(ggplot(temp_data_plot, aes(x=imtype, y=mean, fill=sequence)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.text.x  = element_text(colour = 'black', size=11), 
              axis.text.y  = element_text(colour = 'black', size=11), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size=11),
              legend.title = element_blank(),
              legend.justification=c(1,-2), 
              legend.position=c(1,.2),
              strip.text.x = element_text(colour = 'black', size=13),
              strip.background = element_blank()) +
        facet_grid(.~pref_cat)+
        scale_fill_manual(values=c('white', 'gray'), labels = c('different','same')) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.5), size = .5) +
        geom_bar(aes(fill = sequence), stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black') +
        labs(y = 'signal change (%)') +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=c(0.0,.5,1,1.5,2), expand = c(0,0)))
dev.off()

# smarties plot
temp_data_ips$imtype_names <- factor(temp_data_ips$imtype_names, levels = imtypes_names)
mean_data = aggregate(signal_change ~ pref_cat + imtype_names + sequence, data = temp_data_ips, FUN = mean)
pdf(file=file.path(resdir,'imtypeXprefcatXsequence_smarties.pdf'),
    bg='transparent',
    width=4, 
    height=3.5)
print(ggplot(temp_data_ips, aes(x=sequence, y=signal_change)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.ticks.length=unit(-0.15, "cm"),
              axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
              axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
              strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = 14),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        facet_grid(pref_cat ~ imtype_names) +
        geom_hline(yintercept = 0) +
        geom_violin(colour = NA, fill = 'black', alpha = .25) +
        geom_point(aes(x=sequence, y=signal_change, colour=subject)) +
        geom_line(aes(x = sequence, y=signal_change, colour = subject, group = subject)) +
        geom_point(data = mean_data, aes(x = sequence, y=signal_change, group = pref_cat), colour = 'black', shape = 15, size = 3) +
        geom_path(data = mean_data, aes(x = rep(1:2, 6), y = signal_change)) +
        ylab("signal change (%)") +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=seq(0,2,1), expand = c(0,.2)))
dev.off()

# compare adaptation within preferred category and between image type
temp_data_ip.s = aggregate(signal_change ~ subject + imtype + pref_cat + contrast + imtype_names + pref_cat_names, data = temp_data_ips, FUN = diff) # get adaptation
temp_data_ip.s$signal_change <- -temp_data_ip.s$signal_change
colnames(temp_data_ip.s)[4] <- 'condition'
consX <- c('intact_preferred', 'intact_preferred', 'local_preferred', 
           'intact_nonpreferred', 'intact_nonpreferred', 'local_nonpreferred')
consY <- c('local_preferred', 'global_preferred', 'global_preferred', 
           'local_nonpreferred', 'global_nonpreferred', 'global_nonpreferred')

## pairwise
pwp_unc <- pairwise.t.test(temp_data_ip.s$signal_change, temp_data_ip.s$condition, paired = T, p.adjust.method = "none")
pwp_unc_melt <- melt(pwp_unc$p.value)
pwp_cor <- pairwise.t.test(temp_data_ip.s$signal_change, temp_data_ip.s$condition, paired = T, p.adjust.method = "fdr")
pwp_cor_melt <- melt(pwp_cor$p.value)
df <- data.frame(contrast = 1:length(consX), meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
for (i in 1:length(consX)){
  df$contrast[i] <- paste(consX[i],consY[i],sep=' v ')
  x <- temp_data_ip.s$signal_change[temp_data_ip.s$condition == consX[i]]
  y <- temp_data_ip.s$signal_change[temp_data_ip.s$condition == consY[i]]
  df$meandiff[i] <- mean(x-y)
  test <- t.test(x,y,paired=T)
  df$CI_lower[i] <- test$conf.int[1]
  df$CI_upper[i] <- test$conf.int[2]
  df$t[i] <- test$statistic[[1]]
  df$p[i] <- max(pwp_unc_melt$value[pwp_unc_melt$Var1 == consX[i] & pwp_unc_melt$Var2 == consY[i]], pwp_unc_melt$value[pwp_unc_melt$Var1 == consY[i] & pwp_unc_melt$Var2 == consX[i]], na.rm=T)
  df$p.fdr[i] <- max(pwp_cor_melt$value[pwp_cor_melt$Var1 == consX[i] & pwp_cor_melt$Var2 == consY[i]], pwp_cor_melt$value[pwp_unc_melt$Var1 == consY[i] & pwp_cor_melt$Var2 == consX[i]], na.rm=T)
  df$d[i] <- cohensD(x,y,method='paired')
}
df$p.sig = sig_code(df$p)
df$p.fdr.sig = sig_code(df$p.fdr)

sink(pwFile, append = T)
cat(sprintf('\n# PW contrasts - different > same within category between image type #\n'))
cat(capture.output(df), sep = '\n')
sink() # unsink text file

# plot
data_ips_adapt_plot = aggregate(signal_change ~ imtype + pref_cat, data = temp_data_ip.s, FUN = mean) # get means for each condition
temp = aggregate(signal_change ~ imtype + pref_cat, data = temp_data_ip.s, FUN = std.error) # get standard errors for each condition
data_ips_adapt_plot$se = temp$signal_change # add standard errors to plot data
colnames(data_ips_adapt_plot)[3] = 'mean' # change the column name to mean
data_ips_adapt_plot$mean = abs(data_ips_adapt_plot$mean) # sign flip
data_ips_adapt_plot$imtype = factor(data_ips_adapt_plot$imtype, levels = imtypes)
data_ips_adapt_plot$pref_cat = factor(data_ips_adapt_plot$pref_cat, levels = pref_cats)

pdf(file=file.path(resdir,'adaptation_imtypeXprefcat.pdf'), 
    bg='transparent',
    width=4, 
    height=3.8)
print(ggplot(data_ips_adapt_plot, aes(x=pref_cat, y=mean, fill=imtype)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.text.x  = element_text(colour = 'black', size=11), 
              axis.text.y  = element_text(colour = 'black', size=11), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size=11),
              legend.title = element_blank(),
              legend.justification=c(1,-2), 
              legend.position=c(1,.04),
              strip.text.x = element_text(colour = 'black', size=13),
              strip.background = element_blank()) +
        scale_x_discrete(labels = c("preferred\ncategory", "non-preferred\ncategory")) +
        scale_fill_manual(values=c('white', 'gray66', 'gray33'), labels = c(imtypes_legend)) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.5), size = .5) +
        geom_bar(aes(fill = imtype), stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black') +
        labs(y = 'signal change (%)') +
        coord_cartesian(ylim=c(0, 0.35)) +
        scale_y_continuous(breaks=c(0.0,.1,.2,.3), expand = c(0,0)))
dev.off()

# smarties plot
temp_data_ip.s$imtype_names <- factor(temp_data_ip.s$imtype_names, levels = imtypes_names)
mean_data <- aggregate(data = temp_data_ip.s, signal_change ~ pref_cat + imtype_names, FUN = mean)
pdf(file=file.path(resdir,'adaptation_imtypeXprefcat_smarties.pdf'),
    bg='transparent',
    width=6, 
    height=3.8)
print(ggplot(temp_data_ip.s, aes(x=imtype_names, y=signal_change)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.ticks.length=unit(-0.15, "cm"),
              axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
              axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
              strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = 14),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        facet_grid(. ~ pref_cat) +
        geom_hline(yintercept = 0) +
        geom_violin(colour = NA, fill = 'black', alpha = .25) +
        geom_point(aes(x=imtype_names, y=signal_change, colour=subject)) +
        geom_line(aes(x=imtype_names, y=signal_change, colour = subject, group = subject)) +
        geom_point(data = mean_data, aes(x = imtype_names, y=signal_change), colour = 'black', shape = 15, size = 3) +
        geom_path(data = mean_data, aes(x = rep(1:3, 2), y = signal_change)) +
        ylab("signal change (%)") +
        coord_cartesian(ylim=c(-.3,.65)) +
        scale_y_continuous(breaks=seq(-.2,.6,.2)))
dev.off()

# plot region x imtype x category x sequence, although no 4 way interaction found
data_rips_plot = aggregate(signal_change ~ region + imtype + pref_cat + sequence + pref_cat_names + imtype_names, data = temp_data, FUN = mean) # get means for each condition
temp = aggregate(signal_change ~ region + imtype + pref_cat + sequence, data = temp_data, FUN = std.error) # get standard errors for each condition
data_rips_plot$se = temp$signal_change # add standard errors to plot data
colnames(data_rips_plot)[7] = 'mean' # change the column name to mean
data_rips_plot$imtype = factor(data_rips_plot$imtype, levels = imtypes) # reorder factors
data_rips_plot$pref_cat = factor(data_rips_plot$pref_cat, levels = pref_cats)
pdf(file=file.path(resdir,'imtypeXprefcatXregionXsequence.pdf'), 
    bg='transparent',
    width=6, 
    height=6)
print(ggplot(data_rips_plot, aes(x=imtype, y=mean, fill=sequence)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.text.x  = element_text(colour = 'black', size=11), 
              axis.text.y  = element_text(colour = 'black', size=11), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
              legend.text = element_text(colour = 'black', size=11),
              legend.title = element_blank(),
              legend.justification=c(1,-4.9), 
              legend.position=c(1,0.2),
              strip.text.x = element_text(colour = 'black', size=13),
              strip.text.y = element_text(colour = 'black', size=13),
              strip.background = element_blank()) +
        facet_grid(region ~ pref_cat, labeller=labeller(region = facet_labels)) +
        scale_fill_manual(values=c('white', 'gray'), labels = c("different", "same")) +
        scale_x_discrete(labels = imtypes_names) +
        geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.2, position=position_dodge(.5), size = .5) +
        geom_bar(aes(fill = sequence), stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black') +
        labs(y = 'signal change (%)') +
        coord_cartesian(ylim=ylims) +
        scale_y_continuous(breaks=c(0.0,.5,1,1.5,2), expand = c(0,0)))
dev.off()

# smarties plot
temp_data_smarties <- aggregate(signal_change ~ subject +region + imtype_names + pref_cat_names + sequence, data = temp_data, FUN = mean) # get means across hemi
temp_data_smarties$imtype_names <- factor(temp_data_smarties$imtype_names, levels = imtypes_names)
mean_data <- aggregate(data = temp_data_smarties, signal_change ~ region + imtype_names + pref_cat_names + sequence, FUN = mean)
pdf(file=file.path(resdir,'imtypeXprefcatXregionXsequence_smarties.pdf'),
    bg='transparent',
    width=7, 
    height=6)
print(ggplot(temp_data_smarties, aes(x=sequence, y=signal_change)) + 
        theme(panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(), 
              panel.background = element_blank(), 
              axis.line.x = element_line(colour = 'black', size = .5), 
              axis.line.y = element_line(colour = 'black', size = .5), 
              axis.ticks = element_line(colour = 'black', size = .5),
              axis.ticks.length=unit(-0.15, "cm"),
              axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
              axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
              axis.title.x = element_blank(), 
              axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
              strip.background = element_blank(),
              strip.text = element_text(colour = 'black', size = 14),
              legend.justification=c(1,-2.4), 
              legend.position=c(1,.2)) +
        facet_grid(region ~ pref_cat_names + imtype_names, labeller=labeller(region = facet_labels)) +
        geom_hline(yintercept = 0) +
        geom_violin(colour = NA, fill = 'black', alpha = .25) +
        geom_point(aes(x=sequence, y=signal_change, colour=subject)) +
        geom_line(aes(x=sequence, y=signal_change, colour = subject, group = subject)) +
        geom_point(data = mean_data, aes(x = sequence, y=signal_change, group = pref_cat_names), colour = 'black', shape = 15, size = 3) +
        geom_path(data = mean_data, aes(x = rep(1:2, 12), y = signal_change)) +
        ylab("signal change (%)") +
        coord_cartesian(ylim=c(-.5,2.2)) +
        scale_y_continuous(breaks=seq(0,2,1), expand = c(0,.2)))
dev.off()


# all ROIs - basic analysis --------------------------------------------------------
resdirREGS <- file.path(projDir,'results/eachROI')
dir.create(resdirREGS, showWarnings = F)

for (region in levels(data$region)){
  
  # set up directory and filepaths for results
  resdir <- file.path(projDir,'results/eachROI',region)
  dir.create(resdir, showWarnings = F)
  anovaFile <- file.path(resdir, 'ANOVA.txt')
  
  # get data for this region, run ANOVA and write to file
  temp_data_reg <- droplevels(data[data$region == region,])
  temp_data_reg$condition <- factor(paste(temp_data_reg$imtype, temp_data_reg$category, temp_data_reg$sequence, sep='_'))
  temp_data_reg$contrast <- factor(paste(temp_data_reg$imtype, temp_data_reg$category, sep = '_'))
  cons <- levels(temp_data_reg$contrast)[c(3,5,1,4,6,2)]
  DF <- length(levels(temp_data_reg$subject))-1
  
  temp_ANOVA <- ezANOVA(data = temp_data_reg, dv = signal_change, wid = subject, within = .(imtype, category, sequence), detailed = TRUE)
  temp_ANOVA$ANOVA$pes = temp_ANOVA$ANOVA$SSn/(temp_ANOVA$ANOVA$SSn + temp_ANOVA$ANOVA$SSd) # add effect sizes
  sink(anovaFile, append = F)
  cat('# ANOVA #\n')
  cat(capture.output(temp_ANOVA), sep = '\n')
  sink() # unsink text file
  
  
  ## pairwise
  pwp_unc <- pairwise.t.test(temp_data_reg$signal_change, temp_data_reg$condition, paired = T, p.adjust.method = "none")
  pwp_cor <- pairwise.t.test(temp_data_reg$signal_change, temp_data_reg$condition, paired = T, p.adjust.method = "fdr")
  df <- data.frame(contrast = cons, meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
  for (i in 1:length(cons)){
    different <- temp_data_reg$signal_change[temp_data_reg$contrast == cons[i] & temp_data_reg$sequence == 'diff']
    same <- temp_data_reg$signal_change[temp_data_reg$contrast == cons[i] & temp_data_reg$sequence == 'same']
    df$meandiff[i] <- mean(different-same)
    test <- t.test(different,same,paired=T)
    df$CI_lower[i] <- test$conf.int[1]
    df$CI_upper[i] <- test$conf.int[2]
    df$t[i] <- test$statistic[[1]]
    df$p[i] <- pwp_unc$p.value[colnames(pwp_unc$p.value) == sprintf('%s_diff', cons[i]), rownames(pwp_unc$p.value) == sprintf('%s_same', cons[i])]
    df$p.fdr[i] <- pwp_cor$p.value[colnames(pwp_unc$p.value) == sprintf('%s_diff', cons[i]), rownames(pwp_unc$p.value) == sprintf('%s_same', cons[i])]
    df$d[i] <- cohensD(different,same,method='paired')
  }
  df$p.sig = sig_code(df$p)
  df$p.fdr.sig = sig_code(df$p.fdr)
  
  sink(anovaFile, append = T)
  cat(sprintf('# PAIRWISE - different > same within image type and within category #\n'))
  cat(capture.output(df), sep = '\n')
  sink() # unsink text file
  
  # plot
  temp_data_plot = aggregate(signal_change ~ sequence, data = temp_data_reg, FUN = mean) # get means for each condition
  temp = aggregate(signal_change ~ subject + sequence, data = temp_data_reg, FUN = mean) # collapse across category and imtype, preserving subjects
  temp = aggregate(signal_change ~ sequence, data = temp, FUN = std.error) # get standard error across subjects
  temp_data_plot$se = temp$signal_change # add standard errors to plot data
  seq_labels = c(diff = 'different', same = 'same')
  signif_char = sig_code(temp_ANOVA$ANOVA$`p`[4])
  signif_height = sum(temp_data_plot[1,c(2,3)])+.05
  pdf(file=file.path(resdir,sprintf('barplot.pdf')), 
      bg='transparent',
      width=2, 
      height=3)
  print(ggplot(temp_data_plot, aes(x=sequence, y=signal_change)) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(), 
                axis.line.x = element_line(colour = 'black', size = .5), 
                axis.line.y = element_line(colour = 'black', size = .5), 
                axis.ticks = element_line(colour = 'black', size = .5),
                axis.text.x  = element_text(colour = 'black', size=11), 
                axis.text.y  = element_text(colour = 'black', size=11), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(size=13, margin=margin(0,10,0,0)),
                legend.text = element_text(colour = 'black', size=11),
                legend.title = element_blank()) +
          geom_text(label = signif_char, x = 1.5, y = signif_height, size = 5) +
          geom_path(x = c(1,2), y = c(signif_height-.02, signif_height-.02)) +
          geom_errorbar(aes(ymin=signal_change-se, ymax=signal_change+se), width=.2, position=position_dodge(.5), size = .5) +
          geom_bar(stat = 'identity', width = 0.5, position = position_dodge(width=0.5), colour = 'black', fill = 'gray') +
          labs(y = 'signal change (%)') +
          coord_cartesian(ylim=c(0, 1.5)) +
          scale_x_discrete(labels = seq_labels) +
          scale_y_continuous(breaks=c(0.0,.5,1,1.5,2), expand = c(0,0)))
  dev.off()
  
  # explore overall activity to global, local and intact images
  temp_data_it <- aggregate(signal_change ~ subject + imtype, data = temp_data_reg, FUN = mean)
  
  ## pw contrasts
  pwp_unc <- pairwise.t.test(temp_data_it$signal_change, temp_data_it$imtype, paired = T, p.adjust.method = "none")
  pwp_cor <- pairwise.t.test(temp_data_it$signal_change, temp_data_it$imtype, paired = T, p.adjust.method = "fdr")
  df <- data.frame(contrast = c('intact_local','intact_global','local_global'), meandiff=NA, CI_lower = NA, CI_upper = NA, df = DF, t=NA, p=NA, p.sig='', p.fdr=NA, p.fdr.sig='', d=NA)
  consA <- c('intact','intact', 'local')
  consB <- c('local','global','global')
  for (i in 1:length(consA)){
    A <- temp_data_it$signal_change[temp_data_it$imtype == consA[i]]
    B <- temp_data_it$signal_change[temp_data_it$imtype == consB[i]]
    df$meandiff[i] <- mean(A-B)
    test <- t.test(A,B,paired=T)
    df$CI_lower[i] <- test$conf.int[1]
    df$CI_upper[i] <- test$conf.int[2]
    df$t[i] <- test$statistic[[1]]
    df$p[i] <- pwp_unc$p.value[rownames(pwp_unc$p.value) == consB[i], colnames(pwp_unc$p.value) == consA[i]]
    df$p.fdr[i] <- pwp_cor$p.value[rownames(pwp_cor$p.value) == consB[i], colnames(pwp_cor$p.value) == consA[i]]
    df$d[i] <- cohensD(A,B,method='paired')
  }
  df$p.sig = sig_code(df$p)
  df$p.fdr.sig = sig_code(df$p.fdr)
  
  sink(anovaFile, append = T)
  cat(sprintf('\n# PW - signal change between image type #\n'))
  cat(capture.output(df), sep = '\n')
  sink() # unsink text file
  
  # make face v house stats
  temp_data_c <- aggregate(signal_change ~ subject + category, data = temp_data_reg, FUN = mean)
  test <- t.test(temp_data_c$signal_change[temp_data_c$category == 'face'], temp_data_c$signal_change[temp_data_c$category == 'house'], paired = T)
  effectSize <- cohensD(temp_data_c$signal_change[temp_data_c$category == 'face'], temp_data_c$signal_change[temp_data_c$category == 'house'], method = 'paired') 
  
  sink(anovaFile, append = T)
  cat(sprintf('\n# PW - face v house #\n'))
  cat(capture.output(test), sep = '\n')
  cat(sprintf('\neffect size\n'))
  cat(capture.output(effectSize), sep = '\n')
  sink() # unsink text file
  
  # make smarties plot for V1 figure 5
  temp_data_reg$imtype_names <- revalue(temp_data_reg$imtype, c('global' = 'globally\nscrambled','local' = 'locally\nscrambled'))
  temp_data_reg$imtype_names <- factor(temp_data_reg$imtype_names, levels = imtypes_names)
  mean_data <- aggregate(data = temp_data_reg, signal_change ~ category + imtype_names + sequence, FUN = mean)
  pdf(file=file.path(resdir,'imtypeXcategoryXsequence_smarties.pdf'),
      bg='transparent',
      width=4, 
      height=3.5)
  print(ggplot(temp_data_reg, aes(x=sequence, y=signal_change)) + 
          theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(), 
                panel.background = element_blank(), 
                axis.line.x = element_line(colour = 'black', size = .5), 
                axis.line.y = element_line(colour = 'black', size = .5), 
                axis.ticks = element_line(colour = 'black', size = .5),
                axis.ticks.length=unit(-0.15, "cm"),
                axis.text.x  = element_text(colour = 'black', size=12, margin = margin(8,0,0,0)), 
                axis.text.y  = element_text(colour = 'black', size=12, margin = margin(0,8,0,0)), 
                axis.title.x = element_blank(), 
                axis.title.y = element_text(size=14, margin=margin(0,10,0,0)),
                strip.background = element_blank(),
                strip.text = element_text(colour = 'black', size = 14),
                legend.justification=c(1,-2.4), 
                legend.position=c(1,.2)) +
          facet_grid(category ~ imtype_names) +
          geom_hline(yintercept = 0) +
          geom_violin(colour = NA, fill = 'black', alpha = .25) +
          geom_point(aes(x=sequence, y=signal_change, colour=subject)) +
          geom_line(aes(x = sequence, y=signal_change, colour = subject, group = subject)) +
          geom_point(data = mean_data, aes(x = sequence, y=signal_change, group = category), colour = 'black', shape = 15, size = 3) +
          geom_path(data = mean_data, aes(x = rep(1:2, 6), y = signal_change)) +
          ylab("signal change (%)") +
          coord_cartesian(ylim=c(-.5,3.1)) +
          scale_y_continuous(breaks=seq(0,3,1), expand = c(0,.2)))
  dev.off()
  
}

