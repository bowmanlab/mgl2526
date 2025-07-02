## select the cast you want to use

select.cast <- '20250619_1'
data.directory <- 'Z://public//CTD//'  # may differ for OSX users

## read in some key files

ctd.metadata <- read.csv(paste0(data.directory, 'mgl2506_castmetadata.csv'))
cast.data <- read.csv(paste0(data.directory, select.cast, '_clean.csv'))
bottle.data <- read.csv(paste0(data.directory, 'mgl2506_bottles.csv'))

## some commands to help you understand the data

dim(cast.metadata)
colnames(cast.data)
colnames(bottle.data)

## make a simple profile plot

plot(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg.,
     cast.data$depSM..Depth..salt.water..m.,
     ylim = c(600, 0), # this reverses the y-axis scale
     lty = 1)

## does chlorophyll correlate with oxygen?

plot(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg. ~ cast.data$flECO.AFL..Fluorescence..WET.Labs.ECO.AFL.FL..mg.m.3.)

## neat... two different relationships!  Let's focus on one.

select <- which(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg < 220)

plot(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg.[select] ~ cast.data$flECO.AFL..Fluorescence..WET.Labs.ECO.AFL.FL..mg.m.3.[select])

lm.example <- lm(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg.[select] ~ cast.data$flECO.AFL..Fluorescence..WET.Labs.ECO.AFL.FL..mg.m.3.[select])

abline(lm.example)

summary(lm.example)

## check out new parameter z_1026

plot(ctd.metadata$z_1026_.m.,
     ylab = '1026 isopycnal (m)',
     xlab = 'Cast')

hist(ctd.metadata$z_1026_.m.)

## let's see how much the two oxygen sensors disagree

plot(cast.data$sbox0Mm.Kg..Oxygen..SBE.43..umol.kg.-cast.data$sbox1Mm.Kg..Oxygen..SBE.43..2..umol.kg.,
     type = 'h',
     ylab = 'deltaO2')

#### cast metadata exploration ####

plot(ctd.metadata$cmax_.m. ~ ctd.metadata$comp_depth_.m.,
     ylab = 'Chlorophyll max (m)',
     xlab = '1 % light level (m)') # nice correlation, one outlier

abline(0, 1, lty = 2)

abline(lm(ctd.metadata$cmax_.m. ~ ctd.metadata$comp_depth_.m.))

ctd.metadata[order(ctd.metadata$cmax_.m., decreasing = T),] # identify the outlier, might be interesting

plot(ctd.metadata$maxc_.mg.m.3. ~ ctd.metadata$comp_depth_.m.) # nice correlation!

plot(ctd.metadata$maxc_.mg.m.3. ~ ctd.metadata$ml_depth_.m.) # messy, could be multiple relationships

#### all casts ####

## let's get a little fancier and load all of the casts so we can call them as needed

all.casts <- list()

csv.list <- list.files(path = data.directory,
                     pattern = '*clean.csv',
                     ignore.case = F)

for(f in csv.list){
  temp <- read.csv(paste0(data.directory, select.cast, '_clean.csv'))
  all.casts[f] <- temp
}

## now all the files are stored as objects in the list all.casts and can be called by name.  You can call one like:

working.cast <- all.casts$`20250619_1_clean.csv`
