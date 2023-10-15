
library(survival)
library(msm)
library(coxinterval) 

# Read in the data from http://faculty.washington.edu/jphughes/pubs.html, save as "dat" 

# Add IDs
dat$id <- 1:nrow(dat)

# dtime = 25 means they lived past 24 months, change dstatus to 0
dat$dstatus[dat$dtime==25] <- 0

# Less than 0.19 months from negative test to death was counted as HIV negative immediately 
# prior to death as in Frydman and Szarek 2010
dat$hl[dat$dstatus==1 & dat$hr == 25 & dat$dtime-dat$hl < 0.19]<-
  dat$dtime[dat$dstatus==1 & dat$hr == 25 & dat$dtime-dat$hl < 0.19]
## One person died at time zero, make it time 0.01
dat$dtime[dat$dtime==0] <- 0.01

### Count number of people in each group to match Table II in Frydman and Szarek 2010
table(dat$grp[dat$hr < 25 & dat$dstatus==1]) # tested positive, observed death
table(dat$grp[dat$hr < 25 & dat$dstatus==0]) # tested positive, right cens for death
table(dat$grp[dat$hr == 25 & dat$hl == dat$dtime & dat$dstatus==1]) # tested neg immed prior to death
table(dat$grp[dat$hl == dat$dtime & dat$dstatus==0]) # Right cens for death at last neg HIV test
table(dat$grp[dat$hr == 25 & dat$hl<24 & dat$dstatus==1 & dat$hl< dat$dtime]) # HIV status unknown, died
table(dat$grp[dat$hr == 25 & dat$hl<24 & dat$dstatus==0 & dat$hl< dat$dtime]) # HIV status unknown, censored 

### KM curve and Cox model for HIV free survival using "naive" approach
dat$t2composite <- ifelse(dat$hr<25, dat$hr, dat$dtime)
dat$composite <- ifelse(dat$hr<25, 1, dat$dstatus)

survfit(Surv(t2composite, composite)~grp, data = dat)
coxph(Surv(t2composite, composite)~grp, data = dat)

### Get the data in the right format for the msm package, which implements the parametric approach
dat_msm <- NULL
for (i in 1:nrow(dat)){
  # Start in state 1
  if (dat$dtime[i] > 0)  dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], state = 1, time = 0, obstype = 1))
  
  # Death with negative test immed prior (1->3 transition, so obstype = 2)
  if (dat$hr[i] == 25 & dat$hl[i] == dat$dtime[i] & dat$dstatus[i]==1) dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], state = 3, time = dat$dtime[i], obstype = 2))
  
  # Last negative test for those who didn't have a neg test immed prior to death but had at least one negative test
  if (!(dat$hr[i] == 25 & dat$hl[i] == dat$dtime[i] & dat$dstatus[i]==1) & dat$hl[i]>0 & dat$hl[i]<25) dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], state = 1, time = dat$hl[i], obstype = 1))
  
  # First positive test, if any
  if (dat$hr[i]<25) dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], state = 2, time = dat$hr[i], obstype = 1))

  # Death or cens, no test on the day of death/cens
  if ((dat$dtime[i]!=dat$hl[i] & dat$dtime[i]!=dat$hr[i]) | (dat$dtime[i] == 25 & dat$hr[i]==25 & dat$hl[i] < 25)) dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], 
                                                                                                                                                        state = ifelse(dat$dstatus[i]==1, 3, 4), time = dat$dtime[i], obstype = ifelse(dat$dstatus[i]==1, 3, 1)))
  # Alive and HIV free at 24 months with no tests before then 
  if (dat$dtime[i] == 25 & dat$hl[i] == 25 & dat$hr[i] == 25)  dat_msm <- rbind(dat_msm, data.frame(id = i, grp = dat$grp[i], state = 1, time = dat$dtime[i], obstype = 1))
}

### Check counts of transitions 
statetable.msm(state, id, data=dat_msm[dat_msm$grp == 0,])
# Matches Frydman+Szarek except they had 19 who tested + and were right cens but we have 18 with 2->4 transition (HIV+, cens for death) 
# Difference is likely bc one FF person had a positive test and right censoring on the same day. 
## 1->2 corresponds to tested HIV + (19+12 in Frydman)
## 1->3 corresponds to tested - immed prior to death or died with unknown status (3+24 in Frydman)
## 1->4 corresponds to right cens w unknown status 
## 2->4 corresponds to right cens HIV+ (19 in frydman)

## Parametric approach  with piecewise constant hazards
## Estimating HIV free survival curves separately in the two groups 
piece0 <- msm(state ~ time, subject = id, data = dat_msm[dat_msm$grp == 0,], Q, censor = 4, gen.inits = T,
              pci = c(1.5,3,5, 10, 15), center = F)
piece1 <- msm(state ~ time, subject = id, data = dat_msm[dat_msm$grp == 1,], Q, censor = 4, gen.inits = T,
              pci = c(1.5,3,5, 10, 15), center = F)

# Plot naive KM curve and parametric curves
plot(survfit(Surv(t2composite, composite)~grp, data = dat), col = 'red', lty = 1:2)
lines(seq(from = 0, to =24, length.out = 100), 
      sapply(seq(from = 0, to =24, length.out = 100), function(x) pmatrix.msm(piece0, t=x)[1,1]))
lines(seq(from = 0, to =24, length.out = 100), 
      sapply(seq(from = 0, to =24, length.out = 100), function(x) pmatrix.msm(piece1, t=x)[1,1]),lty = 2)

# Get 2 year HIV-free survival estimates and CIs 
set.seed(897)
pmatrix.msm(piece0, t=24,ci = 'bootstrap', cores = 7, B=250)
set.seed(534)
pmatrix.msm(piece1, t=24,ci = 'bootstrap', cores = 5, B=250)

###### Parametric approach with piecewise constant hazards
## Estimating HR for HIV free survival
msm(state ~ time, subject = id, data = dat_msm, Q, censor = 4, gen.inits = T,
             covariates = ~grp, center = F, pci = c(1.5,3,5, 10, 15))
## Constrain the effect to be the same on both events 
msm(state ~ time, subject = id, data = dat_msm, Q, censor = 4, gen.inits = T,
    covariates = ~grp, center = F, pci = c(1.5,3,5, 10, 15), 
    constraint = list(grp = c(1,1,2)))

##### Put data into the format expected by the coxdual function 
sieve <- NULL
for (i in (1:nrow(dat))){
  # Tested positive
  if (dat$hr[i] < 25) sieve <- rbind(sieve, data.frame(id = dat$id[i], grp =dat$grp[i],
                                                       from = c(0,0,1), to = c(1,2,2), 
                                                       start = c(0,0,dat$hr[i]), stop = c(NA,NA,dat$dtime[i]),
                                                       status = c(0, 0, dat$dstatus[i])))
  
  # Tested negative at the same time as death/censoring
  if (dat$hl[i] == dat$dtime[i]) sieve <- rbind(sieve, data.frame(id = dat$id[i], grp =dat$grp[i],
                                                                  from = c(0,0), to = c(1,2), 
                                                                  start = c(0,0), stop = dat$dtime[i],
                                                                  status = c(0, dat$dstatus[i])))
  # HIV status at time of death/cens unknown
  if (dat$hl[i] <=24 & dat$hl[i]< dat$dtime[i] & dat$hr[i]==25) sieve <- 
      rbind(sieve, data.frame(id = dat$id[i], grp =dat$grp[i],
                              from = c(0,0,NA), to = c(1,2,2), 
                              start = c(0,0,NA), stop = c(rep(max(dat$hl[i],0),2), dat$dtime[i]),
                              status = c(0, 0, dat$dstatus[i])))
}

##### Fit the model with separate effects for each event
fit <- coxdual(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
               + I(grp * (to == 1)) + I(grp * (from %in% 0 & to == 2))
               + I(grp * (from %in% c(NA, 1) & to == 2)), data = sieve)

##### Constrain the two coefficients to be equal 
fit.constr <- coxdual(Surv(start, stop, status) ~ cluster(id) + trans(from, to)
                      + I(grp * (to == 1 |(from %in% 0 & to == 2)))
                      + I(grp * (from %in% c(NA, 1) & to == 2)), data = sieve)
