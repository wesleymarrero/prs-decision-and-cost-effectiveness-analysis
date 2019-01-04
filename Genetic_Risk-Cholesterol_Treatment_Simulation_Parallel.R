# ******************************************************************
# Understanding the Role of Genetic Risk in ASCVD Treatment Planing
# ******************************************************************

# Setup -------------------------------------------------------------------

#Clearing environment
rm(list = ls()[!(ls() %in% c())])

#Selecting home directory
if(Sys.info()["nodename"]=="IOE-TESLA"){ #selecting path by computer name
  home_dir = "~/Genetics CEA"
  setwd(home_dir)
}else{
  home_dir = file.path(Sys.getenv("USERPROFILE"),"/Box Sync/University of Michigan/Research/Current Work/Genetic Risk/R")
  setwd(home_dir)
}

#Loading packages (and installing if necessary)
if(!("data.table" %in% installed.packages()[,"Package"])) install.packages("data.table"); library(data.table) #large data frame manipulation
if(!("doParallel" %in% installed.packages()[,"Package"])) install.packages("doParallel"); library(doParallel) #Parallel computation
if(!("reshape2" %in% installed.packages()[,"Package"])) install.packages("reshape2"); library(reshape2) #Data frames for plots
if(!("ggplot2" %in% installed.packages()[,"Package"])) install.packages("ggplot2"); library(ggplot2) #Plots
if(!("scales" %in% installed.packages()[,"Package"])) install.packages("scales"); library(scales) #% and $ in plots
if(!("ggridges" %in% installed.packages()[,"Package"])) install.packages("ggridges"); library(ggridges) #ridge plots plots

#Selecting number of cores for parallel processing
if(Sys.info()["nodename"]=="IOE-TESLA"){
  cores = 35
}else{
  cores = detectCores()-1 #Number of cores for parallel processing
}

#Generating directories to store results
setwd(home_dir)
dir_names = c()

dir_names[1] = "Merged Replication Results"
if(!dir.exists(dir_names[1])){
  dir.create(dir_names[1])
}

#Multiplier of population to take into account the sampling ratio
multiplier = 1000 #2009-2014 NHANES modified sampling weight

#Generating directory to store plots
setwd(home_dir)
dir_name = "Figures"
if(!dir.exists(dir_name)){
  dir.create(dir_name)
}

line_size = 0.6 #size of lines
point_size_scatter = 0.5 #size of points in scatter plots
point_size_line = 3 #size of points in line plots
width = 5.5 #width of plots
width_scatter = 3.5 #width of scatter plots
height = 3.5 #height of plots
res = 300 #plot resolution

#Defining colors for density plots
cols =  colorRampPalette(c("#000099", "#00FEFF", "#45FE4F", 
                           "#FCFF00", "#FF9400", "#FF3100"))(256)

#Function to get density of points in two diemnsions
#Downloaded from: https://slowkow.com/notes/ggplot2-color-by-density/
get_density = function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Loading Model's Data ----------------------------------------------------

#Continuous NHANES Data
ptdata = fread("Continuous NHANES Synthetic Forecasted Dataset.csv") #2009-2014 continuous NHANES dataset

#Loading Mortality Rates
load("DeathData.RData")

#Loading statin discontinuation and restart probabilities
load("DiscontinuationRestartingRates.RData")

# Initializing Parameters -------------------------------------------------

#Simulation Parameters
reps = 1000 #Number of patient trajectory simulation replications
M = 10 #Number of time periods to simulate
discount = 0.97 #discount factor for QoL and costs

#Identification of states
alive = 1:6 #states at which the patient is alive
dead = 7:10 #states at which the patient is dead
zetachd = c(2,4,5) #states where chd risk is higher
zetastroke = c(3,4,6) #states where stroke risk is higher

#Treatment parameters
allpolicies = 1:3 #index for possible treatment options
allmeds = allpolicies-1 #notreatment and statin intensity level
numalltrt = length(allpolicies) #Total number of treatment options
mod_threshold = 0.1 #risk threshold for moderate intensity treatment in PCE policy
mod_genethreshold = 0.0847 #risk threshold for moderate intensity treatment in GenePCE policy
arr_mod_genethreshold = mod_genethreshold-(mod_genethreshold*(1-(37/40*0.21))) #ARR threshold as per RRR in Collins et al. (2016)
normalization_factor = c(0.944,0.889,0.889,0.889,1.389) #Normalization factors to adjust for effect of genetic information

#Transition probability parameters
N = 10 #Number of stages
numhealth = 10 #Number of health states
events = 2 #Number of events considered in model

#Risk parameters
time_tr = 10 #time period for risk calculations (treatment policies)
time_tp = 1 #time period for risk calculations (transitions in simulation)

zetaval = c(3,3) #multiplier for CHD and stroke, respectively
zeta = matrix(1,numhealth,events) #creating matrix
zeta[zetachd,1] = zetaval[1]
zeta[zetastroke,2] = zetaval[2]
zeta[dead,] = 0

#Testing parameters
test_year = 1 #last year at which genetics can be tested
genetic_cost = 200 #Cost of genetic testing
genetic_lb = 0.0768 #minimum ASCVD risk level to perform genetic testing (1 SD - base case)
genetic_ub = 0.1363 #maximum ASCVD risk level to perform genetic testing (1 SD - base case)
genetic_lb05 = 0.0871 #minimum ASCVD risk level to perform genetic testing (0.5 SD - sensitivity)
genetic_ub05 = 0.1162 #maximum ASCVD risk level to perform genetic testing (0.5 SD - sensitivity)
genetic_lb2 = 0.0624 #minimum ASCVD risk level to perform genetic testing (2 SD - sensitivity)
genetic_ub2 = 0.1911 #maximum ASCVD risk level to perform genetic testing (2 SD - sensitivity)

#Cost parameters
costvector = c(0,144,450) #costs of no treatment, moderate intensity, and high intensity statin ($/year)
cost_chd = 67155 #annual cost of a CHD event
cost_chd_hist = 4499 #annual cost of CHD history
cost_chd_fatal = 18634 #annual cost of a fatal CHD event
cost_stroke = 22143 #annual cost of a stroke event
cost_stroke_hist = 7100 #annual cost of stroke history
cost_stroke_fatal = 11495 #annual cost of a stroke CHD event

#Reward parameters
QoL = c(1,0.90,0.90,0.81,0.88,0.67,0,0,0,0) #QoL weights per state
QoLterm = c(1,0.90,0.90,0.81,0.90,0.90,0,0,0,0) #terminal condition QoL

disut = 0.001 #treatment disutility (assumed to be 0.001*statin intensity)
trtharm = disut*allmeds #total disutility of the treatment using trt disutility

#Standardize mortality rates (Pandya 2015, Smolina 2012, and Dennis 1993)
mortality_rates = list("Male_one_CHD" = c(1, 1/1.6, 1/2.3, (1/1.6)*(1/2.3), 1/1.6, 1/2.3, 0, 0, 0, 0), #Males <2 CHD events
                       "Female_one_CHD" = c(1, 1/2.1, 1/2.3, (1/2.1)*(1/2.3), 1/2.1, 1/2.3, 0, 0, 0, 0), #Females <2 CHD events
                       "Male_multiple_CHD" = c(1, 1/3.4, 1/2.3, (1/3.4)*(1/2.3), 1/3.4, 1/2.3, 0, 0, 0, 0), #Males >=2 CHD events
                       "Female_multiple_CHD" = c(1, 1/2.5, 1/2.3, (1/2.5)*(1/2.3), 1/2.5, 1/2.3, 0, 0, 0, 0) #Females >=2 CHD events
)

#Base case parameters
pce_correct = 0 #indicator for PCE being the true risk (0=GenePCE is the true risk)
adjustment = 0 #indicator of adjustment in RRR using normalization factors based on the PRS (0=no normalization)
or_sd = 1.67 #OR increase per SD of polygenic score
ascvd_rate = 0 #assuming ASCVD event rate is correct
dis_mult = 1 #no change in disutility due to age in main analysis
disc = 0 #indicator for dicontinuation and restart rates (0=no discontinuation)
adherence = 0 #assuming perfect adherence to medications in main analysis
age_subgroup = 0 #indicator of genetic testing based on age (<50,>=50) not risk range (0=no subgroup)
gender_subgroup = 0 #indicator of genetic testing based on gender (male,female) not risk range(0=no subgroup)
mortality_mult = 1 #assuming mortality rates are correct
prs_incorrect = 0 #indicator for Khera's PRS being incorrect while Abraham's PRS was correct (0=Khera's PRS is correct)

# PCE Risk Function -------------------------------------------------------
arisk = function(type,sex,race,age,sbp,smk,tc,hdl,diab,trt,time){

  # ASCVD risk calculator (2013 ACC/AHA Guideline)
  # inputs: type of risk to calculate (1=CHD, 2=stroke), sex (1=male, 0=female),
  # race (1=white,0=black), age, SBP, smoking status (1=smoker, 0=nonsmoker), total
  # cholesterol, HDL, diabetes status (1=diabetic, 0=nondiabetic), trt
  # (1=BP reported is on treatment, 0=BP reported is untreated), time.
  # outputs: likelihood of CHD or stroke in the next "time" years

  if (sex==1){ #male
    if (race==1){ #white
      b_age=12.344
      b_age2=0
      b_tc=11.853
      b_age_tc=-2.664
      b_hdl=-7.990
      b_age_hdl=1.769

      if (trt==1){ #SBP is treated SBP
        b_sbp=1.797
        b_age_sbp=0
      }else{ #SBP is untreated SBP
        b_sbp=1.764
        b_age_sbp=0
      }

      b_smk=7.837
      b_age_smk=-1.795
      b_diab=0.658
      meanz=61.18

      if (time==1){
        basesurv=0.99358
      }else if (time==5){
        basesurv=0.96254
      }else if (time == 10){
        basesurv=0.9144
      }else{
        stop(cat(time, "is an improper time length for risk calculation"))
      }
    }else{ #black
      b_age=2.469
      b_age2=0
      b_tc=0.302
      b_age_tc=0
      b_hdl=-0.307
      b_age_hdl=0

      if (trt==1){ #SBP is treated SBP
        b_sbp=1.916
        b_age_sbp=0
      }else{ #SBP is untreated SBP
        b_sbp=1.809
        b_age_sbp=0
      }

      b_smk=0.549
      b_age_smk=0
      b_diab=0.645
      meanz=19.54

      if (time==1){
        basesurv=0.99066
      }else if (time==5){
        basesurv=0.95726
      }else if (time==10){
        basesurv=0.8954
      }else{
        stop(cat(time, "is an improper time length for risk calculation"))
      }
    }
  }else{ #female
    if (race==1){ #white
      b_age=-29.799
      b_age2=4.884
      b_tc=13.540
      b_age_tc=-3.114
      b_hdl=-13.578
      b_age_hdl=3.149

      if (trt==1){ #SBP is treated SBP
        b_sbp=2.019
        b_age_sbp=0
      }else{ #SBP is untreated SBP
        b_sbp=1.957
        b_age_sbp=0
      }

      b_smk=7.574
      b_age_smk=-1.665
      b_diab=0.661
      meanz=-29.18

      if (time==1){
        basesurv=0.99828
      }else if (time==5){
        basesurv=0.98898
      }else if (time==10){
        basesurv=0.9665
      }else{
        stop(cat(time, "is an improper time length for risk calculation"))
      }
    }else{ #black
      b_age=17.114
      b_age2=0
      b_tc=0.940
      b_age_tc=0
      b_hdl=-18.920
      b_age_hdl=4.475

      if (trt==1){ #SBP is treated SBP
        b_sbp=29.291
        b_age_sbp=-6.432
      }else{ #SBP is untreated SBP
        b_sbp=27.820
        b_age_sbp=-6.087
      }

      b_smk=0.691
      b_age_smk=0
      b_diab=0.874
      meanz=86.61

      if (time==1){
        basesurv=0.99834
      }else if (time==5){
        basesurv=0.98194
      }else if (time==10){
        basesurv=0.9533
      }else{
        stop(cat(time, "is an improper time length for risk calculation"))
      }
    }
  }

  #proportion of ascvd assumed to be CHD or stroke, respectively
  eventprop = c(0.6,0.4)

  indivz = b_age*log(age) + b_age2*(log(age))^2 + b_tc*log(tc) + b_age_tc*log(age)*log(tc) + b_hdl*log(hdl) + b_age_hdl*log(age)*log(hdl) +
    b_sbp*log(sbp) + b_age_sbp*log(age)*log(sbp) + b_smk*smk + b_age_smk*log(age)*smk + b_diab*diab

  risk = eventprop[type]*(1-basesurv^(exp(indivz-meanz)))

}

# PRS Function ------------------------------------------------------------
gor = function(rand=NULL,or){

  #Function to calculate the odds for the incidence of CHD due to genetic factors.
  #These functions assume that the polygenic risk was
  #standardized to zero mean and unit variance and normally distributed

  #Use or = 1.67 for OR per SD from Khera et al. (2017)
  #Use or = 1.38 for OR per SD Abraham et al. (2016).

  #Generating random number to represent polygenic score
  if (missing(rand)){
    polyscore = rnorm(1)
  }else{
    set.seed(rand)
    polyscore = rnorm(1)
  }

  #Identifying genetic risk score quintile using normal distribution
  q1 = qnorm(1/5)
  q2 = qnorm(2/5)
  q3 = qnorm(3/5)
  q4 = qnorm(4/5)

  #Selecting quantile based on polygenic score random number
  if (polyscore < q1){
    q = 1
  }else if (q1 <= polyscore & polyscore < q2){
    q = 2
  }else if (q2 <= polyscore & polyscore < q3){
    q = 3
  }else if (q3 <= polyscore & polyscore < q4){
    q = 4
  }else if (q4 <= polyscore){
    q = 5
  }

  #Using genetic score to estimate CHD odds assuming OR = 1.38 or OR = 1.67 per SD (Khera et al. 2017)
  prs = exp(polyscore*log(or))

  poly = list("prs" = prs,"quintile" = q,"ps"=polyscore)

  return(poly)

}

# Post-treatment Risk Calculation Function --------------------------------
posttrtrisk = function(statins,pretrtrisk,adherence,grsquintile=1,normalizing=c(1,1,1,1,1)){

  #Assumes relative risk reductions from Collins et al. (2016)
  #This function has the option to add a normalizing factors
  #depending on the quintile of the PRS distribution

  if (statins==1){ #No treatment
    risk = pretrtrisk
  }else if(statins==2){ #Moderate intensity
    RR = 1-(37/40*0.21*(1+adherence)*normalizing[grsquintile]) #assuming that adherence only impacts the relative reductions
    risk = RR*pretrtrisk
  }else if(statins==3){ #High intensity
    RR = 1-(53/40*0.21*(1+adherence)*normalizing[grsquintile]) #assuming that adherence only impacts the relative reductions
    risk = RR*pretrtrisk
  }else{
    stop(cat(statins, "is not a valid treatment option"))
  }


  return(risk)
}

# Medications to Treatment Conversion Function ----------------------------
meds2trt = function(nummeds){

  # Converting number of medications to treatment intensity
  if (nummeds==0){
    treatment=1
  }else if(nummeds==1){
    treatment=2
  }else if(nummeds==2){
    treatment=3
  }else{
    stop(cat(nummeds," is not a medication option"))
  }

  return(treatment)
}

# Transition Probabilities Function ---------------------------------------
TP = function(periodrisk, chddeath, strokedeath, alldeath, adherence, grsquintile, adjusted=adjustment){
  
  #Calculating probability of health states transitions given ASCVD risk, fatality likelihoods,
  #and probability of mortality due to non-ASCVD events
  
  #Inputs
  ##periodrisk: 1-year risk of CHD and stroke
  ##chddeath: likelihood of death given a CHD events
  ##strokedeath: likelihood of death given a stroke events
  ##alldeath: likelihood of death due to non-ASCVD events
  ##adherence: probability of perfect adherence
  ##grsquintile: quintile of PRS (for adjusted threshold policies)
  ##adjusted: indicator of RRR adjustment using normalization factors (0=no, 1=yes; for adjusted threshold policies)
  
  #Line for debuggin purposes
  # periodrisk = c.risk.tp; h=1
  
  risk = array(NaN,dim = c(numhealth,N,events,numalltrt)) #stores post treatment risks
  ptrans = array(0,dim = c(numhealth,numhealth,N,numalltrt)) #state transition probabilities--default of 0, to reduce coding/computations
  
  #Incorporating normalizing factor to adjust for effect of genetic information in RRR
  if(adjusted==1){
    nf = normalization_factor
  }else{
    nf = c(1,1,1,1,1)
  }
  
  for (n in 1:N){ #number of decisions remaining (backward for induction purposes)
    t = N-n+1 #time (forward for purposes of input data)
    for (h in 1:numhealth){ #each health state
      for (j in 1:numalltrt){ #each treatment
        for (k in 1:events){ #each event type
          #Calculating post-treatment risks
          if(adjusted==1){ #run sensitivity analysis
            if(k==1){ #Genetic information was only available for CHD
              risk[h,t,k,j] = posttrtrisk(j,periodrisk[h,t,k],adherence,grsquintile,nf)
            }else{ #Stroke (no normalization)
              risk[h,t,k,j] = posttrtrisk(j,periodrisk[h,t,k],adherence)
            }
          }else{ #no normalization
            risk[h,t,k,j] = posttrtrisk(j,periodrisk[h,t,k],adherence)
          }
        }
        
        #Health state transition probabilities: allows for both CHD and stroke in same period
        #Let Dead state dominate the transition to all others
        if (h==7 | h==8 | h==9 | h==10){ #Dead
          ptrans[h,10,t,j] = 1 #must stay dead  [and go to the overall dead state]
        }else{ #alternate denotes the state that is default state if neither a CHD event, stroke, nor death [from CHD, stroke or other] occurs
          if (h==4){  #History of CHD and Stroke
            alternate=4
          }else if (h==6 | h==3){ #Stroke or History of Stroke
            alternate=3
          }else if (h==5 | h==2){ #CHD Event or History of CHD
            alternate=2
          }else{ #Healthy
            alternate=1
          }
          quit = 0
          while (quit==0){ #compute transition probabilities, using a "break" command if you've exceeded 1 (this never happens!)
            ptrans[h,9,t,j] = min(1,strokedeath[t]*risk[h,t,2,j]) #likelihood of death from stroke
            cumulprob = ptrans[h,9,t,j]
            
            ptrans[h,8,t,j] = min(1,chddeath[t]*risk[h,t,1,j]) #likelihood of death from CHD event
            if (cumulprob+ptrans[h,8,t,j]>=1){ #check for invalid probabilities
              ptrans[h,8,t,j]=1-cumulprob
              break} #all other probabilities should be left as 0 [as initialized before FOR loop]
            cumulprob = cumulprob+ptrans[h,8,t,j]
            
            ptrans[h,7,t,j] = min(1,alldeath[t]) #likelihood of death from non CVD cause
            if (cumulprob+ptrans[h,7,t,j]>=1){ #check for invalid probabilities
              ptrans[h,7,t,j] = 1-cumulprob
              break} #all other probabilities should be left as 0 [as initialized before FOR loop]
            cumulprob = cumulprob+ptrans[h,7,t,j]
            
            ptrans[h,6,t,j] = min(1,(1-strokedeath[t])*risk[h,t,2,j]) #likelihood of having stroke and surviving
            if (cumulprob+ptrans[h,6,t,j]>=1){ #check for invalid probabilities
              ptrans[h,6,t,j] = 1-cumulprob
              break} #all other probabilities should be left as 0 [as initialized before FOR loop]
            cumulprob = cumulprob+ptrans[h,6,t,j]
            
            ptrans[h,5,t,j] = min(1,(1-chddeath[t])*risk[h,t,1,j]) #likelihood of having a CHD event and surviving
            if (cumulprob+ptrans[h,5,t,j]>=1){ #check for invalid probabilities
              ptrans[h,5,t,j] = 1-cumulprob
              break} #all other probabilities should be left as 0 [as initialized before FOR loop]
            cumulprob = cumulprob+ptrans[h,5,t,j]
            
            ptrans[h,alternate,t,j] = 1-cumulprob #otherwise, you go to the alternate state
            break #computed all probabilities, now quit
          }
        }
        
      } #[j]
    } #[h]
  } #[n]
  
  transitions = list("posttrtrisk"=risk,"ptrans" = ptrans)

  return(transitions)
}

# Risk Correction Function ------------------------------------------------
riskcorrection = function(patientdata,or_sd=or_sd){
  
  #This function ensures that the overall event rates is the same for PCE risk and GenePCE risk,
  #accounting for asymmetry in the distribution of the PRS.
  
  #Making sure patient's data is in data table format
  patientdata = setDT(patientdata)
  initial_risk = patientdata[!duplicated(ID),]
  
  #10-year risk for CHD calculation
  initial_risk[,chdrisk := arisk(1,sex,race,age,sbp,smk,tc,hdl,diab,0,time_tr),by=ID]
  
  #Simulating genetic test (only genetic risk for CHD parameters was available)
  #Genetic risk using Khera et al. or Abraham et al.
  initial_risk[,c("prs","quintile","ps"):=gor(rand = ID,or=or_sd),by=ID]
  
  #Converting clinical risk to odds
  initial_risk[,chdodds:=chdrisk/(1-chdrisk),by=ID]
  
  #Combining clinical and genetic odds
  initial_risk[,combinedchdodds:=chdodds*prs,by=ID]
  
  #Calculating mean odds
  average_clinical_odds = mean(initial_risk[,chdodds])
  average_combined_odds = mean(initial_risk[,combinedchdodds])
  
  #Correcting PRS
  correction_factor = average_clinical_odds/average_combined_odds
  initial_risk[,prs:=prs*correction_factor]
  
  patientdata$ps = rep(initial_risk$ps,each=10)
  patientdata$prs = rep(initial_risk$prs,each=10)
  patientdata$gquintile = rep(initial_risk$quintile,each=10)
  
  return(patientdata)
}

# Risk Calculation Function -----------------------------------------------
risksimulation = function(id){
  
  # #Line for debugging purposes
  # id = 140;
  
  #Storing patient data
  patientdata = ptdata[ID == id,] #extract all of the patient's data from the larger data matrix
  agemat = patientdata$age[1:N] #saving age of the patient in simulation
  diab = patientdata$diab[1] #saving diabetes status of each patient in simulation
  women = ifelse(patientdata$sex[1]==0,1,0) #saving sex of each patient in simulation
  black = ifelse(patientdata$race[1]==0,1,0) #saving race of each patient in simulation
  ldl = ifelse(patientdata$ldl[1:N]>160,1,0) #saving LDL level of each patient in simulation
  smoke = patientdata$smk[1] #saving smoking status of each patient in simulation
  sbp = patientdata$sbp[1:N] #saving sbp of each patient in simulation
  dbp = (patientdata)$dbp[1:N] #saving dbp of each patient in simulation
  
  #Storing risk calculations
  clinicalrisk.tr = array(0,dim=c(numhealth,N,events)) #10-y CHD and stroke risk (for threshold policies)
  periodrisk.tr = array(0,dim=c(numhealth,N,events)) #10-y risk after scaling
  c.risk.tr = array(0,dim=c(numhealth,N,events)) #10-y CHD combined risk information
  
  clinicalrisk.tp = array(0,dim=c(numhealth,N,events)) #1-y CHD and stroke risk (for transition probabilities)
  periodrisk.tp = array(0,dim=c(numhealth,N,events)) #1-y risk after scaling
  c.risk.tp = array(0,dim=c(numhealth,N,events)) #1-y CHD combined risk information
  
  
  for (h in 1:numhealth){ #each state
    for (t in 1:N){ #each age
      
      #Changing scaling factor according to the age of the patient
      if(patientdata$age[t]>=60){
        zeta_sim = zeta
        zeta_sim[zetastroke,2] = 2
      }else{
        zeta_sim = zeta
      }
      
      for (k in 1:events){ #each event type
        #10-year ASCVD risk calculation (for threshold policies)
        a.risk.tr = arisk(k,patientdata$sex[t],patientdata$race[t],patientdata$age[t],
                          patientdata$sbp[t],patientdata$smk[t],patientdata$tc[t],patientdata$hdl[t],
                          patientdata$diab[t],0,time_tr)
        
        clinicalrisk.tr[h,t,k] = a.risk.tr
        
        #1-year ASCVD risk calculation (for transition probabilities)
        a.risk.tp = arisk(k,patientdata$sex[t],patientdata$race[t],patientdata$age[t],
                          patientdata$sbp[t],patientdata$smk[t],patientdata$tc[t],patientdata$hdl[t],
                          patientdata$diab[t],0,time_tp)
        
        clinicalrisk.tp[h,t,k] = a.risk.tp
        
        if(zeta_sim[h,k] > 1){
          #Scaling odds of 10-year risks
          periododds = clinicalrisk.tr[h,t,k]/(1-clinicalrisk.tr[h,t,k])
          periododds = zeta_sim[h,k]*periododds
          periodrisk.tr[h,t,k] = periododds/(1+periododds)
          
          #Scaling odds of 1-year risks
          periododds = clinicalrisk.tp[h,t,k]/(1-clinicalrisk.tp[h,t,k])
          periododds = zeta_sim[h,k]*periododds
          periodrisk.tp[h,t,k] = periododds/(1+periododds)
        }else if(zeta_sim[h,k]==0){ #set risk to 0
          periodrisk.tr[h,t,k] = 0
          periodrisk.tp[h,t,k] = 0
        }else{ #no scale
          periodrisk.tr[h,t,k] = clinicalrisk.tr[h,t,k]
          periodrisk.tp[h,t,k] = clinicalrisk.tp[h,t,k]
        }
        
        
        if(prs_incorrect==1){ #making decisions according to the wrong prs (sensitivity analysis)
          #if CHD is being considered (only genetic risk for CHD parameters was available)
          if(k==1){
            #Combining clinical and genetic risks in 10-year risks
            odds.tr = periodrisk.tr[h,t,k]/(1-periodrisk.tr[h,t,k])
            c.odds = odds.tr*patientdata$prs_incorrect[1]
            c.risk.tr[h,t,k] = c.odds/(1+c.odds)
            
            #Combining clinical and genetic risks in 1-year risks
            odds.tp = periodrisk.tp[h,t,k]/(1-periodrisk.tp[h,t,k])
            c.odds = odds.tp*patientdata$prs[1]
            c.risk.tp[h,t,k] = c.odds/(1+c.odds)
            
          }else{ #no genetic risk for stroke was available (k=2)
            c.risk.tr[h,t,k] = periodrisk.tr[h,t,k]
            c.risk.tp[h,t,k] = periodrisk.tp[h,t,k]
          }
        }else{ #correct prs (base case)
          #if CHD is being considered (only genetic risk for CHD parameters was available)
          if(k==1){
            #Combining clinical and genetic risks in 10-year risks
            odds.tr = periodrisk.tr[h,t,k]/(1-periodrisk.tr[h,t,k])
            c.odds = odds.tr*patientdata$prs[1]
            c.risk.tr[h,t,k] = c.odds/(1+c.odds)
            
            #Combining clinical and genetic risks in 1-year risks
            odds.tp = periodrisk.tp[h,t,k]/(1-periodrisk.tp[h,t,k])
            c.odds = odds.tp*patientdata$prs[1]
            c.risk.tp[h,t,k] = c.odds/(1+c.odds)
            
          }else{ #no genetic risk for stroke was available (k=2)
            c.risk.tr[h,t,k] = periodrisk.tr[h,t,k]
            c.risk.tp[h,t,k] = periodrisk.tp[h,t,k]
          }
        }
        
        
      } #events [k]
    } #age [t]
  } #states [h]
  
  #CVD event rate multiplier (for sensitivity analysis)
  periodrisk.tp = periodrisk.tp*(1+ascvd_rate)
  periodrisk.tr = periodrisk.tr*(1+ascvd_rate)
  c.risk.tp = c.risk.tp*(1+ascvd_rate)
  c.risk.tr = c.risk.tr*(1+ascvd_rate)
  
  #Ensuring validity of risk estimates
  periodrisk.tp = pmin(periodrisk.tp,1) #prevent risk from exceeding 1
  periodrisk.tp = pmax(periodrisk.tp,0) #prevent risk from going below 0
  
  periodrisk.tr = pmin(periodrisk.tr,1) #prevent risk from exceeding 1
  periodrisk.tr = pmax(periodrisk.tr,0) #prevent risk from going below 0
  
  c.risk.tp = pmin(c.risk.tp,1) #prevent risk from exceeding 1
  c.risk.tp = pmax(c.risk.tp,0) #prevent risk from going below 0
  
  c.risk.tr = pmin(c.risk.tr,1) #prevent risk from exceeding 1
  c.risk.tr = pmax(c.risk.tr,0) #prevent risk from going below 0
  
  riskresults = list("id"=id,"age"=agemat,"diab"=diab,"women"=women,"black"=black,"ldl"=ldl,"smoke"=smoke,"sbp"=sbp,"dbp"=dbp,
                     "ps"=patientdata$ps[1],"gquintile"=patientdata$gquintile[1],"prs"=patientdata$prs[1],
                     "periodrisk.tr"=periodrisk.tr,"c.risk.tr"=c.risk.tr,
                     "periodrisk.tp"=periodrisk.tp,"c.risk.tp"=c.risk.tp)
  
  return(riskresults)
}

# Transition Probabilities Calculation Function ---------------------------
tpsimulation = function(id){
  
  # #Line for debugging purposes
  # id = 19;
  
  #Extracting results from risk simulation to environment
  list2env(risk_results[[id]],.GlobalEnv)
  
  #Loading death rates information
  if(women==0){ #male
    sexcol = 2 #column in deathdata corresponding to male
  }else{
    sexcol = 3 #column in deathdata corresponding to female
  }
  
  #Loading death rates
  chddeath = chddeathdata[is.element(chddeathdata[,1],age),sexcol]
  strokedeath = strokedeathdata[is.element(strokedeathdata[,1],age),sexcol]
  alldeath = alldeathdata[is.element(alldeathdata[,1],age),sexcol]
  
  #Life expectancy
  healthy_lifexp = lifedata[is.element(lifedata[,1],age[N]+1),sexcol]
  
  #Determine transition probabilities
  if(pce_correct==1){
    transitions = TP(periodrisk.tp, chddeath, strokedeath, alldeath, adherence)
  }else if(adjustment==1){
    transitions = TP(c.risk.tp, chddeath, strokedeath, alldeath, adherence, gquintile, adjustment)
  }else{
    transitions = TP(c.risk.tp, chddeath, strokedeath, alldeath, adherence)
  }
  
  tpresults = list("id"=id,"age"=age,"diab"=diab,"women"=women,"black"=black,"ldl"=ldl,"smoke"=smoke,"sbp"=sbp,"dbp"=dbp,
                   "periodrisk.tr"=periodrisk.tr,"c.risk.tr"=c.risk.tr,"prs"=prs,"gquintile"=gquintile,
                   "healthy_lifexp"=healthy_lifexp,"transitions"=transitions)
  
  return(tpresults)  
}

# Patient's Trajectory Simulation Function --------------------------------
ptsimulation = function(id){
  
  #Extracting results from treatment simulation to environment
  list2env(tp_results[[id]],.GlobalEnv)
  
  #Arrays to store results
  alive.notrt = alive10.notrt = 0 #years alive and indicator of 10 years alive under no treatment
  risk.notrt.tr = risk.notrt.gtr = risk.tr = risk.gtr = array(NA, dim = c(M)) #store total risk
  states.notrt = states.tr = states.gtr = array(NA,dim = c(M)) #keeping track of state transitions
  qalys.notrt = qalys.tr = qalys.gtr = array(0, dim = c(M)) #QALYs  
  chd.nonfatal.notrt = chd.nonfatal.tr = chd.nonfatal.gtr = array(0, dim = c(M)) #number of non fatal chd events  
  stroke.nonfatal.notrt = stroke.nonfatal.tr = stroke.nonfatal.gtr = array(0, dim = c(M)) #number of non fatal stroke events  
  chd.fatal.notrt = chd.fatal.tr = chd.fatal.gtr = array(0, dim = c(M)) #number of fatal chd events  
  stroke.fatal.notrt = stroke.fatal.tr = stroke.fatal.gtr = array(0, dim = c(M)) #number of fatal stroke events  
  cost.notrt = cost.tr = cost.gtr = array(NA, dim = c(M)) #total cost under 
  cost.chd.notrt = cost.chd.tr = cost.chd.gtr = array(0, dim = c(M)) #cost of chd events  
  cost.stroke.notrt = cost.stroke.tr = cost.stroke.gtr = array(0, dim = c(M)) #cost of stroke events  
  
  policy.tr = policy.gtr = array(NA, dim = c(M)) #policies in number of meds
  cost.meds.tr = cost.meds.gtr = array(NA, dim = c(M)) #treatment cost  
  d.tr = d.gtr = array(0, dim = c(M)) #Tracking periods at which patients discontinue treatment
  r.tr = r.gtr = array(0, dim = c(M)) #Tracking periods at which patients restart treatment
  
  cost.gt.gtr = array(0, dim = c(M)) #cost of genetic testing 
  
  if(adjustment==1){ #run sesitivity analysis
    #Calculating ARR of moderate intensity treatment
    pretrtrisk = periodrisk.tr[,,1]+periodrisk.tr[,,2] #combining pre-treatment risk of CHD and stroke
    comprisk = array(0,dim = c(numhealth,N,events)) #Post-treatment risk at each state, age and event (ASCVD component)
    posttrtrisk = array(0,dim = c(numhealth,N)) #Post-treatment risk at each state, and age
    arr = array(0,dim = c(numhealth,N)) #Absolute risk reduction at each state, and age
    for (h in seq(numhealth)){ #each health state
      for (t in seq(N)){ #each age
        for (k in seq(events)){ #each event
          #Calculating post-treatment risks after moderate intensity statins
          if(k==1){ #Genetic information was only available for CHD
            comprisk[h,t,k] = posttrtrisk(2,periodrisk.tr[h,t,k],adherence,gquintile,normalization_factor)
          }else{ #Stroke (no normalization)
            comprisk[h,t,k] = posttrtrisk(2,periodrisk.tr[h,t,k],adherence)
          }
        }
        #Combining post-treatment risks of ASCVD components
        posttrtrisk = comprisk[,,1] + comprisk[,,2]
        
        #Calculating ARR
        arr[h,t] = pretrtrisk[h,t] - posttrtrisk[h,t]
      }
    }
  }
  
  #Simulate trajectory of patient over planning horizon
  h.notrt = h.tr = h.gtr = 1 #assume every patient is healthy at first (from NHANES dataset)
  gt.notrt = gt.gtr = 0 #creating indicator for genetic testing
  gt.gtr.cost = 0 #creating indicator for genetic testing cost
  for (t in 1:M){ #each age
    
    #Storing 10-year risk information under no treatment
    risk.notrt.tr[t] = sum(periodrisk.tr[h.notrt,t,]) 
    
    #Storing 10-year risk information under the PCE policy
    risk.tr[t] = sum(periodrisk.tr[h.tr,t,]) 
    
    #Storing 10-year risk information under the GenePCE policy
    risk.gtr[t] = sum(periodrisk.tr[h.gtr,t,])
    
    #Indicators of genetics playing a role in risks
    if(age_subgroup==1){ #Subgroup <50
      if(gt.notrt==0 & t<=test_year & age[t]<50){
        gt.notrt = 1
        gt.gtr = 1
      }
    }else if(age_subgroup==2){ #Subgroup >=50
      if(gt.notrt==0 & t<=test_year & age[t]>=50){
        gt.notrt = 1
        gt.gtr = 1
      }
    }else if(gender_subgroup==1){ #Male subgroup
      if(gt.notrt==0 & t<=test_year & women==0){
        gt.notrt = 1
        gt.gtr = 1
      }
    }else if(gender_subgroup==2){ #Female subgroup
      if(gt.notrt==0 & t<=test_year & women==1){
        gt.notrt = 1
        gt.gtr = 1
      }
    }else{ #base case
      #Indicator of genetics playing a role in untreated risk
      if(gt.notrt==0 & t<=test_year & genetic_lb<=risk.notrt.tr[t] & risk.notrt.tr[t]<=genetic_ub){
        gt.notrt = 1
      }
      
      #Indicator of genetics playing a role in decision making under the GenePCE policy
      if(gt.gtr==0 & t<=test_year & genetic_lb<=risk.gtr[t] & risk.gtr[t]<=genetic_ub){
        gt.gtr = 1
      }
    }
    
    #Storing 10-year risk information under no treatment (with genetics if applies)
    if(gt.notrt==1){
      risk.notrt.gtr[t] = sum(c.risk.tr[h.notrt,t,])
    }else{
      risk.notrt.gtr[t] = risk.notrt.tr[t]
    }
    
    #Modifiying with genetic information (if applies)
    if(gt.gtr==1){
      risk.gtr[t] = sum(c.risk.tr[h.gtr,t,])
    }
    
    
    #Simulating state transitions using all available information
    h.next = array(NA, dim = c(numhealth,numalltrt)) #Store current states in simulation replicate
    for(h in 1:numhealth){ #loop over possible current states
      set.seed(id+r+t+h)
      u = runif(1) #random uniform number (0,1)
      ptrans = apply(transitions$ptrans[h,,t,],2,cumsum) #cummulative sums by treatment intensity over transition probabilities
      for(j in 1:numalltrt){ #loop over treatment intensities
        for(hh in 1:numhealth){ #loop over possible next states
          if(hh==1){ #transition to healthy
            if(u<=ptrans[hh,j]){
              h.next[h,j] = hh
            }
          }else{ #transition to other states
            if(u>ptrans[(hh-1),j] & u<=ptrans[hh,j]){
              h.next[h,j] = hh
            }
          }
        } #[hh]
      } #[j]
    } #[h]
    
    
    #Doubling disutility if age > 70 (for sensitivity analysis)
    if(age[t]>70){
      trtharm_sim = trtharm*dis_mult
    }else{
      trtharm_sim = trtharm
    }
    
    ##No treatment policies
    #Simulating state from transition probabilities
    h.notrt = h.next[h.notrt,1]
    
    #Adjusting state for history of both events
    if(is.element(5,states.notrt)){
      if(states.notrt[(t-1)]==6 & h.notrt==3){
        h.notrt = 4
      }
    }
    if(is.element(6,states.notrt)){
      if(states.notrt[(t-1)]==5 & h.notrt==2){
        h.notrt = 4
      }
    }
    
    #Storing current state
    states.notrt[t] = h.notrt
    
    #Computing number of years alive for no treatment policy
    if(h.notrt<7){
      alive.notrt = alive.notrt + 1
    }
    
    #Calculating expected life and QALYs based on events and gender
    if(t==M){
      if(is.element(5,states.notrt[duplicated(states.notrt)])){ #History of 2 or more CHD events
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[3]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[4]]
        }
      }else{
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[1]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[2]]
        }
      }
      
      Lterm = healthy_lifexp*SMR*mortality_mult #Terminal life years
      Vterm = Lterm*QoLterm #Terminal QALYs
      
    }
    
    #Computing number and cost of CHD events for no treatment policy
    if(h.notrt==2|h.notrt==4){
      if(t==M){ #Cost of history of CHD at time M plus cost for life expectancy
        cost.chd.notrt[t] = cost_chd_hist*discount^(t) + cost_chd_hist*Lterm[h.notrt]*discount^(M+1)
      }else{ #Cost of history of CHD
        cost.chd.notrt[t] = cost_chd_hist*discount^(t)
      }
    }
    if(h.notrt==5){
      chd.nonfatal.notrt[t] = 1
      if(t==M){ #Cost of CHD at time M plus cost for life expectancy
        cost.chd.notrt[t] = cost_chd*discount^(t) + cost_chd_hist*Lterm[h.notrt]*discount^(M+1)
      }else{ #Cost of CHD
        cost.chd.notrt[t] = cost_chd*discount^(t)
      }
    }
    if(h.notrt==8){ #Cost of fatal CHD
      chd.fatal.notrt[t] = 1
      cost.chd.notrt[t] = cost_chd_fatal*discount^(t)
    }
    
    #Computing number and cost of stroke events for no treatment policy
    if(h.notrt==3|h.notrt==4){
      if(t==M){ #Cost of history of stroke at time M plus cost for life expectancy
        cost.stroke.notrt[t] = cost_stroke_hist*discount^(t) + cost_stroke_hist*Lterm[h.notrt]*discount^(M+1)
      }else{ #Cost of history of stroke
        cost.stroke.notrt[t] = cost_stroke_hist*discount^(t)
      }
    }
    if(h.notrt==6){
      stroke.nonfatal.notrt[t] = 1
      if(t==M){ #Cost of stroke at time M plus cost for life expectancy
        cost.stroke.notrt[t] = cost_stroke*discount^(t) + cost_stroke_hist*Lterm[h.notrt]*discount^(M+1)
      }else{ #Cost of stroke
        cost.stroke.notrt[t] = cost_stroke*discount^(t)
      }
    }
    if(h.notrt==9){ #Cost of fatal stroke
      stroke.fatal.notrt[t] = 1
      cost.stroke.notrt[t] = cost_stroke_fatal*discount^(t)
    }
    
    #Computing total cost using no treatment policy
    cost.notrt[t] = cost.chd.notrt[t] + cost.stroke.notrt[t]
    
    #Computing QALYs for no treatment policy
    if (t==M){
      qalys.notrt[t] = max(0,(QoL[h.notrt]-trtharm_sim[1])*discount^(t)) + Vterm[h.notrt]*discount^(M+1)
    }else{
      qalys.notrt[t] = max(0,(QoL[h.notrt]-trtharm_sim[1])*discount^(t))
    }
    
    
    ##PCE policies
    #Determine PCE policies
    if(risk.tr[t]>=mod_threshold){
      policy.tr[t] = 1 
    }else{
      policy.tr[t] = 0
    }
    
    #Adjusting to high intensity statins if the patient has history of ASCVD events
    if(h.tr>=2 & h.tr<=6){
      policy.tr[t] = 2
    }
    
    if(disc==1){ #Run only for sensitivity analysis
      #Simulating statin discontinuation
      #Discontinuation can only happen after the first year of the planning horizon and if patient is receiving treatment
      if(t>1 & max(r.tr[1:(t-1)])>=max(d.tr[1:(t-1)]) & policy.tr[t]>0){
        set.seed(id+r+t) #pseudo-random number sequence
        u = runif(1) #random uniform number (0,1)
        index = t-tail(which(r.tr[1:(t-1)]==max(r.tr[1:(t-1)])),1)
        if(u<=dis_res_prob$DiscontinuationRate[index]){
          d.tr[t] = t
        }
      }
      
      #Simulating statin restarting
      if(t>2 & max(d.tr[1:(t-1)])>max(r.tr[1:(t-1)]) & h.tr<7){ #Restarting can only happen after the second year of the planning horizon and if patient is alive
        set.seed(id+r+t+1) #pseudo-random number sequence
        u = runif(1) #random uniform number (0,1)
        index = t-which.max(d.tr[1:(t-1)])
        if(u<=dis_res_prob$RestartingRate[index]){
          r.tr[t] = t
        }
      }
    }
    
    #Discontinuing treatment
    if(max(d.tr)>max(r.tr)){
      policy.tr[t] = 0
    }
    
    #Storing treatment cost
    if(max(d.tr)>max(r.tr)){
      cost.meds.tr[t] = 0 #statin costs
    }else{
      cost.meds.tr[t] = costvector[meds2trt(policy.tr[t])]*discount^(t) #statin costs
    }
    
    #Simulating state from transition probabilities
    h.tr = h.next[h.tr,meds2trt(policy.tr[t])]
    
    #Adjusting state for history of both events
    if(is.element(5,states.tr)){
      if(states.tr[(t-1)]==6 & h.tr==3){
        h.tr = 4
      }
    } 
    if(is.element(6,states.tr)){
      if(states.tr[(t-1)]==5 & h.tr==2){
        h.tr = 4
      }
    }
    
    #Storing current state
    states.tr[t] = h.tr
    
    #Calculating expected life and QALYs based on events and gender
    if(t==M){
      if(is.element(5,states.tr[duplicated(states.tr)])){ #History of 2 or more CHD events
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[3]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[4]]
        }
      }else{
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[1]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[2]]
        }
      }
      Lterm = healthy_lifexp*SMR*mortality_mult #Terminal life years
      Vterm = Lterm*QoLterm #Terminal QALYs
      
    }
    
    #Computing number and cost of CHD events for PCE policy
    if(h.tr==2|h.tr==4){
      if(t==M){ #Cost of history of CHD at time M plus cost for life expectancy
        cost.chd.tr[t] = cost_chd_hist*discount^(t) + cost_chd_hist*Lterm[h.tr]*discount^(M+1)
      }else{ #Cost of history of CHD
        cost.chd.tr[t] = cost_chd_hist*discount^(t)
      }
    }
    if(h.tr==5){
      chd.nonfatal.tr[t] = 1
      if(t==M){ #Cost of CHD at time M plus cost for life expectancy
        cost.chd.tr[t] = cost_chd*discount^(t) + cost_chd_hist*Lterm[h.tr]*discount^(M+1)
      }else{ #Cost of CHD
        cost.chd.tr[t] = cost_chd*discount^(t)
      }
    }
    if(h.tr==8){ #Cost of fatal CHD
      chd.fatal.tr[t] = 1
      cost.chd.tr[t] = cost_chd_fatal*discount^(t)
    }
    
    #Computing number and cost of stroke events for PCE policy
    if(h.tr==3|h.tr==4){
      if(t==M){ #Cost of history of stroke at time M plus cost for life expectancy
        cost.stroke.tr[t] = cost_stroke_hist*discount^(t) + cost_stroke_hist*Lterm[h.tr]*discount^(M+1)
      }else{ #Cost of history of stroke
        cost.stroke.tr[t] = cost_stroke_hist*discount^(t)
      }
    }
    if(h.tr==6){
      stroke.nonfatal.tr[t] = 1
      if(t==M){ #Cost of stroke at time M plus cost for life expectancy
        cost.stroke.tr[t] = cost_stroke*discount^(t) + cost_stroke_hist*Lterm[h.tr]*discount^(M+1)
      }else{ #Cost of stroke
        cost.stroke.tr[t] = cost_stroke*discount^(t)
      }
    }
    if(h.tr==9){ #Cost of fatal stroke
      stroke.fatal.tr[t] = 1
      cost.stroke.tr[t] = cost_stroke_fatal*discount^(t)
    }
    
    #Computing total cost using PCE policy
    cost.tr[t] = cost.chd.tr[t] + cost.stroke.tr[t] + cost.meds.tr[t]
    
    #Computing QALYs for PCE policy
    if (t==M){
      qalys.tr[t] = max(0,(QoL[h.tr]-trtharm_sim[meds2trt(policy.tr[t])])*discount^(t)) + Vterm[h.tr]*discount^(M+1)
    }else{
      qalys.tr[t] = max(0,(QoL[h.tr]-trtharm_sim[meds2trt(policy.tr[t])])*discount^(t))
    }
    
    ##GenePCE policies
    if(gt.gtr==1){
      if(adjustment==1){#Determine adjusted GenePCE policies
        if(arr[h.gtr,t]>=arr_mod_genethreshold){
          policy.gtr[t] = 1 
        }else{
          policy.gtr[t] = 0
        }
      }else{#Determine GenePCE policies
        if(risk.gtr[t]>=mod_genethreshold){
          policy.gtr[t] = 1 
        }else{
          policy.gtr[t] = 0
        }
      }
    }else{
      if(risk.gtr[t]>=mod_threshold){
        policy.gtr[t] = 1 
      }else{
        policy.gtr[t] = 0
      }
    }
    
    
    #Adjusting to high intensity statins if the patient has history of ASCVD events
    if(h.gtr>=2 & h.gtr<=6){
      policy.gtr[t] = 2
    }
    
    if(disc==1){ #Run only for sensitivity analysis
      #Simulating statin discontinuation
      #Discontinuation can only happen after the first year of the planning horizon and if patient is receiving treatment
      if(t>1 & max(r.gtr[1:(t-1)])>=max(d.gtr[1:(t-1)]) & policy.gtr[t]>0){
        set.seed(id+r+t) #pseudo-random number sequence
        u = runif(1) #random uniform number (0,1)
        index = t-tail(which(r.gtr[1:(t-1)]==max(r.gtr[1:(t-1)])),1)
        if(u<=dis_res_prob$DiscontinuationRate[index]){
          d.gtr[t] = t
        }
      }
      
      #Simulating statin restarting
      if(t>2 & max(d.gtr[1:(t-1)])>max(r.gtr[1:(t-1)]) & h.gtr<7){ #Restarting can only happen after the second year of the planning horizon and if patient is alive
        set.seed(id+r+t+1) #pseudo-random number sequence
        u = runif(1) #random uniform number (0,1)
        index = t-which.max(d.gtr[1:(t-1)])
        if(u<=dis_res_prob$RestartingRate[index]){
          r.gtr[t] = t
        }
      }
    }
    
    #Discontinuing treatment
    if(max(d.gtr)>max(r.gtr)){
      policy.gtr[t] = 0
    }
    
    #Storing treatment cost
    if(max(d.gtr)>max(r.gtr)){
      cost.meds.gtr[t] = 0 #statin costs
    }else{
      cost.meds.gtr[t] = costvector[meds2trt(policy.gtr[t])]*discount^(t)
    }
    
    #Storing genetic testing cost
    if(gt.gtr==1 & gt.gtr.cost==0){
      cost.gt.gtr[t] = genetic_cost*discount^(t) #storing genetic cost
      gt.gtr.cost = 1 #Indicator of genetic cost 
    }
    
    #Simulating state from transition probabilities
    h.gtr = h.next[h.gtr,meds2trt(policy.gtr[t])]
    
    #Adjusting state for history of both events
    if(is.element(5,states.gtr)){
      if(states.gtr[(t-1)]==6 & h.gtr==3){
        h.gtr = 4
      }
    }
    if(is.element(6,states.gtr)){
      if(states.gtr[(t-1)]==5 & h.gtr==2){
        h.gtr = 4
      }
    }
    
    #Storing current state
    states.gtr[t] = h.gtr
    
    #Calculating expected life and QALYs based on events and gender
    if(t==M){
      if(is.element(5,states.gtr[duplicated(states.gtr)])){ #History of 2 or more CHD events
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[3]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[4]]
        }
      }else{
        if(women==0){ #Male mortality rates
          SMR = mortality_rates[[1]]
        }else{ #Female mortality rates
          SMR = mortality_rates[[2]]
        }
      }
      
      Lterm = healthy_lifexp*SMR*mortality_mult #Terminal life years
      Vterm = Lterm*QoLterm #Terminal QALYs
      
    }
    
    #Computing number and cost of CHD events
    if(h.gtr==2|h.gtr==4){
      if(t==M){ #Cost of history of CHD at time M plus cost for life expectancy
        cost.chd.gtr[t] = cost_chd_hist*discount^(t) + cost_chd_hist*Lterm[h.gtr]*discount^(M+1)
      }else{ #Cost of history of CHD
        cost.chd.gtr[t] = cost_chd_hist*discount^(t)
      }
    }
    if(h.gtr==5){
      chd.nonfatal.gtr[t] = 1
      if(t==M){ #Cost of CHD at time M plus cost for life expectancy
        cost.chd.gtr[t] = cost_chd*discount^(t) + cost_chd_hist*Lterm[h.gtr]*discount^(M+1)
      }else{ #Cost of CHD
        cost.chd.gtr[t] = cost_chd*discount^(t)
      }
    }
    if(h.gtr==8){ #Cost of fatal CHD
      chd.fatal.gtr[t] = 1
      cost.chd.gtr[t] = cost_chd_fatal*discount^(t)
    }
    
    #Computing number and cost of stroke events
    if(h.gtr==3|h.gtr==4){
      if(t==M){ #Cost of history of stroke at time M plus cost for life expectancy
        cost.stroke.gtr[t] = cost_stroke_hist*discount^(t) + cost_stroke_hist*Lterm[h.gtr]*discount^(M+1)
      }else{ #Cost of history of stroke
        cost.stroke.gtr[t] = cost_stroke_hist*discount^(t)
      }
    }
    if(h.gtr==6){
      stroke.nonfatal.gtr[t] = 1
      if(t==M){ #Cost of stroke at time M plus cost for life expectancy
        cost.stroke.gtr[t] = cost_stroke*discount^(t) + cost_stroke_hist*Lterm[h.gtr]*discount^(M+1)
      }else{ #Cost of stroke
        cost.stroke.gtr[t] = cost_stroke*discount^(t)
      }
    }
    if(h.gtr==9){ #Cost of fatal stroke
      stroke.fatal.gtr[t] = 1
      cost.stroke.gtr[t] = cost_stroke_fatal*discount^(t)
    }
    
    #Computing total cost
    cost.gtr[t] = cost.chd.gtr[t] + cost.stroke.gtr[t] + cost.meds.gtr[t] + cost.gt.gtr[t]
    
    #Computing QALYs
    if (t==M){
      qalys.gtr[t] = max(0,(QoL[h.gtr]-trtharm_sim[meds2trt(policy.gtr[t])])*discount^(t)) + Vterm[h.gtr]*discount^(M+1)
    }else{
      qalys.gtr[t] = max(0,(QoL[h.gtr]-trtharm_sim[meds2trt(policy.gtr[t])])*discount^(t))
    }
    
  } #[t]
  
  #Indicator of 10-year survival
  if (h.notrt<7){
    alive10.notrt = 1
  }
  
  saveinfo = paste0("Replication ", r," of trajectory simulation for patient ",id,".RData")
  
  represults = list(#"age"=age[1:M],"women"=women, #Include for calibration purposes
    
    "risk_genetic"=as.matrix(prs),
    
    "alive_notrt"=as.matrix(alive.notrt),"alive10_notrt"=as.matrix(alive10.notrt),
    "risk_notrt_threshold"=risk.notrt.tr,
    "risk_notrt_genethreshold"=risk.notrt.gtr,
    "states_notrt"=states.notrt,"qalys_notrt"=qalys.notrt,"cost_notrt"=cost.notrt,
    "chd_nonfatal_notrt"=chd.nonfatal.notrt,"chd_fatal_notrt"=chd.fatal.notrt,"cost_chd_notrt"=cost.chd.notrt,
    "stroke_nonfatal_notrt"=stroke.nonfatal.notrt,"stroke_fatal_notrt"=stroke.fatal.notrt,"cost_stroke_notrt"=cost.stroke.notrt,
    
    "risk_threshold"=risk.tr,
    "states_threshold"=states.tr,
    "policy_threshold"=policy.tr,
    "qalys_threshold"=qalys.tr,"cost_threshold"=cost.tr,"cost_meds_threshold"=cost.meds.tr,
    "chd_nonfatal_threshold"=chd.nonfatal.tr,"chd_fatal_threshold"=chd.fatal.tr,"cost_chd_threshold"=cost.chd.tr,
    "stroke_nonfatal_threshold"=stroke.nonfatal.tr,"stroke_fatal_threshold"=stroke.fatal.tr,"cost_stroke_threshold"=cost.stroke.tr,
    
    "risk_genethreshold"=risk.gtr,
    "states_genethreshold"=states.gtr,
    "policy_genethreshold"=policy.gtr,
    "qalys_genethreshold"=qalys.gtr,"cost_genethreshold"=cost.gtr,"cost_meds_genethreshold"=cost.meds.gtr,"cost_gt_genethreshold"=cost.gt.gtr,
    "chd_nonfatal_genethreshold"=chd.nonfatal.gtr,"chd_fatal_genethreshold"=chd.fatal.gtr,"cost_chd_genethreshold"=cost.chd.gtr,
    "stroke_nonfatal_genethreshold"=stroke.nonfatal.gtr,"stroke_fatal_genethreshold"=stroke.fatal.gtr,"cost_stroke_genethreshold"=cost.stroke.gtr
  )
  
  if(age_subgroup!=0|gender_subgroup!=0){
    #Idicator for initial age>=50
    age_ind = ifelse(age[1]>=50,1,0)
    
    #Adding age and gender indicators to results
    represults = c("age"=age_ind,"women"=women,represults)
  }
  
  return(represults)
}

# Subsetting Patients According to Testing Range --------------------------

#Making sure patient's data is in data table format
ptdata = setDT(ptdata)
patient_risk = ptdata[!duplicated(ID),]

#10-year risk for ASCVD events calculation (CHD + stroke)
patient_risk[,risk := (arisk(1,sex,race,age,sbp,smk,tc,hdl,diab,0,time_tr)+arisk(2,sex,race,age,sbp,smk,tc,hdl,diab,0,time_tr)),by=ID]

#Identifying patients in different genetic testing ranges
SD1 = patient_risk[genetic_lb<=risk & risk<=genetic_ub,ID] #1 SD - base case
SD05 = patient_risk[genetic_lb05<=risk & risk<=genetic_ub05,ID] #0.5 SD - sensitivity analysis
SD2 = patient_risk[genetic_lb2<=risk & risk<=genetic_ub2,ID]  #2 SD - sensitivity analysis

setwd(home_dir)
save(SD05,SD1,SD2,file = "Index of Patients in Testing Ranges.RData",compress = F)
rm(patient_risk);gc()

#Subsetting data for base case
ptdata = ptdata[is.element(ID,SD1),]

# Calculating Risks -------------------------------------------------------

print(paste0("Simulating and Correcting PRS ",Sys.time()))
ptdata = riskcorrection(ptdata,or_sd)

#Running parallel risk simulations
print(paste0("Exporting Data for Risk Calculations ",Sys.time()))
cl = makeCluster(cores)
clusterExport(cl, varlist=c("ptdata",
                            "time_tp","time_tr",
                            "numhealth","N","events",
                            "arisk","zeta","zetachd","zetastroke","time",
                            "prs_incorrect","ascvd_rate",
                            "data.table","risksimulation"))
print(paste("Running Risk Simulations ",Sys.time(),sep = ""))
risk_results = parLapply(cl,unique(ptdata$ID),function(i){
  risksimulation(id=i)
})

#Extracting patient characteristics at the beginning of study
clusterExport(cl, varlist=c("risk_results"))
pt_sim0 = parLapply(cl,seq(risk_results),function(x){

  riskresults = risk_results[[x]]

  #selecting relevant fields
  riskresults[setdiff(names(riskresults),c("id","women","black","ldl","smoke","diab",
                                           "age","sbp","dbp","ps",
                                           "prs","periodrisk.tr","c.risk.tr"))]=NULL

  #Storing just first year
  riskresults$age = riskresults$age[1]
  riskresults$ldl = riskresults$ldl[1]
  riskresults$sbp = riskresults$sbp[1]
  riskresults$dbp = riskresults$dbp[1]
  riskresults$periodrisk.tr = sum(riskresults$periodrisk.tr[1,1,])
  riskresults$c.risk.tr = sum(riskresults$c.risk.tr[1,1,])

  return(riskresults)

})
stopCluster(cl)

#Converting list into data table
pt_sim0 = rbindlist(pt_sim0)

setwd(home_dir)
save(risk_results,file = "Risk Simulation Results.RData",compress = F)
fwrite(pt_sim0,paste0("Merged Risk Results for ",nrow(pt_sim0)," patients.csv"))

print(paste0("Risk Simulation Done ",Sys.time()))

# Calculating Transition Probabilities ------------------------------------

#Running parallel transition probabilities calculation
print(paste0("Exporting Data for Transition Probabilities ",Sys.time()))
cl = makeCluster(cores)
clusterExport(cl, varlist=c("risk_results","chddeathdata","strokedeathdata",
                            "alldeathdata","lifedata","adherence",
                            "numhealth","N","events","numalltrt",
                            "pce_correct","adjustment",
                            "TP","posttrtrisk","tpsimulation"))
print(paste0("Running Transition Probabilities Calculations ",Sys.time()))
tp_results = parLapply(cl,seq(risk_results),function(i){
  tpsimulation(id=i)
})
stopCluster(cl)
print(paste0("Transition Probabilities Calculations Done ",Sys.time()))
rm(risk_results);gc()

setwd(home_dir)
save(tp_results,file = "Transition Probability Results.RData",compress = F)

# Simulating Patients Health Trajectories ---------------------------------

#Loading transition probabilities
setwd(home_dir)
load("Transition Probability Results.RData")

# #Line for debugging purposes
# tp_results = tp_results[1:5]; reps = 3; cores = 2

#Running parallel health trajectory simulations
print(paste0("Exporting Data for Trajectory Simulations ",Sys.time()))
cl = makeCluster(cores)
clusterExport(cl, varlist=c("tp_results","reps","M","discount",
                            "test_year","genetic_lb","genetic_ub",
                            "mod_threshold","mod_genethreshold","arr_mod_genethreshold",
                            "numhealth","N","events","numalltrt",
                            "mortality_rates","mortality_mult","QoL","QoLterm",
                            "trtharm","dis_mult","adherence",
                            "adjustment","normalization_factor",
                            "disc","dis_res_prob",
                            "age_subgroup","gender_subgroup",
                            "costvector","genetic_cost",
                            "cost_chd","cost_chd_hist","cost_chd_fatal",
                            "cost_stroke","cost_stroke_hist","cost_stroke_fatal",
                            "posttrtrisk","meds2trt","ptsimulation"))

for(r in seq(reps)){
  print(paste0("Running Replication ",r," of Health Trajectory Simulations ",Sys.time()))
  clusterExport(cl, varlist=c("r"))
  htsimulation = parSapply(cl,seq(tp_results),function(i){
    ptsimulation(id = i)
  })

  #Formatting data into matrices
  ptresults = apply(htsimulation,1,function(x) do.call(rbind,x))
  rm(htsimulation);gc()

  #Counting ASCVD events
  ##CHD
  ptresults$chd_notrt = ptresults$chd_nonfatal_notrt+ptresults$chd_fatal_notrt
  ptresults$chd_threshold = ptresults$chd_nonfatal_threshold+ptresults$chd_fatal_threshold
  ptresults$chd_genethreshold = ptresults$chd_nonfatal_genethreshold+ptresults$chd_fatal_genethreshold

  ##Stroke
  ptresults$stroke_notrt = ptresults$stroke_nonfatal_notrt+ptresults$stroke_fatal_notrt
  ptresults$stroke_threshold = ptresults$stroke_nonfatal_threshold+ptresults$stroke_fatal_threshold
  ptresults$stroke_genethreshold = ptresults$stroke_nonfatal_genethreshold+ptresults$stroke_fatal_genethreshold

  ##Total
  ptresults$ascvd_notrt = ptresults$chd_notrt + ptresults$stroke_notrt
  ptresults$ascvd_threshold = ptresults$chd_threshold + ptresults$stroke_threshold
  ptresults$ascvd_genethreshold = ptresults$chd_genethreshold + ptresults$stroke_genethreshold

  #Counting treatment intensities
  ptresults$modtrt_threshold = ifelse(ptresults$policy_threshold==1,1,0)
  ptresults$hightrt_threshold = ifelse(ptresults$policy_threshold==2,1,0)

  ptresults$modtrt_genethreshold = ifelse(ptresults$policy_genethreshold==1,1,0)
  ptresults$hightrt_genethreshold = ifelse(ptresults$policy_genethreshold==2,1,0)

  #Saving merged data
  setwd(file.path(home_dir,dir_names[1]))
  save(ptresults,file =  paste0("Replication ", r," of trajectory simulation for ",length(tp_results)," patients.RData"),compress = F)
  # load(paste0("Replication ", r," of trajectory simulation for ",length(tp_results)," patients.RData"))

  #Adding over patients and time
  ptresults = lapply(ptresults,sum)

  #Merging replication results
  if (r==1){
    pt_sim = ptresults
    rm(ptresults);gc()
  }else if (r>1){
    pt_sim = clusterMap(cl,rbind, pt_sim, ptresults,SIMPLIFY = FALSE)
    rm(ptresults);gc()
  }
} #[r]
stopCluster(cl)

setwd(home_dir)
pt_sim = lapply(pt_sim,as.numeric)
save(pt_sim,file = paste0("Health Trajectory Simulation Results ",length(tp_results)," Patients Using ",reps," Replications.RData"),compress = F)
rm(tp_results);gc()

# Population Risk Scores Figures ------------------------------------------

#Loading data
setwd(home_dir)
load("Merged Risk Results for 93552 patients.RData")

# #Complete whole population data
# df = data.frame(PCE=pt_sim0_allpop$periodrisk.tr,GenePCE=pt_sim0_allpop$c.risk.tr)
# df$density = get_density(df$PCE, df$GenePCE)

#Reduced whole population data
digits = 5
index = !duplicated(round(pt_sim0_allpop$periodrisk.tr,digits))
pce = round(pt_sim0_allpop$periodrisk.tr,digits)[index]
genepce = round(pt_sim0_allpop$c.risk.tr,digits)[index]
df = data.frame("PCE"=pce,"GenePCE"=genepce)
df$density = get_density(df$PCE, df$GenePCE)


#Plotting whole population figure
setwd(file.path(home_dir,"Figures"))
jpeg("PopulationDensity.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(df) +
        geom_point(aes(PCE, GenePCE, color = density),size=point_size_scatter)+
        scale_color_gradientn(colors = cols,breaks = c(min(df$density),max(df$density)),
                              labels = c("Low","High"))+
        geom_text(x=0, y=0.8, label="A")+
        # labs(x = "PCE Risk",y = "GenePCE Risk",color="Population Size")+
        labs(x = "",y = "",color="")+
        guides(color = guide_colorbar(title.position="top",ticks=F))+
        scale_x_continuous(labels=percent,limits = c(0,0.8), breaks = seq(0,0.8,0.2)) +
        scale_y_continuous(labels=percent,limits = c(0,0.8), breaks = seq(0,0.8,0.2)) +
        geom_vline(xintercept = genetic_lb, linetype=5,
                   color = "green",size=line_size)+
        geom_vline(xintercept = genetic_ub, linetype=5,
                   color = "green",size=line_size)+
        geom_vline(xintercept = mod_threshold, linetype=2,
                   color = "blue",size=line_size)+
        geom_hline(yintercept=mod_genethreshold, linetype=2,
                   color = "blue",size=line_size)+
        # coord_fixed(ratio = 1) +
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "gray95"),
              panel.background = element_rect(fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="none", legend.box = "horizontal",
              panel.border = element_rect(colour = "black", fill=NA),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black")))
dev.off()

#Re-loading data
setwd(home_dir)
load("Merged Risk Results for 93552 patients.RData")

#Extracting information of patients with genetic test performed (and limiting range of GenePCE)
tested = pt_sim0_allpop[which(genetic_lb<=pt_sim0_allpop$periodrisk.tr&pt_sim0_allpop$periodrisk.tr<=genetic_ub),]
                                # &genetic_lb<=pt_sim0_allpop$c.risk.tr&pt_sim0_allpop$c.risk.tr<=genetic_ub),]

#Complete tested population data
df1 = data.frame(PCE=tested$periodrisk.tr,GenePCE=tested$c.risk.tr)
df1$density = get_density(df1$PCE, df1$GenePCE)

# #Reduced tested population data
# digits = 5
# index = !duplicated(round(tested$periodrisk.tr,digits))
# pce = round(tested$periodrisk.tr,digits)[index]
# genepce = round(tested$c.risk.tr,digits)[index]
# df1 = data.frame("PCE"=pce,"GenePCE"=genepce)
# df1$density = get_density(df1$PCE, df1$GenePCE)


#Plotting tested population figure
setwd(file.path(home_dir,"Figures"))
jpeg("TestedDensity.jpeg", width = width, height = height, units = 'in', res = res)

print(ggplot(df1) +
        geom_point(aes(PCE, GenePCE, color = density),size=point_size_scatter)+
        scale_color_gradientn(colors = cols,breaks = c(min(df1$density),max(df1$density)),labels = c("Low","High"))+
        geom_text(x=0.07, y=0.3, label="B")+
        labs(x = "",y = "",color="")+
        # labs(x = "PCE Risk",y = "GenePCE Risk",color="Density of Risk Score")+
        guides(color = guide_colorbar(title.position="top",ticks=F))+
        scale_x_continuous(labels=percent,limits = c(0.07,0.14),breaks = seq(0.08,0.14,0.02)) +
        scale_y_continuous(labels=percent,limits = c(0.04,0.3),breaks = seq(0,0.3,0.05)) +
        geom_vline(xintercept = mod_threshold, linetype=2,
                   color = "blue",size=line_size)+
        geom_hline(yintercept=mod_genethreshold, linetype=2,
                   color = "blue",size=line_size)+
        # coord_fixed(ratio = 1) +
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(color = "gray95"),
              panel.background = element_rect(fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="none", legend.box = "horizontal",
              panel.border = element_rect(colour = "black", fill=NA),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black")))
dev.off()

# Population Risk Over Time Figures ---------------------------------------

#Loading data
setwd(file.path(home_dir,dir_names[1]))
load("Replication 1 of trajectory simulation for 16174 patients.RData")

##PCE risk scores
#Renaming columns
risk_data = ptresults$risk_notrt_threshold
colnames(risk_data) = 1:ncol(risk_data)
risk_data = melt(risk_data)
names(risk_data) = c("obs","year","risk")
risk_data$year = factor(risk_data$year)

#Excluding deaths
risk_data = risk_data[-which(risk_data$risk==0),]

#Plotting densities
setwd(file.path(home_dir,"Figures"))
jpeg("PCEOverTime.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data = risk_data, aes(x = risk, y = rev(year),fill=rev(year))) +
        geom_density_ridges(scale = 0.9) +
        # geom_density_ridges(stat = "binline", bins = 20, scale = 0.9)+
        # geom_density_ridges(stat="density",scale = 0.9) +
        scale_y_discrete(labels = rev(levels(risk_data$year)))+
        scale_x_continuous(labels = percent,limits = c(0.05,0.7),breaks = seq(0.05,0.7,0.1))+
        labs(x = "PCE Risk",y = "Year of Study")+
        scale_fill_grey()+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.major.y = element_line(color = "gray95"),
              panel.background = element_rect(fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="none", legend.box = "horizontal",
              panel.border = element_rect(colour = "black", fill=NA),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black")))
dev.off()

##GenePCE risk scores
#Renaming columns
risk_data = ptresults$risk_notrt_genethreshold
colnames(risk_data) = 1:ncol(risk_data)
risk_data = melt(risk_data)
names(risk_data) = c("obs","year","risk")
risk_data$year = factor(risk_data$year)

#Excluding deaths
risk_data = risk_data[-which(risk_data$risk==0),]

#Plotting densities
setwd(file.path(home_dir,"Figures"))
jpeg("GenePCEOverTime.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data = risk_data, aes(x = risk, y = rev(year),fill=rev(year))) +
        geom_density_ridges(scale = 0.9) +
        scale_y_discrete(labels = rev(levels(risk_data$year)))+
        scale_x_continuous(labels = percent,limits = c(0.05,0.7),breaks = seq(0.05,0.7,0.1))+
        labs(x = "GenePCE Risk",y = "Year of Study")+
        scale_fill_grey()+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_blank(),
              panel.grid.major.y = element_line(color = "gray95"),
              panel.background = element_rect(fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="none", legend.box = "horizontal",
              panel.border = element_rect(colour = "black", fill=NA),
              axis.text.x=element_text(colour="black"),
              axis.text.y=element_text(colour="black")))
dev.off()


# Table 3 -----------------------------------------------------------------

#Table 3: Description of patients at first year of simulation

#Loading year 0 data
setwd(home_dir)
pt_sim0 = fread("Merged Risk Results at Year 0.csv")

pt_sim0_genetrt = list(
  "below_treatment" = pt_sim0[c.risk.tr<mod_threshold],
  "above_treatment" = pt_sim0[mod_threshold<=c.risk.tr]
)

pt_sim0_trt = list(
  "below_treatment" = pt_sim0[periodrisk.tr<mod_threshold],
  "above_treatment" = pt_sim0[mod_threshold<=periodrisk.tr]
)

#Merging survival results
setwd(file.path(home_dir,dir_names[1])) #setting working directory
tmp = list.files() #loading directory files

cl = makeCluster(cores)

clusterExport(cl, varlist=c("mod_threshold","genetic_lb","genetic_ub"))

tmp1 = parLapply(cl,tmp,function(x){load(x);

  #selecting relevant fields
  ptresults[setdiff(names(ptresults),c("risk_notrt_threshold","risk_notrt_genethreshold","alive10_notrt","alive_notrt"))]=NULL

  #Storing just first year
  ptresults$risk_notrt_threshold = ptresults$risk_notrt_threshold[,1]
  ptresults$risk_notrt_genethreshold = ptresults$risk_notrt_genethreshold[,1]
  ptresults$alive10_notrt = ptresults$alive10_notrt
  ptresults$alive_notrt = ptresults$alive_notrt

  ptresults})

#Separating list by replication number
tmp2 = parLapply(cl,tmp1,data.frame)

#Diving results according to risk groups
pt_sim0_survival_trt = list("below_treatment" = parSapply(cl,tmp2,function(x) colMeans(x[which(genetic_lb<=x$risk_notrt_threshold & x$risk_notrt_threshold<mod_threshold),
                                                                            c("alive_notrt","alive10_notrt")])),
                        "above_treatment" = parSapply(cl,tmp2,function(x) colMeans(x[which(mod_threshold<=x$risk_notrt_threshold & x$risk_notrt_threshold<=genetic_ub),
                                                                            c("alive_notrt","alive10_notrt")])))
stopCluster(cl);rm(tmp,tmp1,tmp2)

pt_sim0_survival_trt = lapply(pt_sim0_survival_trt,function(x) as.data.frame(t(x)))


#Limiting output to 4 digits
options(digits = 4)

#Creating table
t3 = list()
for(i in seq(pt_sim0_trt)){

  t3[[i]] = rbind(
    #Counts
    (multiplier/1e06)*length(pt_sim0_trt[[i]]$id),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$women==1)),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$black==1)),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$ldl==1)),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$smoke==1)),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$diab==1)),

    #Means
    mean(pt_sim0_trt[[i]]$age),
    mean(pt_sim0_trt[[i]]$sbp),
    mean(pt_sim0_trt[[i]]$dbp),
    # mean(pt_sim0_trt[[i]]$prs),
    mean(pt_sim0_survival_trt[[i]]$alive10_notrt),
    mean(pt_sim0_survival_trt[[i]]$alive_notrt),

    #Divider
    NA,

    #Percentages
    length(pt_sim0_trt[[i]]$id)/length(pt_sim0$id),
    length(which(pt_sim0_trt[[i]]$women==1))/length(pt_sim0_trt[[i]]$id),
    length(which(pt_sim0_trt[[i]]$black==1))/length(pt_sim0_trt[[i]]$id),
    length(which(pt_sim0_trt[[i]]$ldl==1))/length(pt_sim0_trt[[i]]$id),
    length(which(pt_sim0_trt[[i]]$smoke==1))/length(pt_sim0_trt[[i]]$id),
    length(which(pt_sim0_trt[[i]]$diab==1))/length(pt_sim0_trt[[i]]$id),

    #SD
    sd(pt_sim0_trt[[i]]$age),
    sd(pt_sim0_trt[[i]]$sbp),
    sd(pt_sim0_trt[[i]]$dbp),
    # sd(pt_sim0_trt[[i]]$prs),
    sd(pt_sim0_survival_trt[[i]]$alive10_notrt),
    sd(pt_sim0_survival_trt[[i]]$alive_notrt)
  )
}
table3 = do.call(cbind,t3)
noquote(formatC(table3,digits = 4,format = "f"))

# Supplemental Table 11 ---------------------------------------------------


#Supplemental Table 11: Reclassification of population with genetic testing using 10% threshold

#Loading year 0 data
setwd(home_dir)
pt_sim0 = fread("Merged Risk Results at Year 0.csv")

pt_sim0_genetrt = list(
  "below_treatment" = pt_sim0[c.risk.tr<mod_threshold],
  "above_treatment" = pt_sim0[mod_threshold<=c.risk.tr]
)

pt_sim0_trt = list(
  "below_treatment" = pt_sim0[periodrisk.tr<mod_threshold],
  "above_treatment" = pt_sim0[mod_threshold<=periodrisk.tr]
)

#Creating table
st11 = list()
for(i in seq(pt_sim0_trt)){

  st11[[i]] = cbind(
    #Counts
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$c.risk.tr<mod_threshold)),
    (multiplier/1e06)*length(which(pt_sim0_trt[[i]]$c.risk.tr>=mod_threshold)),

    #Percentages
    length(which(pt_sim0_trt[[i]]$c.risk.tr<mod_threshold))/length(pt_sim0$c.risk.tr),
    length(which(pt_sim0_trt[[i]]$c.risk.tr>=mod_threshold))/length(pt_sim0$c.risk.tr)
  )
}
sup_table11 = do.call(rbind,st11)
noquote(formatC(sup_table11,digits = 4,format = "f"))

# Table 4 -----------------------------------------------------------------

#Table 4: Health Outcomes for treatment policies

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")
load("Health Trajectory Simulation Results 16174 Patients Using 1000 Replications.RData")

###Remove!
pt_sim = lapply(pt_sim,function(x) head(x,1000))

#Creating table
parameters = c("mean","sd")
table4 = c()
for(j in seq(parameters)){

  if(j==1){#run for mean calculations (more robust)
    tmp = lapply(pt_sim,mean)
  }else{
    tmp = pt_sim
  }

  policies = unique(gsub(".*_","",names(pt_sim)))
  policies = policies[-match(c("notrt","genetic"),policies)]
  t4 = list()
  for(i in seq(policies)){
    t4[[i]] = rbind(
      (multiplier/1e06)*tmp[[paste0("modtrt_",policies[i])]],
      (multiplier/1e06)*tmp[[paste0("hightrt_",policies[i])]],

      (multiplier/1e00)*(tmp[[paste0("chd_","notrt")]]-tmp[[paste0("chd_",policies[i])]]),
      (multiplier/1e00)*(tmp[[paste0("stroke_","notrt")]]-tmp[[paste0("stroke_",policies[i])]]),

      100*((tmp[[paste0("ascvd_","notrt")]]-tmp[[paste0("ascvd_",policies[i])]])/length(SD1)),
      100*((tmp[[paste0("qalys_",policies[i])]]-tmp[[paste0("qalys_","notrt")]])/length(SD1)),

      (multiplier/1e06)*(tmp[[paste0("qalys_",policies[i])]]-tmp[[paste0("qalys_","notrt")]]),

      (multiplier/1e09)*(tmp[[paste0("cost_",policies[i])]]-tmp[[paste0("cost_","notrt")]]),

      (tmp[[paste0("cost_",policies[i])]]-tmp[[paste0("cost_","notrt")]])/
        (tmp[[paste0("qalys_",policies[i])]]-tmp[[paste0("qalys_","notrt")]]),

      (tmp[[paste0("cost_",policies[i])]]-tmp[[paste0("cost_","threshold")]])/
        (tmp[[paste0("qalys_",policies[i])]]-tmp[[paste0("qalys_","threshold")]])
    )
  }

  if(j==2){#run for standard deviation calculations (allowing variability)
    t4 = lapply(t4,function(x) apply(x,1,sd))
  }

  table4 = cbind(table4,do.call(cbind,t4))
} #[j]
print(noquote(formatC(table4,digits = 4,format = "f")))
rm(tmp)

# Table 5 -----------------------------------------------------------------

#Table 5: Comparison of incremetal gains

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")

#Merging simulation results
comparisons = c(NA,-1,0,1)
pt_sim_comp = list()

cl = makeCluster(cores)
for(c in seq(comparisons)){
  clusterExport(cl, varlist=c("c","comparisons","genetic_lb","genetic_ub","home_dir","dir_names","SD1"))

  tmp1 = parSapply(cl,seq(reps),function(r){

    setwd(file.path(home_dir,dir_names[1]))
    load(paste0("Replication ",r," of trajectory simulation for ",length(SD1)," patients.RData"))

    #Selecting agressiveness groups for summary (more/less agressive for at least 3 years)
    if(!is.na(comparisons[c])){
      if (comparisons[c]==-1){ #GeneACSVD less agressive
        genethreshold_la = which(rowSums(ifelse(ptresults$policy_genethreshold<ptresults$policy_threshold & ptresults$policy_threshold==1,1,0))>=3)
        ptresults = lapply(ptresults,function (x) x[genethreshold_la,])
      }else if (comparisons[c]==0){ #GeneACSVD equally agressive
        genethreshold_la = which(rowSums(ifelse(ptresults$policy_genethreshold<ptresults$policy_threshold & ptresults$policy_threshold==1,1,0))>=3)
        genethreshold_ma = which(rowSums(ifelse(ptresults$policy_genethreshold==1 & ptresults$policy_genethreshold>ptresults$policy_threshold,1,0))>=3)

        genethreshold_ea = which(rowSums(ifelse(ptresults$policy_genethreshold==1 & ptresults$policy_genethreshold==ptresults$policy_threshold,1,0))>=3)
        genethreshold_ea = genethreshold_ea[!is.element(genethreshold_ea,intersect(genethreshold_ea,genethreshold_la))
                                                &!is.element(genethreshold_ea,intersect(genethreshold_ea,genethreshold_ma))]
        ptresults = lapply(ptresults,function (x) x[genethreshold_ea,])
      }else if (comparisons[c]==1){ #GeneACSVD more agressive
        genethreshold_ma = which(rowSums(ifelse(ptresults$policy_genethreshold==1 & ptresults$policy_genethreshold>ptresults$policy_threshold,1,0))>=3)
        ptresults = lapply(ptresults,function (x) x[genethreshold_ma,])
      }
    }

    agv = comparisons[c] #Aggressiveness indicator
    number = nrow(ptresults$policy_threshold) #Total population after subset

    ptresults = lapply(ptresults,sum)
    ptresults = c("agressiveness"=agv,"number"=number,ptresults)

    return(ptresults)

    },USE.NAMES = F)

  tmp2 = parApply(cl,tmp1,1,as.numeric)

  tmp3 = data.frame(tmp2)

  pt_sim_comp[[c]] = tmp3
}
stopCluster(cl)
names(pt_sim_comp) = c("no_comparison","genepce_less_agressive","genepce_equally_agressive","genepce_more_agressive")

#Creating table
parameters = c("mean","sd")
table5 = c()
for(j in seq(parameters)){

  if(j==1){#run for mean claculations (more robust)
    pt_sim_comp_t = lapply(pt_sim_comp,function(x) as.data.frame(t(colMeans(x))))
  }else{
    pt_sim_comp_t = pt_sim_comp
  }

  t5 = list()
  for(i in seq(comparisons)){
    t5[[i]] = rbind(
      (multiplier/1e06)*pt_sim_comp_t[[i]]$number,

      (multiplier/1e00)*(pt_sim_comp_t[[i]]$chd_threshold-pt_sim_comp_t[[i]]$chd_genethreshold),
      (multiplier/1e00)*(pt_sim_comp_t[[i]]$stroke_threshold-pt_sim_comp_t[[i]]$stroke_genethreshold),

      100*(pt_sim_comp_t[[i]]$ascvd_threshold-pt_sim_comp_t[[i]]$ascvd_genethreshold)/(pt_sim_comp_t[[i]]$number),
      100*(pt_sim_comp_t[[i]]$qalys_genethreshold-pt_sim_comp_t[[i]]$qalys_threshold)/(pt_sim_comp_t[[i]]$number),

      (multiplier/1e06)*(pt_sim_comp_t[[i]]$qalys_genethreshold-pt_sim_comp_t[[i]]$qalys_threshold),

      (multiplier/1e09)*(pt_sim_comp_t[[i]]$cost_genethreshold-pt_sim_comp_t[[i]]$cost_threshold),

      (pt_sim_comp_t[[i]]$cost_genethreshold-pt_sim_comp_t[[i]]$cost_threshold)/
        (pt_sim_comp_t[[i]]$qalys_genethreshold-pt_sim_comp_t[[i]]$qalys_threshold)
    )
  }

  if(j==2){#run for standard deviation calculations (allowing variability)
    t5 = lapply(t5,function(x) apply(x,1,sd))
  }
  table5 = cbind(table5,do.call(cbind,t5))
} #[j]
print(noquote(formatC(table5,digits = 4,format = "f")))
rm(tmp1,tmp2,tmp3)

# Testing Range and Cost Sensitivity Analyses -----------------------------

#Sensitivity analysis parameters
gcost = c(0,50,100,200,400,800,2000)
testing_ranges = list(c(0.0624,0.1911),c(0.0768,0.1363),c(0.0871,0.1162),c(0.0768,0.1),c(0.0871,0.1))
modgenethreshold = c(0.0894,0.0847,0.0801) #previous 2SD treatment threshold: 0.09164

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")

#Establishing number of replications for sensitivity analysis
reps = 500

# #Line for debugging purposes
# SD1 = SD1[1:2]; SD2 = SD2[1:2]; SD05 = SD05[1:2]; reps = 3; cores = 2

dir_names[2] = "Genetic Cost Sensitivity Analysis Merged Replications Results"
if(!dir.exists(dir_names[2])){
  dir.create(dir_names[2])
}

dir_names[3] = "Genetic Cost Sensitivity Analysis Final Results"
if(!dir.exists(dir_names[3])){
  dir.create(dir_names[3])
}

#Parallel computing parameters
cl = makeCluster(cores)
##Parameters for risk calculations
clusterExport(cl, varlist=c("time_tp","time_tr",
                            "numhealth","N","events",
                            "arisk","zeta","zetachd","zetastroke","time",
                            "prs_incorrect","ascvd_rate",
                            "data.table","risksimulation"))

##Parameters for transition probability calculations
clusterExport(cl, varlist=c("chddeathdata","strokedeathdata",
                            "alldeathdata","lifedata","adherence",
                            "numhealth","N","events","numalltrt",
                            "pce_correct","adjustment",
                            "TP","posttrtrisk","tpsimulation"))

##Parameters for trajectory simulations
clusterExport(cl, varlist=c("reps","M","discount",
                            "test_year",
                            "mod_threshold","arr_mod_genethreshold",
                            "numhealth","N","events","numalltrt",
                            "mortality_rates","mortality_mult","QoL","QoLterm",
                            "trtharm","dis_mult","adherence",
                            "adjustment","normalization_factor",
                            "disc","dis_res_prob",
                            "age_subgroup","gender_subgroup",
                            "costvector","genetic_cost",
                            "cost_chd","cost_chd_hist","cost_chd_fatal",
                            "cost_stroke","cost_stroke_hist","cost_stroke_fatal",
                            "posttrtrisk","meds2trt","ptsimulation"))

for(tr in seq(testing_ranges)[1]){

  #Establishing genetic testing range
  genetic_lb = testing_ranges[[tr]][1]
  genetic_ub = testing_ranges[[tr]][2]

  print(paste0("Running range ",genetic_lb," to ",genetic_ub," ",Sys.time()))

  #Loading data
  setwd(home_dir)
  ptdata = fread("Continuous NHANES Synthetic Forecasted Dataset.csv") #2009-2014 continuous NHANES dataset

  #Subsetting data according to testing range (and adjusting GenePCE treatment threshold)
  if(genetic_lb==testing_ranges[[1]][1]){
    ptdata = ptdata[is.element(ID,SD2),]
    mod_genethreshold = modgenethreshold[1]
  }else if(genetic_lb==testing_ranges[[2]][1]){
    ptdata = ptdata[is.element(ID,SD1),]
    mod_genethreshold = modgenethreshold[2]
  }else if(genetic_lb==testing_ranges[[3]][1]){
    ptdata = ptdata[is.element(ID,SD05),]
    mod_genethreshold = modgenethreshold[3]
  }

  print(paste0("Simulating and Correcting PRS ",Sys.time()))
  ptdata = riskcorrection(ptdata,or_sd)

  print(paste("Running Risk Simulations ",Sys.time(),sep = ""))

  #Running parallel risk simulations
  clusterExport(cl, varlist=c("ptdata"))
  risk_results = parLapply(cl,unique(ptdata$ID),function(i){
    risksimulation(id=i)
  })

  print(paste0("Running Transition Probabilities Calculations ",Sys.time()))

  #Running parallel transition probabilities calculation
  clusterExport(cl, varlist=c("risk_results"))
  tp_results = parLapply(cl,seq(risk_results),function(i){
    tpsimulation(id=i)
  })
  rm(risk_results);gc()

  #Running parallel health trajectory simulations
  clusterExport(cl, varlist=c("tr","tp_results","mod_genethreshold","genetic_lb","genetic_ub"))
  for(r in seq(reps)){
    print(paste0("Running Replication ",r," of Health Trajectory Simulations ",Sys.time()))
    clusterExport(cl, varlist=c("r"))
    htsimulation = parSapply(cl,seq(tp_results),function(i){
      ptsimulation(id = i)
    })

    #Formatting data into matrices
    ptresults = apply(htsimulation,1,function(x) do.call(rbind,x))
    rm(htsimulation);gc()

    #Counting ASCVD events
    ##CHD
    ptresults$chd_notrt = ptresults$chd_nonfatal_notrt+ptresults$chd_fatal_notrt
    ptresults$chd_threshold = ptresults$chd_nonfatal_threshold+ptresults$chd_fatal_threshold
    ptresults$chd_genethreshold = ptresults$chd_nonfatal_genethreshold+ptresults$chd_fatal_genethreshold

    ##Stroke
    ptresults$stroke_notrt = ptresults$stroke_nonfatal_notrt+ptresults$stroke_fatal_notrt
    ptresults$stroke_threshold = ptresults$stroke_nonfatal_threshold+ptresults$stroke_fatal_threshold
    ptresults$stroke_genethreshold = ptresults$stroke_nonfatal_genethreshold+ptresults$stroke_fatal_genethreshold

    ##Total
    ptresults$ascvd_notrt = ptresults$chd_notrt + ptresults$stroke_notrt
    ptresults$ascvd_threshold = ptresults$chd_threshold + ptresults$stroke_threshold
    ptresults$ascvd_genethreshold = ptresults$chd_genethreshold + ptresults$stroke_genethreshold

    #Counting treatment intensities
    ptresults$modtrt_threshold = ifelse(ptresults$policy_threshold==1,1,0)
    ptresults$hightrt_threshold = ifelse(ptresults$policy_threshold==2,1,0)

    ptresults$modtrt_genethreshold = ifelse(ptresults$policy_genethreshold==1,1,0)
    ptresults$hightrt_genethreshold = ifelse(ptresults$policy_genethreshold==2,1,0)

    #Saving merged data
    setwd(file.path(home_dir,dir_names[2]))
    save(ptresults,file =  paste0("Replication ", r," of trajectory simulation for ",length(tp_results)," patients.rdata"),compress = F)

  } #[r]

  print(paste0("Genetic cost sensitivity analysis ",Sys.time()))

  #Merged simulation results
  setwd(file.path(home_dir,dir_names[2]))
  tmp = list.files() #Loading directory files

  #Genetic cost sensitivity analysis
  for(c in seq(gcost)){

    setwd(file.path(home_dir,dir_names[2]))
    clusterExport(cl, varlist=c("c","gcost","discount","multiplier","home_dir","dir_names"))

    tmp1 = parSapply(cl,tmp,function(x){

      setwd(file.path(home_dir,dir_names[2]))
      load(x)

      #Updating genetic testing cost
      ptresults$cost_gt_genethreshold = ifelse(ptresults$cost_gt_genethreshold>0,gcost[c],0)
      ptresults$cost_gt_genethreshold = sapply(seq(ncol(ptresults$cost_gt_genethreshold)),function(x) ptresults$cost_gt_genethreshold[,x]*discount^x) #discounting costs

      #Updating total cost
      ptresults$cost_genethreshold = ptresults$cost_gt_genethreshold+ptresults$cost_meds_genethreshold+
        ptresults$cost_chd_genethreshold+ptresults$cost_stroke_genethreshold

      #Total number genetic tests performed
      number = (multiplier/1e06)*length(which(genetic_lb<=ptresults$risk_notrt_threshold[,1]&ptresults$risk_notrt_threshold[,1]<=genetic_ub))

      #Adding over patients and time
      ptresults = lapply(ptresults,sum)

      #Including number of patients tested
      ptresults = c("number"=number,ptresults)

      return(ptresults)

      },USE.NAMES = F)
    tmp2 = parApply(cl,tmp1,1,cbind)
    pt_sim = parLapply(cl,tmp2,as.numeric)
    rm(tmp1,tmp2)

    setwd(file.path(home_dir,dir_names[3]))
    save(pt_sim,file = paste0("Results for ",length(tp_results)," patients using ",reps," replications - $",gcost[c],
                              " as testing cost ",genetic_lb," to ",genetic_ub," testing range.RData"))

  } #[c]

  #Deleting temporary files
  setwd(file.path(home_dir,dir_names[2]))
  file.remove(tmp)

} #[tr]
stopCluster(cl)

#Returning genetic testing range to original risk range
genetic_lb = testing_ranges[[2]][1]
genetic_ub = testing_ranges[[2]][2]
genetic_cost = gcost[4]
mod_genethreshold = modgenethreshold[2]

# Testing Range and Genetic Cost Figures ----------------------------------

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")

#Setting working directory
setwd(file.path(home_dir,"Genetic Cost Sensitivity Analysis Final Results"))
testing_ranges = list(c(0.0624,0.1911),c(0.0768,0.1363),c(0.0871,0.1162),c(0.0768,0.1),c(0.0871,0.1))
gcost = c(0,50,100,200,400,800,2000)

#Genetic cost sensitivity analysis
gc_sens = sapply(testing_ranges,function(z){
  temp = sapply(gcost,function (y){

    #Subsetting data according to testing range
    if(z[1]==testing_ranges[[1]][1]){
      patients = length(SD2)
    }else if(z[1]==testing_ranges[[2]][1]){
      patients = length(SD1)
    }else if(z[1]==testing_ranges[[3]][1]){
      patients = length(SD05)
    }

    #Loading treatment information for patient i
    load(paste0("Results for ",patients," patients using ",reps," replications - $",y," as testing cost ",z[1]," to ",z[2]," testing range.RData"))

    #Calsulating averages
    pt_sim = lapply(pt_sim,mean)


    pt_sim$ascvd_averted_threshold = (pt_sim$ascvd_notrt-pt_sim$ascvd_threshold) #ascvd events averted under PCE policy
    pt_sim$ascvd_averted_genethreshold = (pt_sim$ascvd_notrt-pt_sim$ascvd_genethreshold) #ascvd events averted under GenePCE policy
    pt_sim$marginal_ascvd_averted_genethreshold = (pt_sim$ascvd_threshold-pt_sim$ascvd_genethreshold) #marginal ascvd events averted under GenePCE policy

    pt_sim$qalys_saved_threshold = (pt_sim$qalys_threshold-pt_sim$qalys_notrt) #qalys saved under PCE policy
    pt_sim$qalys_saved_genethreshold = (pt_sim$qalys_genethreshold-pt_sim$qalys_notrt) #qalys saved under GenePCE policy
    pt_sim$marginal_qalys_saved_genethreshold = (pt_sim$qalys_genethreshold-pt_sim$qalys_threshold) #marginal qalys saved under GenePCE policy

    pt_sim$cost_incurred_threshold = (pt_sim$cost_threshold-pt_sim$cost_notrt) #cost under PCE policy
    pt_sim$cost_incurred_genethreshold = (pt_sim$cost_genethreshold-pt_sim$cost_notrt) #cost under GenePCE policy
    pt_sim$marginal_cost_incurred_genethreshold = (pt_sim$cost_genethreshold-pt_sim$cost_threshold) #marginal cost under GenePCE policy

    pt_sim$cost_qaly_threshold = pt_sim$cost_incurred_threshold/pt_sim$qalys_saved_threshold #cost/qaly saved under PCE policy
    pt_sim$cost_qaly_genethreshold = pt_sim$cost_incurred_genethreshold/pt_sim$qalys_saved_genethreshold #cost/qaly saved under GenePCE policy
    pt_sim$marginal_cost_qaly_genethreshold = pt_sim$marginal_cost_incurred_genethreshold/pt_sim$marginal_qalys_saved_genethreshold #marginal cost/qaly saved under GenePCE policy

    return(pt_sim)
    },USE.NAMES = FALSE)
  temp = apply(temp,1,cbind)
  temp = lapply(temp,as.numeric)

})
gc_sens = apply(gc_sens,1,cbind)
gc_sens = lapply(gc_sens,function(z) t(apply(z,1,unlist)))
for(i in seq(gc_sens)){
  rownames(gc_sens[[i]]) = c("2 SD","1 SD","0.5 SD","Lower 1 SD","Lower 0.5 SD")
  colnames(gc_sens[[i]]) = gcost
}

#Plotting figure cost/qaly saved
df1 = gc_sens$marginal_cost_qaly_genethreshold
rownames(df1) = c("2 SD","1 SD","0.5 SD","Lower 1 SD","Lower 0.5 SD")

df1 = melt(df1)
colnames(df1) = c("risk_range","labels","data")

setwd(file.path(home_dir,"Figures"))
jpeg("TestingRangeandCost.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data=df1, aes(x=labels, y=data, group=risk_range,shape=risk_range)) +
        geom_line(aes(linetype=risk_range,color=risk_range),size=line_size)+
        geom_point(aes(color=risk_range),size=point_size_line)+
        scale_shape_manual(values=c(19,15,17,4,3))+
        scale_x_continuous(labels=dollar,limits = c(0,800),breaks = seq(0,800,100)) +
        scale_y_continuous(labels=dollar,limits = c(-10e03,350e03),breaks = seq(0,350e03,50e03)) +
        labs(x = "Cost of Genetic Testing (2018 USD)",y = "ICER of GenePCE\n (2018 USD/QALY Saved)",
             linetype="Genetic Testing\nScenario",color="Genetic Testing\nScenario",shape="Genetic Testing\nScenario")+
        guides(linetype = guide_legend(title.position="top", title.hjust = 0.5),
               color = guide_legend(title.position="top", title.hjust = 0.5),
               shape = guide_legend(title.position="top", title.hjust = 0.5))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "gray95"),
              panel.background = element_rect(fill=NA),
              panel.border = element_rect(colour = "black", fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="bottom", legend.box = "horizontal",
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")

        ))
dev.off()

#Plotting cost incurred
df1 = gc_sens$marginal_cost_incurred_genethreshold*multiplier/1e09
rownames(df1) = c("2 SD","1 SD","0.5 SD","Lower 1 SD","Lower 0.5 SD")
df1 = melt(df1)
colnames(df1) = c("risk_range","labels","data")

setwd(file.path(home_dir,"Figures"))
jpeg("CostIncurred.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data=df1, aes(x=labels, y=data, group=risk_range,shape=risk_range)) +
        geom_line(aes(linetype=risk_range,color=risk_range),size=line_size)+
        geom_point(aes(color=risk_range),size=point_size_line)+
        scale_shape_manual(values=c(19,15,17,3,4))+
        labs(x = "Cost of Genetic Testing (2018 USD)",y = "Cost of the GenePCE Strategy\n (2018 Billion USD)")+
        scale_x_continuous(labels=dollar,breaks = seq(0,2000,500)) +
        scale_y_continuous(labels=dollar,limits = c(-5,60),breaks = seq(0,60,20)) +
        labs(linetype="Genetic Testing\nRange",color="Genetic Testing\nRange",shape="Genetic Testing\nRange")+
        # labs(linetype="",color="",shape="")+
        labs(linetype="Genetic Testing Range",color="Genetic Testing Range",shape="Genetic Testing Range")+
        guides(linetype = guide_legend(title.position="top", title.hjust = 0.5),
               color = guide_legend(title.position="top", title.hjust = 0.5),
               shape = guide_legend(title.position="top", title.hjust = 0.5))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "gray95"),
              panel.background = element_rect(fill=NA),
              panel.border = element_rect(colour = "black", fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="bottom", legend.box = "horizontal",
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")

        ))
dev.off()

#Plotting ascvd events averted
df1 = cbind("GenePCE"=(gc_sens$ascvd_averted_genethreshold[,1]-gc_sens$ascvd_averted_threshold[,1])*multiplier,
            "PCE"=gc_sens$ascvd_averted_threshold[,1]*multiplier)
df1 = melt(df1)
names(df1)=c("risk_range","strategy","value")

setwd(file.path(home_dir,"Figures"))
jpeg("ascvdtestingranges.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data=df1, aes(x=risk_range, y=value,fill=strategy)) +
        geom_bar(stat="identity",width=0.5)+
        labs(x = "Genetic Testing Range",y = "ASCVD Events Prevented",fill="Treatment Strategy")+
        scale_y_continuous(labels=comma) +
        guides(fill = guide_legend(title.position="top", title.hjust = 0.5))+
        scale_fill_grey(start = 0.5)+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "gray95"),
              panel.background = element_rect(fill=NA),
              panel.border = element_rect(colour = "black", fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="bottom", legend.box = "horizontal",
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")

        ))
dev.off()

#Plotting qalys saved
df1 = cbind("GenePCE"=(gc_sens$qalys_saved_genethreshold[,1]-gc_sens$qalys_saved_threshold[,1])*multiplier/1e06,
            "PCE"=gc_sens$qalys_saved_threshold[,1]*multiplier/1e06)
df1 = melt(df1)
names(df1)=c("risk_range","strategy","value")

setwd(file.path(home_dir,"Figures"))
jpeg("qalystestingranges.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data=df1, aes(x=risk_range, y=value,fill=strategy)) +
        geom_bar(stat="identity",width=0.5)+
        labs(x = "Genetic Testing Range",y = "QALYs Saved (in millions)",fill="Treatment Strategy")+
        scale_y_continuous(labels=comma,limits = c(0,4),breaks = seq(0,4,1)) +
        # scale_y_continuous(labels=comma,limits = c(2,2.15),breaks = seq(2,2.15,0.05),oob = rescale_none) +
        scale_fill_grey(start = 0.5)+
        guides(fill = guide_legend(title.position="top", title.hjust = 0.5))+
        theme_bw()+
        theme(panel.grid.minor = element_blank(),panel.grid.major = element_line(colour = "gray95"),
              panel.background = element_rect(fill=NA),
              panel.border = element_rect(colour = "black", fill=NA),
              plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),
              legend.position="bottom", legend.box = "horizontal",
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")

        ))
dev.off()

# Other Sensitivity Analyses ----------------------------------------------

#Loading index of patients in testing ranges and data
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")
ptdata = fread("Continuous NHANES Synthetic Forecasted Dataset.csv")
ptdata = ptdata[is.element(ID,SD1),] #subsetting for 1 SD testing range

#Establishing number of replications for sensitivity analysis
reps = 500

#Generating directory to store results
dir_names[4] = "Other Sensitivity Analysis Merged Replications Results"
if(!dir.exists(dir_names[4])){
  dir.create(dir_names[4])
}

dir_names[5] = "Other Sensitivity Analysis Final Results"
if(!dir.exists(dir_names[5])){
  dir.create(dir_names[5])
}

#Sensitivity analysis parameters
sens_list = list(list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(1,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,1,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,0.09422,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,0.1,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,1.38,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,2.00,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,c(0,50,50),ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector/2,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,-0.3,trtharm[3],dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,0.004,dis_mult,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],2,disc,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,1,mortality_mult,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,0.75,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,0.5,prs_incorrect),
                 list(pce_correct,adjustment,mod_genethreshold,or_sd,costvector,ascvd_rate,trtharm[3],dis_mult,disc,mortality_mult,1)
                 )

# #Line for debugging purposes
# ptdata = ptdata[is.element(ID,SD1[1:2]),]; reps = 3; cores = 2; sens_list = sens_list[1:2]

#Parallel computing parameters
cl = makeCluster(cores)
##Parameters for risk calculations
clusterExport(cl, varlist=c("time_tp","time_tr",
                            "numhealth","N","events",
                            "arisk","zeta","zetachd","zetastroke","time",
                            "prs_incorrect","ascvd_rate",
                            "data.table","risksimulation"))

##Parameters for transition probability calculations
clusterExport(cl, varlist=c("chddeathdata","strokedeathdata",
                            "alldeathdata","lifedata","adherence",
                            "numhealth","N","events","numalltrt",
                            "pce_correct","adjustment",
                            "TP","posttrtrisk","tpsimulation"))

##Parameters for trajectory simulations
clusterExport(cl, varlist=c("reps","M","discount",
                            "test_year","genetic_lb","genetic_ub",
                            "mod_threshold","mod_genethreshold","arr_mod_genethreshold",
                            "numhealth","N","events","numalltrt",
                            "mortality_rates","mortality_mult","QoL","QoLterm",
                            "trtharm","dis_mult","adherence",
                            "adjustment","normalization_factor",
                            "disc","dis_res_prob",
                            "age_subgroup","gender_subgroup",
                            "costvector","genetic_cost",
                            "cost_chd","cost_chd_hist","cost_chd_fatal",
                            "cost_stroke","cost_stroke_hist","cost_stroke_fatal",
                            "posttrtrisk","meds2trt","ptsimulation"))

for(sa in seq(sens_list)){

  print(paste0("Running scenario ",sa," ",Sys.time()))

  pce_correct = sens_list[[sa]][[1]] #indicator for PCE being the true risk (0=GenePCE is the true risk;1=PCE is the true risk)
  adjustment = sens_list[[sa]][[2]] #indicator of adjustment in RRR using normalization factors based on the PRS (1=run sentivity analysis; 0=base case)
  mod_genethreshold = sens_list[[sa]][[3]]
  or_sd = sens_list[[sa]][[4]] #OR increase per SD of polygenic score
  costvector = sens_list[[sa]][[5]] #costs of no treatment, moderate intensity, and high intensity statin ($/year)
  ascvd_rate = sens_list[[sa]][[6]] #multiplier for sensitivity analysis in CVD rates
  trtharm[3] = sens_list[[sa]][[7]] #high intensity statin treatment disutilty
  dis_mult = sens_list[[sa]][[8]] #multiplier of age for sensitivity analysis in disutility
  disc = sens_list[[sa]][[9]] #indicator for dicontinuation and restart rates (1=run sentivity analysis; 0=base case)
  mortality_mult = sens_list[[sa]][[10]] #multiplier for sensitivity analysis on mortality factors
  prs_incorrect = sens_list[[sa]][[11]] #indicator for sensitivity analysis on making decisions according to Khera's prs while Abraham's prs was true


  #Adjusting GenePCE treatment threshold according to scenario
  if(or_sd==1.38){
    mod_genethreshold = 0.09302
  }else if(or_sd==2){
    mod_genethreshold = 0.07677
  }else if(ascvd_rate==-0.3){
    mod_genethreshold = 0.09075
  }else{
    mod_genethreshold = sens_list[[sa]][[3]]
  }

  #Exporting sensitivity analysis parameters
  clusterExport(cl, varlist=c("pce_correct","adjustment","mod_genethreshold",
                              "or_sd","costvector","prs_incorrect","ascvd_rate",
                              "trtharm","dis_mult","disc","mortality_mult"))

  print(paste0("Simulating and Correcting PRS ",Sys.time()))
  ptdata = riskcorrection(ptdata,or_sd)
  
  if(prs_incorrect==1){ #Run only for "wrong PRS" scenario
    names(ptdata)[match("prs",names(ptdata))] = "prs_incorrect"
    ptdata = riskcorrection(ptdata,1.38)
  }

  print(paste("Running Risk Simulations ",Sys.time(),sep = ""))

  #Running parallel risk simulations
  clusterExport(cl, varlist=c("ptdata"))
  risk_results = parLapply(cl,unique(ptdata$ID),function(i){
    risksimulation(id=i)
  })

  print(paste0("Running Transition Probabilities Calculations ",Sys.time()))

  #Running parallel transition probabilities calculation
  clusterExport(cl, varlist=c("risk_results"))
  tp_results = parLapply(cl,seq(risk_results),function(i){
    tpsimulation(id=i)
  })
  rm(risk_results);gc()

  print(paste0("Running Health Trajectory Simulations ",Sys.time()))

  #Running parallel health trajectory simulations
  clusterExport(cl, varlist=c("tp_results"))
  for(r in seq(reps)){
    clusterExport(cl, varlist=c("r"))
    htsimulation = parSapply(cl,seq(tp_results),function(i){
      ptsimulation(id = i)
    })

    #Formatting data into matrices
    ptresults = apply(htsimulation,1,function(x) do.call(rbind,x))
    rm(htsimulation);gc()

    #Counting ASCVD events
    ##CHD
    ptresults$chd_notrt = ptresults$chd_nonfatal_notrt+ptresults$chd_fatal_notrt
    ptresults$chd_threshold = ptresults$chd_nonfatal_threshold+ptresults$chd_fatal_threshold
    ptresults$chd_genethreshold = ptresults$chd_nonfatal_genethreshold+ptresults$chd_fatal_genethreshold

    ##Stroke
    ptresults$stroke_notrt = ptresults$stroke_nonfatal_notrt+ptresults$stroke_fatal_notrt
    ptresults$stroke_threshold = ptresults$stroke_nonfatal_threshold+ptresults$stroke_fatal_threshold
    ptresults$stroke_genethreshold = ptresults$stroke_nonfatal_genethreshold+ptresults$stroke_fatal_genethreshold

    ##Total
    ptresults$ascvd_notrt = ptresults$chd_notrt + ptresults$stroke_notrt
    ptresults$ascvd_threshold = ptresults$chd_threshold + ptresults$stroke_threshold
    ptresults$ascvd_genethreshold = ptresults$chd_genethreshold + ptresults$stroke_genethreshold

    #Counting treatment intensities
    ptresults$modtrt_threshold = ifelse(ptresults$policy_threshold==1,1,0)
    ptresults$hightrt_threshold = ifelse(ptresults$policy_threshold==2,1,0)

    ptresults$modtrt_genethreshold = ifelse(ptresults$policy_genethreshold==1,1,0)
    ptresults$hightrt_genethreshold = ifelse(ptresults$policy_genethreshold==2,1,0)

    #Saving merged data
    setwd(file.path(home_dir,dir_names[4]))
    save(ptresults,file =  paste0("Replication ", r," of trajectory simulation for ",length(tp_results)," patients.RData"),compress = F)

    #Adding over patients and time
    ptresults = parLapply(cl,ptresults,sum)

    #Merging replication results
    if (r==1){
      pt_sim = ptresults
      rm(ptresults);gc()
    }else if (r>1){
      pt_sim = clusterMap(cl,rbind, pt_sim, ptresults,SIMPLIFY = FALSE)
      rm(ptresults);gc()
    }
  } #[r]

  print(paste0("Adding over Time and Patients ",Sys.time()))

  setwd(file.path(home_dir,dir_names[4]))
  tmp = list.files()

  clusterExport(cl, varlist=c("home_dir","dir_names"))
  tmp1 = parSapply(cl,tmp,function(x){

    #Loading data
    setwd(file.path(home_dir,dir_names[4]))
    load(x)

    #Adding over patients and time
    ptresults = lapply(ptresults,sum)

    return(ptresults)
    },USE.NAMES = F)
  tmp2 = parApply(cl,tmp1,1,cbind)
  pt_sim = parLapply(cl,tmp2,as.numeric)

  #Deleting temporary files
  setwd(file.path(home_dir,dir_names[4]))
  file.remove(tmp)
  rm(tmp,tmp1,tmp2)

  setwd(file.path(home_dir,dir_names[5]))
  save(pt_sim,file = paste0("Results for ",length(tp_results)," patients using ",reps," replications - Scenario ",sa,".RData"))
  rm(tp_results,pt_sim);gc()

} #[sa]
stopCluster(cl)

#Returning values to base case
pce_correct = sens_list[[1]][[1]] #Indicator for PCE being the true risk (0=GenePCE is the true risk;1=PCE is the true risk)
adjustment = sens_list[[1]][[2]] #Indicator of adjustment in RRR using normalization factors based on the PRS (1=run sentivity analysis; 0=base case)
mod_genethreshold = sens_list[[1]][[3]]
or_sd = sens_list[[1]][[4]] #OR increase per SD of polygenic score
costvector = sens_list[[1]][[5]] #costs of no treatment, moderate intensity, and high intensity statin ($/year)
ascvd_rate = sens_list[[1]][[6]] #multiplier for sensitivity analysis in CVD rates
trtharm[3] = sens_list[[1]][[7]] #high intensity statin treatment disutilty
dis_mult = sens_list[[1]][[8]] #multiplier of age for sensitivity analysis in disutility
disc = sens_list[[1]][[9]] #Indicator for dicontinuation and restart rates (1=run sentivity analysis; 0=base case)
mortality_mult = sens_list[[1]][[10]] #multiplier for sensitivity analysis on mortality factors
prs_incorrect = sens_list[[1]][[11]] #indicator for sensitivity analysis on making decisions according to Khera's prs while Abraham's prs was true (1=run sentivity analysis; 0=base case)

# Table 6 -----------------------------------------------------------------

#Table 6: Sensitivity analyses

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")

setwd(file.path(home_dir,"Other Sensitivity Analysis Final Results"))
tmp = list.files() #Loading directory files

#Creating table
parameters = c("mean","sd")
tmp3 = list()
table6 = c()
for(j in seq(parameters)[1]){
  tmp2 = c()
  for(i in seq(tmp)){

    load(paste0("Results for ",length(SD1)," patients using ",reps," replications - Scenario ",i,".RData"))

    if(j==1){#run for mean calculations (more robust)
      pt_sim = lapply(pt_sim,mean)
    }

    tmp1 = list("ascvd_events"=(pt_sim$ascvd_threshold-pt_sim$ascvd_genethreshold)*multiplier,
              "cost_qaly"=(pt_sim$cost_genethreshold-pt_sim$cost_threshold)/
                (pt_sim$qalys_genethreshold-pt_sim$qalys_threshold),
              "modtrt_threshold"=pt_sim$modtrt_threshold*multiplier/1e06,
              "modtrt_genethreshold"=pt_sim$modtrt_genethreshold*multiplier/1e06,
              "hightrt_threshold"=pt_sim$hightrt_threshold*multiplier/1e06,
              "hightrt_genethreshold"=pt_sim$hightrt_genethreshold*multiplier/1e06)

    if(j==2){#run for standard deviation calculations (allowing variability)
      tmp1 = lapply(tmp1,sd)
    }

    tmp2 = rbind(tmp2,do.call(cbind,tmp1))
  } #[i]
  tmp3[[j]] = tmp2


} #[j]
table6 = do.call(cbind,tmp3)
print(noquote(formatC(table6,digits = 4,format = "f")))
rm(tmp,tmp1,tmp2,tmp3)

#Verfication of results
setwd(file.path(home_dir,"Other Sensitivity Analysis Final Results"))
tmp = list.files() #Loading directory files
for(i in seq(tmp)){

  load(paste0("Results for ",length(SD1)," patients using ",reps," replications - Scenario ",i,".RData"))
  cat("Scenario",i,"\n")
  pt_sim = lapply(pt_sim,mean)
  pt_sim$qalys_genethreshold*multiplier/1e06
  pt_sim$qalys_threshold*multiplier/1e06
  pt_sim$cost_genethreshold*multiplier/1e09
  pt_sim$cost_threshold*multiplier/1e09
  print(c((pt_sim$qalys_genethreshold*multiplier/1e06-pt_sim$qalys_threshold*multiplier/1e06),
  (pt_sim$cost_genethreshold*multiplier/1e09-pt_sim$cost_threshold*multiplier/1e09)))

} #[i]


# Subgroup Analysis -------------------------------------------------------

#Loading transition probabilities
setwd(home_dir)
load("Transition Probability Results.RData")

#Establishing number of replications for sensitivity analysis
reps = 500

dir_names[7] = "Merged Subgroup Replication Results"
if(!dir.exists(dir_names[7])){
  dir.create(dir_names[7])
}

dir_names[8] = "Subgroup Analysis Final Results"
if(!dir.exists(dir_names[8])){
  dir.create(dir_names[8])
}

#Subgroup analysis parameters
sg_list = list(c(age_subgroup,gender_subgroup), #base case
               c(1,gender_subgroup), #age <50
               c(2,gender_subgroup), #age >=50
               c(age_subgroup,1), #males
               c(age_subgroup,2)) #females

# #Line for debugging purposes
# tp_results = tp_results[1:5]; reps = 3; cores = 2; sg_list = sg_list[1:2]

#Parallel computing parameters
cl = makeCluster(cores)
clusterExport(cl, varlist=c("tp_results","reps","M","discount",
                            "test_year","genetic_lb","genetic_ub",
                            "mod_threshold","mod_genethreshold","arr_mod_genethreshold",
                            "numhealth","N","events","numalltrt",
                            "mortality_rates","QoL","QoLterm",
                            "trtharm","dis_mult","adherence",
                            "adjustment","normalization_factor",
                            "disc","dis_res_prob",
                            "age_subgroup","gender_subgroup",
                            "costvector","genetic_cost",
                            "cost_chd","cost_chd_hist","cost_chd_fatal",
                            "cost_stroke","cost_stroke_hist","cost_stroke_fatal",
                            "posttrtrisk","meds2trt","ptsimulation"))

for(sg in seq(sg_list)){

  print(paste0("Running scenario ",sg," ",Sys.time()))

  age_subgroup = sg_list[[sg]][1] #Indicator for age subgroup analysis (0=base case;1=age<50;2=age>=50)
  gender_subgroup = sg_list[[sg]][2] #Indicator for gender subgroup analysis (0=base case;1=males;2=females)

  #Running parallel health trajectory simulations
  for(r in seq(reps)){
    print(paste0("Running Replication ",r," of Health Trajectory Simulations ",Sys.time()))
    clusterExport(cl, varlist=c("sg","r","age_subgroup","gender_subgroup"))
    htsimulation = parSapply(cl,seq(tp_results),function(i){
      ptsimulation(id = i)
    })

    #Formatting data into matrices
    ptresults = apply(htsimulation,1,function(x) do.call(rbind,x))
    rm(htsimulation);gc()

    #Counting ASCVD events
    ##CHD
    ptresults$chd_notrt = ptresults$chd_nonfatal_notrt+ptresults$chd_fatal_notrt
    ptresults$chd_threshold = ptresults$chd_nonfatal_threshold+ptresults$chd_fatal_threshold
    ptresults$chd_genethreshold = ptresults$chd_nonfatal_genethreshold+ptresults$chd_fatal_genethreshold

    ##Stroke
    ptresults$stroke_notrt = ptresults$stroke_nonfatal_notrt+ptresults$stroke_fatal_notrt
    ptresults$stroke_threshold = ptresults$stroke_nonfatal_threshold+ptresults$stroke_fatal_threshold
    ptresults$stroke_genethreshold = ptresults$stroke_nonfatal_genethreshold+ptresults$stroke_fatal_genethreshold

    ##Total
    ptresults$ascvd_notrt = ptresults$chd_notrt + ptresults$stroke_notrt
    ptresults$ascvd_threshold = ptresults$chd_threshold + ptresults$stroke_threshold
    ptresults$ascvd_genethreshold = ptresults$chd_genethreshold + ptresults$stroke_genethreshold

    #Counting treatment intensities
    ptresults$modtrt_threshold = ifelse(ptresults$policy_threshold==1,1,0)
    ptresults$hightrt_threshold = ifelse(ptresults$policy_threshold==2,1,0)

    ptresults$modtrt_genethreshold = ifelse(ptresults$policy_genethreshold==1,1,0)
    ptresults$hightrt_genethreshold = ifelse(ptresults$policy_genethreshold==2,1,0)

    #Saving merged data
    setwd(file.path(home_dir,dir_names[7]))
    save(ptresults,file =  paste0("Replication ", r," of trajectory simulation for ",length(tp_results)," patients.rdata"),compress = F)

  } #[r]


  print(paste0("Adding over Time and Patients ",Sys.time()))

  setwd(file.path(home_dir,dir_names[7]))
  tmp = list.files()

  clusterExport(cl, varlist=c("home_dir","dir_names"))
  tmp1 = parSapply(cl,tmp,function(x){

    #Loading data
    setwd(file.path(home_dir,dir_names[7]))
    load(x)

    #Total number genetic tests performed
    if(sg==1){#base case (PCE risk range)
      number = length(which(genetic_lb<=ptresults$risk_notrt_threshold[,1]&ptresults$risk_notrt_threshold[,1]<=genetic_ub))
    }else if(sg==2){ #below 50 years old
      number = length(which(ptresults$age==0))
    }else if(sg==3){ #at least 50 years old
      number = length(which(ptresults$age==1))
    }else if(sg==4){ #males
      number = length(which(ptresults$women==0))
    }else if(sg==5){ #females
      number = length(which(ptresults$women==1))
    }

    #Adding over patients and time
    ptresults = lapply(ptresults,sum)

    #Including number of patients tested
    ptresults = c("number"=number,ptresults)

    ptresults},USE.NAMES = F)
  tmp2 = parApply(cl,tmp1,1,cbind)
  pt_sim = parLapply(cl,tmp2,as.numeric)

  #Deleting temporary files
  setwd(file.path(home_dir,dir_names[7]))
  file.remove(tmp)
  rm(tmp,tmp1,tmp2)

  setwd(file.path(home_dir,dir_names[8]))
  save(pt_sim,file = paste0("Results for ",length(tp_results)," patients using ",reps," replications - subgroup ",sg,".RData"))

} #[sg]
stopCluster(cl)

#Returning values to base case
age_subgroup = sg_list[[1]][1] #Indicator for age subgroup analysis (0=base case;1=age<50;2=age>=50)
gender_subgroup = sg_list[[1]][2] #Indicator for gender subgroup analysis (0=base case;1=males;2=females)

# Supplemental Table 13 ---------------------------------------------------


#Supplemental Table 13: Subgroup analysis

#Loading index of patients in testing ranges
setwd(home_dir)
load("Index of Patients in Testing Ranges.RData")

setwd(file.path(home_dir,"Subgroup Analysis Final Results"))
tmp = list.files()

#Creating table
parameters = c("mean","sd")
tmp2 = list()
sup_table13 = c()
for(j in seq(parameters)){
  tmp1 = c()
  for(i in seq(tmp)){

    load(paste0("Results for ",length(SD1)," patients using ",reps," replications - subgroup ",i,".RData"))

    if(j==1){#run for mean calculations (more robust)
      pt_sim = lapply(pt_sim,mean)
    }

    policies = c("genethreshold")
    st13 = list()
    for(p in seq(policies)){
      st13[[p]] = rbind(
        (multiplier/1e06)*pt_sim[["number"]],

        (multiplier/1e06)*(pt_sim[[paste0("modtrt_",policies[p])]]),
        (multiplier/1e06)*(pt_sim[[paste0("hightrt_",policies[p])]]),

        (multiplier/1e00)*(pt_sim[[paste0("chd_","notrt")]]-pt_sim[[paste0("chd_",policies[p])]]),
        (multiplier/1e00)*(pt_sim[[paste0("stroke_","notrt")]]-pt_sim[[paste0("stroke_",policies[p])]]),

        100*((pt_sim[[paste0("ascvd_","notrt")]]-pt_sim[[paste0("ascvd_",policies[p])]])/pt_sim[["number"]]),
        100*((pt_sim[[paste0("qalys_",policies[p])]]-pt_sim[[paste0("qalys_","notrt")]])/pt_sim[["number"]]),

        (multiplier/1e06)*(pt_sim[[paste0("qalys_",policies[p])]]-pt_sim[[paste0("qalys_","notrt")]]),

        (multiplier/1e09)*(pt_sim[[paste0("cost_",policies[p])]]-pt_sim[[paste0("cost_","notrt")]]),

        (pt_sim[[paste0("cost_",policies[p])]]-pt_sim[[paste0("cost_","notrt")]])/
          (pt_sim[[paste0("qalys_",policies[p])]]-pt_sim[[paste0("qalys_","notrt")]]),

        (pt_sim[[paste0("cost_",policies[p])]]-pt_sim[[paste0("cost_","threshold")]])/
          (pt_sim[[paste0("qalys_",policies[p])]]-pt_sim[[paste0("qalys_","threshold")]])
      )
    }

    if(j==2){#run for standard deviation calculations (allowing variability)
      st13 = lapply(st13,function(x) apply(x,1,sd))
    }

    tmp1 = cbind(tmp1,do.call(cbind,st13))
  }
  tmp2[[j]] = tmp1
} #[j]
sup_table13 = do.call(cbind,tmp2)
print(noquote(formatC(sup_table13,digits = 4,format = "f")))
rm(tmp,tmp1,tmp2)

# Number of Replications Plot - Convergence -------------------------------

setwd(home_dir)
load("Health Trajectory Simulation Results 16174 Patients Using 2000 Replications.RData")

marginal_cost_qaly = c()
for(r in seq(length(pt_sim$cost_genethreshold))){
  set.seed(r)
  index = sample(length(pt_sim$cost_genethreshold),r,replace = F)

  #Mean
  marginal_cost_qaly[r] = mean(pt_sim$cost_genethreshold[index]-pt_sim$cost_threshold[index])/
    mean(pt_sim$qalys_genethreshold[index]-pt_sim$qalys_threshold[index])

}

df = data.frame(reps = seq(length(pt_sim$cost_genethreshold))[1:1500],m = marginal_cost_qaly[1:1500])

#Creating plot
setwd(file.path(home_dir,"Figures"))
jpeg("ConvergenceAnalysis.jpeg", width = width, height = height, units = 'in', res = res)
print(ggplot(data=df, aes(x=reps, y=m, group=1)) +
        geom_line(size=line_size,color="black")+
        scale_x_continuous(breaks = seq(0,2e03,5e02))+
        scale_y_continuous(labels=dollar,breaks = seq(20e03,100e03,2e04),limits = c(20e03,100e03))+
        labs(x = "Number of Replications",y = "ICER")+
        geom_hline(yintercept=marginal_cost_qaly[2000], linetype=1,
                   color = "red",size=line_size)+
        theme_bw() +
        theme(panel.grid.minor = element_blank(),
              panel.grid.major = element_line(colour = "gray95"),
              panel.background = element_rect(fill=NA),
              legend.position="none",
              panel.border = element_rect(colour = "black", fill=NA),
              axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black")))
dev.off()


# Selection based on SE ---------------------------------------------------

#Loading data
setwd(home_dir)
load("Health Trajectory Simulation Results 16174 Patients Using 2000 Replications.RData")

N = 10 #Number of resamples
marginal_cost_qaly_exp = c()

for(n in seq(N)){
  marginal_cost_qaly = c()
  for(r in seq(length(pt_sim$qalys_genethreshold))){
    set.seed(r+n)
    index = sample(length(pt_sim$qalys_genethreshold),r,replace = T)

    marginal_cost_qaly[r] = mean(pt_sim$cost_genethreshold[index]-pt_sim$cost_threshold[index])/
      mean(pt_sim$qalys_genethreshold[index]-pt_sim$qalys_threshold[index])

  }
  marginal_cost_qaly_exp = cbind(marginal_cost_qaly_exp,marginal_cost_qaly)
}

ses = apply(marginal_cost_qaly_exp,1,function(x) sd(x))
ses[500]
ses[1000]
