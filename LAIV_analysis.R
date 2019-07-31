library(data.table)
library(ggplot2)
library(cowplot)
library(HDInterval)

# Directory where this script and analysis data files are kept
root = "~/Dropbox/LAIV/Analysis/";

# Plot setup
theme_set(theme_cowplot(font_size = 7, font_family = "Helvetica", line_size = 0.25))

# Age distributions for England and Wales are derived from ONS survey data.
age.dist = fread(paste0(root, "ONS_UKpop.csv"));
age.dist = age.dist[country == "E" | country == "W", lapply(.SD, sum), by = "Age", .SDcols = 6:21];

# Convert Fleming age-stratified measures to Cromer age-stratified measures
# (0-4, 5-17, 18-49, 50-64, 65+) to (<6m, 6m-4, 5-14, 15-44, 45-64, 65+)
f2c = function(fm)
{
    # We use the distribution from 2001 as that is the midpoint of Fleming's study.
    popsize = age.dist$population_2001;
    rates = c(rep(fm[1],  5),  # 0-4
              rep(fm[2], 13),  # 5-17
              rep(fm[3], 32),  # 18-49
              rep(fm[4], 15),  # 50-64
              rep(fm[5], 26)); # 65+

    # Reweight based on age groups used in this analysis
    c(weighted.mean(rates[(0 :0)+1],  popsize[(0 :0)+1]),  # <6m
      weighted.mean(rates[(0 :4)+1],  popsize[(0 :4)+1]),  # 6m-4
      weighted.mean(rates[(5 :14)+1], popsize[(5 :14)+1]), # 5-14
      weighted.mean(rates[(15:44)+1], popsize[(15:44)+1]), # 15-44
      weighted.mean(rates[(45:64)+1], popsize[(45:64)+1]), # 45-64
      weighted.mean(rates[(65:90)+1], popsize[(65:90)+1])) # 65+
}

# Draw from triangular distribution specified according to mode and HDI (l,r)
# containing proportion p of density, clamped to range [minimum, maximum]
rtriang = function(n, mode, l, r, p = 0.95, minimum = -Inf, maximum = Inf)
{
    # bounds of triangular distribution: a = left vertex, b = right vertex, c = peak
    h = sqrt(1 - p);
    a = (h*mode - l) / (h - 1);
    b = (r - h*mode) / (1 - h);
    c = mode;
    
    # choose variates from triangular distribution
    u = runif(n);
    x = ifelse(u < (c - a) / (b - a), a + sqrt(u * (b - a) * (c - a)), b - sqrt((1 - u) * (b - a) * (b - c)));
    return (pmin(maximum, pmax(minimum, x)));
}

# Draw from success probability distribution according to number of successes in
# a series of Bernoulli trials. This assumes a uniform prior on the success probability.
rprop = function(n, successes, trials)
{
    return (rbeta(n, successes + 1, trials - successes + 1));
}


######################################################################
#          PART I. IMPACT OF LAIV ON ANTIBIOTIC PRESCRIBING          #
######################################################################

# age groups over which to run the analysis
AGs = c("0-6m", "6m-4y", "5-14y", "15-44y", "45-64y", "65y+");
age.groups = factor(AGs, levels = AGs);
n.age.groups = length(age.groups);
pop = age.dist$population_2015;
group.sizes = c(pop[1] / 2,
                pop[1] / 2 + sum(pop[2:5]),
                sum(pop[6:15]),
                sum(pop[16:45]),
                sum(pop[46:65]),
                sum(pop[66:length(pop)]));

# Function to run analysis.
# n.samples:      number of Monte Carlo samples for analysis.
# gp.scenario:    GP visits scenario. 1 = Hayward (FluWatch), 2 = Fleming, 3 = Cromer
# uptake.scenario LAIV uptake scenario (paediatric uptake). Must be 0, .15, .3, .5, or .7
# abx.scenario    Antibiotic prescribing scenario. 1 = Pouwels, 2 = Fleming, 3 = 100%
LAIV.Impact = function(n.samples, gp.scenario, uptake.scenario, abx.scenario)
{
    # data table with main analysis components
    analysis = data.table(nsamp = rep(1:n.samples, each = n.age.groups),
        age.group = rep(age.groups, n.samples));
    
    # Uncertainty for estimates in GP visits attributable to influenza from Fleming et al, 2016.
    # These can be read from the paper's Fig. 2, which shows 95% confidence intervals for
    # influenza-attributable GP visits. Using the width of the error bars on this plot, extracted
    # using WebPlotDigitizer, we can estimate the mean relative standard deviation for this
    # uncertainty distribution.
    ff2 = fread(paste0(root, "fleming2016_fig2.txt"));
    ff2[, influ.mean := influ.A.mid + influ.B.mid];
    ff2[, influ.sd := sqrt((influ.A.2sd / 1.96)^2 + (influ.B.2sd / 1.96)^2)];
    ff2[, influ.rel.sd := influ.sd / influ.mean];
    fleming.rel.sd = ff2[, mean(influ.rel.sd)];
    
    
    # I (1) GP visits attributable to influenza
    # -----------------------------------------
    
    # Vector of standard normal variates for uncertainty calculations
    gp.variation = rep(rnorm(n.samples), each = n.age.groups);
    
    if (gp.scenario == 1) {
        # Rates (per person season) of PCR confirmed illness: Hayward et al., 2014
        fl00_04 = (rtriang(n.samples, 7.76, 2.07, 29.12, 0.95, 0) +
            rtriang(n.samples, 0, 0, 9.46, 0.95, 0) +
            rtriang(n.samples, 6.32, 1.70, 23.43, 0.95, 0) +
            rtriang(n.samples, 9.27, 5.06, 16.97, 0.95, 0) +
            rtriang(n.samples, 14.31, 7.06, 28.99, 0.95, 0)) / 500;
        
        fl05_15 = (rtriang(n.samples, 4.43, 1.57, 12.52, 0.95, 0) +
            rtriang(n.samples, 5.94, 2.33, 15.16, 0.95, 0) +
            rtriang(n.samples, 25.78, 15.23, 43.62, 0.95, 0) +
            rtriang(n.samples, 9.21, 6.50, 13.05, 0.95, 0) +
            rtriang(n.samples, 12.99, 8.09, 20.84, 0.95, 0)) / 500;
        
        fl16_44 = (rtriang(n.samples, 3.51, 1.32, 9.31, 0.95, 0) +
            rtriang(n.samples, 1.53, 0.54, 4.35, 0.95, 0) +
            rtriang(n.samples, 7.42, 3.95, 13.91, 0.95, 0) +
            rtriang(n.samples, 3.38, 2.02, 5.67, 0.95, 0) +
            rtriang(n.samples, 13.23, 8.37,  20.90, 0.95, 0)) / 500;
        
        fl45_64 = (rtriang(n.samples, 1.45, 0.46, 4.54, 0.95, 0) +
            rtriang(n.samples, 3.16, 1.48, 6.74, 0.95, 0) +
            rtriang(n.samples, 6.86, 4.18, 11.26, 0.95, 0) +
            rtriang(n.samples, 1.97, 1.19, 3.24, 0.95, 0) +
            rtriang(n.samples, 4.79, 2.71, 8.47, 0.95, 0)) / 500;
        
        fl65_pl = (rtriang(n.samples, 2.26, 0.56, 9.12, 0.95, 0) +
            rtriang(n.samples, 0, 0, 4.10, 0.95, 0) +
            rtriang(n.samples, 1.41, 0.37, 5.40, 0.95, 0) +
            rtriang(n.samples, 0.18, 0.03, 1.21, 0.95, 0) +
            rtriang(n.samples, 4.62, 2.11, 10.10, 0.95, 0)) / 500;
        
        # Probability of GP visit given PCR confirmed illness. Hayward et al (2014) provide an age-
        # stratified estimate of this but the numbers are low. So we use the rate of GP visits given
        # PCR-confirmed illness in 5-64 year olds as a "base value" (proportion 12/82), inflating 
        # the estimates for 0-4 year olds and 65+ year olds by 12% from this base value. This
        # results in an overall rate which is consistent with that measured by Hayward et al.
        gpfl = rprop(n.samples, 12, 82);
        gp00_04 = gpfl + 0.12;
        gp05_15 = gpfl;
        gp16_44 = gpfl;
        gp45_64 = gpfl;
        gp65_pl = gpfl + 0.12;
        
        # Save intermediate steps
        analysis$flu.rate = c(rbind(fl00_04, fl00_04, fl05_15, fl16_44, fl45_64, fl65_pl));
        analysis$gp.visit.prob.given.flu = c(rbind(gp00_04, gp00_04, gp05_15, gp16_44, gp45_64, gp65_pl));
        
        # Combine: GP consultation rate as rate of disease times probability of consulting
        analysis$gp.visit.rate = c(rbind(fl00_04 * gp00_04 * 1000,   # <6m
                                         fl00_04 * gp00_04 * 1000,   # 6m-4
                                         fl05_15 * gp05_15 * 1000,   # 5-14
                                         fl16_44 * gp16_44 * 1000,   # 15-44
                                         fl45_64 * gp45_64 * 1000,   # 45-64
                                         fl65_pl * gp65_pl * 1000)); # 65+
    } else if (gp.scenario == 2) {
        # GP consultations attributable to influenza: Fleming et al., 2016
        analysis$gp.visit.rate = rep(f2c(c(2436 + 539, 1217 + 996, 899 + 288, 1076 + 179, 1140 + 76)), n.samples);
        # Add variation due to uncertainty: Fleming et al., 2016
        analysis$gp.visit.rate = analysis$gp.visit.rate * (1 + gp.variation * fleming.rel.sd);
        # Adjust to be per 1000 rather than per 100,000
        analysis$gp.visit.rate = analysis$gp.visit.rate / 100;
    } else if (gp.scenario == 3) {
        # GP consultations attributable to influenza: Cromer et al., 2014
        analysis$gp.visit.rate = rep(c(7361, 6090, 3875, 1878, 1829, 582), n.samples);
        # Add variation due to uncertainty: Cromer et al., 2014
        analysis$gp.visit.rate = analysis$gp.visit.rate + gp.variation * c(304, 170, 107, 35, 34, 26) / 1.96;
        # Adjust to be per 1000 rather than per 100,000
        analysis$gp.visit.rate = analysis$gp.visit.rate / 100;
    }
    
    
    # I (2) Impact of LAIV on GP visits
    # ---------------------------------
    
    # Load LAIV impact as data table.
    # Source: PHE influenza model with paediatric LAIV impact (Baguelin et al., 2013)
    load(paste0(root, "2-16yo.rda"));
    influenza = as.data.table(prediction_2013_df);
    setkey(influenza, np);
    
    # Classify by age group to combine high risk and low risk groups
    influenza[Group %like% "0 -",     age.group := age.groups[1]];
    influenza[Group %like% "4 y",     age.group := age.groups[2]];
    influenza[Group %like% "14 y",    age.group := age.groups[3]];
    influenza[Group %like% "44 y",    age.group := age.groups[4]];
    influenza[Group %like% "64 y",    age.group := age.groups[5]];
    influenza[Group %like% "65\\+ y", age.group := age.groups[6]];
    setorder(influenza, np, age.group);
    
    # Calculate impact of LAIV
    pre.LAIV  = influenza[uptake == 0,               .(value = sum(value)), by = .(np, age.group)];
    post.LAIV = influenza[uptake == uptake.scenario, .(value = sum(value)), by = .(np, age.group)];
    analysis$LAIV.impact = rep((pre.LAIV$value - post.LAIV$value) / pre.LAIV$value, length.out = n.samples * n.age.groups);
    
    
    # I (3) Antibiotic prescriptions resulting from GP visits attributable to influenza
    # ---------------------------------------------------------------------------------
    
    if (abx.scenario == 1) {
        # Estimates based on Pouwels et al's (2018) analysis of prescribing rates for different
        # conditions in English primary care. From Hayward et al. 2014, 88.1% of influenza-attributable
        # consultations are for ILI while the rest are for acute respiratory illness without fever.
        # Pouwels et al. 2018 report that 29% of ILI consultations and 48% of acute cough consultations
        # result in a prescription for antibiotics within 30 days. We use cough for ARI-F.
        
        # In order to age-stratify the results, we assume prescribing in under-5s is around 20% less,
        # and in over-45s is around 20% more, than prescribing in 5-44yos, in keeping with Fleming et
        # al (2016). This is also consistent with results from Meier et al (Eur J Clin Microbiol Infect
        # Dis (2000) 19:834–842) and Pitman et al (Vaccine 31 (2013) 927–942).
        abx_skew = rnorm(n.samples, 0.2, 0.05);
        abx_factor = cbind(1 - abx_skew, 1 - abx_skew, 1 - abx_skew, 1, 1 + abx_skew, 1 + abx_skew);
        abx_norm = rowSums(t(t(abx_factor) * group.sizes)) / sum(group.sizes);
        abx_factor = abx_factor / abx_norm;
        
        # Final result
        analysis$abx.gp.rate = (0.881 * 0.29 + 0.119 * 0.48) * c(t(abx_factor));
    } else if (abx.scenario == 2) {
        # Divide antibiotic rate by consultation rate to get rate of prescribing given consultation.
        # GP consultation rates attributable to influenza: Fleming et al. 2016, Table 2
        gpr = f2c(c(2436 + 539, 1217 + 996, 899 + 288, 1076 + 179, 1140 +  76));
    
        # Antibiotic prescription rates attributable to influenza: same source
        rxr = f2c(c(1429 + 348, 699 + 602, 612 + 210, 911 + 148, 980 +  62));
    
        # Fleming et al. do not report uncertainty for antibiotic prescription rates (just seasonal
        # fluctuations). So we assume the same variability as for consultations.
        analysis$abx.gp.rate = (rxr / gpr) * rep(rnorm(n.samples, 1, fleming.rel.sd), each = n.age.groups);
    } else if (abx.scenario == 3) {
        # Assume prescription rate follows that of acute rhinosinusitis in England, 0.90 (Pouwels et al 2018)
        abx_skew = rnorm(n.samples, 0.2, 0.05);
        abx_factor = cbind(1 - abx_skew, 1 - abx_skew, 1 - abx_skew, 1, 1 + abx_skew, 1 + abx_skew);
        abx_norm = rowSums(t(t(abx_factor) * group.sizes)) / sum(group.sizes);
        abx_factor = abx_factor / abx_norm;
        
        # Final result
        analysis$abx.gp.rate = 0.9 * c(t(abx_factor));
    }
    
    # Overall calculation of LAIV impact on antibiotic prescription incidence
    analysis[, abx.inc.reduction := gp.visit.rate * LAIV.impact * abx.gp.rate];
    
    return (analysis)
}


# I (4) Results and figures
# -------------------------

# Get "headline" reduction in antibiotic prescription incidence
CrudeAbxReduction = function(analysis)
{
    cbind(analysis[, mean(abx.inc.reduction), by = age.group], group.sizes)[, weighted.mean(V1, group.sizes)]
}

# Get table of age-stratified estimates for LAIV impact analysis
AgeStratified = function(analysis, gs)
{
    # Overall...
    analysis$pop.size = rep(gs, nrow(analysis) / n.age.groups);
    overall = analysis[, lapply(.SD, function(x) weighted.mean(x, pop.size)),
                       by = nsamp, .SDcols = gp.visit.rate:abx.inc.reduction];
    
    # Results...
    res = rbind(
        cbind(measure = "Mean",
              analysis[, lapply(.SD, mean), by = age.group, .SDcols = gp.visit.rate:abx.inc.reduction]),
        cbind(measure = c("Lo95HDI", "Hi95HDI"),
              analysis[, lapply(.SD, hdi), by = age.group, .SDcols = gp.visit.rate:abx.inc.reduction]),
        cbind(measure = "Mean", age.group = "Overall",
              overall[, lapply(.SD, mean), .SDcols = gp.visit.rate:abx.inc.reduction]),
        cbind(measure = c("Lo95HDI", "Hi95HDI"), age.group = "Overall",
              overall[, lapply(.SD, hdi), .SDcols = gp.visit.rate:abx.inc.reduction]));
    
    res$measure = factor(res$measure, levels = c("Mean", "Lo95HDI", "Hi95HDI"));
    res$age.group = factor(res$age.group, levels = c(as.character(age.groups), "Overall"));
    res = res[order(age.group, measure)];
    return (res)
}

# Printable table
PrintAgeStratified = function(agestrat)
{
    table1 = paste0("\t", paste0(names(agestrat)[3:ncol(agestrat)], collapse = "\t"), "\n");
    prec = 3;
    for (r in seq(1, nrow(agestrat) - 2, by = 3))
    {
        row = agestrat$age.group[r];
        for (co in 3:ncol(agestrat)) {
            row = paste0(row, "\t", signif(agestrat[r, ..co], prec),
                " (", signif(agestrat[r + 1, ..co], prec),
                " – ", signif(agestrat[r + 2, ..co], prec), ")");
        }
        table1 = paste0(table1, row, "\n");
    }
    cat(table1)
}

# Run analysis: n.samples, gp.scenario, uptake.scenario, abx.scenario
n.samples = 1000000;
analysis = LAIV.Impact(n.samples, 2, 0.5, 3);

# Get results
agestrat = AgeStratified(analysis, group.sizes);
PrintAgeStratified(agestrat)

# normal:         5.320942
# GP scenario 1:  4.3081028
# GP scenario 3:  7.712143
# uptake 0.3:     3.515785
# uptake 0.7:     6.666193
# abx scenario 1: 2.254137
# abx scenario 3: 6.997558
# all low:        1.2355660
# all high:       12.992868

# Direct protection "RCT" results
analysis[, abx.rate := gp.visit.rate * abx.gp.rate];
analysis[, vacc.efficacy := rep(c(0, 0.56, 0.56, 0.56, 0, 0), n.samples)];
analysis[, direct.reduction := abx.rate * vacc.efficacy];
analysis[, mean(direct.reduction), by = "age.group"];
setcolorder(analysis, c(1,2,3,4,5,7,8,9,6));
agestrat = AgeStratified(analysis, c(0, p2_4, p5_14, p15_15, 0, 0))
PrintAgeStratified(agestrat);


######################################################################
#        PART II. IMPACT OF PRESCRIBING ON RESISTANT DISEASE         #
######################################################################

# II (1) Loading data
# -------------------

# load.drugs: loads drug consumption for a subset of years, sectors (AC, HC, or ACHC), and ATC codes.
load.drugs = function(years, sectors, atcs)
{
    drugs = fread(paste0(root, "drugs.txt"));
    drugs = drugs[year %in% years & sector %in% sectors & atc %in% atcs];
    return (dcast(drugs, country ~ atc + sector + year, value.var = "ddd.per.thousand"))
}

# load burden: loads resistant disease burdens normalised according to norm (none, pop for
# population by country, or BSI for total bloodstream infections by each species). BSI.agg
# should be either mean or max. Note that this loads the number of cases from original
# files, so the calculated incidences are not standardised for different age distributions.
load.burden = function(norm, BSI.agg = max)
{
    burden = fread(paste0(root, "burden.txt"));

    if (norm == "none") {
        # no normalization
        burden$norm = 1;
    } else if (norm == "pop") { 
        # normalize disease burden relative to population in thousands
        pop.cov = fread(paste0(root, "population.and.coverage.txt"));
        burden$norm = pop.cov$pop.2015.eurostat[match(burden$country, pop.cov$country)] / 1000;
    } else if (norm == "BSI") {
        # normalise disease burden relative to total number of bloodstream infections
        # caused by the species in question (whether resistant or sensitive)
        pop.cov = fread(paste0(root, "population.and.coverage.txt"));
        isolates = fread(paste0(root, "resistance_2015.csv"), na.strings = "-");
        isolates = isolates[Indicator == "Total tested isolates", c("Population", "RegionName", "NumValue")];
        
        # Collate total tested isolates and coverage for all strains
        spp = c("Acinetobacter/ColRACI/CRACI/MDRACI",
                "Enterococcus faecalis/VRE",
                "Enterococcus faecium/VRE",
                "Escherichia/ColREC/CREC/3GCREC",
                "Klebsiella/ColRKP/CRKP/3GCRKP",
                "Pseudomonas/ColRPA/CRPA/MDRPA",
                "Staphylococcus/MRSA",
                "Streptococcus/PRSP/PMRSP");
        bsi = NULL;
        for (s in spp) {
            args = strsplit(s, "/")[[1]];
            species = args[1];
            for (strain in args[2:length(args)]) {
                bsi.uncorr = isolates[Population %like% species,
                                      .(strain = ..strain, tested = BSI.agg(NumValue)),
                                      by = .(country = RegionName)];
                bsi.correction = cbind.data.frame(country = pop.cov$country, coverage = pop.cov[[paste0("coverage.", strain)]]);
                bsi = rbind(bsi, merge(bsi.uncorr, bsi.correction, by = "country"));
            }
        }

        # Sum together both Enterococcus species
        bsi = rbind(bsi[strain != "VRE"],
                    bsi[strain == "VRE", .(tested = sum(tested)), by = .(country, strain, coverage)][, c(1,2,4,3)]);
        bsi[, norm := tested / (coverage / 100)];
        
        burden = merge(burden, bsi[, .(country, strain, norm)], by = c("country", "strain"))
    } else {
        stop("norm must be one of: none, pop, BSI.");
    }
    
    # perform normalisation
    burden[, 3:11] = burden[, 3:11] / burden$norm;
    
    # Greece's unreported S. pneumoniae burdens are given as 0 in this data set, so remove these
    burden[country == "Greece" & strain %like% "^PM?RSP$", 3:11] = NA;
    
    return (burden)
}

# load.data: loads both drug consumption and burden data.
load.data = function(years, sectors, atcs, norm, BSI.agg = max)
{
    d = load.drugs(years, sectors, atcs);
    b = load.burden(norm, BSI.agg);
    return (merge(b, d, by = "country"))
}

# II (2) Predictive model
# -----------------------

# predict: predict the impact of reducing prescribing on outcomes. 
# data is a data.table from the function load.data;
# countries is the vector of countries;
# strains is a vector of strains to predict the impact upon;
# outcomes is a vector of outcomes (e.g. DALY.mid);
# predictors is a vector of predictors (e.g. J01_AC_2015);
# projections is a vector of the projected change in each predictor;
predict = function(data, countries, strains, outcomes, predictors, projections, pop.factor)
{
    impacts = NULL;
    for (outcome in outcomes) {
        for (str in strains) {
            # Build model
            subdata = data[strain == str, c(..outcome, ..predictors)];
            model = lm(paste(outcome, "~ ."), data = subdata);
            
            # Make predictions
            for (co in countries) {
                change = 0;
                for (p in 1:length(predictors)) {
                    change = change + projections[p] * model$coeff[p + 1] * 
                        data[country == co & strain == str, norm];
                }
                impact = data.table(country = co, strain = str, outcome = outcome, 
                                    predictors = paste(predictors, collapse = "+"),
                                    orig = data[country == co & strain == str][[outcome]] * data[country == co & strain == str, norm] * pop.factor,
                                    change = change * pop.factor,
                                    final = data[country == co & strain == str, ..outcome][[outcome]] * 
                                        data[country == co & strain == str, norm] * pop.factor + change * pop.factor);
                impacts = rbind(impacts, impact);
            }
        }
    }
    return (impacts)
}

predict0 = function(data, countries, strains, outcomes, predictors, projections, pop.factor)
{
    impacts = NULL;
    for (outcome in outcomes) {
        for (str in strains) {
            # Build model, removing predictors with negative coefficients
            pred = predictors;
            proj = projections;
            repeat {
                subdata = data[strain == str, c(..outcome, ..pred)];
                model = lm(paste(outcome, "~ ."), data = subdata);
                negatives = model$coeff[-1] < 0;
                if (any(negatives)) {
                    pred = pred[!negatives];
                    proj = proj[!negatives];
                } else {
                    break;
                }
            }
            
            # Make predictions
            for (co in countries) {
                change = 0;
                for (p in 1:length(pred)) {
                    change = change + proj[p] * model$coeff[p + 1] * 
                        data[country == co & strain == str, norm];
                }
                impact = data.table(country = co, strain = str, outcome = outcome, 
                                    predictors = paste(pred, collapse = "+"), change = change * pop.factor);
                impacts = rbind(impacts, impact);
            }
        }
    }
    return (impacts)
}

# II (3) Analysis
# ---------------

# overall reduction for main analysis: 5.322 (3.731-6.985)
delta = -7 * 5.322 / 365;
delta = -7 * 3.731 / 365;
delta = -7 * 6.985 / 365;

reg = load.data(2015, "AC", "J01", "BSI", max);

impacts = predict(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  "J01_AC_2015", delta, sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]

imp = impacts[outcome == "DALY.mid"]
imp[strain %like% "EC$", species := "E. coli"]
imp[strain %like% "SA$", species := "S. aureus"]
imp[strain %like% "KP$", species := "K. pneumoniae"]
imp[strain %like% "E$", species := "Enterococcus spp."]
imp[strain %like% "PA$", species := "P. aeruginosa"]
imp[strain %like% "SP$", species := "S. pneumoniae"]
imp[strain %like% "ACI$", species := "Acinetobacter spp."]

imp2 = imp[, .(DALYs.orig = sum(orig), DALYs.final = sum(final), DALYs.change = sum(change)), by = .(species)]
imp2[, label := round(DALYs.change, 0)]
imp2[, label := paste0(ifelse(label > 0, "+", "-"), abs(label))]

imp2$species = factor(imp2$species, levels = imp2$species[order(imp2$DALYs.orig)])
ggplot(imp2) + 
    geom_col(aes(x = species, y = DALYs.orig), fill = "#ff9999") + 
    geom_col(aes(x = species, y = DALYs.final), fill = "#9999ff") + 
    geom_text(aes(x = species, y = DALYs.orig + 300, label = label), hjust = 0, size = 2) +
    ylim(0, 35000) +
    labs(x = NULL, y = "Resistance-attributable DALYs") +
    coord_flip()
ggsave("~/Dropbox/LAIV/Analysis/byspecies.pdf", width = 8.7, height = 6, units = "cm")

imp[, label := round(change, 0)]
imp[, label := paste0(ifelse(label > 0, "+", "-"), abs(label))]
imp$strain = factor(imp$strain, levels = imp$strain[order(imp$orig)])
ggplot(imp) + 
    geom_col(aes(x = strain, y = orig/1000), fill = "#ff9999") + 
    geom_col(aes(x = strain, y = final/1000), fill = "#9999ff") + 
    geom_text(aes(x = strain, y = orig/1000 + 0.2, label = label), hjust = 0, size = 3) +
    ylim(0, 30) +
    labs(x = NULL, y = "Disability-adjusted life years (thousands)") +
    coord_flip()
ggsave("~/Dropbox/LAIV/Analysis/bystrain.pdf", width = 8, height = 8, units = "cm")

reg = load.data(2015, "AC", c("J01A", "J01C", "J01F"), "BSI", max)
impacts = predict(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  c("J01A_AC_2015", "J01C_AC_2015", "J01F_AC_2015"), delta * c(0.0620, 0.7545, 0.1835), sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]

reg = load.data(2015, "AC", c("J01AA", "J01CA", "J01CE", "J01FA"), "BSI", max)
impacts = predict(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  c("J01AA_AC_2015", "J01CA_AC_2015", "J01CE_AC_2015", "J01FA_AC_2015"), delta * c(0.0620, 0.4752, 0.2793, 0.1835), sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]




reg = load.data(2015, "AC", "J01", "BSI", max)
impacts = predict0(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  "J01_AC_2015", delta, sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]

reg = load.data(2015, "AC", c("J01A", "J01C", "J01F"), "BSI", max)
impacts = predict0(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  c("J01A_AC_2015", "J01C_AC_2015", "J01F_AC_2015"), delta * c(0.0620, 0.7545, 0.1835), sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]

reg = load.data(2015, "AC", c("J01AA", "J01CA", "J01CE", "J01FA"), "BSI", max)
impacts = predict0(reg, "United Kingdom", unique(reg$strain), c("DALY.mid", "cases.mid", "deaths.mid"),
                  c("J01AA_AC_2015", "J01CA_AC_2015", "J01CE_AC_2015", "J01FA_AC_2015"), delta * c(0.0620, 0.4752, 0.2793, 0.1835), sum(pop) / 65128861)
impacts
impacts[, .(change = sum(change)), by = outcome]




######################################################################
#                            MISCELLANEOUS                           #
######################################################################

# Consumption vs burdens illustrative scatter plot
dbc = load.data(2015, "AC", "J01", "pop")[strain == "total"];
dbc = melt(dbc, measure.vars = patterns("lo", "mid", "hi"), value.name = c("lo", "mid", "hi"))
dbc[variable == 1]$variable = "Attributable DALYs";
dbc[variable == 2]$variable = "Resistant cases";
dbc[variable == 3]$variable = "Attributable deaths";
dbc$variable = factor(dbc$variable, levels = c("Resistant cases", "Attributable DALYs", "Attributable deaths"));

ggplot(dbc, aes(x = J01_AC_2015)) +
  geom_smooth(method = "lm", aes(y = mid), size = 0.5, colour = "#7788ff", fill = "#cccccc") + 
  geom_linerange(aes(ymin = lo, ymax = hi), size = 0.3) +
  geom_segment(aes(x = J01_AC_2015 - 0.1, xend = J01_AC_2015 + 0.1, y = mid, yend = mid), size = 0.3) +
  geom_point(data = dbc[country == "United Kingdom"], aes(y = mid), shape = 1, size = 2, colour = "red") +
  facet_wrap(~variable, ncol = 1, scales = "free") +
  labs(x = "Primary care antibiotic consumption (DDD per thousand person-days)", y = "Incidence (per 1000 person-years)") +
  theme(strip.background = element_blank(), legend.position = "none")
ggsave("~/Dropbox/LAIV/Analysis/regression.pdf", width = 8.7, height = 8, units = "cm");

# Economic calculations
# 1415 in USD 2016 / 101.38% (Bureau of Labour Statistics CPI) -> USD 2014
1415 / 1.0138
# = 1395.739 USD 2014 x 79 (GBR) / 130 (USD) health care PPP 2014 (https://www.oecd.org/health/health-systems/
# International-Comparisons-of-Health-Prices-and-Volumes-New-Findings.pdf) x 1.647701 GBP
1395.739 * 79 / 130 * 1 / 1.647701
# = 514.7656 GBP (2014) x 1.0099 -> 2015 GBP
514.7656 * 1.0099
# = 519.8618 GBP 2015

# Total price of 432 resistant infections in 2015 pounds
519.8618 * 431.62048

# LAIV program cost 2008 GBP to 2015 GBP, millions
205*1.2035

# PLOT: Uncertainty analysis
tornado = fread("Factor,Lo,Mid,Hi,Range
Consultation rate,520,642,931,411
Vaccine uptake,424,642,804,380
Prescription rate,272,642,642,370
Statistical model,414,642,762,348");
torfactors = c("Consultation\nrate", "Vaccine\nuptake", "Prescription\nrate", "Statistical\nmodel");
tornado$Factor = factor(torfactors, levels = rev(torfactors));

ggplot(tornado) + 
    geom_col(aes(x = Factor, y = Lo - Mid), fill = "#ff9999") + 
    geom_col(aes(x = Factor, y = Hi - Mid), fill = "#9999ff") +
    geom_hline(aes(yintercept = 0), colour = "#880000", linetype = "dashed") +
    scale_y_continuous(breaks = c(-300, -200, -100, 0, 100, 200, 300) - 42, labels = c(342, 442, 542, 642, 742, 842, 942) - 42) +
    labs(x = NULL, y = "DALYs averted") +
    coord_flip()
ggsave("~/Dropbox/LAIV/Analysis/sensitivity.pdf", width = 8.7, height = 6, units = "cm")

# PLOT: Secular trend in antibiotic prescribing
usage = fread("~/Documents/UK/England.BNF05.txt")
usage = usage[BNF.CODE %like% "^0501", sum(items), by = PERIOD]
usage[, year := floor(PERIOD/100)]
usage[, month := PERIOD - year * 100]
usage[, time := year + (month-1)/12]
ggplot(usage[time >= 2012]) + geom_line(aes(x = time, y = V1/1000)) + geom_smooth(aes(x = time, y = V1/1000), method = "lm") + labs(x = "Year", y = "Thousands of antibiotic prescriptions")
