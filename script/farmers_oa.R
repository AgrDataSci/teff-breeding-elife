# Coerce farmers evaluation to rankings
library("PlackettLuce")
library("janitor")
library("tidyr")
library("gosset")

# read the data
load("data/phenotypic.data.EtNAM.Rdata")
load("data/BLUP.values.EtNAM.Rdata")

# farmers evaluation on the genotypes
farmnam <- farmN
# blups agronomic metrics 
blups <- blups.met.nam

names(blups) <- make_clean_names(names(blups))

rm(blups.farm.nam, blups.met.nam, 
   h2.farm.nam, h2.met.nam, metN, farmN)

names(blups)[1] <- "geno"

# retain the averaged blups
blups <- blups[, c("geno", "bm", "db", "df", "dh", "gy",
                   "ph","spl","sps", "tgw")]

str(blups)

blups$geno <- as.character(blups$geno)

# deal with the farm nam data
head(farmnam)

names(farmnam) <- make_clean_names(names(farmnam))

names(farmnam)[names(farmnam)=="id"] <- "geno"

names(farmnam)[6:15]

names(farmnam)[6:15] <- paste0("farmer_", names(farmnam)[6:15])

head(farmnam)

farmnam[1:5] <- lapply(farmnam[1:5], as.character)

# reshape the table with farmers as rows
dat <- pivot_longer(data = farmnam, 
                    cols = starts_with("farmer_"),
                    names_to = "farmer",
                    names_prefix = "farmer_",
                    values_to = "oa",
                    values_drop_na = TRUE)

head(dat)

# fix genotype names
dat$geno[grepl("local_check", dat$geno)] <- "localcheck"

# add ril and family as separated colunms
dat <- dat %>%
  separate(geno, c("family","ril"), sep = "_", remove = FALSE)

unique(dat$family)

f <- as.integer(dat$family)

dat$family[!is.na(f)] <- paste0("N", f[!is.na(f)])

unique(dat$family)

# get an id which is the combination of farmer , repetition and location
dat$id2 <- as.factor(paste(dat$farmer, dat$rep, dat$location, sep = "_"))

dat$id <- as.integer(dat$id2)

dat$gender <- ifelse(grepl("f_", dat$farmer), "F", "M")

# test sample
#dat <- dat[grepl("^N|RF", dat$geno), ]

dat <- dat[!is.na(dat$geno), ]

items <- unique(c(dat$geno, blups$geno))

items <- items[items %in% dat$geno & items %in% blups$geno]

dat <- dat[dat$geno %in% items, ]

blups <- blups[blups$geno %in% items, ]

# make pairwise comparisons
cc <- gosset:::.combn2(items, 2)

dat <- split(dat, dat$id)

pair_oa <- lapply(dat, function(y) {
  # get the rankings as pair comparisons
  # ties are not considered and will be NA's
  pair <- apply(cc, 2, function(x){
    
    # take the first item in the comparison
    i <- x[1]
    # and the second one
    j <- x[2]
    
    # combine the rankings for these two items
    # with i as first and j as the second colunm
    p <- c(y$oa[y$geno == i], y$oa[y$geno == j])
    
    if  (length(p) != 2) {
      p <- c(0, 0)
    }
    
    # if i is higher than j, add 1, this means that i beats j
    # if i is lower than j, add -1, this means that j beats i
    # if none of these options, add NA
    p <- ifelse(p[1] > p[2], 1, ifelse(p[1] < p[2] , -1, NA))
    
  })
  
})

pair_oa

pair <- do.call("rbind", pair_oa)
dat <- do.call("rbind", dat)

id <- names(pair_oa)

# convert this matrix into a paircomp object
pair <- psychotools::paircomp(pair, labels = as.character(items))

# pairwises into grouped rankings 
G <- as.grouped_rankings(pair)

# dataset with rank features 
rank_features <- as.data.frame(dat[!duplicated(dat$id), c("gender","location")])
rank_features[1:2] <- lapply(rank_features[1:2], as.factor)

# dataset with the genotypes features 
geno_features <- blups
rownames(geno_features) <- 1:nrow(geno_features)

# the formula to analyse the genotype features
form <- formula(paste("~", paste(names(blups[-c(1)]), collapse = " + ")))

form

pld <- cbind(G, rank_features)

str(pld)

pl1 <- pltree(G ~ gender + location,
              data = pld,
              alpha = 0.05,
              minsize = 10,
              verbose = TRUE)


pl1

plot(pl1)

top_items(pl1)

worth_map(pl1)

wr <- worst_regret(pl1)

pl2 <- pltree(G ~ gender + location, 
              worth = form,
              data = list(rank_features,
                          geno_features),
              minsize = 20,
              verbose = TRUE)

summary(pl2)

save(pl1, pl2, file = "output/pladmm_model.rda")






