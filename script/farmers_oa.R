# Coerce farmers evaluation to rankings
library("PlackettLuce")
library("janitor")
library("tidyverse")
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

# deal with the farmnam data
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
dat %>%
  separate(geno, c("family","ril"), sep = "_", remove = FALSE) ->
  dat

dat

unique(dat$family)

f <- as.integer(dat$family)

dat$family[!is.na(f)] <- paste0("N", f[!is.na(f)])

unique(dat$family)

# do the same with blups
blups %>% 
  separate(geno, c("family","ril"), sep = "_", remove = FALSE) ->
  blups

f <- as.integer(blups$family)

blups$family[!is.na(f)] <- paste0("N", f[!is.na(f)])

unique(blups$family)

# get an id which is the combination of farmer , repetition and location
dat$id2 <- as.factor(paste(dat$farmer, dat$rep, dat$location, sep = "_"))

dat$id <- as.integer(dat$id2)

dat$gender <- ifelse(grepl("f_", dat$farmer), "F", "M")

length(unique(dat$id))

#.....................................................
#.....................................................
# Sample the data ####
# PlackettLuce cannot work with so many combinations 
# so we sample the data to deal with this issue
dat <- dat[grepl("^N|RF", dat$family), ]
dat <- dat[!is.na(dat$geno), ]

unique(dat$family)

# retain only the items (genotypes) that are present in 
# both datasets
items <- unique(c(dat$geno, blups$geno))

items <- items[items %in% dat$geno & items %in% blups$geno]

length(items)

# filter the data to keep only the items that are present in both
dat <- dat[dat$geno %in% items, ]

blups <- blups[blups$geno %in% items, ]

# make the sample by family
datS <- split(dat, dat$family)

RF <- which(names(datS) %in% "RF")

# samplesize <- 25
# datS[-RF] <- lapply(datS[-RF], function(x){
#   set.seed(4032022)
#   # take a sample
#   genoS <- sample(unique(x$geno), samplesize)
#   # filter to keep sampled items
#   x <- x[x$geno %in% genoS, ]
#   # rename items to avoid computation issues in PlackettLuce
#   x$geno <- paste0(x$family, "-" , as.integer(as.factor(x$geno)))
#   x
# })
# datS <- do.call("rbind", datS)

# take the average per family and assess by family
dat %>% 
  group_by(id2, family) %>% 
  summarise(oa = mean(oa)) %>% 
  mutate(geno = family,
         id = id2) %>% 
  ungroup() ->
  datS

datS

length(unique(datS$geno))

itemsS <- unique(datS$geno)

# make pairwise comparisons
cc <- gosset:::.combn2(itemsS, 2)

ncol(cc)
datS <- split(datS, datS$id)

pair_oa <- lapply(datS, function(y) {
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
datS <- do.call("rbind", datS)

# convert this matrix into a paircomp object
pair <- psychotools::paircomp(pair, labels = as.character(itemsS))

# pairwises into grouped rankings
G <- as.grouped_rankings(pair)

mod <- PlackettLuce(G)

summary(mod, ref = "RF")

# dataset with rank features 
unique(datS$id)
datS$gender <- ifelse(grepl("f_", datS$id), "F", "M")
datS$location <- ifelse(grepl("_adet", datS$id), "Adet",
                        ifelse(grepl("_kulumsa", datS$id), "Kulumsa",
                               "Geregera")) 
 
rank_features <- as.data.frame(datS[!duplicated(datS$id), c("gender","location")])
rank_features[1:2] <- lapply(rank_features[1:2], as.factor)

pld <- cbind(G, rank_features)

str(pld)

pl1 <- pltree(G ~ gender + location,
              data = pld,
              alpha = 0.01,
              minsize = 20,
              verbose = TRUE)


pl1

plot(pl1)

top_items(pl1)

worth_map(pl1)

reliability(pl1, ref = "RF")

coef(pl1, log = F)

# dataset with the genotypes features 
geno_features <- blups
names(geno_features)
geno_features <- geno_features[,-which(grepl("geno|ril", names(geno_features)))]
# group blups by family
geno_features <- aggregate(. ~ family, data = geno_features, mean)
names(geno_features)[1] <- "geno"

# the formula to analyse the genotype features
f <- formula(paste("~", paste(names(geno_features[-c(1)]), collapse = " + ")))

f

pl2 <- pltree(G ~ gender + location, 
              worth = f,
              data = list(pld,
                          geno_features),
              minsize = 20,
              verbose = TRUE)

summary(pl2)

coef(pl2)

save(pl1, pl2, file = "output/pladmm_model.rda")


plot(pl1)

plot(pl2)

class(pl2)

