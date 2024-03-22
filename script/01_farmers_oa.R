# Assess farmers evaluation of genotypes 
# using PlackettLuce model
library("PlackettLuce")
library("multcompView")
library("janitor")
library("tidyverse")
library("gosset")
library("gtools")
library("lme4")

pladmm_coeffs <- function(object, ...) {
  
  # Extract ids from terminal nodes
  node_id <- partykit::nodeids(object, terminal = TRUE)
  
  # get models from each node
  nodes <- list()
  for (i in seq_along(node_id)) {
    obj_i <- object[[ node_id[i] ]]$node$info$object
    coefs <- coef(obj_i)
    coefficients <- matrix(NA, nrow = length(coefs), ncol = 4L, 
                           dimnames = list(names(coefs), c("Estimate", "Std. Error", 
                                                           "z value", "Pr(>|z|)")))
    coefficients[, 1L] <- coefs
    se <- sqrt(diag(vcov(obj_i)))
    coefficients[names(se), 2L] <- se
    coefficients[, 3L] <- coefficients[, 1L]/coefficients[, 2L]
    coefficients[, 4L] <- 2L * pnorm(-abs(coefficients[, 3L]))
    
    coefficients <- as.data.frame(coefficients)
    
    coefficients[, 5] <- gtools::stars.pval(coefficients[, 4])
    
    coefficients[, 4] <- formatC(coefficients[, 4], format = "e", digits = 2)
    
    coefficients[, 6] <- node_id[i]
    
    coefficients[, 7] <- rownames(coefficients)
    
    rownames(coefficients) <- 1:nrow(coefficients)
    
    coefficients <- coefficients[,c(6, 7, 1:5)]
    
    names(coefficients)[1] <- "Node"
    
    names(coefficients)[c(2, 7)] <- ""
    
    nodes[[i]] <- coefficients
    
  }
  
  result <- do.call("rbind", nodes)
  
  rownames(result) <- 1:nrow(result)
  
  result
  
}

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

# get the combination of farmer , repetition and location as an id
dat$id2 <- as.factor(paste(dat$farmer, dat$rep, dat$location, sep = "_"))

dat$id <- as.integer(dat$id2)

dat$gender <- ifelse(grepl("f_", dat$farmer), "F", "M")

length(unique(dat$id))

#.....................................................
#.....................................................
# Combine the data by families ####
# PlackettLuce cannot work with so many combinations 
# so we combine the data to deal with this issue
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

# split by ids (farmer)
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

itemsS_levels <- c("RF", "N1", "N3", "N5","N8", "N10", "N16", "N19",
                   "N32", "N36", "N45", "N46", "N51")

# convert this matrix into a paircomp object
pair <- psychotools::paircomp(pair, labels = itemsS_levels)

# pairwises into grouped rankings
G <- as.grouped_rankings(pair)

G

mod <- PlackettLuce(G)

reference <- "RF"

summary(mod, ref = reference)

re <- reliability(mod, ref = reference)

# dataset with rank features 
unique(datS$id)
datS$gender <- ifelse(grepl("f_", datS$id), "F", "M")
datS$location <- ifelse(grepl("_adet", datS$id), "Adet",
                        ifelse(grepl("_kulumsa", datS$id), "Kulumsa",
                               "Geregera")) 
 
rank_features <- as.data.frame(datS[!duplicated(datS$id), c("gender","location")])
names(rank_features) <- c("Gender", "Location")
rank_features$Gender <- ifelse(rank_features$Gender == "M", "Man", "Woman")
rank_features[1:2] <- lapply(rank_features[1:2], as.factor)

pld <- cbind(G, rank_features)

str(pld)

# empty model (no covariates)
pl <- pltree(G ~ 1, 
             data = pld,
             alpha = 0.01,
             minsize = 10,
             verbose = TRUE)

# model with location and gender 
pl1 <- pltree(G ~ Location + Gender,
              data = pld,
              alpha = 0.01,
              gamma = TRUE,
              minsize = 10,
              verbose = TRUE)

plot(pl1, ref = reference)

# model with gender only
pl2 <- pltree(G ~ Gender,
               data = pld,
               alpha = 0.01,
               minsize = 10,
               verbose = TRUE)

deviance(pl)
deviance(pl1)
deviance(pl2)

# ......................................
# ......................................
# ......................................
# PL model with item covariates ####
# dataset with the genotypes features 
geno_features <- blups
names(geno_features)
geno_features <- geno_features[,-which(grepl("geno|ril", names(geno_features)))]


# # group blups by family
geno_features <- blups
names(geno_features)
geno_features <- geno_features[,-which(grepl("geno|ril", names(geno_features)))]

# group blups by family
genotype_features <- aggregate(. ~ family, data = geno_features, mean)
names(genotype_features)[1] <- "geno"

# the formula to analyse the genotype features
f <- formula(paste("~", paste(names(genotype_features[-c(1)]), collapse = " + ")))

f

# head(geno_features)
# 
# geno_features$family <- factor(geno_features$family, 
#                                levels = union("RF", unique(geno_features$family)))
# 
# genotype_features <- data.frame(genotype = unique(geno_features$family))
# 
# metrics <- names(geno_features)[-1]
# 
# for(i in seq_along(metrics)){
#   
#   form <- paste0(metrics[i], " ~ family")
#   
#   form <- as.formula(form)
#   
#   mod_i <- lm(form, data = geno_features)
#   
#   coeffs <- as.vector(coef(mod_i))
#   
#   coeffs[1] <- 0
#   
#   genotype_features <- cbind(genotype_features, coeffs)
# 
# }
# 
# names(genotype_features)[-1] <- metrics
# 
# genotype_features
# 
# # the formula to analyse the genotype features
# f <- formula(paste("~", paste(metrics, collapse = " + ")))
# 
# f

# PLadmm with location and gender
lvls <- c("RF", "N1", "N3", "N5", "N8",
          "N10", "N16", "N19", "N32", "N36",
          "N45", "N45", "N51")

pl3 <- pltree(G ~ Gender + Location, 
              worth = f,
              data = list(pld,
                          genotype_features),
              minsize = 10,
              verbose = TRUE)

summary(pl3)
plot(pl3)
coef(pl3)

summary(pl3)

nodes <- predict(pl3, type = "node")

node_ids <- sort(unique(nodes))

models <- list()
nobs <- integer()

for (i in seq_along(node_ids)) {
  x <- pld[nodes == node_ids[i], ]
  nobs <- cbind(nobs, nrow(x))
  x <- PlackettLuce(x$G)
  models[[i]] <- x
}

branch <- gosset:::build_tree_branches(pl3)
node <- gosset:::build_tree_nodes(models,
                                  ref = reference,
                                  log = TRUE,
                                  multcomp = FALSE,
                                  node.ids = node_ids,
                                  n.obs = nobs,
                                  levels = rev(itemsS_levels))

tree <- branch / node

tree <-
tree +
theme(axis.text.x = element_text(size = 12, angle = 0))


tree

pl3_coef <- pladmm_coeffs(pl3)

write.csv(pl3_coef, file = "output/pltree-farmers-choices.csv", row.names = FALSE)

ggsave("output/pltree-farmers-choices.png",
       width = 28,
       height = 25,
       units = "cm",
       dpi = 600)

ggsave("output/pltree-farmers-choices.svg",
       width = 28,
       height = 25,
       units = "cm",
       dpi = 600)

capture.output(summary(pl3),
               file = "output/pltree-farmers-choices.txt")

# write outputs

save(pl1, pl2, pl3, file = "output/PL_models.rda")


# # ...........................................
# # ...........................................
# # ...........................................
# # Run within families to get the performance of each genotype
# # PL models by family 
# families <- itemsS_levels[-1]
# 
# family_models <- list()
# 
# for(i in seq_along(families)) {
#   dat_i <- dat[dat$family == families[i], ]
#   
#   # remove the parental line
#   dat_i <- dat_i[!dat_i$geno %in% families[i], ]
#   
#   items_i <- unique(dat_i$geno)
#   
#   # make pairwise comparisons
#   cc <- gosset:::.combn2(items_i, 2)
#   
#   ncol(cc)
#   
#   # split by ids (farmer)
#   dat_i <- split(dat_i, dat_i$id)
#   
#   pair_oa <- lapply(dat_i, function(y) {
#     # get the rankings as pair comparisons
#     # ties are not considered and will be NA's
#     pair <- apply(cc, 2, function(x){
#       
#       # take the first item in the comparison
#       i <- x[1]
#       # and the second one
#       j <- x[2]
#       
#       # combine the rankings for these two items
#       # with i as first and j as the second colunm
#       p <- c(y$oa[y$geno == i], y$oa[y$geno == j])
#       
#       if  (length(p) != 2) {
#         p <- c(0, 0)
#       }
#       
#       # if i is higher than j, add 1, this means that i beats j
#       # if i is lower than j, add -1, this means that j beats i
#       # if none of these options, add NA
#       p <- ifelse(p[1] > p[2], 1, ifelse(p[1] < p[2] , -1, NA))
#       
#     })
#     
#   })
#   
#   pair_oa
#   
#   pair <- do.call("rbind", pair_oa)
#   
#   # convert this matrix into a paircomp object
#   pair <- psychotools::paircomp(pair, labels = items_i)
#   
#   # pairwises into grouped rankings
#   G <- as.grouped_rankings(pair)
#   
#   mod <- PlackettLuce(G)
#   
#   family_models[[i]] <- mod
#   
# }
# 
# # write (update) outputs
# capture.output(print(pl1),
#                cat("\n\n"),
#                summary(pl1), 
#                cat("\n\n\n\n\n\n\n\n"),
#                print(pl2),
#                cat("\n\n"),
#                summary(pl2),
#                cat("\n\n\n\n\n\n\n\n"),
#                print(pl3),
#                cat("\n\n"),
#                summary(pl3),
#                cat("\n\n\n\n\n\n\n\n"),
#                print(pl4),
#                cat("\n\n"),
#                summary(pl4),
#                file = "output/PL_models.text")
# 
# save(pl1, pl2, pl3, pl4, family_models, file = "output/PL_models.rda")
# 
# # get the top items in each family 
# best_families <- data.frame()
# for (i in seq_along(families)) {
#   
#   coef_i <- coef(family_models[[i]], log = FALSE)
#   
#   b_i <- data.frame(family = families[i],
#                     geno = names(rev(sort(coef_i))[1:20]),
#                     worth = as.vector(rev(sort(coef_i))[1:20]),
#                     rank = 1:20)
#   
#   best_families <- rbind(best_families, b_i)
#     
# }
# 
# write.csv(best_families, file = "output/top20_by_family.csv", row.names = FALSE)
# 
# 
# 
# 
