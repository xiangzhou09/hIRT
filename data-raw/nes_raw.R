# rm(list = ls())
#
# library("pryr")
# library("dplyr")
# library("splines")
# library("ggplot2")
# library("devtools")
#
# load("C:/Users/Xiang/Dropbox/hIRT/code/nes8012.RData")
#
# nes_v1 <- nes_v1 %>% data.frame(nvalid_econ = rowSums(!is.na(dplyr::select(., health_ins7:FS_assistblacks3))),
#                                 nvalid_civil = rowSums(!is.na(dplyr::select(., negro_chan3:blacks_deserve_more5))),
#                                 nvalid_moral = rowSums(!is.na(dplyr::select(., women_role7:abort4))),
#                                 nvalid_foreign = rowSums(!is.na(dplyr::select(., urss_coop7:FS_space3)))) %>%
#   group_by(year) %>% mutate(nvalid_econ_yr_level = max(nvalid_econ),
#                             nvalid_civil_yr_level = max(nvalid_civil),
#                             nvalid_moral_yr_level = max(nvalid_moral),
#                             nvalid_foreign_yr_level = max(nvalid_foreign)) %>%
#   ungroup()
#
# nes_econ <- nes_v1 %>%
#   filter(nvalid_econ >=1, nvalid_econ_yr_level>=3) %>%
#   dplyr::select(year, gender, party, educ, health_ins7:FS_assistblacks3)
#
# nes_civil <- nes_v1 %>%
#   filter(nvalid_civil >=1, nvalid_civil_yr_level>=3) %>%
#   dplyr::select(year, gender,  party, educ, negro_chan3:blacks_deserve_more5)
#
# nes_moral <- nes_v1 %>%
#   filter(nvalid_moral >=1, nvalid_moral_yr_level>=3) %>%
#   dplyr::select(year, gender,  party, educ, women_role7:abort4)
#
# nes_foreign <- nes_v1 %>%
#   filter(nvalid_foreign >=1, nvalid_foreign_yr_level>=3) %>%
#   dplyr::select(year, gender,  party, educ, urss_coop7:FS_space3)
#
# nes_econ2012 <- subset(nes_econ, year == 2012) %>%
#   dplyr::select(-year, -FS_aids3, -FS_aidcollege3, -FS_homeless3,
#                 -FS_foodstamps3, -FS_assistblacks3) %>%
#   mutate(party = factor(party, levels = 1:3,
#                         labels = c("Democrat", "independent",
#                                    "Republican")))
#
# nes_racial2008 <- subset(nes_civil, year == 2008) %>%
#   dplyr::select(gender, party, educ,
#                 hard_blacks5, no_favor_blacks5,
#                 blacks_try_harder5, blacks_deserve_more5) %>%
#   mutate(party = factor(party, levels = 1:3,
#                         labels = c("Democrat", "independent",
#                                    "Republican")))
#
# summary(nes_racial2008)
#
# use_data(nes_racial2008, overwrite = TRUE)
