# Written in Tidy mostly for convenience,
# on slow(er) machines and/or larger corpora probably better in data.table
library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpointdensity)
library(collostructions)
library(zipfR)

set.seed(667)
"%!in%" <- function(x, y) !("%in%"(x, y))

# {{{ Basic table

# --- Data compiled and extracted with the help of cwb tools
# token count
corpsize <- 542341719

# import
# --- cwb-scan-corpus COCA-S lemma+0 pos+0 ?pos+0=/n.*/ > "coca_nouns"
df_raw <- read.delim("data/coca_nouns",
  col.names = c("counts", "type", "pos"))

# remove types containing special characters except hivens
df_clean <- df_raw[grepl("^[a-zA-Z]+\\-?[a-zA-Z]+$", df_raw$type), ]

# find types attracted to proper noun tags for later filtering
df_clean$np <- ifelse(grepl("^np", df_clean$pos), yes = "np", no = "nn")
df_proper <- aggregate(counts ~ type + np, data = df_clean, sum) %>%
  pivot_wider(id_cols = type, names_from = np, values_from = counts)
df_proper[is.na(df_proper)] <- 0

df_proper <- mutate(df_proper, freq = np + nn, rel_np = np / freq)

# FIXME: using collostruction::collex. more of a dirty hack
proper_clx <- collex(as.data.frame(df_proper[c("type", "np", "freq")]),
       corpsize, str.dir = TRUE)

proper_nouns <- proper_clx %>% filter(STR.DIR > 0) %>% pull(COLLEX)

# aggregate plural and singular pos
df_clean$pl <- ifelse(grepl("^n.{1,2}2", df_clean$pos), yes = "pl", no = "sg")
df_aggr <- aggregate(counts ~ type + pl, data = df_clean, sum) %>%
  pivot_wider(id_cols = type, names_from = pl, values_from = counts)
df_aggr[is.na(df_aggr)] <- 0

# freqs
df_aggr <- mutate(df_aggr, freq = pl + sg, rel_pl = pl / freq)

# fix column type
nums <- c("pl", "sg", "freq")
df_aggr[nums] <- lapply(df_aggr[nums], as.numeric)

# }}}

# {{{ Add dispersion
# import dispersion value (Gries' dp) and add to table
# data extracted with dispersion scripts on same data set
# https://github.com/alex-raw/dispersion

df <- read.delim("data/coca_dispersion_nn", header = TRUE, quote = NULL,
  col.names = c("type", "counts", "dp")) %>%
  left_join(df_aggr, disp, by = "type") %>%
  select(!counts)

# }}}

# {{{ Sg/Pl: Collex. and associations
# functions to add association measures; FIXME: still hacky and slow

clx <- function(x) {
  y <- as.data.frame(df[c("type", x, "freq")])
  collex(y, corpsize, str.dir = TRUE)
}

add_cxn <- function(path, name) {
  freqs <- read.delim(path, header = FALSE,
    col.names = c(name, "type"))
  freqs["type"] <- sub("\\s.*", "", freqs[, "type"])  # remove ids from types
  clx <- left_join(df, freqs, by = "type")
  clx[is.na(clx)] <- 0
  clx
}

add_ass <- function(x, clx) {
  colname1 <- paste0("as_", x)
  colname2 <- paste0("sig_", x)
  names(clx)[names(clx) == "COLLEX"] <- "type"
  names(clx)[names(clx) == "STR.DIR"] <- colname1
  names(clx)[names(clx) == "SIGNIF"] <- colname2
  left_join(df, clx[c("type", colname1, colname2)])
}

# }}}

# {{{ Modifying collocation

# --- data queried on the same cwb version using cqp
# [lemma = "a"] [] [word = "of" %c] [pos != "at1.*" & pos != "y.*"]? \
# [pos = "nn1"] [pos != "n.*"]
df <- add_cxn("data/x_of", "x_of")

# definite article
df <- add_cxn("data/det", "det")

# [lemma = "pair"] [word = "of" %c] [pos != "nn.*" & pos != "y.*"]* \
# [pos = "nn.*"] [pos != "n.*"] within s
df <- add_cxn("data/piece_of", "piece_of")
df <- add_cxn("data/lot_of", "lot_of")
df <- add_cxn("data/pair_of", "pair_of")

# [pos = "dd(_.*)?"] [pos != "nn.*"]? [pos = "n.{1,2}1"] [pos != "nn.*"]
df <- add_cxn("data/some_sg", "some_sg")
df <- add_cxn("data/some_pl", "some_pl")

# combine features
df <- add_cxn("data/count_sg", "count_sg")
df <- add_cxn("data/mass", "mass")

# }}}

# {{{ tables

# plurale tantum
fig_plurale_tb <- df %>%
  filter(dp < .5959, rel_pl > .97,
    type %!in% c("police", "personnel", "cattle", "fridays", "thursdays")) %>%
  select(type, pl, sg, freq) %>%
  arrange(desc(freq)) %>%
  head(20)

fig_singulare_tb <- df %>%
  filter(dp < .5959, sg / freq > .97 & sg / freq < 1,
    type %!in% c("saturdays", "us", "work", "fridays", "thursdays", "sundays")) %>%
    select(type, pl, sg, freq) %>%
    arrange(desc(freq)) %>%
    head(20)

# }}}

# {{{ *Pair of* collostruction analysis

# add association
pair_collex <- clx("pair_of")
pair_collex$COLLEX <- reorder(pair_collex$COLLEX, pair_collex$STR.DIR)
df <- add_ass("pair_of", pair_collex)

# coding top 50 for lexical fields
fields <- list(
  footwear = c("shoe", "boot", "sock", "sneakers", "ski", "heel", "slipper",
               "sandal", "slack", "skate", "loafer", "pump", "flip-flop",
               "moccasin", "stocking", "snowshoe", "cleat"),
  handwear = c("glove", "handcuff", "mitten"),
  eyewear = c("binoculars", "glass", "sunglasses", "eyeglass", "spectacle",
              "goggle", "bifocals"),
  earwear = c("earring", "headphone", "earplug", "shade"),
  legwear = c("jeans", "pant", "shorts", "underwear", "panties", "trousers",
              "underpants", "sweatpants", "tights", "pantyhose", "overall",
              "pajamas", "legging", "khakis", "brief", "boxer", "breech",
              "dungaree", "sweat", "suspenders", "coveralls", "undershorts"),
  tool = c("scissors", "pliers", "heddles", "tweezers", "tongs", "shear",
           "forceps"),
  `body part` = c("eye", "hand", "leg", "wing", "foot", "arm", "lip", "ear",
                  "lung", "breast"),
  other = c("dumbbell",  "chromosome", "door", "headlight", "chopsticks",
            "dice", "eagle", "variate",  "chair", "vase", "figure", "hose",
            "ticket", "throw", "touchdown", "particle", "speaker", "opposite",
            "pistol", "ace", "chip"),
  animate = c("twin", "man", "woman", "star", "lover", "bird", "student", "owl",
              "mule")
)

# FIXME: attrocious! I'm embarrassed. Will be fixed soon
collex_raw <- pair_collex %>% head(50) %>%
  mutate(field = ifelse(COLLEX %in% fields[["footwear"]], "footwear", ""),
         field = ifelse(COLLEX %in% fields[["handwear"]], "handwear", field),
         field = ifelse(COLLEX %in% fields[["eyewear"]], "eyewear", field),
         field = ifelse(COLLEX %in% fields[["earwear"]], "earwear", field),
         field = ifelse(COLLEX %in% fields[["legwear"]], "legwear", field),
         field = ifelse(COLLEX %in% fields[["tool"]], "tool", field),
         field = ifelse(COLLEX %in% fields[["body part"]], "body part", field),
         field = ifelse(COLLEX %in% fields[["other"]], "other", field),
         field = ifelse(COLLEX %in% fields[["animate"]], "animate", field),
         field_redux = ifelse(field %in% c("handwear", "footwear"), "hand/footwear", field),
         field_redux = ifelse(field %in% c("earwear", "eyewear"), "ear/eyewear", field_redux)
  )

fig_collex <- collex_raw %>%
  ggplot(aes(x = STR.DIR, y = COLLEX)) +
    scale_x_log10(name = "Association strength (logl)") + ylab("type") +
    geom_point()

fig_collex_col <- collex_raw %>%
  ggplot(aes(x = STR.DIR, y = COLLEX, color = field_redux)) +
    scale_x_log10(name = "Association strength (logl)") + ylab("type") +
    scale_color_discrete(name = "semantic field") +
    geom_point()

# }}}

# {{{ *Pair of* on pl continuum
# tupelize

# FIXME: attrocious! I'm embarrassed. Will be fixed soon
df_pair <- df %>%
  mutate(field = ifelse(type %in% fields[["footwear"]], "footwear", ""),
           field = ifelse(type %in% fields[["handwear"]], "handwear", field),
           field = ifelse(type %in% fields[["eyewear"]], "eyewear", field),
           field = ifelse(type %in% fields[["earwear"]], "earwear", field),
           field = ifelse(type %in% fields[["legwear"]], "legwear", field),
           field = ifelse(type %in% fields[["tool"]], "tool", field),
           field = ifelse(type %in% fields[["body part"]], "body part", field),
           field = ifelse(type %in% fields[["other"]], "other", field),
           field = ifelse(type %in% fields[["animate"]], "animate", field),
           field_redux = ifelse(field %in% c("handwear", "footwear"), "hand/footwear", field),
           field_redux = ifelse(field %in% c("earwear", "eyewear"), "ear/eyewear", field_redux)
  ) %>%
  filter(dp < .599, dp > 0, as_pair_of > 0) %>%
  mutate(rel = pair_of / freq,
         ass = dense_rank(desc(as_pair_of)),
         labels = ifelse(ass <= 100, type, ""),
         dp_norm = (dp - min(dp)) / (max(dp) - min(dp)))

pair_plain <- df_pair %>% ggplot(aes(x = rel_pl, y = rel)) +
  scale_y_continuous(name = "a pair of") +
  scale_x_continuous(name = "plural")

fig_pair_plain <- pair_plain + geom_point()

fig_pair_col <- pair_plain +
    geom_point(aes(color = as_pair_of,
             alpha = pair_of,
             size = 1 - dp_norm,
                   )) +
    scale_y_continuous(name = "a pair of") +
    scale_color_gradient2(name = "association (logl)",
                           trans = "log10") +
    scale_alpha(name = "frequency in construction\n(log)",
                trans = "log10", range = c(.5, 1)) +
    scale_size_continuous(name = "dispersion\n1 - dp_norm", trans = "sqrt")
    ggtitle("tiea")


fig_pair_cline <- fig_pair_col +
    geom_smooth(show.legend = FALSE, se = FALSE) +
    geom_text_repel(aes(label = labels), alpha = 2, size = 4)

fig_pair_field <- df_pair %>%
  filter(type %!in% c("tongs", "shade", "breech"), field_redux != "") %>%
  ggplot(aes(x = rel_pl, y = rel)) +
  scale_y_continuous(name = "a pair of") +
  scale_x_continuous(name = "plural") +
    geom_point(aes(color = field_redux,
             size = 1 - dp_norm,
                   )) +
    scale_y_continuous(name = "a pair of") +
    scale_color_discrete(name = "semantic field") +
    scale_alpha(name = "frequency in construction\n(log)",
                trans = "log10", range = c(.5, 1)) +
    scale_size_continuous(name = "dispersion\n1 - dp_norm", trans = "sqrt") +
    geom_smooth(show.legend = FALSE, se = FALSE) +
    geom_text_repel(aes(label = labels), alpha = 2, size = 3) +
    ggtitle("Top 100 attracted collexemes of (a) pair(s) of")

# mass nouns

df_final <- df %>%
  mutate(rel_nmass = (mass + det + x_of) / freq,
         rel_ncount = (count_sg + pl) / freq,
         labels = ifelse(freq > 135000, type, ""),
         ) %>%
  filter(type %!in% proper_nouns, dp < .598 & dp > 0, freq > 100)

fig_mass_raw <- df_final %>% ggplot(aes(x = rel_ncount, y = rel_nmass)) +
    scale_y_continuous(name = "uncountable or neutral determiner") +
    scale_x_continuous(name = "-s or countable determiner")

fig_mass <- fig_mass_raw + geom_pointdensity(adjust = .04, alpha = .3) +
  scale_color_viridis_c()

fig_mass_text <- fig_mass_raw + geom_density_2d(color = "black", alpha = .2) +
  geom_text(aes(label = labels), size = 4)


# {{{ Productivity

growth <- function(x) {
  tfl(x) %>% tfl2spc() %>% lnre("fzm", ., exact = TRUE) %>% lnre.vgc(N = 1:N(.))
}

lnre_sg <- tfl(df$sg) %>% tfl2spc() %>% lnre("fzm", ., exact = TRUE)
vgc_sg <- lnre.vgc(lnre_sg, N = 1:N(lnre_sg))

lnre_pl <- tfl(df$pl) %>% tfl2spc() %>% lnre("fzm", ., exact = TRUE)
vgc_pl <- lnre.vgc(lnre_pl, N(vgc_sg))

fig_piece <- growth(df$piece_of)
fig_pair <- growth(df$pair_of)
fig_lot <- growth(df$lot_of)

fig_some_sg <- growth(df$some_sg)
fig_some_pl <- growth(df$some_pl)

# }}}

# export figures
figures <- ls()[grepl("fig_", ls())]
save(list = figures, file = "data.RData")

# rm(df_aggr, df_clean, df_raw, disp, nums)
# gc()
