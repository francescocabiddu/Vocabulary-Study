# Vocabulary analysis - 15/06/19
# Francesco Cabiddu, cabiddu@hotmail.it

# Shortcuts sections (see bottom page http://tiny.cc/hau68y):
# Collapse one section: mac (Alt+Cmd+L), win (Alt+L)
# Expand one section: mac (Shift+Alt+Cmd+L), win (Shift+Alt+L)
# Collapse all sections: mac (Alt+Cmd+O), win (Alt+O)
# Expand all sections: mac (Shift+Alt+Cmd+O), win (Shift+Alt+O)

# load libraries
lib <- c("magrittr", "tidyverse",
         "beepr", "data.table", "fastmatch",
         "mailR", "stringdist")
lapply(lib, require, character.only = TRUE)
rm(lib)

#### homemade funs ####
match <- fmatch # drop-in replacement; see http://tiny.cc/9clb8y

phonemic_match <- function(df, ref_var_phon, ref_var_word, match_type) {
  # function that matches a column of orthographic words with their correspondent
  # phonetic transcription, working out compounds using single reference words
  # or using already prepared reference compounds
  if (match_type == "single") {
    manch_tok_trans <- df %>% filter(!str_detect(phon, "NA"))
    
    manch_tok_not_trans <- df %>%
      filter(str_detect(phon, "NA")) %>%
      mutate(phon = word %>%
               str_split("\\+") %>%
               sapply(function(x) {
                 ref_var_phon[match(x, ref_var_word)] %>%
                   paste(collapse = "_")
               }))
    
    rbind(manch_tok_trans, manch_tok_not_trans)
  } else {
    manch_tok_trans <- df %>% filter(!str_detect(phon, "NA"))
    
    manch_tok_not_trans <- df %>%
      filter(str_detect(phon, "NA")) %>%
      mutate(phon = word %>%
               (function(x) {
                 ref_var_phon[match(x, ref_var_word)]
               })) %>%
      mutate(phon = phon %>%
               (function(x) {
                 x[which(is.na(x))] <- "NA"
               }))
    
    
    rbind(manch_tok_trans, manch_tok_not_trans)
  }
} 

paste5 <- function(..., sep = " ", collapse = NULL, na.rm = F) {
  # paste escaping NAs; see http://tiny.cc/yspb8y
  if (na.rm == F)
    paste(..., sep = sep, collapse = collapse)
  else
    if (na.rm == T) {
      paste.na <- function(x, sep) {
        x <- gsub("^\\s+|\\s+$", "", x)
        ret <- paste(na.omit(x), collapse = sep)
        is.na(ret) <- ret == ""
        return(ret)
      }
      df <- data.frame(..., stringsAsFactors = F)
      ret <- apply(df, 1, FUN = function(x) paste.na(x, sep))
      
      if (is.null(collapse))
        ret
      else {
        paste.na(ret, sep = collapse)
      }
    }
}

notify_me <- function(subj = "R: Come back to me", pass, recipient) {
  # send an email to notify a step has finished
  send.mail(from = "fcabiddur@gmail.com", # change sender
            to = recipient, # change recipient email address
            subject= subj,
            body = "DO NOT REPLY",
            smtp = list(host.name = "smtp.gmail.com", port = 465, 
                        user.name="fcabiddur@gmail.com", # change sender
                        passwd=pass, # password sender email address
                        ssl=TRUE),
            authenticate = TRUE,
            send = TRUE)
}

tok_to_typ <- function(df) {
  # it takes unique types by stage not present in previous stages
  df %>%
    select(-gram, -word, -hour, -half_hour) %>%
    group_by(set, id, stage) %>%
    distinct(phon, .keep_all = TRUE) %>% # take unique types in each stage
    ungroup %>%
    filter(!phon %in% manch_typ_plu$phon) %>% # filter plural nouns out
    group_by(set, id, phon) %>%
    filter(stage == min(stage)) %>% # take only first presentation of a type
    ungroup %>%
    arrange(set, id, stage, phon) 
}

rmse <- function(error) {
  sqrt(mean(error^2))
}

assign_tertiles <- function(var, df_ref_var, range, subset_phonemic) {
  var %>%
    (function(x) {
      # tertiles based on word types present in manchester corpus
      # or in manchester corpus subset depending on length/frequency range
      tertiles <- manch_typ %>%
        distinct(phon, .keep_all = TRUE)
      
      if (subset_phonemic == TRUE) {
        tertiles %<>% filter(phonemic_len %in% range)
      } else if (subset_phonemic == FALSE) {
        tertiles %<>% filter(freq_ter %in% range)
      }
      
      tertiles %<>%
        (function(y) {
          y[[df_ref_var]] %>%
            quantile(probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
        })
      
      ifelse(x < tertiles[2], 1,
             ifelse(x >= tertiles[2] & x < tertiles[3], 2,
                    3))
      
    })
}

convert_phonemes <- function(df) {
  # create a column of converted phonemes to single characters for each word
  df %>%
    mutate(phon_converted = phon %>%
             str_split("_") %>%
             sapply(function(i) {
               paste0(phonemes_converted[i], collapse="")
             }))
}

compute_neighbours <- function(target_converted, refts_converted, reft_var_converted) {
  ref_words <- refts_converted
  ref_words %<>%
  {.[which(. != target_converted)]} # exclude identical words
  
  stringdist(target_converted, ref_words, method = "lv") %>% # apply levenstein distance
  {which(. == 1)} %>% # select only distance == 1
  {ref_words[.]} %>%
  {.[nchar(.) != (nchar(target_converted) == 1)]} %>% # exclude monophonemic words if target is monophonemic
    length
  # delete line above and add code below to see actual neighbours
  #{which(spok_bnc_typ[[reft_var_converted]] %in% .)} %>%
  #{spok_bnc_typ$phon[.]}
}

adj.poss <- function(x) {
  # all possible adjacent combinations of phonemes in a word, including the word itself
  n <- length(x)
  if(n == 1L) return(NA) # if a word is monophonemic return NA
  idx <- expand.grid(start = 1L:n, len = 2L:(n))
  idx$end <- idx$start + idx$len - 1L
  idx <- idx[idx$end <= n, ]
  Map(function(start, end) x[start:end], idx$start, idx$end)
}

add_sequences <- function(df) {
  # add 2to5 phoneme sequences for each word in a dataframe
  df %>%
    mutate(seqs = phon %>%
             str_split("_") %>%
             sapply(function(word) {
               word %>%
                 adj.poss %>%
                 sapply(paste5, collapse ="_", na.rm = TRUE) %>%
                 unique
             })) %>%
    mutate(seq2 = sapply(seqs, str_subset, 
                         "^[A-Z]+_[A-Z]+$"),
           seq3 = sapply(seqs, str_subset, 
                         "^[A-Z]+_[A-Z]+_[A-Z]+$"),
           seq4 = sapply(seqs, str_subset, 
                         "^[A-Z]+_[A-Z]+_[A-Z]+_[A-Z]+$"),
           seq5 = sapply(seqs, str_subset, 
                         "^[A-Z]+_[A-Z]+_[A-Z]+_[A-Z]+_[A-Z]+$"))
  
  
}

count_seqs_types <- function(df, df_count, seq_count_name, seq_count) {
  # for each type in a dataframe
  # compute the proportion of types sharing seq (2 to target length) with target
  df %>%
    mutate(!!seq_count_name := sapply(1:nrow(.), function(i) {
      lapply(df_count %>% 
               filter(id == df$id[i]) %>% # by child
               filter(stage %in% 0:df$stage[i]) %>% # up to a target stage
               filter(phon != df$phon[i]) %>% # exclude target from reference df
               {.[[seq_count]]},
             function(x) {
               c(x, df[[seq_count]][[i]])
             }) %>%
        # duplicated returns the positional value for each duplicated element of x
        # modify here if the matched word needs to be counted once for each shared sound sequence in the target
        lapply(duplicated) %>%
        sapply(function(x) {
          any(x > 0) 
        }) %>%
        {sum(.) / length(.)}
    })) %>%
    (function(x) {
      x[[seq_count_name]][sapply(x[[seq_count]], is_empty)] <- list(NA)
      x
    })
}

count_seqs_tokens <- function(df, df_count, seq_count_name, seq_count) {
  # for each type in a dataframe
  # compute the proportion of tokens sharing seq (2 to target length) with target
  df %>%
    (function(DF) {
      proportion_tokens <- sapply(1:nrow(DF), function(i) {
        df_filtered <- df_count %>% 
          filter(id == DF$id[i]) %>%
          filter(stage %in% 0:DF$stage[i]) %>%
          filter(phon != DF$phon[i]) 
        
        df_freq <- df_filtered %>%
          group_by(phon) %>%
          summarise(N_word = n()) %>%
          ungroup
        
        df_filtered <- df_filtered %>%
          left_join(., df_freq, by = "phon") %>%
          distinct(phon, .keep_all = TRUE)
        
        lapply(df_filtered %>%
        {.[[seq_count]]},
        function(x) {
          c(x, DF[[seq_count]][[i]])
        }) %>%
          lapply(duplicated) %>%
          sapply(function(x) {
            any(x > 0) 
          }) %>%
          {. * df_filtered$N_word} %>%
          {sum(.) / sum(df_filtered$N_word)}
      })
      
      DF %>%
        mutate(!!seq_count_name := proportion_tokens) %>%
        (function(x) {
          x[[seq_count_name]][sapply(x[[seq_count]], is_empty)] <- list(NA)
          x
        })
    })
}

#### import and clean data ####
  #### pronunciation dictionaries ####
CMU_DICT <- read_lines("CMU-Lexicon-words-to-phonemes.txt") %>%
{tibble(word_phon = .)} %>%
  mutate(word_phon = word_phon %>%
           str_remove(" $")) %>%
  separate(word_phon, c("word", "phon"), " ") %>%
  arrange(word) %>%
  distinct(word, .keep_all = TRUE) %>%
  mutate(phon = phon %>%
           str_replace_all("UU", "AH")) %>%
  mutate(word = word %>%
           str_replace_all("-", "\\+")) %>%
  (function(x) {
    list(words = x %>% filter(!str_detect(word, "\\+")),
         compounds = x %>% filter(str_detect(word, "\\+")))
  })

# a collection of maternal types is also used, as some frequent orthographic types 
# were manually transcribed; this collection includes compounds created from CMU_DICT single words
# plus other frequent maternal types manually transcribed
MOT_DICT <- read_lines("Manchester-corpus-mot-filt-orth-and-phon.txt") %>%
{tibble(line = .)} %>%
  mutate(line = line %>%
           str_remove_all("^[0-9]+ ") %>%
           str_split("\\) \\(")) %>%
  (function(x) {
    tibble(word_phon = x$line %>% unlist)
  }) %>%
  mutate(word_phon = word_phon %>%
           str_remove_all("^ | $") %>%
           str_remove_all("^\\(|\\)$")) %>%
  separate(word_phon, c("word", "phon"), " ") %>%
  mutate(phon = phon %>%
           str_remove_all("0|1|2")) %>%
  distinct(word, .keep_all = TRUE) %>%
  arrange(word) %>%
  filter(phon != "UNKNOWN") %>%
  (function(x) {
    list(words = x %>% filter(!str_detect(word, "\\+")),
         compounds = x %>% filter(str_detect(word, "\\+")))
  })

  #### manchester corpus ####
# import orthographic corpus to 
# (1) list plural nouns, 
# (2) extract mothers and children's utterances
MANCHESTER <- read_lines("manchester_corpus.txt")

# remove punctuation/not useful markers
manch_tok <- MANCHESTER %>% 
{tibble(line = .)} %>%
  mutate(line = line %>%
           str_remove(" $")) %>%
  separate(line, c("set_utt", "gram", "id_hour_halfhour"), "\\) ") %>%
  separate(set_utt, c("set", "utt"), ": ") %>%
  mutate(gram = gram %>%
           str_remove("^\\(%MOR: "),
         set = set %>%
           str_remove("^\\(\\(\\*"),
         utt = utt %>%
           str_remove(" [.!?]{1}$"),
         gram = gram %>%
           str_remove(" [.!?]{1}$")) %>%
  # utterances with comments within utt excluded (not regular comment pattern!) (~500utts)
  filter(!str_detect(utt, "[!]+")) %>% 
  filter(!str_detect(utt, "=")) %>% 
  mutate(utt = utt %>%
           str_remove_all("-0[^ ]+|^0[^ ]+ ")) %>%
  mutate(utt = utt %>%
           str_remove_all(" 0[^ ]+")) %>%
  filter(!utt %in% c("0", ".")) %>%
  mutate(utt = utt %>%
           str_replace("0[^ ]+ ", " ")) %>%
  mutate(utt = utt %>%
           str_remove(" \\..*$")) %>%
  mutate(utt = utt %>%
           str_remove("^ | $"),
         gram = gram %>%
           str_remove("^ | $")) %>%
  mutate(utt = utt %>%
           str_remove_all(" &[^ ]+")) %>%
  mutate(utt = utt %>%
           str_remove_all("&[^ ]+ ")) %>%
  mutate(utt = utt %>%
           str_remove("^&[^ ]+$")) %>%
  filter(utt != "?") %>%
  mutate(utt = utt %>%
           str_remove("\\?$")) %>%
  mutate(utt = utt %>%
           str_remove("@L$")) %>%
  mutate(utt = utt %>%
           str_replace_all("@L ", "XXX") %>%
           str_replace_all("XXX", " ")) %>%
  mutate(utt = utt %>%
           str_replace_all("@L@", "@")) %>%
  mutate(utt = utt %>%
           str_replace_all("  ", " ")) %>%
  mutate(utt = utt %>%
           str_remove("^ | $")) %>%
  mutate(utt = utt %>%
           str_split(" "),
         gram = gram %>%
           str_split(" ")) %>%
  filter(set %in% c("CHI", "MOT")) # filter for children and mothers

manch_tok %<>%
  mutate(utt = utt %>%
           (function(x) {
             # collection of corrected orthographic words
             extra_lexicon <- read_tsv("extra-lexicon-words.txt", col_names = FALSE) %>%
               rename(word_correction = X1) %>%
               separate(word_correction, c("word", "correction"), " ") %>%
               distinct(word, .keep_all = TRUE) %>%
               arrange(word) %>%
               filter(!str_detect(correction, "_")) %>% # exclude wrong corrections (in phonetic)
               (function(x) {
                 df <- x
                 
                 df %>% # exclude wrong corrections in external file!
                   filter(!word %in% (x %>%
                                        filter(!str_detect(correction, "\\+")) %>%
                                        {.$word[which(!.$correction %in% unique(c(MOT_DICT$words$word, 
                                                                                  MOT_DICT$compounds$word,
                                                                                  CMU_DICT$words$word, 
                                                                                  CMU_DICT$compounds$word)))]}))
               })
             
             # replace manchester typos using external file
             lapply(x, function(y) {
               y %>%
                 paste5(
                   extra_lexicon$correction[match(y, extra_lexicon$word)], na.rm = TRUE)
             })
           }) %>%
           lapply(str_remove, "^[^ ]+ "))

# exclude utterances containing words not in dictionaries
manch_tok %<>%
  mutate(utt_split_compounds = utt %>%
           lapply(function(x) {
             str_split(x, "\\+")
           })) %>%
  mutate(utt_MOT_DICT_compounds = utt %>%
           lapply(function(x) {
             MOT_DICT$compounds$phon[match(x, MOT_DICT$compounds$word)]
           }),
         utt_CMU_DICT_compounds = utt %>%
           lapply(function(x) {
             CMU_DICT$compounds$phon[match(x, CMU_DICT$compounds$word)]
           })) %>%
  mutate(utt_MOT_DICT_single = utt_split_compounds %>%
           sapply(function(x) {
             sapply(x, function(y) {
               MOT_DICT$words$phon[match(y, MOT_DICT$words$word)] %>%
                 paste(collapse = "_")
             })
           }),
         utt_CMU_DICT_single = utt_split_compounds %>%
           sapply(function(x) {
             sapply(x, function(y) {
               CMU_DICT$words$phon[match(y, CMU_DICT$words$word)] %>%
                 paste(collapse = "_")
             })
           })) %>%
  (function(df) {
    for (i in seq_along(df$set)) {
      df$dict_check[i] <- rbind(df$utt_MOT_DICT_compounds[[i]],
                                df$utt_MOT_DICT_single[[i]],
                                df$utt_CMU_DICT_compounds[[i]],
                                df$utt_CMU_DICT_single[[i]]
      ) %>% 
        t %>% 
        paste5(na.rm = TRUE) %>%
        str_split(" ") %>%
        sapply(function(x) {
          str_detect(x, "NA") %>%
          {sum(./length(.))}
        }) %>%
        {. == 1} %>% 
        sum
    }
    
    df
  })

manch_tok %<>%
  filter(dict_check == 0) 

# check for mismatch between utterance and grammatical variable
manch_tok %<>%
  mutate(utt_lengths = sapply(utt, length),
         gram_lenghts = sapply(gram, length)) %>%
  mutate(equal_lengths = utt_lengths == gram_lenghts)

manch_tok %<>%
  (function(x) { # filter mismatch for not matching gram, assign NA to gram, reunite
    manch_tok_mismatch <- x %>%
      filter(equal_lengths == FALSE)
    
    for (i in seq_along(manch_tok_mismatch$set)) {
      manch_tok_mismatch$gram[[i]] <- rep(NA, length(manch_tok_mismatch$utt[[i]]))
    }
    
    x %>%
      filter(equal_lengths == TRUE) %>%
      rbind(
        manch_tok_mismatch
      )
  })

manch_tok %<>%
  mutate(id = id_hour_halfhour %>%
           str_extract("^[a-zA-Z]+") %>%
           tolower,
         hour = id_hour_halfhour %>%
           str_extract("[0-9]{2}") %>%
           as.numeric,
         half_hour = id_hour_halfhour %>%
           str_extract("[ab]{1}\\.") %>%
           str_remove("\\.")) %>%
  select(-id_hour_halfhour)

# create 20 stages with approximately same number of utterances
manch_tok %<>%
  group_by(set, id) %>%
  mutate(n_utt = n()/20,
         utt_cum = 1) %>%
  mutate(utt_cum = cumsum(utt_cum)) %>%
  ungroup %>% 
  mutate(stage = ifelse(utt_cum < n_utt, 1,
                        ifelse(utt_cum >= n_utt & utt_cum < n_utt*2, 2,
                               ifelse(utt_cum >= n_utt*2 & utt_cum < n_utt*3, 3,
                                      ifelse(utt_cum >= n_utt*3 & utt_cum < n_utt*4, 4,
                                             ifelse(utt_cum >= n_utt*4 & utt_cum < n_utt*5, 5,
                                                    ifelse(utt_cum >= n_utt*5 & utt_cum < n_utt*6, 6,
                                                           ifelse(utt_cum >= n_utt*6 & utt_cum < n_utt*7, 7,
                                                                  ifelse(utt_cum >= n_utt*7 & utt_cum < n_utt*8, 8,
                                                                         ifelse(utt_cum >= n_utt*8 & utt_cum < n_utt*9, 9,
                                                                                ifelse(utt_cum >= n_utt*9 & utt_cum < n_utt*10, 10,
                                                                                       ifelse(utt_cum >= n_utt*10 & utt_cum < n_utt*11, 11,
                                                                                              ifelse(utt_cum >= n_utt*11 & utt_cum < n_utt*12, 12,
                                                                                                     ifelse(utt_cum >= n_utt*12 & utt_cum < n_utt*13, 13,
                                                                                                            ifelse(utt_cum >= n_utt*13 & utt_cum < n_utt*14, 14,
                                                                                                                   ifelse(utt_cum >= n_utt*14 & utt_cum < n_utt*15, 15,
                                                                                                                          ifelse(utt_cum >= n_utt*15 & utt_cum < n_utt*16, 16,
                                                                                                                                 ifelse(utt_cum >= n_utt*16 & utt_cum < n_utt*17, 17,
                                                                                                                                        ifelse(utt_cum >= n_utt*17 & utt_cum < n_utt*18, 18,
                                                                                                                                               ifelse(utt_cum >= n_utt*18 & utt_cum < n_utt*19, 19, 20
                                                                                                                                               )))))))))))))))))))) %>%
  select(-n_utt, -utt_cum)

# dataset in long format
manch_tok %<>%
  apply(., 1, function(x) {
    x %>%
      (function(x) {
        tibble(set = x$set,
               id = x$id,
               stage = x$stage,
               word = unlist(x$utt), 
               gram = unlist(x$gram),
               hour = x$hour,
               half_hour = x$half_hour)
      })
  }) %>%
  rbindlist

# convert othographic tokens to MOT_DICT and (secondly) CMU_DICT phonetics
manch_tok %<>%
  mutate(phon = "NA") %>%
  phonemic_match(MOT_DICT$compounds$phon, MOT_DICT$compounds$word, "compound") %>%
  phonemic_match(MOT_DICT$words$phon, MOT_DICT$word$word, "single") %>%
  phonemic_match(CMU_DICT$words$phon, CMU_DICT$words$word, "single") %>%
  phonemic_match(CMU_DICT$compounds$phon, CMU_DICT$compounds$word, "compound") %>%
  select(set:stage, word:gram, phon, hour, half_hour) %>%
  filter(!str_detect(phon, "NA"))

# for each phonemic word type, work out its most frequent gram.
# NOTE. The model output is only in phonemic form.
# Therefore, this grammatical summary considers phonemic types (not orthographic),
# Later, this will allow to assign gram to phonemic model types and exclude plurals
# using the same criterion for every set (mothers, children, model).
manch_typ_gram <- manch_tok %>%
  group_by(phon) %>% 
  summarise(gram = list(table(gram) %>% 
                          sort(decreasing = TRUE)) %>%
              lapply(function(x) {
                names(x[1])
              }) %>%
              (function(x) {
                x[unlist(lapply(x , is.null))] <- NA
                x 
              })) %>%
  ungroup %>%
  mutate(gram = unlist(gram)) 

  #### plurals ####
# save table of plural nouns
manch_typ_plu <- manch_typ_gram %>% 
  filter(str_detect(gram, "^N:+.*PL$|^N\\|+.*PL$")) 

#### Study 1 ####
  #### manchester (types by stage) ####
# for each stage take word types (excluding types present in previous stages)
manch_typ <- tok_to_typ(manch_tok)

  #### model (types by stage) ####
MODEL <- read_lines("Model.txt")

mod_typ <- MODEL %>%
{tibble(line = .)} %>%
  filter(line %>%
           str_detect("vocab-learnt")) %>%
  separate(line, c("id_stage", "phon"), " vocab-learnt ") %>%
  separate(id_stage, c("id", "stage"), " ") %>%
  mutate(id = id %>% tolower,
         stage = stage %>% as.numeric) %>% 
  mutate(n_types = phon %>%
           str_extract("^[0-9]+ ") %>%
           as.numeric) %>%
  mutate(phon = phon %>%
           str_remove("^[0-9]+ ") %>%
           str_split(" ")) %>%
  apply(., 1, function(x) {
    x %>%
      (function(x) {
        tibble(id = x$id,
               stage = x$stage,
               phon = unlist(x$phon))
      })
  }) %>%
  rbindlist %>%
  filter(!phon %in% manch_typ_plu$phon) %>%
  mutate(phon = phon %>%
           str_replace_all("UU", "AH")) %>%
  mutate(id = id %>%
           str_replace("warren", "warr") %>%
           str_replace("dominic", "domin")) %>%
  group_by(id) %>%
  distinct(phon, .keep_all = TRUE) %>%
  ungroup

  #### mothers random samples ####
# create 10 maternal samples randomly sampling the same number of children's tokens
mot_random <- list(tok = list(),
                   typ = list())

for (i in 1:10) {
  set.seed(i)
  
  mot_random$tok[[paste("df", i, sep = "_")]] <-
    manch_tok %>%
    filter(set == "MOT") %>%
    (function(x) {
      x %>%
        left_join(., manch_tok %>%
                    filter(set == "CHI") %>%
                    group_by(id, stage) %>% 
                    summarise(chi_n_tok = n()) %>%
                    ungroup,
                  by = c("id", "stage"))
    }) %>%
    (function(x) {
      matched_x <- x[0,]
      
      while(nrow(x) > 0) {
        x %<>%
          group_by(id, stage) %>%
          sample_frac(0.99, replace = FALSE) %>%
          mutate(mot_n_tok = n()) %>%
          ungroup() %>%
          # small range given (less then 1% of tokens at a stage) to increase the match
          mutate(token_comparison = mot_n_tok - chi_n_tok <= 20)
        
        matched_x <- rbind(matched_x, x %>% filter(token_comparison == TRUE))
        
        x %<>%
          filter(token_comparison == FALSE)
      }
      
      matched_x %>%
        select(set:half_hour)
    })
  
  mot_random$typ[[paste("df", i, sep = "_")]] <-
    mot_random$tok[[paste("df", i, sep = "_")]] %>%
    tok_to_typ
  
} ; rm(i)
  #### syllabic length ####
vowels <- c("AA", "AE", "AH", "AO", "AW", "AY", "EH", "ER", "EY", "IH", "IY", "OW", "OY", "UH", "UW") 

manch_typ %<>%
  mutate(len = phon %>%
           str_split("_") %>%
           sapply(function(x) {
             (x %in% vowels) %>%
               sum
           }),
         phonemic_len = phon %>% # add phonemic length as well for control plots later
           str_split("_") %>%
           sapply(length)) 

mod_typ %<>%
  mutate(len = phon %>%
           str_split("_") %>%
           sapply(function(x) {
             (x %in% vowels) %>%
               sum
           }),
         phonemic_len = phon %>%
           str_split("_") %>%
           sapply(length)) 

mot_random$typ %<>%
  lapply(function(df) {
    df %<>%
      mutate(len = phon %>%
               str_split("_") %>%
               sapply(function(x) {
                 (x %in% vowels) %>%
                   sum
               }),
             phonemic_len = phon %>%
               str_split("_") %>%
               sapply(length)) 
  })

  #### frequency (Spoken BNC) ####
SPOK_BNC <- read_tsv("SpokenBNC.txt") %>%
  na.omit

SPOK_BNC %<>%
  # unite split cases (e.g., have, n't) that are considered together in Manchester corpus and Model
  (function(x) {
    w2 <- x$word %>% str_which("^'|^n't$")
    w1 <- w2 - 1
    
    all_pos <- c(w1, w2) %>% sort
    
    x[-all_pos, ] %>%
      rbind(
        tibble(c5 = matrix(c(x$c5[w1], x$c5[w2]), ncol = 2) %>%
        {split(., seq(nrow(.)))},
        hw = matrix(c(x$hw[w1], x$hw[w2]), ncol = 2) %>%
        {split(., seq(nrow(.)))},
        pos = matrix(c(x$pos[w1], x$pos[w2]), ncol = 2) %>%
        {split(., seq(nrow(.)))},
        w1 = x$word[w1], 
        w2 = x$word[w2]) %>%
          unite("word", c("w1", "w2"), sep = "")
      )
  })

SPOK_BNC %<>%
  mutate(word = word %>%
           str_replace_all("'", "@") %>%
           toupper) 

# convert spoken bnc to phonemic
SPOK_BNC %<>%
  mutate(word = word %>%
           str_replace_all("-", "\\+")) %>%
  left_join(., SPOK_BNC %>%
              distinct(word) %>%
              mutate(word = word %>%
                       str_replace_all("-", "\\+")) %>% 
              arrange(word) %>% 
              mutate(phon = "NA") %>%
              phonemic_match(MOT_DICT$compounds$phon, MOT_DICT$compounds$word, "compound") %>%
              phonemic_match(MOT_DICT$words$phon, MOT_DICT$word$word, "single") %>%
              phonemic_match(CMU_DICT$words$phon, CMU_DICT$words$word, "single") %>%
              phonemic_match(CMU_DICT$compounds$phon, CMU_DICT$compounds$word, "compound") %>% 
              filter(!str_detect(phon, "NA")), by = "word") %>%
  na.omit 

# create spoken bnc typ adding frequency columns
spok_bnc_typ <- SPOK_BNC %>%
  group_by(phon) %>%
  summarise(freq = n() / 1000000,
            # freq + 1 before taking log10 to avoid negative numbers, as in IPhOD formula
            freq_log10 = log10((n() / 1000000) + 1)) %>%
  ungroup

# delete SPOK_BNC to get space in ws
rm(SPOK_BNC)

# assign frequency to children, model, mothers
# words not in spoken bnc have no frequency value therefore are excluded
manch_typ %<>%
  left_join(., spok_bnc_typ %>%
              select(-freq_log10), by = "phon")

mod_typ %<>%
  left_join(., spok_bnc_typ %>%
              select(-freq_log10), by = "phon")

mot_random$typ %<>%
  lapply(function(df) {
    df %<>%
      left_join(., spok_bnc_typ %>%
                  select(-freq_log10), by = "phon")
  })

# assign frequency tertiles
manch_typ %<>%
  mutate(freq_ter = freq %>%
           assign_tertiles("freq", 1:100, TRUE),
         freq_ter_26 = freq %>%
           assign_tertiles("freq", 2:6, TRUE)) %>%
  group_by(phonemic_len) %>%
  mutate(freq_ter_each_len = freq %>%
           assign_tertiles("freq", phonemic_len, TRUE)) %>%
  ungroup

mod_typ %<>%
  mutate(freq_ter = freq %>%
           assign_tertiles("freq", 1:100, TRUE),
         freq_ter_26 = freq %>%
           assign_tertiles("freq", 2:6, TRUE)) %>%
  group_by(phonemic_len) %>%
  mutate(freq_ter_each_len = freq %>%
           assign_tertiles("freq", phonemic_len, TRUE)) %>%
  ungroup

mot_random$typ %<>%
  lapply(function(df) {
    df %<>%
      mutate(freq_ter = freq %>%
               assign_tertiles("freq", 1:100, TRUE),
             freq_ter_26 = freq %>%
               assign_tertiles("freq", 2:6, TRUE)) %>%
      group_by(phonemic_len) %>%
      mutate(freq_ter_each_len = freq %>%
               assign_tertiles("freq", phonemic_len, TRUE)) %>%
      ungroup
  })

  #### ND (Spoken BNC) ####
# convert phonemes into single characters (for levenstein distance in compute_neighbours fun)
phonemes_converted <- c(letters, LETTERS)[1:39] %>%
  (function(x) {
    names(x) <- manch_tok$phon %>%
      str_split("_") %>%
      unlist %>%
      unique %>%
      sort
    
    x
  })

spok_bnc_typ %<>%
  (function(df) {
    ref_df <- spok_bnc_typ %>%
      convert_phonemes
    
    df_converted <- df %>%
      convert_phonemes
    
    df_converted %>%
      mutate(nd = sapply(phon_converted, function(x) {
        spok_bnc_lengths <- ref_df$phon_converted %>% nchar
        
        abs(spok_bnc_lengths - nchar(x)) %>%
        {compute_neighbours(x, 
                            ref_df$phon_converted[. %in% 0:1], "phon_converted")} 
      })) %>%
      select(-phon_converted)
  })

# assign nd to children, model, mothers
# words not in spoken bnc are excluded to use identical samples across measures (freq, nd, pp and control plots)
manch_typ %<>%
  left_join(., spok_bnc_typ %>%
              select(phon, nd), by = "phon")

mod_typ %<>%
  left_join(., spok_bnc_typ %>%
              select(phon, nd), by = "phon")

mot_random$typ %<>%
  lapply(function(df) {
    df %<>%
      left_join(., spok_bnc_typ %>%
                  select(phon, nd), by = "phon")
  })


# assign nd tertiles
manch_typ %<>%
  mutate(nd_ter = nd %>%
           assign_tertiles("nd", 1:100, TRUE),
         nd_ter_26 = nd %>%
           assign_tertiles("nd", 2:6, TRUE)) %>%
  group_by(phonemic_len) %>%
  mutate(nd_ter_each_len = nd %>%
           assign_tertiles("nd", phonemic_len, TRUE)) %>%
  ungroup %>%
  group_by(freq_ter) %>%
  mutate(nd_ter_each_freq = nd %>%
           assign_tertiles("nd", freq_ter, FALSE)) %>%
  ungroup

mod_typ %<>%
  mutate(nd_ter = nd %>%
           assign_tertiles("nd", 1:100, TRUE),
         nd_ter_26 = nd %>%
           assign_tertiles("nd", 2:6, TRUE)) %>%
  group_by(phonemic_len) %>%
  mutate(nd_ter_each_len = nd %>%
           assign_tertiles("nd", phonemic_len, TRUE)) %>%
  ungroup %>%
  group_by(freq_ter) %>%
  mutate(nd_ter_each_freq = nd %>%
           assign_tertiles("nd", freq_ter, FALSE)) %>%
  ungroup

mot_random$typ %<>%
  lapply(function(df) {
    df %<>%
      mutate(nd_ter = nd %>%
               assign_tertiles("nd", 1:100, TRUE),
             nd_ter_26 = nd %>%
               assign_tertiles("nd", 2:6, TRUE)) %>%
      group_by(phonemic_len) %>%
      mutate(nd_ter_each_len = nd %>%
               assign_tertiles("nd", phonemic_len, TRUE)) %>%
      ungroup %>%
      group_by(freq_ter) %>%
      mutate(nd_ter_each_freq = nd %>%
               assign_tertiles("nd", freq_ter, FALSE)) %>%
      ungroup
  })

#### Study 2 ####
  #### noun sets ####
# extract maternal nouns of 1, 2, 3 syllables
mot_nouns <- manch_typ %>%
  filter(set == "MOT",
         len %in% 1:3) %>%
  select(id:len) %>%
  left_join(., manch_typ_gram, by = "phon") %>%
  filter(str_detect(gram, "^N[|]{1}[A-Z+]+$|N:PROP[|]{1}[A-Z+]+$"))
  
mot_nouns %<>%
  (function(df) {
    chi_nouns <- manch_typ %>%
      filter(set == "CHI",
             len %in% 1:3) %>%
      select(id:len) %>%
      left_join(., manch_typ_gram, by = "phon") %>%
      filter(str_detect(gram, "^N[:|]{1}")) 
    
    # check if A mother's noun was learned by A child
    for (i in seq_along(df$id)) {
      df$learned[i] <- chi_nouns %>%
        filter(id == df$id[i],
               phon == df$phon[i]) %>%
        nrow
    }
    
    # assign stage at which a mot noun was learned by a child
    df_learned <- df %>%
      filter(learned == TRUE) %>%
      (function(x) {
        for (i in seq_along(x$id)) {
          x$stage[i] <- chi_nouns %>%
            filter(id == x$id[i],
                   phon == x$phon[i]) %>%
                   {.$stage}
          
        }
        
        x
      })
    
    # assign new stage to not learned nouns
    # average stage at which a noun was learned (for each id and length)
    df_not_learned <- df %>%
      filter(learned == FALSE) %>%
      select(-stage) %>%
      left_join(., df_learned %>%
                  group_by(id, len) %>%
                  summarise(stage = mean(stage) %>% round(0)) %>%
                  ungroup, by = c("id", "len")) %>%
      select(id, stage, everything())
    
    rbind(
      df_learned,
      df_not_learned
    )
  })

# match learned/not learned sets on maternal frequency up to word stage (for each id and len)
# assign maternal frequency
for (i in seq_along(mot_nouns$id)) {
  mot_nouns$mot_freq[i] <- manch_tok %>%
    filter(set == "MOT",
           id == mot_nouns$id[i],
           stage %in% 0:(mot_nouns$stage[i]),
           phon == mot_nouns$phon[i]) %>%
    nrow
  
} ; rm(i)

# match frequency
set.seed(10011)
mot_nouns %<>%
  group_by(id, len) %>%
  group_split %>% # at least dplyr 0.8.0 needs to be installed to use group_split
  lapply(function(df) {
    df_learned <- df %>%
      ungroup %>%
      filter(learned == 1) %>%
      arrange(desc(mot_freq))
    
    df_not_learned <- df %>%
      ungroup %>%
      filter(learned == 0) %>%
      arrange(desc(mot_freq))
    
    mean_to_match <- df_not_learned$mot_freq %>% mean
    df_matched <- df_learned
    row_deletion <- 1
    range_difference <- sample(seq(-0.05, 0.15, 0.05), 1) # small flexible random range given
    
    while((mean(df_matched$mot_freq) - mean_to_match) > range_difference) {
      df_matched <- df_learned %>%
      {.[-c(1:row_deletion), ]}
      
      row_deletion %<>% {. + 1}
    } 
    
    n_row_df_learned <- nrow(df_learned)
    
    df_learned %<>%
      mutate(freq_match_include = rep(FALSE, row_deletion) %>%
               c(rep(TRUE, n_row_df_learned - row_deletion))) 
    
    df_not_learned %<>%
      mutate(freq_match_include = TRUE)
    
    df_learned %>%
      rbind(df_not_learned)
  }) %>%
  rbindlist

  #### reference vocabularies ####
# create children and maternal reference vocabularies (tokens, types)
# and generate phonological sequences for each word
vocabularies <- list(chi_tok = manch_tok %>% 
                  filter(set == "CHI") %>%
                  select(id, stage, phon) %>%
                  mutate(phonemic_len = phon %>%
                           str_split("_") %>%
                           sapply(length)) %>%
                  # exclude monophonemic words as
                  # they don't contain sequences from 2to5 phonemes
                  filter(phonemic_len != 1) %>% 
                  select(-phonemic_len),
                chi_typ = manch_tok %>% 
                  filter(set == "CHI") %>%
                  select(-gram, -word, -hour, -half_hour) %>%
                  group_by(set, id, stage) %>%
                  distinct(phon, .keep_all = TRUE) %>%
                  ungroup %>%
                  group_by(set, id, phon) %>%
                  filter(stage == min(stage)) %>%
                  ungroup %>%
                  arrange(set, id, stage, phon) %>%
                  select(id, stage, phon) %>%
                  mutate(phonemic_len = phon %>%
                           str_split("_") %>%
                           sapply(length)) %>%
                  filter(phonemic_len != 1) %>%
                  select(-phonemic_len),
                mot_tok = manch_tok %>% 
                  filter(set == "MOT") %>%
                  select(id, stage, phon) %>%
                  mutate(phonemic_len = phon %>%
                           str_split("_") %>%
                           sapply(length)) %>%
                  # exclude monophonemic words as
                  # they don't contain sequences from 2to5 phonemes
                  filter(phonemic_len != 1) %>% 
                  select(-phonemic_len),
                mot_typ = manch_tok %>% 
                  filter(set == "MOT") %>%
                  select(-gram, -word, -hour, -half_hour) %>%
                  group_by(set, id, stage) %>%
                  distinct(phon, .keep_all = TRUE) %>%
                  ungroup %>%
                  group_by(set, id, phon) %>%
                  filter(stage == min(stage)) %>%
                  ungroup %>%
                  arrange(set, id, stage, phon) %>%
                  select(id, stage, phon) %>%
                  mutate(phonemic_len = phon %>%
                           str_split("_") %>%
                           sapply(length)) %>%
                  filter(phonemic_len != 1) %>%
                  select(-phonemic_len))

vocabularies$chi_typ %<>%
  add_sequences

vocabularies$mot_typ %<>%
  add_sequences

vocabularies$chi_tok %<>%
  left_join(., vocabularies$chi_typ %>%
              select(-c(id, stage)) %>%
              distinct(phon, .keep_all = TRUE), by = "phon")

vocabularies$mot_tok %<>%
  left_join(., vocabularies$mot_typ %>%
              select(-c(id, stage)) %>%
              distinct(phon, .keep_all = TRUE), by = "phon")
  
# generate phoneme sequences for mot_nouns as well
mot_nouns %<>%
  add_sequences

  #### proportion of types ####
# compute proportion of types sharing phoneme sequences with each word
mot_nouns %<>%
  # in children's vocabulary (types)
  count_seqs_types(vocabularies$chi_typ, "prop_chi_typ_seq2", "seq2") %>%
  count_seqs_types(vocabularies$chi_typ, "prop_chi_typ_seq3", "seq3") %>%
  count_seqs_types(vocabularies$chi_typ, "prop_chi_typ_seq4", "seq4") %>%
  count_seqs_types(vocabularies$chi_typ, "prop_chi_typ_seq5", "seq5") %>%
  # in mothers' vocabulary (types)
  count_seqs_types(vocabularies$mot_typ, "prop_mot_typ_seq2", "seq2") %>%
  count_seqs_types(vocabularies$mot_typ, "prop_mot_typ_seq3", "seq3") %>%
  count_seqs_types(vocabularies$mot_typ, "prop_mot_typ_seq4", "seq4") %>%
  count_seqs_types(vocabularies$mot_typ, "prop_mot_typ_seq5", "seq5")

  #### proportion of tokens ####
# compute proportion of tokens sharing phoneme sequences with each word
mot_nouns %<>%
  # in children's vocabulary (tokens)
  count_seqs_tokens(vocabularies$chi_tok, "prop_chi_tok_seq2", "seq2") %>%
  count_seqs_tokens(vocabularies$chi_tok, "prop_chi_tok_seq3", "seq3") %>%
  count_seqs_tokens(vocabularies$chi_tok, "prop_chi_tok_seq4", "seq4") %>%
  count_seqs_tokens(vocabularies$chi_tok, "prop_chi_tok_seq5", "seq5") %>%
  # in mothers' vocabulary (tokens)
  count_seqs_tokens(vocabularies$mot_tok, "prop_mot_tok_seq2", "seq2") %>%
  count_seqs_tokens(vocabularies$mot_tok, "prop_mot_tok_seq3", "seq3") %>%
  count_seqs_tokens(vocabularies$mot_tok, "prop_mot_tok_seq4", "seq4") %>%
  count_seqs_tokens(vocabularies$mot_tok, "prop_mot_tok_seq5", "seq5")
