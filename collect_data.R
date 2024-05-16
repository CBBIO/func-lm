library(tidyverse)
library(furrr)
library(glue)
furrr_options(globals = TRUE,  seed = T)
plan(strategy = 'multisession', workers = 3)

prefix="HUMAN"

check_terms <- function(gold,preds){
  
  ontology <- case_match(preds$ONTOLOGY[[1]],
                         "Biological Process"~"BP",
                         "Molecular Function"~"MF",
                         "Cellular Component"~"CC")
  
  # print(gold)
  # print(preds$GOID.pred)
  # print(ontology)
  
  parents   <- unlist(mget(gold,get(paste0("GO",ontology,"PARENTS")),ifnotfound = NA))
  children  <- unlist(mget(gold,get(paste0("GO",ontology,"CHILDREN")),ifnotfound = NA))
  ancestors <- unlist(mget(gold,get(paste0("GO",ontology,"ANCESTOR")),ifnotfound = NA))
  offspring <- unlist(mget(gold,get(paste0("GO",ontology,"OFFSPRING")),ifnotfound = NA))
  return(map_chr(preds$GOID.pred, ~if_else(.x==gold,'HIT',
                                 if_else(.x %in% parents,"PARENT",
                                         if_else(.x %in% children,"CHILD",
                                                 if_else(.x %in% ancestors,"ANCESTOR",
                                                         if_else(.x %in% offspring,"OFFSPRING","NO RELATION")))))
  ))
}



library(GO.db)
select(GO.db,keys = keys(GO.db),columns = c("TERM","ONTOLOGY"))  %>%
  mutate(ONTOLOGY=case_when(ONTOLOGY == 'BP' ~ "Biological Process",
                            ONTOLOGY == 'MF' ~ "Molecular Function",
                            ONTOLOGY == 'CC' ~ "Cellular Component")) -> terms

gold_files<-dir(path='gold',pattern=prefix,full.names = T)
gold<-map_dfr(gold_files,read_delim,delim=';',col_names=c('UNIPROT_ID','GOID')) %>%
  separate_longer_delim(cols = GOID,delim = "|") %>% 
  left_join(terms) %>% filter(!is.na(TERM))

gold<-gold %>% filter(!TERM %in% c("biological_process","molecular_function","cellular_component"))

prediction_files <- dir(path='predictions',pattern=prefix,full.names = T)
names(prediction_files) <- str_match(prediction_files,paste0(prefix,"_(.+).csv"))[,2]

all_pred <- map_dfr(prediction_files,read_delim,delim=';',col_names= c("UNIPROT_ID","ONTOLOGY","GOID","TERM","SCORE"),.id = "method") %>%
  mutate(UNIPROT_ID=str_match(UNIPROT_ID,"\\|(.+)\\|")[,2])


all_pred<-all_pred %>% filter(!TERM %in% c("biological_process","molecular_function","cellular_component","all"))

#all_pred<-bind_rows(all_pred,gold %>% mutate(SCORE=Inf,method='gold') %>% dplyr::select(method,UNIPROT_ID,ONTOLOGY,GOID,TERM,SCORE) ) 
  
big_table <- left_join(gold,all_pred,
                           by = c('UNIPROT_ID','ONTOLOGY'),
                           multiple = "all", relationship = "many-to-many",suffix=c(".gold",".pred")) 


big_table %>%
    nest(data = c(ONTOLOGY,GOID.pred,TERM.pred,method,SCORE)) %>%
  #filter(row_number()<10) %>% 
    mutate(results = future_map2(GOID.gold, data, ~check_terms(.x, .y),
                                 .options = furrr_options(seed = T, packages = "GO.db"),
                                 .progress = TRUE)) %>% 
    unnest(cols = c(data, results)) -> super_big_table


save(super_big_table, file=paste0("super_big_table_rr_",prefix,".Rdata"))
