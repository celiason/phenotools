# work with markup comments on chars

txt <- readLines("/Users/chadeliason/Dropbox/Chad/phenome/data/2015-09-02/modified/final_reorderedJAC+CME.txt")
txt <- paste0(txt, collapse="\n")

# find actions, comments associated with characters
tmp <- str_match_all(txt, regex("\n([A-Z\\?]+)(\\d+)\\.(.*?)\\{(.*?)\\}", multiline=TRUE, dotall=TRUE))
tmp <- tmp[[1]]
todo <- tmp[!is.na(tmp[,2]),2] # todo (e.g., DD - duplicated, XX - cut, KKM - keep and modify...)
charnum <- tmp[!is.na(tmp[,3]),3] # charnum
comments <- tmp[!is.na(tmp[,3]),5] # comments

# charnum
# todo
# comments

write.csv(data.frame(charnum, todo, comments), file = "/Users/chadeliason/Dropbox/Chad/phenome/data/2015-09-02/modified/regex_extracted.csv")

