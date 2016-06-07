# work with markup comments on chars

file <- "/Users/chadeliason/Dropbox/phenome dataset/data/2015-09-02/modified/final_reorderedJAC+CME.txt"

txt <- readLines(file)

txt <- paste0(txt, collapse="\n")

# find actions, comments associated with characters

tmp <- str_match_all(txt, regex("\n([A-Z\\?]+)(\\d+)", multiline=TRUE, dotall=TRUE))

todo <- tmp[[1]][,2]
charnum <- as.numeric(tmp[[1]][,3])

todo

charnum








str_match_all(txt, regex("\n([A-Z\\?]+)(\\d+)\\..*?\\{(.*?)\\}", multiline=TRUE, dotall=TRUE))

# tmp <- str_match_all(txt, regex("\\{(.*?)\\}", multiline=TRUE, dotall=TRUE))
tmp <- tmp[[1]]
todo <- tmp[!is.na(tmp[,2]),2] # todo (e.g., DD - duplicated, XX - cut, KKM - keep and modify...)
charnum <- tmp[!is.na(tmp[,3]),3] # charnum
comments <- tmp[!is.na(tmp[,3]),5] # comments

todo

charnum

comments

# charnum
# todo
# comments

write.csv(data.frame(charnum, todo, comments), file = "/Users/chadeliason/Dropbox/Chad/phenome/data/2015-09-02/modified/regex_extracted.csv")



# length(charnum)

# tmp <- str_match_all(txt, regex("\\{(.*?)\\}", multiline=TRUE, dotall=TRUE))

# txt[str_detect(txt, "\\{.*\\}")]

