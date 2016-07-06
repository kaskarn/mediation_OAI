#### SIMPLE DATA FIXES ON EXISTING VARIABLES GO HERE

##### Data fixes and reformats
# all discrete ints to factor
for(i in colnames(d2)) if(is.integer(d2[[i]])) if(length(unique(d2[[i]])) < 10) d2[[i]] <- factor(d2[[i]])

#smoking categories:
for(i in grep("smoker$", colnames(d2))) levels(d2[[i]]) <- c("Never", "Current", "Former", "Never") #Second "Never" is occasional, not regular ever smokers)

#infer missing age using time between T's:
ageloc <- match("v00age", colnames(d2))
for(i in 1:length(times[-1])) d2[is.na(d2[ageloc+i]),ageloc+i] <- d2[is.na(d2[ageloc+i]),ageloc+i-1] + tyears[i+1] - tyears[i]

#categorize education and change referent
d2$v00edcv %<>% factor(labels = c("<HS", "HS", "S.C.", "Coll", "S.G.", "Grad"))
d2$v00edcv %<>% relevel(ref = last(levels(.)))