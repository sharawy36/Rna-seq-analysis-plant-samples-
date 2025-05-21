 count <- read.delim("C:/Users/pc/Downloads/count.txt", comment.char="#")
count$Chr = NULL 
count$Start = NULL
count$End = NULL
count$Strand =NULL
count$Length =NULL
colnames(count)
names(count)[names(count)=="X.home.m.sharawy.yarb.aligned.sample1.bam"] = "sample_1"
names(count)[names(count)=="X.home.m.sharawy.yarb.aligned.sample2.bam"] = "sample_2"
names(count)[names(count)=="X.home.m.sharawy.yarb.aligned.sample3.bam"] = "sample_3"
names(count)[names(count)=="X.home.m.sharawy.yarb.aligned.sample4.bam"] = "sample_4"
dim(count)
mode(count)
class(count)
df_cleaned <- count[apply(count, 1, function(row) all(row != 0)), ]
varRow = apply(count, 1, var, na.rm = T )
constRow <- (varRow == 0 | is.na(varRow))
df_clean =count[ !constRow, ]
dim(df_clean)
View(df_clean)
rownames(df_clean) = df_clean$Geneid
df_clean$Geneid = NULL
df = as.matrix(df_clean)
genes=row.names(df)
df=apply(df,2,as.integer)
row.names(df)=genes
View(df)
after_normalized = log2(df+1)
par(mfrow  = c(1, 2))
p1 = plot(density(df), main = "Before Transformation")
p2 = plot(density(after_normalized), main = "AFter Normalization")
##box plot before normalization for features 
options(repr.plot.width=10, repr.plot.height=8)
par(mar = c(8, 5, 2, 2))

boxplot(t(df), 
        main = "Boxplot Before Normalization", 
        horizontal = T, 
        names = rownames(df), 
        las = 2,
        col ="lightgreen"
)

##box plot after  normalization 
options(repr.plot.width=10, repr.plot.height=8)
par(mar = c(8, 5, 2, 2))

boxplot(t(after_normalized), 
        main = "Boxplot after Normalization", 
        horizontal = T, 
        names = rownames(after_normalized), 
        las = 2,
        col = "red"
)
###visualization for samples 
options(repr.plot.width=10, repr.plot.height=8)
par(mar = c(8, 5, 2, 2))

boxplot(df, 
        main = "Boxplot Before Normalization", 
        horizontal = T, 
        names = colnames(df), 
        las = 2,
        col ="lightgreen"
)

##box plot after  normalization 
options(repr.plot.width=10, repr.plot.height=8)
par(mar = c(8, 5, 2, 2))

boxplot(after_normalized, 
        main = "Boxplot after Normalization", 
        horizontal = T, 
        names = colnames(after_normalized), 
        las = 2,
        col = "red"
)
#######3scaling 
mat_scaled= scale(t(after_normalized), center = T, scale = T)
### Density plot for Raw, normalizd and scaled Data
options(repr.plot.width = 15, repr.plot.height = 6)
par(mar = c(8, 5, 2, 2), mfrow = c(1, 3), cex.axis = 1.5)



#plotting before Normalization
plot(
  density(apply(df, 1, mean, na.rm= T)), 
  main = "Before Log2", 
  xlab = "Mean Expression", 
  ylab = "Density", 
  col = "blue", 
  lwd = 2
)

# plotting after Normalization with log2
plot(
  density(apply(after_normalized, 1, mean, na.rm= T)), 
  main = "After Log2", 
  xlab = "Mean Expression", 
  ylab = "Density", 
  col = "red", 
  lwd = 2
)
# plotting After Z score scaling

plot(
  density(apply(mat_scaled, 2, mean, na.rm= T)), 
  main = "After Log2", 
  xlab = "Mean Expression", 
  ylab = "Density", 
  col = "lightgreen", 
  lwd = 2
)
##################### box plot for samples 
# Set plot options
options(repr.plot.width = 10, repr.plot.height = 8)

# Set layout and margins
par(mar = c(8, 5, 2, 2), mfrow = c(1, 3))

# Create boxplot for the  samples before log2 transformation
boxplot(df,
        main = "Before log2",
        horizontal = TRUE,
        names = colnames(df),
        las = 2,
        col = "red")

# Create boxplot for t samples after log2 transformation
boxplot(after_normalized,
        main = "After log2",
        horizontal = TRUE,
        names = colnames(after_normalized),
        las = 2,
        col = "lightgreen")

# Create boxplot for the  after log2 transformation and scaling
boxplot(t(mat_scaled),
        main = "After log2 + scaled",
        horizontal = TRUE,
        names = rownames(mat_scaled),
        las = 2,
        col = "blue")

