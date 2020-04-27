# library imports
library(tidyverse)
library(scales)
library(limma)
library(edgeR)
library(psych)

# ============== CV function ===================================================
CV <- function(df) {
    # Computes CVs of data frame rows
        # df - data frame, 
        # returns vector of CVs (%)
    ave <- rowMeans(df)    # compute averages
    sd <- apply(df, 1, sd) # compute standard deviations
    cv <- 100 * sd / ave   # compute CVs in percent (last thing gets returned)
}

# =========== Boxplot with median label ========================================
labeled_boxplot <- function(df, ylim, title) {
    # Makes a box plot with the median value labeled
        # df - data frame with data to compute CVs of
        # ylim - upper limit for y-axis
        # title - plot title
    cv = CV(df)
    boxplot(cv, ylim = c(0, ylim), notch = TRUE, main = title)
    text(x = 0.65, y = boxplot.stats(cv)$stats[3], 
         labels = round(boxplot.stats(cv)$stats[3], 1))
}

# ================== TMM normalization from DGEList object =====================
apply_tmm_factors <- function(y, color = NULL, plot = TRUE) {
    # computes the tmm normalized data from the DGEList object
        # y - DGEList object
        # returns a dataframe with normalized intensities
    
    # compute and print "Sample loading" normalization factors
    lib_facs <- mean(y$samples$lib.size) / y$samples$lib.size
    cat("\nLibrary size factors:\n", 
        sprintf("%-5s -> %f\n", colnames(y$counts), lib_facs))
    
    # compute and print TMM normalization factors
    tmm_facs <- 1/y$samples$norm.factors
    cat("\nTrimmed mean of M-values (TMM) factors:\n", 
        sprintf("%-5s -> %f\n", colnames(y$counts), tmm_facs))
    
    # compute and print the final correction factors
    norm_facs <- lib_facs * tmm_facs
    cat("\nCombined (lib size and TMM) normalization factors:\n", 
        sprintf("%-5s -> %f\n", colnames(y$counts), norm_facs))

    # compute the normalized data as a new data frame
    tmt_tmm <- as.data.frame(sweep(y$counts, 2, norm_facs, FUN = "*"))
    colnames(tmt_tmm) <- str_c(colnames(y$counts), "_tmm")
    
    # visualize results and return data frame
    if(plot == TRUE) {
        boxplot(log10(tmt_tmm), col = color, notch = TRUE, main = "TMM Normalized data")
    }
    tmt_tmm
}

# ================= reformat edgeR test results ================================
collect_results <- function(df, tt, x, xlab, y, ylab) {
    # Computes new columns and extracts some columns to make results frame
        # df - data in data.frame
        # tt - top tags table from edgeR test
        # x - columns for first condition
        # xlab - label for x
        # y - columns for second condition
        # ylab - label for y
        # returns a new dataframe
    
    # condition average vectors
    ave_x <- rowMeans(df[x])
    ave_y <- rowMeans(df[y])
    
    # FC, direction, candidates
    fc <- ifelse(ave_y > ave_x, (ave_y / ave_x), (-1 * ave_x / ave_y))
    direction <- ifelse(ave_y > ave_x, "up", "down")
    candidate <- cut(tt$FDR, breaks = c(-Inf, 0.01, 0.05, 0.10, 1.0), 
                     labels = c("high", "med", "low", "no"))
    
    # make data frame
    temp <- cbind(df[c(x, y)], data.frame(logFC = tt$logFC, FC = fc, 
                                          PValue = tt$PValue, FDR = tt$FDR, 
                                          ave_x = ave_x, ave_y = ave_y, 
                                          direction = direction, candidate = candidate, 
                                          Acc = tt$genes)) 
    
    # fix column headers for averages
    names(temp)[names(temp) %in% c("ave_x", "ave_y")]  <- str_c("ave_", c(xlab, ylab))    
    
    temp # return the data frame
}

# =============== p-value plots ================================================
pvalue_plots <- function(results, ylim, title) {
    # Makes p-value distribution plots
        # results - results data frame
        # ylim - ymax for expanded view
        # title - plot title
    p_plot <- ggplot(results, aes(PValue)) + 
        geom_histogram(bins = 100, fill = "white", color = "black") +
        geom_hline(yintercept = mean(hist(results$PValue, breaks = 100, 
                                     plot = FALSE)$counts[26:100]))

    # we will need an expanded plot
    p1 <- p_plot + ggtitle(str_c(title, " p-value distribution"))
    p2 <- p_plot + coord_cartesian(xlim = c(0, 1.0), ylim = c(0, ylim)) + ggtitle("p-values expanded")
    grid.arrange(p1, p2, nrow = 2) # from gridExtra package
}

# ============= log2 fold-change distributions =================================
log2FC_plots <- function(results, range, title) {
    # Makes faceted log2FC plots by candidate
        # results - results data frame
        # range - plus/minus log2 x-axis limits
        # title - plot title
    ggplot(results, aes(x = logFC, fill = candidate)) +
        geom_histogram(binwidth=0.1, color = "black") +
        facet_wrap(~candidate) +
        ggtitle(title) + 
        coord_cartesian(xlim = c(-range, range))
}

# ========== Setup for MA and volcano plots ====================================
transform <- function(results, x, y) {
    # Make data frame with some transformed columns
        # results - results data frame
        # x - columns for x condition
        # y - columns for y condition
        # return new data frame
    df <- data.frame(log10((results[x] + results[y])/2), 
                     log2(results[y] / results[x]), 
                     results$candidate,
                     -log10(results$FDR))
    colnames(df) <- c("A", "M", "candidate", "P")
    
    df # return the data frame
}

# ========== MA plots using ggplot =============================================
MA_plots <- function(results, x, y, title) {
    # makes MA-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots 
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # 2-fold change lines
    ma_lines <- list(geom_hline(yintercept = 0.0, color = "black"),
                     geom_hline(yintercept = 1.0, color = "black", linetype = "dotted"),
                     geom_hline(yintercept = -1.0, color = "black", linetype = "dotted"))

    # make main MA plot
    ma <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("logFC (", y, "/", x, ")")) +
        scale_x_continuous("Ave_intensity") +
        ggtitle(title) + 
        ma_lines
    
    # make separate MA plots
    ma_facet <- ggplot(temp, aes(x = A, y = M)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scale_y_continuous(paste0("log2 FC (", y, "/", x, ")")) +
        scale_x_continuous("log10 Ave_intensity") +
        ma_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)"))

    # make the plots visible
    print(ma)
    print(ma_facet)
}    

# ========== Scatter plots using ggplot ========================================
scatter_plots <- function(results, x, y, title) {
    # makes scatter-plot DE candidate ggplots
        # results - data frame with edgeR results and some condition average columns
        # x - string for x-axis column
        # y - string for y-axis column
        # title - title string to use in plots
        # returns a list of plots
    
    # 2-fold change lines
    scatter_lines <- list(geom_abline(intercept = 0.0, slope = 1.0, color = "black"),
                          geom_abline(intercept = 0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          geom_abline(intercept = -0.301, slope = 1.0, color = "black", linetype = "dotted"),
                          scale_y_log10(),
                          scale_x_log10())

    # make main scatter plot
    scatter <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        ggtitle(title) + 
        scatter_lines

    # make separate scatter plots
    scatter_facet <- ggplot(results, aes_string(x, y)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        scatter_lines +
        facet_wrap(~ candidate) +
        ggtitle(str_c(title, " (separated)")) 

    # make the plots visible
    print(scatter)
    print(scatter_facet)
}

# ========== Volcano plots using ggplot ========================================
volcano_plot <- function(results, x, y, title) {
    # makes a volcano plot
        # results - a data frame with edgeR results
        # x - string for the x-axis column
        # y - string for y-axis column
        # title - plot title string
    
    # uses transformed data
    temp <- transform(results, x, y)
    
    # build the plot
    ggplot(temp, aes(x = M, y = P)) +
        geom_point(aes(color = candidate, shape = candidate)) +
        xlab("log2 FC") +
        ylab("-log10 FDR") +
        ggtitle(str_c(title, " Volcano Plot"))
}

# ============== individual protein expression plots ===========================
# function to extract the identifier part of the accesssion
get_identifier <- function(accession) {
    identifier <- str_split(accession, "\\|", simplify = TRUE)
    identifier[,3]
}

set_plot_dimensions <- function(width_choice, height_choice) {
    options(repr.plot.width=width_choice, repr.plot.height=height_choice)
}

plot_top_tags <- function(results, nleft, nright, top_tags) {
    # results should have data first, then test results (two condition summary table)
    # nleft, nright are number of data points in each condition
    # top_tags is number of up and number of down top DE candidates to plot
    # get top ipregulated
    up <- results %>% 
        filter(logFC >= 0) %>%
        arrange(FDR)
    up <- up[1:top_tags, ]
    
    # get top down regulated
    down <- results %>% 
        filter(logFC < 0) %>%
        arrange(FDR)
    down <- down[1:top_tags, ]
    
    # pack them
    proteins <- rbind(up, down)
        
    color = c(rep("red", nleft), rep("blue", nright))
    for (row_num in 1:nrow(proteins)) {
        row <- proteins[row_num, ]
        vec <- as.vector(unlist(row[1:(nleft + nright)]))
        names(vec) <- colnames(row[1:(nleft + nright)])
        title <- str_c(get_identifier(row$Acc), ", int: ", scientific(mean(vec), 2), 
                       ", FDR: ", scientific(row$FDR, digits = 3), 
                       ", FC: ", round(row$FC, digits = 1),
                       ", ", row$candidate)
        barplot(vec, col = color, main = title,
                cex.main = 1.0, cex.names = 0.7, cex.lab = 0.7)
    }    
}

# load the IRS-normalized data and check the table
data_import <- read_tsv("labeled_grouped_protein_summary_TMT_9_IRS_normalized.txt", guess_max = 5250)

# "Filter" column flags contams and decoys
# "Missing" column flags proteins without reporter ion intensities (full sets missing)
# the table from pandas is sorted so the rows we want come first
data_all <- filter(data_import, is.na(Filter), is.na(Missing))
data_sl <- data_all %>% select(., contains("SLNorm_")) %>% 
  select(., -contains("_Unused")) %>% 
  select(., -contains("_Pool"))
data_irs <- data_all %>% select(., contains("IRSNorm_")) %>% 
  select(., -contains("_Unused")) %>% 
  select(., -contains("_Pool"))
data_pool <- data_all %>% select(., contains("_Pool_"))

# save a few columns for the results table
all_results <- data_all %>% select(., ProtGroup, Counter, Accession, Description, starts_with("PSMs_Used"))

# save gene names for edgeR so we can double check that results line up
accessions <- data_all$Accession

# see how many rows of data we have
length(accessions)

# all categories of metaplastic breast cancer tissue
mbc_sl <- select(data_sl, contains("_C"), contains("_SP"), contains("_SQ"))
mbc_irs <- select(data_irs, contains("_C"), contains("_SP"), contains("_SQ"))

# triple negative breast cancer tissue
tn_sl <- select(data_sl, contains("_TN"))
tn_irs <- select(data_irs, contains("_TN"))

# Normal tissue
n_sl <- select(data_sl, contains("_N"))
n_irs <- select(data_irs, contains("_N"))

# collect the biological replicate channels
bio_sl <- cbind(n_sl, tn_sl, mbc_sl)
bio_irs <- cbind(n_irs, tn_irs, mbc_irs)

# set a color vector for plots
color <- c(rep("red", 7), rep("blue", 6), rep("green", 14))

# boxplots
boxplot(log10(bio_sl), col = color, notch = TRUE, main = "SL Normalized data")
boxplot(log10(bio_irs), col = color, notch = TRUE, main = "IRS Normalized data")

# clustering of the replicates before IRS
plotMDS(log2(bio_sl), col = color, main = "Replicates before IRS")

# clustering after IRS
plotMDS(log2(bio_irs), col = color, main = "Replicates after IRS")

# the normal tissue samples (includes the one mis-classified sample)
pairs.panels(log10(n_sl), lm = TRUE, main = "Normal, SLNorm")
pairs.panels(log10(n_irs), lm = TRUE, main = "Normal, IRSNorm")

# the triple negative samples
pairs.panels(log10(tn_sl), lm = TRUE, main = "Triple Negative, SLNorm")
pairs.panels(log10(tn_irs), lm = TRUE, main = "Triple Negative, IRSNorm")

# the MBC C samples
pairs.panels(log10(mbc_sl[1:4]), lm = TRUE, main = "Metaplastic C, SLNorm")
pairs.panels(log10(mbc_irs[1:4]), lm = TRUE, main = "Metaplastic C, IRSNorm")

# the MBC Sp samples
pairs.panels(log10(mbc_sl[5:10]), lm = TRUE, main = "Metaplastic Sp, SLNorm")
pairs.panels(log10(mbc_irs[5:10]), lm = TRUE, main = "Metaplastic Sp, IRSNorm")

# the MBC Sq samples
pairs.panels(log10(mbc_sl[11:14]), lm = TRUE, main = "Metaplastic Sq, SLNorm")
pairs.panels(log10(mbc_irs[11:14]), lm = TRUE, main = "Metaplastic Sq, IRSNorm")

# check the pooled standards
pairs.panels(log10(data_pool[1:3]), lm = TRUE, main = "Pools, SLNorm")
pairs.panels(log10(data_pool[4:6]), lm = TRUE, main = "Pools, IRSNorm")

# put groups together into a single data frame
tmt_sl <- bio_sl
tmt_irs <- bio_irs

# define the positions of the groups
N <- 1:7
TN <- 8:13
MBC <- 14:27

# set some colors by condition
group <- c(rep("N", 7), rep("TN", 6), rep("MBC", 14))

# get the biological sample data into a DGEList object
y <- DGEList(counts = tmt_irs, group = group, genes = accessions)

# run TMM normalization (also includes a library size factor)
y <- calcNormFactors(y)

tmt_tmm <- apply_tmm_factors(y, color = color)

# check the clustering
plotMDS(y, col = color, main = "all samples after TMM")

# put CVs in data frames to simplify plots and summaries
cv_sl <- data.frame(N = CV(bio_sl[N]), TN = CV(bio_sl[TN]), MBC = CV(bio_sl[MBC])) 
cv_irs <- data.frame(N = CV(bio_irs[N]), TN = CV(bio_irs[TN]), MBC = CV(bio_irs[MBC]))
cv_tmm <- data.frame(N = CV(tmt_tmm[N]), TN = CV(tmt_tmm[TN]), MBC = CV(tmt_tmm[MBC]))

# see what the median CV values are
medians <- apply(cv_sl, 2, FUN = median)
print("SLNorm median CVs by condition (%)")
round(medians, 2)

medians <- apply(cv_irs, 2, FUN = median)
print("IRSNorm median CVs by condition (%)")
round(medians, 2)

medians <- apply(cv_tmm, 2, FUN = median)
print("Final median CVs by condition (%)")
round(medians, 2)

# see what the CV distibutions look like
# need long form for ggplot
long_cv_sl <- gather(cv_sl, key = "group", value = "cv")

# traditional boxplots
ggplot(long_cv_sl, aes(x = group, y = cv, fill = group)) +
  geom_boxplot(notch = TRUE) +
  ggtitle("SL CV distributions")

# density plots
ggplot(long_cv_sl, aes(x = cv, color = group)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 200)) +
  ggtitle("SL CV distributions")

# need long form for ggplot
long_cv_tmm <- gather(cv_tmm, key = "group", value = "cv") 

# traditional boxplots
ggplot(long_cv_tmm, aes(x = group, y = cv, fill = group)) +
  geom_boxplot(notch = TRUE) +
  ggtitle("Final CV distributions")

# density plots
ggplot(long_cv_tmm, aes(x = cv, color = group)) +
  geom_density() +
  coord_cartesian(xlim = c(0, 200)) +
  ggtitle("Final CV distributions")

# compute dispersions and plot BCV
y <- estimateDisp(y)
plotBCV(y, main = "BCV plot of IRS/TMM normalized data")

# the exact test object has columns like fold-change, CPM, and p-values
et <- exactTest(y, pair = c("N", "TN"))

# this counts up, down, and unchanged genes (proteins) at 10% FDR
summary(decideTestsDGE(et, p.value = 0.10))

# the topTags function adds the BH FDR values to an exactTest data frame 
# make sure we do not change the row order (the sort.by parameter)!
topTags(et)$table
tt <- topTags(et, n = Inf, sort.by = "none")

# make an MD plot (like MA plot)
plotMD(et, p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# check the p-value distribution
ggplot(tt$table, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(et$table$PValue, breaks = 100, 
                                    plot = FALSE)$counts[26:100])) +
  ggtitle("N versus TN p-value distribution")

# get the results summary
results <- collect_results(tmt_tmm, tt$table, N, "N", TN, "TN")

# make column names unique by adding comparison (for the accumulated frame)
results_temp  <- results
colnames(results_temp) <- str_c(colnames(results), "_N_TN")

# accumulate the testing results
all_results <- cbind(all_results, results_temp)

# see how many candidates by category
results %>% count(candidate)

# plot log2 fold-changes by category
ggplot(results, aes(x = logFC, fill = candidate)) +
  geom_histogram(binwidth=0.1, color = "black") +
  facet_wrap(~candidate) +
  coord_cartesian(xlim = c(-4, 4)) +
  ggtitle("N vs TN logFC distributions by candidate")

# make MA plots
MA_plots(results, "ave_N", "ave_TN", "N vs TN")

# make scatter plots
scatter_plots(results, "ave_N", "ave_TN", "N vs TN")

# make a volcano plot
volcano_plot(results, "ave_N", "ave_TN", "N vs TN")

# look at the top 20 candidates in each direction (up in TN, then down in TN)
set_plot_dimensions(6, 3.5)
plot_top_tags(results, 7, 6, 20)
set_plot_dimensions(7, 7)

# the exact test object has columns like fold-change, CPM, and p-values
et <- exactTest(y, pair = c("N", "MBC"))

# this counts up, down, and unchanged genes (proteins) at 10% FDR
summary(decideTestsDGE(et, p.value = 0.10))

# the topTags function adds the BH FDR values to an exactTest data frame 
# make sure we do not change the row order (the sort.by parameter)!
topTags(et)$table
tt <- topTags(et, n = Inf, sort.by = "none")

# make an MD plot (like MA plot)
plotMD(et, p.value = 0.10)
abline(h = c(-1, 1), col = "black") # 2-fold change lines

# check the p-value distribution
ggplot(tt$table, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(et$table$PValue, breaks = 100, 
                                    plot = FALSE)$counts[26:100])) +
  ggtitle("N vs MBC p-value distribution")

# get the results summary
results <- collect_results(tmt_tmm, tt$table, N, "N", MBC, "MBC")

# make column names unique by adding comparison
results_temp  <- results
colnames(results_temp) <- str_c(colnames(results), "_N_MBC")

# accumulate the testing results
all_results <- cbind(all_results, results_temp)

# see how many candidates by category
results %>% count(candidate)

# plot log2 fold-changes by category
ggplot(results, aes(x = logFC, fill = candidate)) +
  geom_histogram(binwidth=0.1, color = "black") +
  facet_wrap(~candidate) +
  coord_cartesian(xlim = c(-4, 4)) +
  ggtitle("N vs MBC logFC distributions by candidate")

# make MA plots
MA_plots(results, "ave_N", "ave_MBC", "N vs MBC")

# make scatter plots
scatter_plots(results,  "ave_N", "ave_MBC", "N vs MBC")

# make a volcano plot
volcano_plot(results,  "ave_N", "ave_MBC", "N vs MBC")

# look at the top 20 candidates (up in MBC, then down in MBC)
set_plot_dimensions(6, 3.5)
plot_top_tags(results, 7, 14, 20)
set_plot_dimensions(7, 7)

# the exact test object has columns like fold-change, CPM, and p-values
et <- exactTest(y, pair = c("TN", "MBC"))

# this counts up, down, and unchanged genes (proteins) at 10% FDR
summary(decideTestsDGE(et, p.value = 0.10))

# the topTags function adds the BH FDR values to an exactTest data frame 
# make sure we do not change the row order (the sort.by parameter)!
topTags(et)$table
tt <- topTags(et, n = Inf, sort.by = "none")

# make an MD plot (like MA plot)
plotMD(et, p.value = 0.10)
abline(h = c(-1, 1), col = "black")

# check the p-value distribution
ggplot(tt$table, aes(PValue)) + 
  geom_histogram(bins = 100, fill = "white", color = "black") + 
  geom_hline(yintercept = mean(hist(et$table$PValue, breaks = 100, 
                                    plot = FALSE)$counts[26:100])) +
  ggtitle("TN vs MBC p-value distribution")

# get the results summary
results <- collect_results(tmt_tmm, tt$table, TN, "TN", MBC, "MBC")

# make column names unique by adding comparison
results_temp  <- results
colnames(results_temp) <- str_c(colnames(results), "_TN_MBC")

# accumulate the testing results
all_results <- cbind(all_results, results_temp)

# see how many candidates by category
results %>% count(candidate)

# plot log2 fold-changes by category
ggplot(results, aes(x = logFC, fill = candidate)) +
  geom_histogram(binwidth=0.1, color = "black") +
  facet_wrap(~candidate) +
  coord_cartesian(xlim = c(-4, 4)) +
  ggtitle("TN vs MBC logFC distributions by candidate")

# make MA plots
MA_plots(results, "ave_TN", "ave_MBC", "TN vs MBC")

# make scatter plots
scatter_plots(results, "ave_TN", "ave_MBC", "TN vs MBC")

# make a volcano plot
volcano_plot(results, "ave_TN", "ave_MBC", "TN vs MBC")

# look at the top 10 candidates (up in MBC, then down in MBC)
set_plot_dimensions(6, 3.5)
plot_top_tags(results, 6, 14, 10)
set_plot_dimensions(7, 7)

# write the results to disk
write.table(all_results, "three_tissues_results.txt", sep = "\t",
           row.names = FALSE, na =  " ")

# log the session details
sessionInfo()


