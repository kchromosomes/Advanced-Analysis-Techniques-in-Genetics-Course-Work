 set.seed(123)
> n_ind = 30
> n_loci = 10
> geno = matrix(sample(0:2, n_ind * n_loci, replace = true, prob = c(0.4, 0.4,0.2)), nrow + n_ind, ncol = n_loci)
Error: object 'true' not found

> geno = matrix(sample(0:2, n_ind * n_loci, replace = TRUE, prob = c(0.4, 0.4,0.2)), nrow + n_ind, ncol = n_loci)
Error in nrow + n_ind : non-numeric argument to binary operator

> geno = matrix(sample(0:2, n_ind * n_loci, replace = TRUE, prob = c(0.4, 0.4,0.2)), nrow = n_ind, ncol = n_loci)
> rownames(geno) = paste0("ind", 1:n_ind)
> colnames(geno) = paste0("locus", 1:n_loci)
> geno[1:5, 1:5]
     locus1 locus2 locus3 locus4 locus5
ind1      0      0      1      2      0
ind2      1      1      1      1      0
ind3      0      0      2      0      0
ind4      0      2      1      0      1
ind5      0      1      1      1      2
> alt_freq = colMeans(geno) / 2
> ref_freq = 1 - alt_freq
> data.frame(locus = names(alt_freq), alt_freq = round(alt_freq, 3), ref_freq = round(ref_freq, 3))
          locus alt_freq ref_freq
locus1   locus1    0.433    0.567
locus2   locus2    0.367    0.633
locus3   locus3    0.400    0.600
locus4   locus4    0.467    0.533
locus5   locus5    0.383    0.617
locus6   locus6    0.383    0.617
locus7   locus7    0.450    0.550
locus8   locus8    0.450    0.550
locus9   locus9    0.467    0.533
locus10 locus10    0.483    0.517
> #part 2: H-W equilibrium
> set.seed(42)
> p_1 = 0.8
> q_1 = 1-p_1
> pop_1 = data.frame(population = "Pop_1")
> genotype = samples (c("AA", "Aa", "aa") size = 20, replace = TRUE, prob = c(p_1^2), 2 * p_1 * q_1, q_1^2)) head(pop_1)
Error: unexpected symbol in "genotype = samples (c("AA", "Aa", "aa") size"

> set.seed(42) p_1 = 0.8 q_1 = 1 - p_1
Error: unexpected symbol in "set.seed(42) p_1"

> set.seed(42)
> p_1 = 0.80
> q_1 = 1 - p_1
> pop_1 = data.frame(population = "Pop_1"), genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_1^2, 2 * p_1 * q_1, q_1^2))
Error: unexpected ',' in "pop_1 = data.frame(population = "Pop_1"),"

> pop_1 = data.frame(population = "Pop_1"), genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_1^2, 2 * p_1 * q_1, q_1^2)))
Error: unexpected ',' in "pop_1 = data.frame(population = "Pop_1"),"

> pop_1 = data.frame(population = "Pop_1") genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_1^2, 2 * p_1 * q_1, q_1^2)))
Error: unexpected symbol in "pop_1 = data.frame(population = "Pop_1") genotype"

> # yay, syntax
> set.seed(42)
> 
> # Population 1: 20 individuals; allele A frequency p = 0.80
> p_1 <- 0.80  
> q_1 <- 1 - p_1
> 
> pop_1 <- data.frame(population = "Pop_1",
+                     # sample function for choosing genotypes
+                     genotype = sample(c("AA", "Aa", "aa"), # possible genotype values
+                                       size = 20,                           # no. of individuals to simulate
+                                       replace = TRUE,                      # genotypes can  occur more than once
+                                       # p2 + 2pq + q2
+                                       prob = c(p_1^2, 2 * p_1 * q_1, q_1^2)) # probability of each genotype
+ )
> 
> head(pop_1) # peek into the pop_1 dataframe
  population genotype
1      Pop_1       Aa
2      Pop_1       Aa
3      Pop_1       AA
4      Pop_1       Aa
5      Pop_1       Aa
6      Pop_1       AA
> #population 2..
> p_2 = 0.30
> q_2 = 1 - p_2
> pop_2 <- data.frame(population = Pop_2, genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_2^2, 2 * p_2 * q_2, q_2^2)))
Error: object 'Pop_2' not found

> pop_2 <- data.frame(population = pop_2, genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_2^2, 2 * p_2 * q_2, q_2^2)))
Error: object 'pop_2' not found

> pop_2 <- data.frame(population = "Pop_2", genotype = sample(c("AA", "Aa", "aa"), size = 20, replace = TRUE, prob = c(p_2^2, 2 * p_2 * q_2, q_2^2)))
> head(pop_2)
  population genotype
1      Pop_2       Aa
2      Pop_2       aa
3      Pop_2       AA
4      Pop_2       AA
5      Pop_2       aa
6      Pop_2       Aa
> #calculating genotype freq...
> geno_count_pop1 = table(pop_1$genotype)
> geno_count_pop2 = table(pop_2$genotype)
> #conv. genotype counts to frequencies
> geno_prop_pop1 = prop.table(table(pop_1$genotype))
> geno_prop_pop2 = prop.table(table(pop_2$genotype))
> #combining indvs from both populations
> pooled_population = rbind(pop_1, pop_2)
> geno_count_pooled = table(pooled_population$genotype)
> geno_prop_pooled = prop.table(table(pooled_population$genotype))
> #building comparision table...
> geno_summary = data.frame(genotype = c("AA", "Aa", "aa"), count_pop1 = as.integer(geno_count_pop1[c("AA", "Aa", "aa")]), freq_pop1 = as.numeric(geno_prop_pop1[c("AA", "Aa", "aa")]), count_pop2 = as.integer(geno_prop_pop2[c("AA", "Aa", "aa")]), freq_pop2 = as.numeric(geno_prop_pop2[c("AA", "Aa", "aa")]), count_pooled = as.integer(geno_count_pooled[c("AA", "Aa", "aa")]), freq_pooled = as.numeric(geno_prop_pooled[c("AA", "Aa", "aa")])
+ geno_summary
Error: unexpected symbol in:
"= as.integer(geno_prop_pop2[c("AA", "Aa", "aa")]), freq_pop2 = as.numeric(geno_prop_pop2[c("AA", "Aa", "aa")]), count_pooled = as.integer(geno_count_pooled[c("AA", "Aa", "aa")]), freq_pooled =
geno_summary"

> 
> 
> geno_summary <- data.frame(
+     genotype = c("AA", "Aa", "aa"),
+     count_pop1 = as.integer(geno_count_pop1[c("AA", "Aa", "aa")]),
+     freq_pop1  = as.numeric(geno_prop_pop1[c("AA", "Aa", "aa")]),
+     count_pop2 = as.integer(geno_count_pop2[c("AA", "Aa", "aa")]),
+     freq_pop2  = as.numeric(geno_prop_pop2[c("AA", "Aa", "aa")]),
+     count_pooled = as.integer(geno_count_pooled[c("AA", "Aa", "aa")]),
+     freq_pooled  = as.numeric(geno_prop_pooled[c("AA", "Aa", "aa")])
+ )
> 
> geno_summary
  genotype count_pop1 freq_pop1 count_pop2 freq_pop2 count_pooled freq_pooled
1       AA          9      0.45          2       0.1           11       0.275
2       Aa         10      0.50         10       0.5           20       0.500
3       aa          1      0.05          8       0.4            9       0.225
> sum(geno_prop_pop1)
[1] 1
> sum(geno_prop_pop2)
[1] 1
> sum(geno_prop_pooled)
[1] 1
> 
> 
> 
> #calc allele frequencies
> 
> counts_1 = geno_count_pop1
> n_1 = NROW(pop_1)
> n_A_1 = 2 * counts_1["AA"] + counts_1["Aa"]
> n_a_1 = 2 * counts_1["aa"] + counts_1["Aa"]
> 
> #allele counts to frequencies
> p_A_1 = n_A_1 / (2 * N_1)
Error: object 'N_1' not found

> p_A_1 = n_A_1 / ( 2 * n_1)
> q_A_1 = n_a_1 / ( 2* n_1)
> p_A_1
 AA 
0.7 
> q_A_1
 aa 
0.3 
> p_A_1 + q_A_1
AA 
 1 
> 
> 
> #Exercise 1
> #Estimate allele freq. for Population 2
> 
> #(count genotype for pop2)
> counts_2 = geno_count_pop2
> 
> #(amnt of indvs)
> N_2 = nrow(pop_2)
> 
> #counting a" & a' alleles
> n_A_2 = 2 * counts_2["AA"] + counts_2["Aa"]
> n_a_2 = 2 * counts_2["aa"] + counts_2["Aa"]
> 
> #convert allele counts to allele freqs
> p_A_2 = n_A_2 / ( 2 * n_2)
Error: object 'n_2' not found

> p_A_2 - n_A_2 / ( 2 * N_2)
Error: object 'p_A_2' not found

> p_A_2 <- n_A_2 / ( 2 * N_2)
> #guess this command actually needs "<-" instead of "="
> q_A_2 <- n_a_2 / (2 * N_2)
> p_A_2
  AA 
0.35 
> q_A_1
 aa 
0.3 
> oops
Error: object 'oops' not found

> q_A_2
  aa 
0.65 
> p_A_2 + q_A_2
AA 
 1 
> #...misnamed q as A instead of a
> 
> q_A_2 = "q_a_2"
> p_A_2 + q_a_2
Error: object 'q_a_2' not found

Error: unexpected ':' in ":"

> 
> 
> View(pop_2)
> View(pooled_population)
> ....
Error: object '....' not found

> 
> 
> allele_count_pop1 = table(pop_1$allele)
> allele_count_pop2 = table(pop_2$allele)
> 
> allele_prop_pop1 = prop.table(table(pop_1$allele))
> head(pop_2$allele)
NULL
> agh
Error: object 'agh' not found

> #revisit binding populations
> 
> 
> 
> 
> #PART III
> #Obsv. & Expected Heterozygosity
> 
> Ho = colMeans( geno == 1)
> He = 2 * ref_freq * alt_freq
> 
> het_table <- data.frame (locus = names (Ho), Ho = round(Ho, 3), He = round(He, 3), diff = round (Ho - He, 3))
> het_table
          locus    Ho    He   diff
locus1   locus1 0.467 0.491 -0.024
locus2   locus2 0.333 0.464 -0.131
locus3   locus3 0.467 0.480 -0.013
locus4   locus4 0.533 0.498  0.036
locus5   locus5 0.433 0.473 -0.039
locus6   locus6 0.433 0.473 -0.039
locus7   locus7 0.500 0.495  0.005
locus8   locus8 0.433 0.495 -0.062
locus9   locus9 0.333 0.498 -0.164
locus10 locus10 0.300 0.499 -0.199
> 
> 
> 
> #PCA for population structure
> 
> pca_sim <- prcomp(geno, scale. = TRUE)
> plot(pca_sim$x[, 1], pca_sim$x[, 2], xlab = "PC1", ylab = "PC2", main = "Simulated genotype PCA")
Error in plot.new() : figure margins too large

> #(paste)
> 
> pca_sim <- prcomp(geno, scale. = TRUE)
> plot(pca_sim$x[, 1], pca_sim$x[, 2],
+      xlab = "PC1", ylab = "PC2", main = "Simulated genotype PCA")
Error in plot.new() : figure margins too large

> ??plot.new
> 
> #unsure of how to adjust margins
> #could be lit. pane size..
> pca_sim <- prcomp(geno, scale. = TRUE)
> > plot(pca_sim$x[, 1], pca_sim$x[, 2],
Error: unexpected '>' in ">"

> pca_sim <- prcomp(geno, scale. = TRUE) plot(pca_sim$x[, 1], pca_sim$x[, 2],
Error: unexpected symbol in "pca_sim <- prcomp(geno, scale. = TRUE) plot"

> pca_sim <- prcomp(geno, scale = TRUE) 
> plot(pca_sim$x[, 1], pcs_sim$x[, 2], xlab = "PC1", ylab = "PC2", main = "Simulated genotype PCA")
Error: object 'pcs_sim' not found

> plot(pca_sim$x[, 1], pca_sim$x[, 2], xlab = "PC1", ylab = "PC2", main = "Simulated genotype PCA")
Error in plot.new() : figure margins too large

> ??pch
> ?pch
> pch
Error: object 'pch' not found

> plot(pca_sim$x[, 1], pca_sim$x[, 2], pch = 19, xlab = "PC1", ylab = "PC2", main = "Simulated genotype PCA")
> #
> #Needed to add pch (plotting character value) to determine point size/type? Solved margin issue. 
> #For future reference: value changes size or sahpe; 8 can result in a star while "" can let you set the symbol for the plot points
> 
> 
> #Working with Red-Spruce PC1/PC2 data
> traits <- read.table("data/raw/red_spruce_fitness_traits.txt", header = TRUE, sep = "\t")
> 
> 
> traits$Region <- factor(traits$Region)
> plot(trats$PC1, traits$PC2), col = as.integer(traits$Region), pch = 19, xlab = "PC1", ylab = "PC2", main = "Genetic PCA Colored by Region"
Error: unexpected ',' in "plot(trats$PC1, traits$PC2),"

>  plot(traits$PC1, traits$PC2, col = as.integer(traits$Region), pch = 19, xlab = "PC1", ylab = "PC2", main = "Genetic PCA Colored by Region")
> legend("topright", legend = levels(traits$Region), col = 1:nlevels(traits$Region), pch = 19)
> 
> 
> #Notes:"In population genetics, PC1/PC2 values from a genetic marker PCA usually captures the strongest axes of population structure - often correlated with geography (isolation by distance) or major environmental gradients"


> 
> View(pooled_population)
> View(geno)
> View(geno)
> View(pop_1)
> View(pop_1)
> View(traits)
> View(pop_2)
> View(geno_summary)
> View(pooled_population)
> View(pca_sim)
> View(het_table)
> View(geno_summary)
> View(pca_sim)
> pca_sim[["x"]]
              PC1         PC2         PC3          PC4          PC5         PC6         PC7        PC8
ind1  -0.54152354 -1.94120659 -2.13694385  1.035913421 -0.792827620  0.14780062 -0.13319993 -0.7536204
ind2   1.22532949  0.33995242  0.15019549  1.165531683 -0.467449626  0.26624691  0.61754151  0.1463899
ind3  -1.44329264 -1.69530672  0.88219697  0.062731060 -0.001022344  2.20835119  1.21746187 -0.5322827
ind4  -1.03161694 -0.03278653  1.26367351 -0.521862619 -0.675947277  1.54747539 -0.90976128  0.5453483
ind5   1.76671762  1.91042163 -1.17044089  0.170150126  0.750369059  1.25186692 -0.09231409  0.2312667
ind6  -2.21141244 -0.20838866 -0.69791074  1.579247139  0.204551410  0.78473978  0.10494498 -0.5012024
ind7  -0.70123342 -0.20483856 -0.18502790  1.378113340  1.836808487  0.09879450  0.06047906 -0.1865367
ind8  -0.38307771 -1.53800752 -0.91025362 -0.544025476  0.537970812 -0.42667667 -1.39868714  0.1906098
ind9  -2.45587600 -0.39911927 -1.05249417 -0.316384302  0.154800968  0.35032338 -0.71448759 -0.2250345
ind10 -1.30757737  2.71509014 -0.77907230  0.406008005 -1.704780040  0.19406814  1.19890760 -0.1294024
ind11 -0.19125049  2.48760789 -1.69930662  0.888178100  1.379363519 -0.74095489  0.91478826  0.5140267
ind12  1.00850931 -0.54474259 -0.68402770 -0.353026329 -1.341344363  0.88935739 -0.05761962  0.2764375
ind13  1.22144786 -1.23813752 -1.81478976 -1.852666359  0.134704771  0.27511481  1.17217414  1.2395170
ind14 -1.68542510  0.95851316  0.61309220 -0.978697507  0.586893938  0.07994359 -1.83941786  0.5642012
ind15  1.51426444  0.47543099 -1.23936278 -0.770876573  0.583759350 -0.78075479 -0.95830693 -0.7934078
ind16  2.36345031 -0.59971441  1.97158503  1.403437778  0.941044800  0.23794501  0.26332125 -0.7120314
ind17  0.09385855  1.85551870  2.53586686 -0.647540370 -0.529311178 -0.84396906 -0.03320234 -0.2079016
ind18  0.62116007 -0.81892998  0.55249103 -0.137989295 -0.972590307 -0.31789016  0.55553833 -1.1352935
ind19  0.29712734  0.44960452  0.58614308 -1.630142256  0.018728292 -0.28735292  0.84598281  1.4201804
ind20 -1.00931717 -1.24550855  2.10135248 -0.006614578  2.468268860 -0.04471651  0.45090443  0.5007426
ind21 -1.34140302  0.93796072  1.22918166 -0.007395484 -0.305858880 -1.54839259  0.58572036 -0.6170983
ind22 -0.62968584  0.45216514 -0.46648418  2.776729652 -0.931343324 -0.82350334 -0.84883068  1.4198388
ind23  2.13267184 -0.97667388 -0.76307034 -0.256751521  0.338953052 -0.32570520 -0.26318521 -0.5168112
ind24 -0.66987426  0.95908778  0.31402726 -1.479100737 -1.215576868  1.10425671 -0.36425418 -0.9600651
ind25  0.22649233  1.23503831 -0.85493871 -1.606925427  0.623749540 -0.79905753 -0.12692150 -1.1574312
ind26  2.89480280 -0.90772080  0.54140906  0.185710555 -0.995887178  0.26013817 -0.24789136  0.9138486
ind27 -0.47778136 -2.12815011 -0.21948001 -0.267613242 -0.248354484 -2.29660374  0.61068291 -0.3434258
ind28 -1.75034070 -0.67835371  0.04966899 -0.928127543  0.409074521 -0.21502067  0.53086570  0.8185436
ind29  2.16970530  1.14513269  0.59809582  0.534310977  0.920962411  0.87681113 -0.55007837 -0.3880198
ind30  0.29515072 -0.76393867  1.28462415  0.719677784 -1.707710299 -1.12263557 -0.59115512  0.3786136
               PC9        PC10
ind1  -0.470651264 -0.42242045
ind2   1.273285997 -0.15745030
ind3   0.739393928 -0.53147515
ind4  -1.365158205  0.50762761
ind5  -0.743773205 -0.46005547
ind6   0.693478152  1.09509174
ind7  -0.857095750  0.26220066
ind8  -0.003051882 -0.87860929
ind9  -0.070612857  0.45545296
ind10 -1.307020719 -0.66865189
ind11  0.741576868 -0.40612974
ind12  0.862482124  0.17402276
ind13 -0.531482585  0.14381323
ind14  0.195122124 -0.54026149
ind15  0.406241043  0.49379463
ind16 -0.169794801 -0.31029414

