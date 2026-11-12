#R Module 1
#Assignment (& Exercise) Code


#Create all needed directories & download file into specified folder
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create(here::here("results"), showWarnings = FALSE)
dir.create(here::here("src"), showWarnings = FALSE)

 # you wouldnt download a car....
download.file("https://raw.githubusercontent.com/anoobvinu07/Genomic_assisted_selection/master/data/FitnessTraits_GeneticParameters_RedSpruce.txt",
              destfile = "data/raw/red_spruce_fitness_traits.txt")


#defining traits as the data file, so that the shorthand (traits) refers to the data within the file
traits <- read.delim(here::here("data/raw/red_spruce_fitness_traits.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)

#counting dimensions & checking variables to make sure data is being read properly. 
dim(traits)
names(traits)
head(traits, 5)

#Getting values of data, height specifically

summary(traits$Height)

#Exercises

#E2
str(traits)

#E4
table(traits$Population)
table(traits$Region)


#(Assignment:Computing summary stats for 4 traits broken down by population and region...)
#(E3)
summary((tapply(traits$Survival, traits$Population)), na.rm = TRUE)
summary((tapply(traits$Survival, traits$Region)), na.rm = TRUE)

summary((tapply(traits$Genetic_Diversity, traits$Population)), na.rm = TRUE)
summary((tapply(traits$Genetic_Diversity, traits$Region)), na.rm = TRUE)

summary((tapply(traits$Germination, traits$Population)), na.rm = TRUE)
summary((tapply(traits$Germination, traits$Region)), na.rm = TRUE)

summary((tapply(traits$Fitness, traits$Population)), na.rm = TRUE)
summary((tapply(traits$Fitness, traits$Region)), na.rm = TRUE)

#Producing figures

#Histogram
hist(traits$Survival, main = "Survival Distribution", xlab = "Rate of Survival", breaks = 30)

#box-plot (per chosen variable, out of curiosity)
boxplot(Fitness ~ Region, data = traits, main = "Fitness by Region", ylab = "Fitness", pch = 19)
boxplot(Survival ~ Region, data = traits, main = "Survival by Region", ylab = "Survival", pch = 19)
boxplot(Genetic_Diversity ~ Region, data = traits, main = "Genetic Diversity by Region", ylab = "GenDiverse", pch = 19)
boxplot(Germination ~ Region, data = traits, main = "Germination by Region", ylab = "Germination :)", pch = 19)

#^funky

#Scatter-plot w/trend lines
plot(traits$Germination, traits$Survival, xlab = "Germination", ylab = "Survival", main = "Germination Rates vs Survival Rates", abline(lm(Germination ~ Survival, data = traits), col = "maroon", lwd = 2))
plot(traits$Genetic_Diversity, traits$Fitness, xlab = "Gen Diversity", ylab = "Fitness", main = "Genetic Diversity vs Fitness Rates", abline(lm(Genetic_Diversity ~ Fitness, data = traits), col = "maroon", lwd = 2))
#more funky graphs. Would these data need to be transformed? How do you pick which is best? Is this where skewness & kurtosis comes into play?

#Scatter-plots for the factors you asked about (Elevation, Fitness or Seed weight)
plot(traits$Elevation, traits$Fitness, xlab = "Elevation", ylab = "Fitness", main = "Elevation vs Fitness Rates", abline(lm(Elevation ~ Fitness, data = traits), col = "maroon", lwd = 2))
plot(traits$SeedWeight, traits$Fitness, xlab = "Seed Weight", ylab = "Fitness", main = "Seed Weight vs Fitness Rates", abline(lm(SeedWeight ~ Fitness, data = traits), col = "maroon", lwd = 2))
plot(traits$Elevation, traits$Fitness, xlab = "Elevation", ylab = "SeedWeight", main = "Elevation vs SeedWeight", abline(lm(Elevation ~ SeedWeight, data = traits), col = "maroon", lwd = 2))
#I'm not sure why the trend line isn't showing for the 1st & 3rd chart

#Trying out a fitted GGplot for untransformed data
#ggplot(data = traits, aes(x = Genetic_Diversity, y = Fitness)) + geom_point() + geom_smooth(method = "loess", se = TRUE, level = 0.95) + scale_color_brewer
#(palette = "Dark2") + theme_classic() + xlab("Genetic Diversity") + ylab("Fitness") + labs(title = "Fitted Fitness & Genetic Diversity")
#...ggplot not loading properly .-.

pop_summary <- aggregate(cbind(Germination, Genetic_Diversity, Fitness, Survival) ~ Population, data = traits, FUN = mean)
pop_summary

#sub-sets & arrangement BUT doesn't work due to missing default.
#sorted_fitness <-pop_summary[order(pop_summary$fitness, )]
#sorted_fitness


#Saving...
write.csv(pop_summary, here::here("results/R_module 1_population_trait_summary.csv"), row.names = FALSE)

# exporting plot as a .png
png(here::here("results/Genetic_Diversity_by_Region.png"), width = 800, height = 600)
boxplot(Genetic_Diversity ~ Region, data = traits)
dev.off()


################################################################################################

#Idea check responses
#What is the difference between vector and a data-frame in R?
  #Vectors are 1 dimensional columns where all elements have to be homogeneous 
  #or the same type of variable whether its numeric or a character. A data-frame is a way 
  #to present this information . Each column acts as a vector, but they have to be 
  #the same size. Data-frames can have heterogeneous data, so numerical and characters. 
  
  
#Why prefer a relative path over an absolute path?
  #This allows others to run your code successfully, 
  #without having to edit file paths or break file links. 
  #Absolute paths are user-centric whereas relative paths can be used by others
  #When they run your code locally on their device or on a shared space such as the cluster. 
  
  
#What does str(traits) tell you that head(traits) doesn't?
  #str() will show the name of the column(s) whereas head() will show what the actual file's rows
  #say. It is useful to check both in case R misreads the file and shows the name of the column incorrectly. 
  #This can help show if data needs to be rearranged within the file itself for R to read it properly.

  
#Why does treating Region as a factor matter for boxplot(Height ~ Region, data = traits)
  #Treating Region as a factor allows the function to work, since it requires numerical data. 
  #I think that keeping Region as a character wouldn't be compatible with this function.
  
#If mean(traits$height) returned NA, what would you check first & how would you fix it?
  #I could review the column to have sums(is.na(x)) count the amount of NAs in that part of the data. to prevent this from interfering with the function again, 
  #I can add rm.na = TRUE to tell R to ignore these values and only work with the ones present
  
#Assignment Answers:

#Which population has the highest and lowest mean fitness?
#Population G18 has the highest mean fitness, at 25.53 & population XWS has the lowest mean fitness value at 2.18!

#Is there a visible relationship between elevation and fitness or seed weight?
#...not sure. Having issue with scatter-plots trend line showing up and ggplot isn't loading correctly on posit cloud







