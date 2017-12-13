library(ANTsR)
library(knitr)
##C57 SR
# Load in Behavior and Imaging Data
behavior <- read.csv("/Users/omega/Documents/Natalie/youngmice/VBA_results/data/C57_data_trials.csv")
labled.set <-read.csv("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/legendsCHASS2symmetric.csv")
labeled.brain.img <-antsImageRead("/Users/omega/Documents/Natalie/SCCAN_tutorial/Eig_zscore/data/MDT_labels_chass_symmetric.nii.gz")
mask <-antsImageRead('/Users/omega/Documents/Natalie/youngmice/VBA_results/data/MDT_mask_e3.nii')
mang_files<-list.files(path = "/Users/omega/Documents/Natalie/youngmice/VBA_results/data", pattern = "T2_to_MDT",full.names = T,recursive = T)
mang_mat <-imagesToMatrix(mang_files,mask)
jac_files<-list.files(path = "/Users/omega/Documents/Natalie/youngmice/VBA_results/data", pattern = "jac_to_MDT",full.names = T,recursive = T)
jac_mat <-imagesToMatrix(jac_files,mask)
suscept_files<-list.files(path = "/Users/omega/Documents/Natalie/youngmice/VBA_results/data", pattern = "X_to_MDT",full.names = T,recursive = T)
suscept_mat <-imagesToMatrix(suscept_files,mask)
la <-mang_mat[,90:dim(mang_mat)[2]]

# Split training and testing


regress.prob1 <- sparseRegression(inmatrix=jac_mat, outcome='d4', 
                                  demog=behavior, mask=mask, sparseness=0.1, 
                                  nvecs=4, its=3, cthresh=250)
regress.prob2 <- sparseRegression(inmatrix=suscept_mat, outcome='d4', 
                                  demog=behavior, mask=mask, sparseness=0.1, 
                                  nvecs=4, its=3, cthresh=350)
regress.prob3 <- sparseRegression(inmatrix=la, outcome='d4', 
                                  demog=behavior, mask=mask, sparseness=0.1, 
                                  nvecs=4, its=3, cthresh=250)

jac.eigmat <- imageListToMatrix(regress.prob1$eigenanatomyimages,mask) #convert sparse eigen regions into a matrix
suscept.eigmat <- imageListToMatrix(regress.prob2$eigenanatomyimages,mask) #convert sparse eigen regions into a matrix
mang.eigmat <- imageListToMatrix(regress.prob3$eigenanatomyimages,mask) #convert sparse eigen regions into a matrix

myprojs_jac <- jac_mat %*% t(jac.eigmat) # project those eigen images onto our original matrix of 6 images
myprojs_suscept <- suscept_mat %*% t(suscept.eigmat) # project those eigen images onto our original matrix of 6 images
myprojs_mang <- mang_mat %*% t(mang.eigmat) # project those eigen images onto our original matrix of 6 images

w <- data.frame(cbind(behavior[,'d4'],myprojs_suscept)) # column combind the behavior wth the projections
names(w) <- c('Distance_4', names(w)[-1]) # insert column names
mylm <- lm('Distance_4 ~ .', data=w) # behavior correlation with the number of projections
distpred <- predict(mylm, newdata=w) # based on the linear model predict the distances for the same day
distlm <-lm(distpred~behavior[,'d4']) # correlate real distances with predicted distances
modsum <-summary(distlm) # produce summary of the linear model
r2 <- modsum$adj.r.squared
my.p <- modsum$coefficients[2,4]
plot(behavior[,'d4'], distpred, xlab = 'Real Total Distance on Day 4', ylab = 'Predicted Total Distance on day 4', 
     main='Predicted vs. Real Total Distance on Day 4') # generate plot
abline(lm(distpred~behavior[,'d4']))
rp = vector('expression',2)
rp[1] = substitute(expression(italic(R)^2 == MYVALUE), 
                   list(MYVALUE = format(r2, digits = 3)))[2]
rp[2] = substitute(expression(italic(p) == MYOTHERVALUE), 
                   list(MYOTHERVALUE = format(my.p, digits = 2)))[2]
legend('topleft', legend = rp, bty = 'n')

smooth.eig.one <-smoothImage(regress.prob1$eigenanatomyimages[[1]], 0.1)
antsImageWrite(smooth.eig.one, "/Users/omega/Documents/Natalie/Thesis/C57/New_Method/suscept_c57_0p1_350_eig1_nm.nii.gz")

smooth.eig.two <-smoothImage(regress.prob1$eigenanatomyimages[[2]], 0.1)
antsImageWrite(smooth.eig.two, "/Users/omega/Documents/Natalie/Thesis/C57/New_Method/suscept_c57_0p1_350_eig2_nm.nii.gz")

smooth.eig.three <-smoothImage(regress.prob1$eigenanatomyimages[[3]], 0.1)
antsImageWrite(smooth.eig.three, "/Users/omega/Documents/Natalie/Thesis/C57/New_Method/suscept_c57_0p1_350_eig3_nm.nii.gz")

smooth.eig.four <-smoothImage(regress.prob1$eigenanatomyimages[[4]], 0.1)
antsImageWrite(smooth.eig.four, "/Users/omega/Documents/Natalie/Thesis/C57/New_Method/suscept_c57_0p1_350_eig4_nm.nii.gz")
