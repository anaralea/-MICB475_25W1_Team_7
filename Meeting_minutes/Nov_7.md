1. Discuss Qiime2 processing
- no need to run alpha rarefaction or any filtering in Qiime, will do it all in R
- Classifier that is listed in the Qiime code covers only V4 region, will need a different classifier
- Need to train our own classifier using the provided primers in the paper
- Added the training code to Qiime code
- Added export code for the phyloseq components

2. Phyloseq Code
- Removed load package code, not good practice to load them everytime
- Added comments to the code for clarity
- Removed shannon entropy is.na code
- Adding code for filtering male and female
    - subset samples, creating two new phyloseq objects that need to be saved to computer
- Removed all of the alpha diversity code, not relevant until we subset
  
3. Next steps
- Need to create rarefied versions of each phyloseq
- Get processing done ASAP
- Want results by the time we have Avril back in to meet the week after next
  
5. Plan actionables for next week/meeting
- A little behind on the processing, will meet next week in person
- Finish running qiime code and export all the files - Joshua
- Upload files to Github
- Will allocate first 3 aims next week
