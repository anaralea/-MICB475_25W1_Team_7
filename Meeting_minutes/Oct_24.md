  1. Discuss new metadata
  - correct for sex of deer (compounding variable) -> sex was not significant according to the paper (but let's still do this)
  - put filtering for sex in aim 1 before the other aims
  - Sample SRR23559504 is missing the sex of that sample (filter that one out)
  - supplemental data from paper has geographical location that doesn't align
  -   different sample ids from ours
  -   Evelyn will see if it is possible to be linked (either nothing or she'll provide new data by Monday)
  -   might not have enough samples for things if we analyse by geographical location
  -   region of origin was only significant for two metrics but not for alpha diversity
  -   some regions are just controls, some are just chronic-wasting disease
    
  2. Discuss QIIME2 processing
  - script for QIIME2 looks good, if there's issues it could be the dataset
  - manifest has 2 columns (forward, absolute pathway and reverse, absolute pathway)
  -   delete reverse-absolute filepath column
  -   remove "forward-" from forward column
  -   import again as single-read
  -   denoise it as single-read
  -   can do this in excel, ask Evelyn for help if needed
  - quality looks really high for demux file - Avril investigated
  -   just trim to majority so 250
    
  3. Review draft timeline for Gantt chart
  - quite ambitious right now, can push things back and move things around
  - specifically manuscript writing can start later
    
  4. Plan actionables for next week/meeting
  - get processing done early next week
  - get started on creating phyloseq object in R
  -   will make processing faster

  5. Hypothesis
  - Write prediction and have literature support
     
  6. Deadlines for the proposal
  - Soft deadline: October 31st (next meeting)
  - Final deadline: November 2nd 11:59PM
