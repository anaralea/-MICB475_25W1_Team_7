1. How to interpret results from analyses:
  - Alpha & beta diversity (aim 1) - Joshua
    - no significance
    - present the abundance-based alpha diversity ones
    - beta diversity we might just want bray curtis (matches well with the Shannon, other ones could be supplemental)
    - keep bray curtis in the manuscript
    
  - Core microbiome & indicator taxa (aim 2) - Clara
    - LN and Brain has huge differences
    - 0 for abundance (should increase because we are getting every rare taxa, change to 0.1% (0.001) or 0.01% (0.0001)) threshold
    - 0.2 prevalence thresholds
    - how many bugs are unique to each condition vs shared
    - males have very different profile
    - female microbiomes are more sensitive
    - keep the ones that are 0.6 or higher or unique to both disease conditions
    - male indicators are all significant so we can keep all of them, for females we should filter
    - core microbiome has a lot of significance but indicator taxa not many
    - keep the venndiagrams for presentation also
      
  - Differential abundance analysis (aim 3) - Theresa
    - keep the female and males separated
    - increase the bar plot axis and remove the ASV names
    - look up the bugs phylum in the literature and see what exactly they do
    - make a progression to see if they're the same things to see if they're further depleted
    - line the graphs up to the progression
    - male has same pattern as females but still a lot
    - keep volcano only for presentation, don't include the bar graphs but have everything in the manuscript
   
    - convert everything to relative abundance and plot it out
    - take the phyloseq object and take the transform sample counts command (makes it relative abundance)
    - use the melt command (psmelt)
    - use that data frame and input that into ggplot
    - make a boxplot of the 3 different indicator groups and see if it matches with DESeq2
    - use the highest DESeq2 ones that are -25 or +25 etc.
   
2. Next steps
  - Do we have enough to talk about from the figures?
  - Functional analysis (aim 4), does this need to be done?
    - could add on to the differential abundance analysis
    - PICRus don't do it on the server because our servers are not powerful enough
    - instead send Avril the location of the files (absolute file paths) we want to run and the code we want to run
    - do all 3 pairwise comparisons like we did for DESeq2, follow instructions in the module (download the file, there's one function broken)
    - don't worry if we don't have it done for the presentation

3. Presentation
  - We are paired with team 8 (machine learning), contact them on the Monday because their stuff is a little complicated
  - for figures don't include everything because other team is presenting our slides
  - just include the main points we want to include in our story
  - could pick and choose some specific genus/genera for the presentation, if the patterns are obvious and easy
  - otherwise keep the general findings
  - Soft deadline: November 28th 11:59PM (next meeting)
  - Due date: November 30th 11:59PM

4. Manuscript
  - Have draft structure and outline for the 5th
  - Draft due date: December 14th 11:59PM
  - Final due date: December 21st 11:59PM

5. Plan actionables for next week/meeting
  - Discuss schedules/availabilities for next week(s)
  - Delegate tasks
  - Discuss goals for next meeting
