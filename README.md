# -MICB475_25W1_Team_7
## October 10, 2025
### Agenda:
  1. Decide on a dataset and project
  - Create a research question and decide on aims
  2. Discuss group roles
  - Individual strengths and weaknesses
  - Discuss who will be making meeting agendas and taking meeting notes
  - Discuss breakup of work in the main project
  3. Discuss layout of GitHub lab notebook

### Meeting Minutes
1. Dataset Comparisons
   - one dataset on deer and one on mouse
   - the mouse, metadata doesn't contain a lot of information, mostly timepoints (weeks post-infection, whether infected or not infected), just alpha and beta diversity that may have already been analyzed.
   - deer model has less metadata but is better, analyzed alpha and beta diversity among prion status (infection status), sex of the deers, and whether it came from farm or freerange (which wasn't annotated in the metadata). pathology, negative deer are negative, positive are either positive in the lymph nodes or the brain

2. Deer dataset
  - not sure what we can infer for the deer, and not sure when deer was infected dependent on where the deer are from, disease progression from lymph nodes -> brain
  - significant difference between farms and freerange (source of the sample), but no way of distinguishing that from metadata
  - deer didn't do functional analysis (optional module 19, similar to aligning ASVs to taxa, but instead aligning them to metabolic pathways)
  - R techniques we could apply: core microbiome, differential abundance (deer did a little bit), indicator taxa
  - if we separate based on source of the sample, would reduce sample size?

3. Project ideas
  - pathology comparison (Wisconsin has mix of farm 1 and 2, there are 132 samples, 100 are freerange (control) from Ohio)
  - farm 2 has 30 samples, farm 1 has 100 (look at how many samples per pathology group)
  - 30 samples start with P, 100 start with W, and then another 100 starting with C (control) -> P vs W pathology comparison
  - not many freerange, freerange has lowest diversity, could include if we have enough data from freerange
  - age ranges from 1-13 years

4. Methodologies
- filter out freerange
- reproduce figure 1 (unweighted unifrac) -> to confirm farm 1 vs farm 2
- separate farms
- diversity (alpha and beta) metrics analysis on each farm
- coremicrobiome
- indicator taxa -> maybe just farm 1 might not have enough samples for farm 2
- differential abundance
- functional analysis -> dependent on how much data we obtain from previous steps

5. Actionables
- read through rubric for proposal
- build draft proposal and figure out which methods fall within which aims for next week to discuss
- qiime processing has to be done before the proposal
- proposal is due in 2 weeks -> due on 26th
- ask for support earlier in week if needed
