# -MICB475_25W1_Team_7
## November 28, 2025
### Agenda:
1. Discuss Functional Analysis
  - Will focus on doing metacyc analysis only
  - Where to find ggpicrust tutorial R file?
  - Which differential abundance analysis to run?
  - Should we generate a heatmap and PCA plot?
  - Should we do an enrichment analysis?
2. Discuss Presentation powerpoint
  - Finalise title
  - Get feedback on background slides
  - Discuss which figures to include
  - Discuss what main points we want talked about for each figure
  - Should we add a slide discussing where we got the dataset from?
  - How pretty do we want the presentation to be? Ana is not inclined towards powerpoint design :(
  - https://docs.google.com/presentation/d/1ZBg2K2fEGqfh-9jn72PzBNw0xwfwf0P_NF1B1iTidC0/edit?usp=sharing
3. Discuss Group 8's presentation
  - When/who to reach out to to get primed on the machine learning aspect
  - Should touch base with each other once we receive the presentation on who will be presenting what part
## November 21, 2025
### Agenda:
1. How to interpret results from analyses (https://docs.google.com/document/d/1hUvy8htjWEY0QiQ4kn2A1pbDCHXjczSSJCuUWtARLBE/edit?usp=sharing):
  - Alpha & beta diversity (aim 1) - Joshua
  - Core microbiome & indicator taxa (aim 2) - Clara
      - Venn diagrams (x2) for female & male
          - While creating this, I set:
            choice of detection (abundance) = 0
            choice of detection threshold (prevalance) = 0.2. 0.3, 0.2
female_neg_ASVs <- core_members(female_neg, detection = 0, prevalence = 0.2)
female_ln_pos_ASVs <- core_members(female_ln_pos, detection = 0, prevalence = 0.3)
female_br_ln_pos_ASVs <- core_members(female_br_ln_pos, detection = 0, prevalence = 0.2)
          - Is this okay?
      - ISA tables in .tsv (x2) for female & male
  - Differential abundance analysis (aim 3) - Theresa
    - Volcano plots (x6) for female & male
    - Bar plots (x6) for female & male
2. Next steps
  - Do we have enough to talk about from the figures?
  - Functional analysis (aim 4), does this need to be done?
3. Presentation
  - We are paired with team 8
  - Soft deadline: November 28th 11:59PM (next meeting)
  - Due date: November 30th 11:59PM
4. Manuscript
  - Have draft structure and outline ready for next meeting
  - Draft due date: December 14th 11:59PM
  - Final due date: December 21st 11:59PM
5. Plan actionables for next week/meeting
  - Discuss schedules/availabilities for next week(s)
  - Delegate tasks
  - Discuss goals for next meeting

## November 14, 2025
### Agenda:
1. Phyloseq Code
- Phyloseq objects are made
- Updated code and saved phyloseq objects can be found in phyloseq folder on Github
- Does everything look good?
2. Next steps / Divide up tasks
- What processing and analysis are we doing next? And who is doing what?
- Alpha & beta diversity (aim 1)
- Core microbiome & indicator taxa (aim 2)
- Differential abundance analysis (aim 3)
- Functional analysis (aim 4)
3. Project Proposal Resubmission
- Go over grading rubric
4. Plan actionables for next week/meeting
- Discuss schedules/availabilities for next week
- Delegate tasks
- Discuss goals for next meeting

## November 7, 2025
### Agenda:
_Note: Please edit/update the Nov_7.md file in the Meeting_minutes folder with today's meeting notes!_
1. Discuss Qiime2 processing
- Metadata was modified based on Keemei feedback
  - sample_name column changed to samplename
- Which classifier to use for taxonomy?
- What kind of filtering needs to be done in QIIME2 (if any)?
  - Mitochondria and Chloroplast, frequency based, metadata based
- Are we doing Alpha rarefaction and diversity metrics in QIIME2?
2. Phyloseq Code
- Go over code that was uploaded
- Does it look good?
- How are we separating male/female?
3. Next steps
- What processing and analysis are we doing next?
- Alpha & beta diversity (aim 1)
- Core microbiome & indicator taxa (aim 2)
- Differential abundance analysis (aim 3)
- Functional analysis (aim 4)
4. Plan actionables for next week/meeting
  - Discuss schedules/availabilities for next week
  - Delegate tasks
  - Discuss goals for next meeting
    
## October 31, 2025
### Agenda:
  [Link to Proposal Here](https://docs.google.com/document/d/1P7zjuJCClrZfTSDmK4jr_p3a6KggbvnTL8H-YwDolqA/edit?usp=sharing)
  1. Discuss QIIME2 processing
  - Issues with artificial quality scores
  - Updates from Avril
  - Next steps
  2. Gantt Chart
  - How to get rid of blue bars and only keep orange?
  - How to select for certain bars only to change colours?
  3. New Dataset Wrangling Section
  - Go over what is on proposal so far
  - What else should be included? What level of detail?
  4. Plan actionables for next week/meeting
  - Discuss schedules/availabilities for next week
  - Discuss goals for next meeting
  5. Deadlines for the proposal - do we need more time or no?
  - Soft deadline: TODAY!
  - Final deadline: November 2nd 11:59PM

## October 24, 2025
### Agenda:
  [Link to Proposal Here](https://docs.google.com/document/d/1P7zjuJCClrZfTSDmK4jr_p3a6KggbvnTL8H-YwDolqA/edit?usp=sharing)
  1. Discuss new metadata
  - Changes to the focus of the project (farm comparison, or sex comparison?)
  - Changes to the aims - to review
  - Sample SRR23559504 is missing the sex of that sample
  2. Discuss QIIME2 processing
  - Go over rough QIIME2 code
  - Troubleshooting difficulties
  3. Review draft timeline for Gantt chart
  - Check feasibility
  4. Plan actionables for next week/meeting
  - Discuss schedules/availabilities for next week
  - Discuss potential to take on more beyond just QIIME2 processing and proposal
  5. Deadlines for the proposal
  - Soft deadline: October 31st (next meeting)
  - Final deadline: November 2nd 11:59PM

## October 17, 2025
### Agenda:
  [Link to Proposal Here](https://docs.google.com/document/d/1P7zjuJCClrZfTSDmK4jr_p3a6KggbvnTL8H-YwDolqA/edit?usp=sharing)
  1. Discuss research question
  - Draft: How do microbial communities vary across prion pathology states in white-tailed deer, and are these differences consistent between Farm 1 and Farm 2?
  2. Discuss aims of the project
  - Aim 0: Reproduce unweighted UniFrac PcoA plots to confirm previously reported differences between Farm 1 and Farm 2
  - Aim 1: Analyze alpha and beta diversity in deer gut microbiomes across pathology states in farm vs free-range environment
  - Aim 2: Identify core microbiome members and indicator taxa that are associated with pathology states
  - Aim 3: Perform differential abundance testing between pathology states within farms
  3. Discuss what approaches belong to which aims and what needs to be done for QIIME2 processing
  - Diversity (Alpha & Beta) analysis
  - Coremicrobiome Analysis
  - Indicator Taxa
  - Differential Abundance Analysis
  - Functional Analysis
  4. Divide the workload for the project proposal & discuss schedules/availability for next week
  5. Deadlines for proposal
  - Soft deadline: October 24th (next meeting)
  - Final deadline: October 26th 11:59PM

    
## October 10, 2025
### Agenda:
  1. Decide on a dataset and project
  - Create a research question and decide on aims
  2. Discuss group roles
  - Individual strengths and weaknesses
  - Discuss who will be making meeting agendas and taking meeting notes
  - Discuss breakup of work in the main project
  3. Discuss layout of GitHub lab notebook
