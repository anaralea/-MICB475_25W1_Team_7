1. Functional Analysis
   - ana has meeting with avril today at 4pm, file needed is not in the server
   - there might not be anything functionally interesting, so will decide to leave in/out after meeting
     
2. Manuscript
   - figures are outlined in Nov_28.md
   - use title from presentation
   - don't focus conclusion on sex-stratifying because could be consequence of sample size
   - forgot the proposal ever existed and aims

   Introduction
   - here's a research gap, here's what we did
   - can have a short background on CWD and prion-disease
   - clearly establish why it's relevant to microbiome studies and why we're doing this
   - what we're trying to investigate and what we found (punchline right at the end)
   - paragraph or so on previous studies (dataset) and what they found (what stuff has been done in the field already)
     
   Methods
   - talk about the techniques used in the figures/tables (those are the sections for the methods)
   - we did this analysis with this package, we filtered this way, used these thresholds, etc. (assume person knows the technique but no need to go into detail)
   - one key thing: comprehensive of all the tools and how we reference them
     
   Discussion:
   - either literature, or focus on biological model
   - results become the outline
   - talking about no significance in alpha/beta (i.e. we didn't find, in the paper they did/didn't, is it consistent with literature)
   - we looked at two levels of severity which the paper didn't (ours is more nuanced)
   - variance in beta diversity is very low, doesn't capture/isn't fully representative
   - core microbiome same thing
   - deseq2 can compare specific bugs (these bugs are associated more with.., are the species more pathogenic in males/females, also say it isn't described in literature)
   - quoting literature (extrapolates) is in the discussion

   Results:
   - just be pure observations
   - just repeat what we found from the figures, talk about it more at the biological level (model)
   - i.e. if the diversity not showing significance, means the communities are not that fundamentally different from each other, etc.
   - i.e. males have more unique microbes than females, etc.
   - don't go over the functionality here, in methods
   - can be short
  
   Limitations:
   - uneven sample sizes (males have significantly less)
   - 12 different regions (geographic origin) that they didn't annotate, they found regional differences that could affect our result (p-value: 0.001)
   - might explain why we're not getting diversity differences, could be a confounding variable

   Future Directions:
    - shortterm and longterm plans
    - functional analysis could shortterm
    - longterm could be analyzing specific bugs
     
3. Github Cleanup / Team Meeting Grade
  - create more sub-directories in code, figures, etc.
  - only keep final version of the code
  - except for indicator taxa (keep the full one also)
  - for data, separate QIIME outputs and our outputs into folders
  - otherwise, that's okay
