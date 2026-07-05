# sceptre3 Slack — gRNA-assignment discussion (channel C05MA6XN5LZ)


**[2023-08-08 18:17] Timothy Barry**
```sum(response_matrix_highmoi_experimental["ENSG00000159023", sceptre_object@grna_assignments$grna_group_idxs$chr1.4081_top_two] &gt;= 1)```

**[2023-08-09 09:55] Eugene Katsevich**
I was envisioning the plots Tim points to as being relevant to <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?pane=issue&amp;itemId=35280140|Write a function to visualize gRNA count distributions across cells> rather than <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?pane=issue&amp;itemId=35280056|Write a function to visualize positive control results>.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?pane=issue&amp;itemId=35280140
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?pane=issue&amp;itemId=35280056

**[2023-08-09 18:42] Timothy Barry**
set the channel description: Here is a link to the project: <https://github.com/orgs/Katsevich-Lab/projects/6>
  - link: https://github.com/orgs/Katsevich-Lab/projects/6

**[2023-08-09 19:56] Eugene Katsevich** _(thread reply)_
Does this mean that if I reassign `b` in `o2`, then `o@b` will also change?

**[2023-08-10 12:26] Timothy Barry**
This might be somewhat lower priority, but I wanted to flag this before I forget. There might be situations in which we want to perform cell QC before gRNA assignment. For example, suppose the user selects a mixture model option to assign gRNAs to cells. If some of the cells are doublets (and thus have extremely high gRNA expression), this could cause the mixture model to fail. Thus, we might need to do some filtering of low-quality cells before assigning gRNAS. (The thresholding and max assignment strategies should be more robust to these issues.)

**[2023-08-10 12:28] Timothy Barry**
Of course, we need a cell QC step after gRNA assignment as well to remove cells with 0 or 2+ gRNAs in low MOI.

**[2023-08-10 12:29] Eugene Katsevich**
I think the gRNA assignment step should this kind of cell QC built into it.

**[2023-08-10 12:31] Eugene Katsevich**
A bigger concern I have is that gRNAs need to be assigned to cells before pairwise QC is done, but as you suggest it would be good to do cellwise QC prior to assigning gRNAs to cells. I wanted to have all of the QC steps in one function. Doing some QC, followed by gRNA assignment, followed by more QC, sounds onerous.

**[2023-08-10 12:33] Eugene Katsevich**
Unless we build in “some filtering of low-quality cells before assigning gRNAs” as part of the gRNA assignment step. Either way, users should be aware that some cell QC may happen at the time of gRNA assignment.

**[2023-08-10 13:36] Timothy Barry**
Got it. I think this may be something that we need to think a bit more carefully about. Presumably the gRNA assignment function would have an argument `grna_expression_cutoff_quantile` (or something to this effect), where cells with a gRNA library size above this threshold would be removed. I imagine that we also might want `grna_expression_cutoff_quantile` as an argument in `run_qc`. So we would need to check for consistency of this argument across `assign_grnas` and `run_qc`, for example.

**[2023-08-10 15:36] Timothy Barry**
I’ve been going back-and-forth with a few clients today (some old, some new). I have a few thoughts based on these conversations (feel free to ignore). Lots of people (in fact probably most people) seem to want to analyze singleton gRNAs instead of gRNA groups. We might consider splitting up specification of the gRNA groups into two pieces. First, within `create_sceptre_object`, the user passes a (required) vector indicating which of the gRNAs are NT gRNAs. Next, within `plan_analysis`, the user passes an optional data frame mapping each _targeting_ gRNA to its gRNA group. If this latter data frame is not passed, we carry out a singleton analysis.

**[2023-08-17 18:30] Timothy Barry**
Hey Joseph, yes, I’d recommend using the example high MOI data in the function documentation of `create_sceptre_object`. The data frame `pc_pairs_highmoi_experimental` contains positive control pairs. If you run the code in the example, you’ll get the negative and positive control results. I hope this helps.

**[2023-09-04 11:24] Timothy Barry**
Hey Joseph, thank you for your PR. I merged the PR into the dev branch and generated the PC plot using the high MOI example data (see below). I think this is a good starting point, but I have a few comments.

1. It looks like the text overlaps in several places (see the title and x-axis labels). It would be good if we could avoid this, e.g. by shortening the names of the titles and axis labels.
2. The idea to include both a violin plot and histogram is interesting. However, the two plots contain the same information, so I’m not sure there’s much of a value-add to having the histogram there. We might consider removing the histogram.
3. The violin plot is a bit hard to read: the PC p-values look like a vertical line and the NC p-values look like a horizontal line. I think we should try a jitter plot (via `geom_jitter()`) instead of a violin plot, which would enable us to see the individual p-values more clearly. You might need to play around with the parameters of the jitter plot to get this looking nice. Color might help as well.
4. I think that the text at the bottom might be somewhat misleading. In general the number of PC pairs does not match the number of NC pairs. Thus, it is not really an apples-to-apples comparison to apply BH to the two sets of p-values and then report the number of discoveries. I think it might be best to remove the text from this plot.
In summary I think a good starting point for revising this plot would be to create a single jitter plot containing the negative control p-values in one column and the positive control p-values in the other. Would you be able to give this a shot (if you think it is reasonable)?

Thanks again for the PR. I am happy that we are up-and-running with you making PRs to the codebase. Also, thanks for including such a detailed comment in your PR.

**[2023-09-04 13:00] Louis Deutsch** _(thread reply)_
from here? <https://gist.github.com/JoFrhwld/2266961>
  - link: https://gist.github.com/JoFrhwld/2266961

**[2023-09-04 15:32] Timothy Barry**
Thanks Joseph, it looks quite good. I merged the PR. I also made a couple small changes (see latest <https://github.com/Katsevich-Lab/sceptre/commit/a6f376d5f8eff5de6dd0abbe57d8321cec45d96f|commit>).
  - link: https://github.com/Katsevich-Lab/sceptre/commit/a6f376d5f8eff5de6dd0abbe57d8321cec45d96f

**[2023-09-05 20:49] Louis Deutsch**
In `plot_grna_count_distributions` if the user doesn't provide `grnas_to_plot`, so a random 4 are used, would it be worth printing a message notifying the user that these were random? So it's clear that there isn't any particular reason for those 4 having been picked? (in addition to the docs saying it too)

**[2023-09-05 21:28] Timothy Barry** _(thread reply)_
It’s an interesting idea, but I think I’d lean against it in the interest of not overwhelming the user with console messages. However, we should consider removing the `set.seed` so that the user sees n random gRNAs each time, thereby enabling the user to cycle through the gRNAs. We might even consider keeping track of which gRNAs have been shown so far so that the user can cycle through all of the gRNAs. We could do this via closure. I think Gene suggested something along these lines.

**[2023-09-05 21:32] Louis Deutsch**
I'm exploring a bar plot as an alternative to the histogram in `plot_grna_count_distributions`. Here's how it currently looks, and a bar plot instead. Thoughts?

**[2023-09-06 17:08] Timothy Barry** _(thread reply)_
If we were to do something along these lines, I think the bin width should increase exponentially, e.g. 1 -&gt; 2 -&gt; 4 -&gt; 8 etc. Furthermore, there would need to be enough bins to accommodate gRNAs with very high counts. Is the hybrid scale not an option?

**[2023-09-07 21:12] Louis Deutsch**
(new thread) Finished a draft of the bar plot approach to `plot_grna_count_distributions`. This is `plot_grna_count_distributions(sceptre_object, threshold=11)` with the low MOI data. If `threshold` is not provided, the single-width bins go up to 10, and then the exponentially growing bins start with 11-12. If `threshold` is provided, the single-width bins go to that value (or to 10, whichever is larger) and the exponentially growing bins start from there. A vertical line is also drawn in this case.

I didn't submit a PR yet but I pushed the code in the branch `update-plot_grna_count_distributions`

**[2023-09-07 21:16] Louis Deutsch** _(thread reply)_
I made the vertical line go above the bar but maybe it should go through it so it's actually at the exact value of `threshold`?

**[2023-09-07 21:45] Timothy Barry** _(thread reply)_
As for `threshold`, we threshold counts greater than or equal to this value. Thus, if the `threshold` is 5, then the vertical line should be drawn between 4 and 5. So I think we’ll want to subtract 1 before plotting the threshold.

**[2023-09-07 22:00] Louis Deutsch** _(thread reply)_
Re: threshold, I was thinking I’d make it unambiguous that it is “less than or equal to” by drawing it at the half integer above it. Would drawing it between 4 and 5 be ambiguous with it not including 5?  

**[2023-09-07 22:11] Timothy Barry** _(thread reply)_
I do think that drawing the threshold between 4 and 5 is unambiguous. We want to communicate to the user that cells to the left of the vertical line are categorized as unperturbed while those to the right of the vertical line are categorized as perturbed. IMO drawing the threshold at 4.5 is confusing: are cells with a count of 4 classified as perturbed or unperturbed?

**[2023-09-07 22:44] Louis Deutsch** _(thread reply)_
Hmm so if the threshold is 11 are we saying that all cells with a gRNA count of up to 11 inclusive are unperturbed?

**[2023-09-07 22:48] Louis Deutsch** _(thread reply)_
Here's with the simpler bar labels. I think sideways is still worth it so they don't overlap. That looks a lot cleaner.

Re: threshold, this has it at 11, so my thinking is to draw the line at 11.5 because the counts &lt;= threshold are {0, ..., 11}, so going between integers makes it clear that 11 is below the threshold and 12 is above, versus drawing it on 11 and then it's not as clear (to me) what happens exactly at 11. Is that what you are saying?

**[2023-09-07 23:15] Timothy Barry** _(thread reply)_
Nice! And that’s close to what I’m saying, I think. A threshold of 11 means that all cells up to 10 inclusive are unperturbed.

**[2023-09-08 10:41] Louis Deutsch** _(thread reply)_
ohh ok gotcha! Does this seem right then? Still threshold=11 but now the bar at (threshold-1) is the last that has a width of 1

**[2023-09-08 15:05] Eugene Katsevich** _(thread reply)_
Note that AnnData is the format used by <https://pertpy.readthedocs.io/en/latest/|pertpy>, who have processed many existing datasets into a common format.
  - link: https://pertpy.readthedocs.io/en/latest/

**[2023-09-08 20:07] Timothy Barry** _(thread reply)_
This is what that plot looks like. We probably should remove the “multiple gRNAs” category for high MOI, but otherwise, I cannot think of any changes.

**[2023-09-09 11:22] Eugene Katsevich** _(thread reply)_
It would be good if there were some functionality to draw vertical lines at the thresholds under consideration.

**[2023-09-09 12:08] Eugene Katsevich** _(thread reply)_
I recall John told us that a standard step in choosing QC thresholds is to plot distributions of various cellwise covariates and figure out where to chop off the tails. We should have some functionality to do this, so our users are not tempted to go through Seurat for QC.

**[2023-09-09 16:46] Eugene Katsevich** _(thread reply)_
I image there could be very useful associated plotting functionality to visualize the p-values for each target across gRNAs.

These tasks are getting at the variable effectiveness of gRNAs, a more principled solution to which Joseph happens to be working on in a parallel project.

**[2023-09-09 16:52] Eugene Katsevich**
Another item that might fall under Joseph’s task of testing the code is to make sure our input checks are bulletproof. The data we are asking users to supply is somewhat complicated, and there are many ways to mess up. (We saw an example of this with Evvie, where she did not order her cells in the same way across gene and gRNA expression matrices. Tim has since fixed this issue.) Ideally, Joseph would try every possible combination of malformed input and update the input checking code where necessary to make sure we have all our bases covered. There may also be situations where we need to choose between throwing an error and a warning. For example, if users specify a gRNA in the gRNA group data frame that is not present in the data, we could throw a warning and simply ignore that gRNA or we could throw an error.

**[2023-09-09 17:50] Eugene Katsevich**
One idea for the documentation is to add a flowchart that reviews the main steps. For example, I made such a <https://katsevich-lab.github.io/simulatr/|flowchart for simulatr>. Another source of inspiration is the <https://user-images.githubusercontent.com/45662603/208345270-31443000-600f-4f46-810f-9432e8ed70e0.png|PerturbDecode flowchart>. I could imagine differing levels of detail for this flowchart, but we might want to indicate which functions are optional and which ones are required, and which functions can be run in parallel versus which functions must be run prior to other functions. Especially now that sceptre entails a multi-step process, I think such a flowchart is particularly important at this stage.
  - link: https://katsevich-lab.github.io/simulatr/
  - link: https://user-images.githubusercontent.com/45662603/208345270-31443000-600f-4f46-810f-9432e8ed70e0.png

**[2023-09-09 18:19] Timothy Barry**
<@U0239E8QPH6>, thanks for these helpful suggestions. Are you able to add relevant tasks to the <https://github.com/orgs/Katsevich-Lab/projects/6/views/1|Github project>? Or update existing tasks to incorporate these suggestions? There is a lot to be done, and I myself will be adding a couple tasks for Joseph in the near future. I think we should try to proceed systematically, with Joseph finishing all the plotting tasks (which are fairly high priority) before moving onto other tasks.

Edit: I went ahead and updated the Github task to incorporate Gene’s suggestions.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1

**[2023-09-09 18:26] Timothy Barry** _(thread reply)_
This is already a <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=35281055|task> that has been assigned to me.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=35281055

**[2023-09-09 18:27] Timothy Barry** _(thread reply)_
See this <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=35400061|task>. I have updated the description to include input checking.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=35400061

**[2023-09-09 18:38] Timothy Barry** _(thread reply)_
Joesph already is working on functionality to downsample points in the QQ-plot for speed purposes. I’ve added your suggestion to the <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=37771516|task> as well.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=37771516

**[2023-09-09 19:17] Timothy Barry**
I thought of a couple additional tasks, which I added to the Github project. I have partially completed both of these tasks.

1. *Computing the set of cis pairs for enhancer-targeting screens.* I have started building functionality to automatically compute the set of cis pairs for the user. This feels like an important component of a complete product to me. See this <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=38140402|task>. <@U0239E8QPH6>, I believe that for genes on the + strand, the TSS is located at `start`, while for genes on the - strand, the TSS is located at `end`. Are you able to confirm? It would be nice if Joseph could take this from here.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=38140402

**[2023-09-09 20:01] Timothy Barry**
2. *Add a `plot_cell_covariates` function*. I opened a task/issue <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=38141104|here>.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=38141104

**[2023-09-09 20:09] Eugene Katsevich** _(thread reply)_
&gt; Joseph already is working on functionality to downsample points in the QQ-plot for speed purposes. I’ve added your suggestion to the <https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=37771516|task> as well.
One issue with downsampling discovery pairs is that the extreme tail (which constitutes a small fraction of all pairs) may largely disappear.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1?sortedBy%5Bdirection%5D=asc&amp;sortedBy%5BcolumnId%5D=Status&amp;pane=issue&amp;itemId=37771516

**[2023-09-15 19:19] Eugene Katsevich**
Folks, I’ve <https://github.com/Katsevich-Lab/sceptre3-project/blob/main/documentation-proposal.pdf|taken a crack> at outlining what the new sceptre website would look like. The main structure would be as follows:
• A “Get started” page with an outline of the sceptre steps (communicated primarily via schematic) and a minimal working example of the code.
• Two “SCEPTRE demonstration” vignettes (aimed at plug-and-play users) that would use SCEPTRE to walk through the analyses of the example low- and high-MOI data without delving too deeply into all of the analysis choices.
• Several “SCEPTRE steps” vignettes (aimed at power users), to break down in much more detail what each step entails, what the options are, what accompanying visualizations there are, and how to interpret them. 
What do you guys think of this structure? Once we reach a consensus, we need to figure out how to divide the relevant tasks. I am unlikely to have much bandwidth to actually create this content, but I would be very happy to look over drafts and actively go back and forth with you via Slack. I can imagine that Joseph can help write some of the vignettes.
  - link: https://github.com/Katsevich-Lab/sceptre3-project/blob/main/documentation-proposal.pdf

**[2023-09-16 13:41] Timothy Barry**
Hi everyone, I am back from the conference. It was great. Here are some thoughts.

Main takeaways:

- The IGVF has proposed a map-perturb-predict framework to improve understanding of noncoding elements, (enhancers, silencers, noncoding GWAS SNPs, and noncoding GWAS indels), which involves the following steps: first, identify (i.e., “map”) noncoding elements that are likely causal; second, experimentally perturb these noncoding elements (as well as genes found to be regulated by these noncoding elements) to assess their effect on the transcriptome, proteome, and chromatin landscape; finally, train deep learning models on the perturbation data to predict the effects of noncoding elements not prioritized for experimental perturbation.
- The main data types used for the perturb step (respectively, map step) will be single-cell CRISPR screens and massively parallel reporter assays (respectively, single-cell RNA-seq and single-cell ATAC-seq).
- About one third of the talks involved some discussion of single-cell CRISPR screens. It seems inevitable that we will see a deluge of perturb-seq/tap-seq-style experiments in the coming few years.
- Jesse Engreitz’s lab plans to generate many CRISPRi perturb-seq/tap-seq datasets within the next year and many CRISPR prime editing perturb-seq/tap-seq datasets within the next two years. Other PIs made similar comments.
- I briefly talked to Jesse and Jay Shendure, and they remain quite positive about sceptre. Jesse has directed two of his PhD students (including Ronghao Zhou) to use sceptre to analyze a new in vivo tap-seq dataset. This dataset was generated via the Parse Biosciences kit. I am communicating with Ronhao about the analysis of these data over email.
- Markus Ramste (Stanford) is using tap-seq to link GWAS variants to genes in muscle cells. Markus tried to use sceptre but found it difficult to use. (Although he is not coming from a bioinformatics/statistics/computing background.) Markus want to give sceptre another shot once we have released version 1.0.
- Both Markus and Andreas Gshwind knew that I am working with Evvie to analyze a TAP-seq dataset. Thus, word seems to be spreading about sceptre at Stanford.
- Andrew Allen presented two posters related to single-cell CRISPR screen analysis. One was about combining p-values across gRNAs that target the same site, and the other was about using mixture models to assign gRNAs to cells. Their mixture model seems to be a latent variable GLM like ours. However (and I say this with humility), the first author was not really able to answer questions about how they are fitting their model (e.g., are they using an EM algorithm)?
- I talked pretty extensively with Andrew about their method for pooling p-values across gRNAs and suggested that he look at our updated GLM-EIV paper.
- Luca and Lucas continue to work on their infrastructure for benchmarking single-cell CRISPR screen analysis methods. I likely will help Lucas integrate the latest version of sceptre into their benchmarking pipeline within the next few months. Lucas also continues work on sceptre competitor perturbo, which he says is faster than sceptre, especially for trans analyses.
- As a side note, we might consider adding an inference option “asymptotic,” which would return the standard, asymptotic GLM p-value based on the score statistic. This option would be extremely fast and should handle trans analyses quite easily. Of course, we would be trading statistical accuracy for speed.

People in attendance:

- Jesse Engreitz
- Jay Shendure
- Lior Pachter
- Thomas Norman
- Luca Pinello
- Mike Love
- Andrew Allen
- Andreas Gschwind
- Jian Ma
- Xihong Lin (virtual)

People not in attendance:

- John Morris
- Neville Sanjana

**[2023-09-16 15:06] Timothy Barry**
I think that the proposal for the sceptre 1.0 documentation is quite good. I am a bit conflicted as to whether we need separate vignettes for low and high MOI analyses. Aside from that I am totally on board.

**[2023-09-16 15:12] Eugene Katsevich**
&gt; I think that the proposal for the sceptre 1.0 documentation is quite good. I am a bit conflicted as to whether we need separate vignettes for low and high MOI analyses. Aside from that I am totally on board.
I’m very glad to hear that! I agree that the low- and high-MOI SCEPTRE demonstration vignettes would probably be somewhat similar. Perhaps that’s ok?
&gt; We also at some point will want to add a vignette about steps to take if the method is miscalibrated.
I think this discussion belongs in the SCEPTRE steps vignette for the calibration check.

**[2023-09-16 15:24] Eugene Katsevich** _(thread reply)_
Also, gRNA-to-cell assignment is actually optional for two reasons: Because we have sensible defaults that will be applied automatically and because users may specify their own gRNA-to-cell assignments. Perhaps the schematic should reflect this.

**[2023-09-16 21:09] Timothy Barry**
Hey <@U0524GR916C>,

I wanted to apologize for two pieces of confusion I created. First, I asked you to add functionality to downsample the pairs in the QQ plots so that the plots would render more quickly. I did not realize that this functionality had already been built into sceptre by Gene. Second, I asked you to handle the case when the number of negative control pairs does not equal the number of discovery pairs. However, I took care of this issue myself and failed to let you know that I had done so. I understand from Gene that you lost several hours of work trying to implement these tasks, and I am very sorry about that.

As leader on this project, it is my job to assign you tasks and manage your contributions to the codebase. I have not done the best job at this. (In spite of this, you’ve done really well to acclimate to the codebase, and you’ve made some great initial pushes.) Moving forward, I think it would be best for you to work on one issue at a time. I have marked one of your tasks on the <https://github.com/orgs/Katsevich-Lab/projects/6|project> as “high priority” and all others as “low priority.” I think it would be best for you to finish the task labeled as high priority before moving onto any of the other tasks. (At that point we can label a different task as high priority.) The task currently labeled as high priority is the one that involves writing a unit test to check `n_nonzero_trt`.

I hope that this change will help prevent problems similar to the one that we ran into this week. Thanks for understanding!
  - link: https://github.com/orgs/Katsevich-Lab/projects/6

**[2023-09-21 20:24] Louis Deutsch**
hi all, I've just pushed some functions for mocking the 4 standard data sets we use. They're in a new branch `update/unit-testing` , and live here <https://github.com/Katsevich-Lab/sceptre/blob/update/unit-testing/tests/testthat/helper-mock-data.R|https://github.com/Katsevich-Lab/sceptre/blob/update/unit-testing/tests/testthat/helper-mock-data.R>.

There are 3 of them and there's an example of their usage at the end. I tried to write pretty comprehensive documentation for them.

The main one I'd like feedback on is the one for mocking UMI count matrices (`make_mock_count_matrix`). Right now this makes fake data with 14 different patterns, with the patterns being at either the cell-level or the gene-level. If we have the patterns at the cell level then they are as follows:
1. cells with constant values: one cell with all 0, one cell with all 1, and one cell with all `big` (just some large value)
2. nearly constant cells: cells with all values constant but one
3. cells with linear sequences of values
4. cells with iid Pois(1) values
An example is below showing how the patterns can either be across cells or across genes. Does this seem like it covers the sort of things that we might run into and need to be able to handle? This is all supposed to be *valid* data that things execute correctly on.
  - link: https://github.com/Katsevich-Lab/sceptre/blob/update/unit-testing/tests/testthat/helper-mock-data.R

**[2023-09-27 15:48] Eugene Katsevich**
&gt; I finished a draft of the Get Started page. I would appreciate feedback when you get the chance.
That’s great news! I will read it and provide you with feedback as soon as I can. It appears the “Get Started” vignette was not built when you pushed to GitHub, but I have built it locally.
&gt; It might be fun to turn the documentation into an ebook (via `Bookdown`). The “Get Started” vignette could be the first chapter, and the subsequent vignettes could be the following chapters.
That’s an interesting idea! Would this e-book be in addition to, or instead of, the documentation on the website? Speaking of e-books, I’m in the process of converting and expanding the lab wiki into an <https://ekatsevi.github.io/statistical-computing/|e-book>. Hopefully I will also eventually teach a graduate class at Penn based on these materials. Anyways, I discovered that <https://quarto.org/|Quarto> has replaced bookdown as the software of choice for creating e-books. I’ve been using it and like it so far.
  - link: https://ekatsevi.github.io/statistical-computing/
  - link: https://quarto.org/

**[2023-09-27 16:40] Timothy Barry**
&gt;  It appears the “Get Started” vignette was not built when you pushed to GitHub, but I have built it locally.
Huh, that’s kind of strange. I’ve attached a pdf version of the vignette so we can make sure that we are looking at the same document.

&gt; Would this e-book be in addition to, or instead of, the documentation on the website?
I am thinking that the ebook would be pretty much identical to the documentation on the website, just in book format.

&gt; I’m in the process of converting and expanding the lab wiki into an <https://ekatsevi.github.io/statistical-computing/|e-book>.
Nice!

&gt; discovered that <https://quarto.org/|Quarto> has replaced bookdown as the software of choice for creating e-books.
That’s a good point. JJ Allaire presented Quarto at Bioc2023. If we do turn the documentation into a book, we should probably use Quarto.
  - link: https://ekatsevi.github.io/statistical-computing/
  - link: https://quarto.org/

**[2023-09-30 16:20] Eugene Katsevich**
Hey <@U0239H5UC9E>, Ziang ran across the issue that `sceptre` requires at least on non-targeting gRNA for discovery analyses based on the complement set. This is a bug, right? I have taken the liberty of <https://github.com/Katsevich-Lab/sceptre/commit/4e07f555f4c14a5a7fe4779f0dc56084171895dc|pushing a fix> to the `dev` branch. I hope this is ok with you.
  - link: https://github.com/Katsevich-Lab/sceptre/commit/4e07f555f4c14a5a7fe4779f0dc56084171895dc

**[2023-10-04 14:03] Timothy Barry**
Hi all, I have completed drafts of the first three articles (i.e., import data, set analysis parameters, and assign gRNAs) on the `dev` branch. It would be great to receive feedback when you have the chance.

**[2023-10-04 14:09] Eugene Katsevich**
<@U0239H5UC9E>: Ziang ran across another issue: When `positive_control_pairs` is not passed to `set_analysis_parameters()`, an error is thrown due to the following line:
```sceptre_object@positive_control_pairs &lt;- positive_control_pairs |&gt; dplyr::mutate(grna_target = as.character(grna_target), response_id = as.character(response_id))```
It can be fixed by setting the default argument value `positive_control_pairs = data.frame(grna_target = character(0), response_id = character(0))` rather than `positive_control_pairs = data.frame()`.

**[2023-10-04 14:28] Eugene Katsevich**
<@U0239H5UC9E>: Ziang and I would like to extract the gRNA group to cell assignments from the SCEPTRE object after calling `assign_grnas()`. I noticed that there were multiple fields in the `sceptre_object` that appeared relevant to gRNA assignments: `initial_grna_assignments_list`, `grna_assignments`, and `grna_assignments_raw`. Could you let me know what each of these fields stores?

**[2023-10-04 14:45] Timothy Barry**
We might consider adding an extractor function for this purpose. (Would that be useful?)

`initial_grna_assignments_list` is the list of gRNA-to-cell assignments. `grna_assignments_raw` is the grouped version of this list. And `grna_assignments` is the grouped and QC’ed version of this list.

**[2023-10-05 14:28] Eugene Katsevich**
&gt; `initial_grna_assignments_list` is the list of gRNA-to-cell assignments. `grna_assignments_raw` is the grouped version of this list. And `grna_assignments` is the grouped and QC’ed version of this list.
Ok, got it.

**[2023-10-07 14:57] Louis Deutsch**
In `import_data` currently if `extra_covariates = NULL` is used then the function will allow `grna_matrix` and `response_matrix` to have different column names. For example, the following runs without error:
```  response_matrix <- matrix(0, 4, 2, dimnames = list(letters[1:4], c("A", "B")))
  grna_matrix <- matrix(1, 6, 2, dimnames = list(letters[5:10], c("A", "C")))
  grna_target_data_frame <- data.frame(grna_id = letters[5:10], grna_target = "z")
  import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low",
    extra_covariates = NULL
  )```
This is because in section 9 of `check_import_data_inputs` we only check if the colnames of grna_matrix and response_matrix are identical if `extra_covariates` has a non-NULL rownames attribute. Should we change this so that if response_matrix and grna_matrix have non-null colnames, they must be identical?

This throws the intended error
```  extra_covariates <- data.frame(batch = rep("b1", 2), row.names = c("A", "B"))
  import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low",
    extra_covariates = extra_covariates
  )```

**[2023-10-07 15:56] Louis Deutsch** _(thread reply)_
similarly if the check is not triggered, such as by extra_covariates not having row names, then having invalid response or grna matrix colnames like `"111"` will also not throw an error

**[2023-10-11 15:49] Louis Deutsch**
<@U0239H5UC9E> (starting a new thread for this) After having merged in `dev` to my unit testing branch here's the up-to-date code that does not throw an error but I think should:

```  response_matrix <- matrix(0, 2, 2, dimnames = list(c("response1", "response2"), c("cell1", "cell2")))
  grna_matrix <- matrix(0, 2, 2, dimnames = list(c("grna1", "grna2"), c("cell1", "cell2_mismatched")))
  grna_target_data_frame <- data.frame(grna_id = c("grna1", "grna2"), grna_target = "target")
  import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low"
  )
  import_data(
    response_matrix = response_matrix,
    grna_matrix = grna_matrix,
    grna_target_data_frame = grna_target_data_frame,
    moi = "low",
    extra_covariates = data.frame(batch = c("b", "b")) # if no row.names then this still runs without error
  )```
I propose we modify section 9 of `check_import_data_inputs`  as follows:
```  barcode_list <- list(response_matrix = colnames(response_matrix), grna_matrix = colnames(grna_matrix), 
                       extra_covariates = rownames(extra_covariates))
  non_nulls <- which(!sapply(barcode_list, is.null))
  # if at least two of these have non-null names, they all must be identical
  if(length(non_nulls) >= 2)  { 
    for(i in 1:(length(non_nulls) - 1)) {
      for(j in (i+1):length(non_nulls)) {
        if(!identical(barcode_list[[i]], barcode_list[[j]])) {
          stop(paste0("You have provided cell barcodes in the `", names(barcode_list)[i], "` and `", names(barcode_list)[j], "`. These cell barcodes must be identical across objects."))
        }
      }
    }
  }
  
  # checking that no barcodes are all numbers
  for(i in non_nulls) {
    if(any(grepl(pattern = "^[0-9]+$", x = barcode_list[[i]]))) {
      stop(paste0("Some barcodes in `", names(barcode_list)[i], "` are all numeric which is not allowed."))
    }
  }```
I'm kind of guessing here on what the intended behavior is when we come across any barcodes that match the regex `^[0-9]+$`. Are we trying to forbid those? That's how I interpreted it. The code as it currently stands would actually just not enter the `identical` checks due to everything being a conjunction in

```check_barcodes_provided(response_cell_barcodes) &&
      check_barcodes_provided(grna_cell_barcodes) &&
      check_barcodes_provided(covariate_cell_barcodes)```
If this all sounds good I'll make a PR for this

**[2023-10-11 17:33] Eugene Katsevich**
Hi <@U0239H5UC9E>, I’ve now read through the Get Started vignette. On the whole, it’s quite well-written! I have taken the liberty of pushing some minor edits myself, including re-styling a few code chunks and bolding a key sentence here and there for skimmers to more easily get the main points. Beyond that, I have a few comments:
• For “The Whole Game,” it would be nice to include at least some indication of what the results look like. To this end, perhaps we can show a subset of the discovery results data frame (see below)?
• It might be good to include a vertical line at the mean cells per gRNA in Section III. Also, is the cells per gRNA plot referring to individual gRNAs or gRNA groups?
• Whether in “Get Started” or in the e-book, we might want to have a section devoted to terminology, kind of like <https://www.archrproject.com/bookdown/a-brief-primer-on-atac-seq-terminology.html|this one>. We have a lot of terms (gRNA, target, response, treatment cells, etc) that would be good to define all in one place for users to be able to refer to. 
• It takes a while to build the vignette. We might want to cache the results somehow to facilitate faster iteration.
• For some reason the “standard sceptre pipeline” figure is not centered despite the `fig.align="center"` code.
  - link: https://www.archrproject.com/bookdown/a-brief-primer-on-atac-seq-terminology.html

**[2023-10-11 18:20] Timothy Barry** _(thread reply)_
To clarify, the intended behavior of the function is as follows: if barcodes have been provided, then we verify the barcodes match across `response_matrix`, `grna_matrix`, and `extra_covariates`.

**[2023-10-11 18:31] Eugene Katsevich**
Hey <@U0239H5UC9E>, I’ve taken the liberty of creating a <https://github.com/Katsevich-Lab/sceptre-book|Quarto Book>. I recommend you switch to working within the latter repo.
  - link: https://github.com/Katsevich-Lab/sceptre-book

**[2023-10-12 21:13] Eugene Katsevich** _(thread reply)_
Please see <https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#discovery-pairs|this section> of the documentation for `sceptre`. Quoting from there, cis pairs are “target-response pairs for which the target and response are located on the same chromosome and in close physical proximity to one another” while trans pairs encompass the “entire set of possible target-response pairs.”
  - link: https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#discovery-pairs

**[2023-10-13 15:25] Eugene Katsevich**
I’ve just pushed some minor cosmetic changes to Import Data. One comment on this chapter is that “Import data from a collection of R objects” should discuss the issue of column names for `response_matrix` and `grna_matrix` (and perhaps also row names for `extra_covariates`). My understanding is that if column names are omitted, users must make sure that the columns of `response_matrix` and `grna_matrix` and the rows of `extra_covariates` are given in a consistent order. We should stress this point, because making a mistake here will completely invalidate the results. On the other hand, if column names are given, then `sceptre` will take care of making sure things are in the right order. I would probably encourage users to include the column names and row names to be on the safe side. Tim: What do you think?

**[2023-10-13 21:44] Eugene Katsevich**
Update on my end: I’ve started the <https://timothy-barry.github.io/sceptre-book/glossary.html|Glossary>. Also, <@U0239H5UC9E>, I think the Acknowledgments are short enough that we might as well add them as the last section of <https://timothy-barry.github.io/sceptre-book/|the first chapter>. What do you think?
  - link: https://timothy-barry.github.io/sceptre-book/glossary.html
  - link: https://timothy-barry.github.io/sceptre-book/

**[2023-10-13 21:58] Timothy Barry**
Also, Gene, note that “gRNA group” has been phased out in favor of “target.”

**[2023-10-14 11:00] Eugene Katsevich** _(thread reply)_
I was motivated to add that term due to the section titled <https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#grna-grouping-strategy|“gRNA grouping strategy”>. So the term “gRNA group” appears to still be present, even if it is much less central than before. What do you think?
  - link: https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#grna-grouping-strategy

**[2023-10-14 12:12] Eugene Katsevich** _(thread reply)_
I do see <https://github.com/timothy-barry/sceptre-book/blob/main/rm_old_files.sh|rm_old_files.sh>. How do you recommend we use this?
  - link: https://github.com/timothy-barry/sceptre-book/blob/main/rm_old_files.sh

**[2023-10-14 12:29] Eugene Katsevich** _(thread reply)_
Hi <@U0239H5UC9E>, I believe the reason the cross-package links do not work is that sceptre’s current <https://katsevich-lab.github.io/sceptre/reference/|pkgdown reference page> does not have the new functions. Once we push the latest `sceptre` version to the main branch, I believe the cross-package links will start working.
  - link: https://katsevich-lab.github.io/sceptre/reference/

**[2023-10-14 14:54] Timothy Barry** _(thread reply)_
Could we perhaps change “gRNA group” to “gRNA grouping strategy”?

**[2023-10-14 16:37] Eugene Katsevich**
Hey folks, I’ve just pushed some changes to the e-book. Most notably, I’ve completed the <https://timothy-barry.github.io/sceptre-book/glossary.html|glossary>, which I think came out quite nicely and I hope our users will appreciate. Besides this, I’ve made the following small changes:
• I have removed the Acknowledgments chapter and added it as the last section of <https://timothy-barry.github.io/sceptre-book/|Welcome> instead, linking to “The whole game” in “Navigating this book”. <@U0239H5UC9E>: It is in fact possible to link to unnumbered chapters. In this case, I did so via `[The whole game](sceptre.qmd)`.
• I have added a preview of what the discovery analysis data frame looks like at the end of the initial code demonstration in <https://timothy-barry.github.io/sceptre-book/sceptre.html|The whole game>.
• I have made a couple notes in <https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects|Import data from a collection of R objects> about making sure that the columns of `response_matrix` and `grna_matrix` are in the same order as each other and as the rows of `extra_covariates`. 
Tim and I discussed that I would also make some comments about the reference genome in `construct_cis_pairs()`, but I chose not to do this because that discussion may change if Tim changes the corresponding functionality to allow users to supply their own reference genomes.
  - link: https://timothy-barry.github.io/sceptre-book/glossary.html
  - link: https://timothy-barry.github.io/sceptre-book/
  - link: https://timothy-barry.github.io/sceptre-book/sceptre.html
  - link: https://timothy-barry.github.io/sceptre-book/import-data.html#import-data-from-a-collection-of-r-objects

**[2023-10-16 12:15] Eugene Katsevich**
A couple comments about Chapter 2 that Tim and I have discussed, but wanted to put here for reference:
• I would like to think some more about “gRNA grouping strategy” and what the API should be. Two distinct things are going on here: (1) whether we want to get p-values for individual gRNAs or for groups of gRNAs sharing a target and (2) if we want to get p-values for gRNA groups, how we integrate information across multiple gRNAs with the same target. We might consider letting users separately specify whether they want a grouped analysis or singleton analysis, and if they want a grouped analysis, what kind of grouped analysis they want (i.e. two arguments instead of one). We also want to leave room in the API for more sophisticated grouped analyses (i.e. beyond Bonferroni) where ineffective gRNAs are weeded out or downweighted in a principled fashion. However, such analyses may have tuning parameters, which poses additional challenges.
• Relating to `formula`: We automatically compute `grna_n_nonzero` and include it for high-MOI analyses. Would it make more sense to use as a covariate the actual number of gRNAs we have determined that cell to contain? Is this annoying because we would need to compute this after the gRNA assignment step?

**[2023-10-16 12:22] Timothy Barry** _(thread reply)_
I like the term “integration.” Might a middle ground be to have a single argument called `grna_integration_strategy` that takes values “singleton” (i.e., no integration), “union,” and “bonf”?

**[2023-10-16 12:25] Eugene Katsevich** _(thread reply)_
I think `grna_integration_strategy` is an excellent name for this argument. I’m happy for you to go ahead and make this change to the code. On my end, I can change the corresponding section of the e-book.

**[2023-10-16 14:22] Timothy Barry**
I have updated `grna_grouping_strategy` to `grna_integration_strategy` throughout the package.

**[2023-10-16 14:24] Timothy Barry**
Now we need to change all `grna_grouping_strategy` references to `grna_integration_strategy` throughout the book. <@U0239E8QPH6>, do you think you could handle doing this in section 2.4 (as you are rewriting this section anyway)? I will handle this change throughout all other parts of the book.

**[2023-10-18 12:23] Eugene Katsevich**
<@U0239H5UC9E>: I’ve read through Chapter 3: Assign gRNAs, making minor changes here and there. On the whole, it looks good. One minor question I have is whether you intended the printing of the `sceptre_object` at the end of the Chapter to involve `sceptre_object_lowmoi_mixture` or `sceptre_object_lowmoi`. You have the former, but then it’s unclear why you ran the line `sceptre_object_lowmoi &lt;- sceptre_object_lowmoi_mixture`, because `sceptre_object_lowmoi` is not used anywhere after that.

**[2023-10-18 12:26] Timothy Barry** _(thread reply)_
Thanks. I printed `sceptre_object_lowmoi_mixture` because the MOI field is meaningful for this object but not for `sceptre_object_lowmoi`. If this is confusing then perhaps we can change it.

**[2023-10-18 13:03] Eugene Katsevich** _(thread reply)_
Aren’t these objects the same after calling `sceptre_object_lowmoi &lt;- sceptre_object_lowmoi_mixture`?

**[2023-10-18 16:45] Eugene Katsevich**
Hey <@U0239H5UC9E>, I’ve found two contradictory statements in Chapter 5:
• “The number of negative control gRNAs per group (by default) is set equal to the average number of gRNAs per discovery target.”
• “By default `calibration_group_size` is set to the median number of gRNAs per discovery target.”
If you let me know which of these is correct, I can change the other.

**[2023-10-18 18:00] Eugene Katsevich**
I’ve finished my pass through the <https://timothy-barry.github.io/sceptre-book/|e-book>. I think it’s about ready for release.
  - link: https://timothy-barry.github.io/sceptre-book/

**[2023-10-19 12:34] Louis Deutsch**
There's a bug in `print`  right now:

```&gt; data("lowmoi_example_data")
&gt; sceptre_object &lt;- import_data(
+ response_matrix = lowmoi_example_data$response_matrix,
+ grna_matrix = lowmoi_example_data$grna_matrix,
+ extra_covariates = lowmoi_example_data$extra_covariates,
+ grna_target_data_frame = lowmoi_example_data$grna_target_data_frame,
+ moi = "low")
&gt; print(sceptre_object)
Error in .local(x, ...) : '...' used in an incorrect context```
looks like it's because there's no `...` in the arguments of the print function now:
```setMethod("print", signature = signature("sceptre_object"), function(x) {
  args &lt;- list(...)```

**[2023-10-19 13:35] Louis Deutsch**
Is there some inconsistency in the use of the generic `plot` right now in the e-book? E.g. in section 3.2 `plot_grna_count_distributions` is explicitly called but then in section 3.4 it's just the generic `plot` . Or maybe that's deliberately using the generic to demo that feature?

**[2023-10-19 13:42] Louis Deutsch**
In the first paragraph of section 3.4 there's a slight inaccuracy: it says
```n_grnas_to_plot is the number of (randomly selected) gRNAs to plot; grnas_to_plot is a vector of names of one or more gRNAs to plot```
but the gRNAs selected are only random if `grnas_to_plot` is kept at its default of `NULL` . In the example specific `grnas_to_plot` are provided so there is no random selection

**[2023-10-19 14:16] Timothy Barry** _(thread reply)_
Hey Joseph, `plot()` never calls the function `plot_grna_count_distributions()`. This is because `plot_grna_count_distributions()` is not associated with the output of any particular step in the pipeline.

**[2023-10-19 14:20] Timothy Barry** _(thread reply)_
If `grnas_to_plot` is specified, then `n_grnas_to_plot` is ignored. Would you be able to clarify this (in a single sentence) in the text?

**[2023-10-19 14:30] Louis Deutsch** _(thread reply)_
How about

`n_grnas_to_plot` is the number of gRNAs to plot; `grnas_to_plot` is a vector of names of one or more gRNAs to plot (if not provided, then the gRNAs in the plot are randomly selected)

?

**[2023-10-19 14:31] Timothy Barry** _(thread reply)_
I think that’s good. I might say, “`n_grnas_to_plot` is the number of gRNAs to plot; `grnas_to_plot` is a vector of names of one or more gRNAs to plot. (If the latter argument is not provided, then the gRNAs in the plot are selected randomly.)

**[2023-10-20 15:16] Eugene Katsevich**
Below is a draft of an email to our users. Note that the link to the vignette does not yet work because we have not yet pushed the pkgdown site. Please let me know what you think.

--------------------

Dear sceptre user,

We are very excited to announce the release of sceptre v0.99.0, which serves as a beta version of sceptre v1. This milestone includes the following developments:

• A reimagined user experience based on a modular, object-oriented workflow.
• Further improvements in computational and memory efficiency.
• A unified interface for low- and high-MOI analyses.
• A suite of plotting functions facilitating visualization of each step in the pipeline.
• Expanded support for gRNA assignment (including a new mixture model approach) and quality control.
• Functionality to run sceptre in parallel on multi-core Mac and Linux machines.
• Interoperability with output from 10X Cell Ranger (available now) and Parse Biosciences (coming soon).
• Helper functions to facilitate the construction of _cis_ and _trans_ discovery sets.
• An <https://timothy-barry.github.io/sceptre-book/|e-book> guiding users through the entire process of analyzing their data using sceptre.
Given these developments, sceptre v1 facilitates the entire analysis pipeline for single-cell CRISPR screens, starting from UMI count data obtained from tools like 10X Cell Ranger. To learn more, please consult the new <https://katsevich-lab.github.io/sceptre/articles/sceptre.html|vignette> and <https://timothy-barry.github.io/sceptre-book/|e-book>. We recommend you switch over to this latest release of sceptre as soon as is convenient for you. Please let us know if you have any questions!

Best wishes,
Tim, Joseph, and Gene
  - link: https://timothy-barry.github.io/sceptre-book/
  - link: https://katsevich-lab.github.io/sceptre/articles/sceptre.html
  - link: https://timothy-barry.github.io/sceptre-book/

**[2023-10-29 13:21] Timothy Barry**
Hey Gene, yes, `dev` and `main` currently are equivalent. (This also can be checked on the <https://github.com/Katsevich-Lab/sceptre/branches|“branches” page> of the Github repo.) I do think we should make all updates to dev and only pull into main when we are ready to issue a new release. As we’ve discussed, each push on the main branch should be accompanied by a Github release and an update to the log file.
  - link: https://github.com/Katsevich-Lab/sceptre/branches

**[2023-10-29 14:59] Eugene Katsevich**
Related to point 3: If we use sceptre’s “read from 10X” functionality in the data package, then this will induce a “circular dependency” between the sceptre and sceptreData packages, which <https://github.com/Bioconductor/Contributions#submitting-related-packages|is discouraged>. However, “if circular dependencies are truly unavoidable…Bioconductor will support a ‘Suggest / Depend’ circular dependency. Most often an accompanied data package will “Suggest” the software package, while the software package will “Depend” on the data package.”
  - link: https://github.com/Bioconductor/Contributions#submitting-related-packages

**[2023-10-30 12:09] Eugene Katsevich**
On the topic of power analyses, we may need to clarify that `run_power_check()` is not intended as a way to actually assess the power of `sceptre`. Indeed, the power of `sceptre` for a discovery pair is influenced by the expression level of the corresponding gene and the effect size of the perturbation on the gene; both expression level and effect size are probably different on the positive control pairs and the discovery pairs. We could imagine adding a function to `sceptre` that actually tries to quantify its power via simulations with a given effect size and expression level (and gRNA frequency). I know this is something Jesse is interested in. This isn’t our top priority right now, but I just wanted to make sure we and our users are aware of the distinction between such a function and the current `run_power_check()`.

**[2023-10-30 15:20] Timothy Barry** _(thread reply)_
This seems reasonable to me. So would the workflow be as follows? (I’ll assume Neville using CRISPRko or CRISPRi).

1. Set the analysis parameters to do a singleton analysis with a left-tailed test.
2. Run a calibration check to verify that the method is well-calibrated.
3. Run the power check. Assume that the gRNAs with large p-values are ineffective.
4. Set the analysis parameters again to do a grouped analysis with a two-tailed test. Exclude all ineffective gRNAs from the discovery pairs. Then proceed as normal.
I think this sort of workflow makes sense. Perhaps we can add a section to the book describing this workflow, although I’ve been holding off on strategies for removing ineffective gRNAs, because I know you and Joseph are working on this problem.

**[2023-10-30 15:27] Timothy Barry** _(thread reply)_
4. Fair enough. Fortunately, the Bioconductor people (including Kasper) seem to be excited enough about sceptre that they might allow us to do some things in a nonstandard way. Perhaps for now we can add a runnable example to each man page and include a pdf of the book inside of the package (assuming that the pdf is small enough to fit).

**[2023-11-04 19:02] Eugene Katsevich**
Folks: I’ve been thinking about “future-proofing” the `sceptre` API. In particular, a few changes that might occur in the not-too-distant future are:
• Adding the saddlepoint approximation as an alternative to resampling (based on my project with Ziang).
• Adding more sophisticated strategies for integrating information across gRNAs (based on my project with Joseph).
• Adding support for other test statistics (distilled versus full score statistic, numerator-only score statistic).
• Adding support for other GLMs beyond NB regression.
*Saddlepoint approximation.* The most relevant argument for the addition of the saddlepoint approximation is `fit_parametric_curve`, because both deal with approximations to the resampling distribution. From this perspective, I think it makes sense to change the logical `fit_parametric_curve` argument to a character `resampling_approximation` argument. Then, `resampling_approximation` can take values `"none"` (corresponding to using the resampling distribution as is), `"skew-normal"` (corresponding to fitting a skew-normal distribution, as currently implemented), or `"saddlepoint"` (corresponding to the saddlepoint approximation we might implement in the future). Making this change would also leave room for still other strategies for approximating the resampling distribution (e.g., suppose we wanted an option implementing the skew-t fit from the original `sceptre`).

*More sophisticated gRNA integration strategies.* I think that adding more options to the existing `grna_integration_strategy` parameter will be sufficient.

*Other test statistics.* We have considered a few different test statistics for `sceptre` in the past. Furtheremore, we currently only have a version of the saddlepoint approximation for the test statistic corresponding to the numerator of the score statistic. To accommodate other test statistics, we might at some point want to add another argument called `test_statistic`. The saddlepoint approximation, at least initially, may not be compatible with all test statistic choices.

*Other GLMs.* One current limitation of SCEPTRE is that it applies only to count-based responses. We should be able to easily extend it to continuous responses, including those obtained from various dimensionality-reduction techniques. This might require us to add another argument called `family`, which could take values such as `negative.binomial`, `gaussian`, etc.

Note that among the four changes proposed above, only one (the first) impacts a currently existing argument; the others involve adding new arguments or new options for existing arguments. Therefore, the decision about the first one is more pressing than the decisions about the other three.

**[2023-11-04 19:18] Eugene Katsevich**
One logistical question is what you think would be the best way of organizing and prioritizing to-do items for `sceptre`. Should we continue using our <https://github.com/orgs/Katsevich-Lab/projects/6/views/1|GitHub project>? I have added the <https://github.com/Katsevich-Lab/sceptre/issues/29|issue> Tim opened about running a standard GLM score test to the project. Other to-do list items I’ve added are to accelerate the compilation of the vignette, to make the links on the webpage to the e-book open new tabs, and to update the sentence “support for Parse Biosciences outputs is coming soon” in the vignette.

On a related note, the Signac package has a <https://github.com/orgs/stuart-lab/projects/1/views/6|“roadmap” GitHub project> that lists a set of features and bug fixes that will accompany planned releases. Would having such a thing be useful for us? Or are we not there yet?

On a tangentially related note, now that `sceptre` is approaching stability, we might want to solicit contributions from our users. Maintaining an organized list of tasks needing to be done would be especially helpful if the list of `sceptre` contributors starts to grow.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1
  - link: https://github.com/Katsevich-Lab/sceptre/issues/29
  - link: https://github.com/orgs/stuart-lab/projects/1/views/6

**[2023-11-04 22:55] Timothy Barry**
Thanks Gene for these thoughtful and interesting suggestions.

I’ll start with a more philosophical point before I address the specifics of Gene’s message. In my opinion we have reached the point where we need to be very strict and disciplined about adding new features to the package. I have tried my best to keep the codebase as simple as possible, but as you guys will see, the codebase is fairly insane. I feel like I am close to my limit with respect to being able to reason about the codebase comfortably. Thus, my (somewhat austere) view is this: unless we expect a prospective feature to confer an enormous value add to the package, we probably should lean toward not adding that feature to the package. Furthermore, in deciding whether to add new features to the package, we should prioritize features that implement functionality that is largely orthogonal to the functionality that sceptre currently implements. In my opinion we need to transition away from “statistical researcher” mode and into “software engineer” when thinking about this package, especially as the package sees increasing use out in the real world by biologists, who are using it to do some pretty serious work. Software engineers in general are concerned about <https://en.wikipedia.org/wiki/Feature_creep|feature creep>, which can lead to <https://en.wikipedia.org/wiki/Software_bloat|software bloat>, which in turn can cause a software package to burst at the seams.

OK, rant over. :sweat_smile:

Gene, with respect to your list, the features that strike me as most urgent are (2) and (4): more sophisticated gRNA integration strategies and GLMs beyond NB regression. I know that your and Joseph’s project (which I am excited about) is still in an early stage, so I do not have strong thoughts about that currently. I recall that you once floated the idea that you and Joseph might build a separate package that links up with sceptre; that seems appealing to me as a solution as well. I agree that should at some point add a `family` argument to allow for different exponential family response distributions (maybe NB and Gaussian to start). This seems doable.

Finally, I am on board with continuing to use the Github project to track issues. I think the roadmap idea is interesting, but it does seem a ways off for us. I’ve not yet heard any user express interest in contributing to the package, but I certainly would be openminded about that. We might politely request that such users attend at least one or two developer meetings.
  - link: https://en.wikipedia.org/wiki/Feature_creep
  - link: https://en.wikipedia.org/wiki/Software_bloat

**[2023-11-05 11:05] Timothy Barry** _(thread reply)_
As a brief addendum, the main reason that the codebase is complex is because we are implementing a wide variety of analysis options for users (e.g., low vs. high MOI, permutations vs. CRT, calibration check vs. discovery analysis vs. PC analysis, singleton vs. grouped gRNA integration, complement set vs. NT cells, etc.). These analysis options unfortunately are not independent of one another and can interact in subtle ways, leading to a combinatorial explosion. :exploding_head:

**[2023-11-18 19:04] Eugene Katsevich**
I have found a few bugs in `sceptre`:
1. An error is thrown when `import_data` is called with a factor `response_names`. We might want to cast factor vectors to character vectors. In general, I see that we are not systematically checking that all the specified arguments are of the types specified by our API. Is this something that would be good to do?
2. An error is thrown when `import_data` is called with a tibble `extra_covariates`. The reason for this is the line `v &lt;- extra_covariates[,extra_covariate_name]` in `check_functions.R`. Tibbles, when subsetted in this way, return tibbles rather than vectors. I believe we can fix this by replacing `extra_covariates[,extra_covariate_name]` by `extra_covariates[[extra_covariate_name]]`, which will pull out the column as desired for both data frames and tibbles.
3. If `assign_grnas()` is called but `run_qc()` is not, then calling `run_power_check()` will result in an error saying that `run_qc()` should be called prior to `run_power_check()` instead of running QC automatically. The reason for this is that `skip_assign_grnas_and_run_qc()` checks if _neither_ of `assign_grnas()` and `run_qc()` has been run; it does not handle the case when just one of these has been run.

**[2023-11-29 22:09] Louis Deutsch**
in case anyone needs it, per Tim's suggestion, I made a fork of `ondisc` set to the last commit that has the functions for loading the `.odm` files on the HPC3. I renamed the package to `ondiscFork` to avoid confusion. <https://github.com/jdeu1023/ondiscFork>
  - link: https://github.com/jdeu1023/ondiscFork

**[2023-12-03 10:28] Eugene Katsevich**
<@U0239H5UC9E>: I recall that Replogle et al used a dual-guide CRISPRi construct to get stronger repression. Shouldn’t this affect how we do our gRNA assignment for those data?

**[2023-12-03 14:46] Timothy Barry**
Reading the paper, these two sentences seem relevant.

&gt;  "To maximize CRISPRi efficacy, we used multiplexed CRISPRi libraries in which each element contains two distinct sgRNAs targeting the same gene."
&gt;  "After excluding cells bearing sgRNAs targeting different genes, which are an expected byproduct of lentiviral recombination between sgRNA cassette or doublet encapsulation during scRNA-seq, we retained &gt;2.5 million high-quality cells with a median coverage of &gt;100 cells per perturbation."
I am not exactly sure about the dual-guide business, but it seems that they removed cells containing multiple gRNAs. However, if a given cell contained two gRNAs _targeting the same site_, then that cell was retained.

This relates to cellwise QC as opposed to gRNA assignment. We are not currently carrying out cellwise QC in this way.

**[2023-12-03 14:51] Eugene Katsevich**
What they did, I believe, is deliver lentiviral vectors bearing two gRNAs targeting the same element. So in a sense, the unit at which perturbations are delivered to cells are not _gRNAs_ but _lentiviral vectors_. Each lentiviral vector has _one target_ but achieves knockdown of that target based on _two gRNAs._ So we could imagine detecting these lentiviral vectors, for example, by summing over the counts of the gRNAs they contain.

**[2023-12-03 14:53] Eugene Katsevich**
I suspect your finding that their MOI is close to two might be caused by their dual-guide vector construction.

**[2024-01-11 16:26] Timothy Barry**
Hey all, good meeting today. Here are a few tasks that I think it might be useful for Joseph to work on in the near future.

1. *Data*
• Install the latest version of `ondisc`. Also, install the latest version of the `nf-pipeline` branch of the `sceptre` package.
• Take a look at the "rd7" dataset of the Replogle paper, which is in the HPC3 directory `data/external/replogle-2022/processed/rd7`. Also take a look at the raw data at `data/external/replogle-2022/raw`.
• Take a look at the scripts used to download and process the Replogle data <https://github.com/Katsevich-Lab/import-replogle-2022|here>.
• You can load the `sceptre_object` associated with the rd7 data as follows:
```library(sceptre)
x &lt;- read_ondisc_backed_sceptre_object("sceptre_object.rds", "gene.odm", "grna.odm")```
• You can take a look at the gene/gRNA expression matrix as follows:
```gene_mat &lt;- get_response_matrix(x)
grna_mat &lt;- get_grna_matrix(x)
gene_mat[10,] # etc```
2. *Paper*
• It would be excellent to really get familiar with the Replogle 2022 <https://pubmed.ncbi.nlm.nih.gov/35688146/|paper>.
*3. Testing*
• To the extent that you can write tests that are based on the `nf-pipeline` branch and are focused on the API of the package, I think it would be good to do so.
4. Nextflow
• It would be good to become proficient with Nextflow, both on the cluster and locally. The <https://www.nextflow.io/docs/latest/index.html|docs> would be a good place to start.
Once my revision is done, I should have a bit more bandwidth to think carefully about how to proceed with the sceptre-3 project!
  - link: https://github.com/Katsevich-Lab/import-replogle-2022
  - link: https://pubmed.ncbi.nlm.nih.gov/35688146/
  - link: https://www.nextflow.io/docs/latest/index.html

**[2024-01-17 12:34] Eugene Katsevich**
I’ve created an <https://github.com/Katsevich-Lab/sceptre/issues/73|issue> and added it to our <https://github.com/orgs/Katsevich-Lab/projects/6/views/1|project>. There are now 19 tasks in the project, whose prioritization we need to update. I can try to take a first pass at this, taking into account our previous discussions. Perhaps we can also discuss a bit at tomorrow’s meeting.
  - link: https://github.com/Katsevich-Lab/sceptre/issues/73
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1

**[2024-01-18 11:39] Eugene Katsevich**
Could you give me write access to <https://github.com/timothy-barry/sceptre-pipeline|sceptre-pipeline>? I’m trying to create issues in our project in the appropriate repos, and I can’t create issues in `sceptre-pipeline` because I do not have write access to that repo.
  - link: https://github.com/timothy-barry/sceptre-pipeline

**[2024-01-18 13:28] Eugene Katsevich**
Folks, I spent some time thinking about how we might organize our efforts in this project, and updating our GitHub project accordingly. Here is what I came up with. Let’s discuss during our meeting this afternoon.

*Thoughts on Bioconductor submission.* In our manuscript, we will be advertising `sceptre` as a fairly mature product that is ready for primetime. Given this, I think we should put in the effort to submit `sceptre` to Bioconductor prior to submission. This will take some extra time, but will make our manuscript more compelling. Furthermore, it will save us time down the road, because Bioconductor submission will be more of a headache the longer we wait to do it.

*Remaining phases of `sceptre3` project.* Similarly to what we discussed last week, we can split up the remaining work into phases:
1. *Core functionality:* The first phase is to finish the core functionality of the `sceptre` Nextflow pipeline, including any changes necessary to the `sceptre` and `ondisc` packages. This phase will be led by Tim.
2. *Software polishing:* The second phase is to polish the pipeline and two R packages. By “polishing”, I mean testing, documentation, critical bug fixes, and preparation for Bioconductor submission. I propose for Joseph to lead the polishing for the `sceptre` package and for Tim to lead the polishing for the `sceptre` pipeline and `ondisc` package.
3. *Application:* The third phase is to apply the `sceptre` package to analyze the Replogle and Gasperini datasets, in order to demonstrate that `sceptre` has excellent computational and statistical performance on large datasets. This phase may overlap with the second, to the extent that these datasets can serve as test cases as we polish our software. In this phase of the project, Tim’s involvement will be less critical. We can postpone determining a more exact division of labor until we get there.
4. *Manuscript writing:* The last phase is to write the manuscript. 
*Organization of GitHub project.* Our GitHub project had become somewhat overgrown, because it included not just what we need to do for the `sceptre3` project but also tasks we should consider down the road as part of the longer-term maintenance and improvement of `sceptre`. For this reason, I split the initial project into two: one called <https://github.com/orgs/Katsevich-Lab/projects/6/views/1|sceptre3> and one called <https://github.com/orgs/Katsevich-Lab/projects/7|post-sceptre3>. For the `sceptre3` project, I have broken up the tasks into the aforementioned phases. I have also introduced the “to do next” status as an indicator of which task a given person should be working on in a given week, if they have completed their previous tasks. Finally, I have added preliminary assignees for each task.
  - link: https://github.com/orgs/Katsevich-Lab/projects/6/views/1
  - link: https://github.com/orgs/Katsevich-Lab/projects/7

**[2024-01-26 18:41] Timothy Barry** _(thread reply)_
Hm. Looking at those plots, the two methods seem to be producing pretty highly concordant results. Maybe the estimates of sceptre are slightly smaller (in magnitude)? I've not tried assessing sceptre's fold change estimates on simulated data.

**[2024-01-27 10:44] Timothy Barry** _(thread reply)_
Finally, in most screens (especially in high MOI), we would expect many more control cells than treatment cells.

**[2024-01-27 11:02] Eugene Katsevich** _(thread reply)_
(Assuming it’s full versus distilled that’s causing the discrepancy. I could definitely see this as being the case.)

**[2024-01-27 11:41] Eugene Katsevich** _(thread reply)_
I had ChatGPT whip up a quick simulation, considering a single covariate that is either correlated or uncorrelated with gRNA presence. The finding was that the two-step approach (i.e the distilled approach) produced a downwardly biased estimate in the correlated case but not the uncorrelated case. In Maddie’s simulation though, she’s permuting the perturbation presence indicator, which means one would expect it to be uncorrelated with the covariates…

**[2024-02-08 23:29] Timothy Barry**
<https://github.com/Katsevich-Lab/sceptre/issues/89|This issue> is now resolved, so Joseph, I believe the tests should work now.
  - link: https://github.com/Katsevich-Lab/sceptre/issues/89

**[2024-02-09 00:26] Timothy Barry**
Hi all, I finished a draft of Part II of the book. The easiest way to access the book would be to pull the `sceptre-book` repo and then to change branches to the `large-scale` branch. You can then view the book by clicking on `docs/assign-grnas.html`, for example. I'd appreciate any feedback.

**[2024-02-10 17:20] Eugene Katsevich**
Hey <@U0239H5UC9E>, I’ve made a first pass through Part II of the book. Overall, it’s quite good! Below are a few comments, most of which are minor:
• Chapter 7: Should we add the requisite branches to the installation instructions of  `sceptre` and `ondisc`?
• Chapter 7: “linux” and “sun grid engine” should be capitalized.
• Chapter 7.4: `list.files()` shows more than the three files you claim are there.
• Chapter 8.1: You say that the numbers of pods are chosen by the user. This may be true indirectly, but in fact the user chooses the _sizes_ of the pods rather than their _number_. It would be useful to clarify.
• Chapter 9.1.1: I assume many operating systems come with a Java installation. Perhaps asks folks to check if they have Java of a minimum version prior to installing it.
• Chapter 9.2.1: Is there a reason why the formula object argument to the pipeline cannot be specified as a string?
• Chapter 9.2.1: There is a typo in `response_n_umis_range_uppper`.

**[2024-02-13 22:06] Louis Deutsch**
What happens to cells that are not determined to express any gRNA, NT including? Are they treated as NT? (assume they are not removed in cell-wise QC)

**[2024-02-13 22:58] Louis Deutsch**
Tim, I've just pushed a new commit to the `api-unit-tests` branch adding a new self-contained example of a possible bug to the `temp-bugs-or-exceptions.R` file. It's something to do with when every cell is determined to be expressing some gRNA, and `positive_control_pairs` seem to matter. The error is from `assign_grnas()` and is:
```Error in (function (cl, name, valueClass)  : 
  c("assignment of an object of class "matrix" is not valid for @'initial_grna_assignment_list' in an object of class "sceptre_object"; is(value, \"list\") is not TRUE", "assignment of an object of class "array" is not valid for @'initial_grna_assignment_list' in an object of class "sceptre_object"; is(value, \"list\") is not TRUE")```
(edit) I've pushed another commit with another bug reproduced in that file (more of an unhandled exception). `run_qc()` crashes with the error
```Error in seq.default(start[i], stop[i]) : 'to' must be a finite number```
when no grnas were assigned (eg if no UMI counts exceeded the threshold). (This is the crash I was trying to replicate last week and wasn't finding).

**[2024-02-13 23:09] Louis Deutsch** _(thread reply)_
never mind I've answered this. I forgot that `moi = "high"` overrides `control_group = "nt_cells"` so I was wondering why those cells were always being counted as part of the control group

**[2024-02-14 16:22] Timothy Barry** _(thread reply)_
Thanks Joseph. This is exactly the kind of careful work that really helps to improve the robustness of the package!

Is the idea that there are three (possible) bugs corresponding to the functions `bug_threshold_assign_grnas()`, `bug_error_in_assign_grnas_involving_matrices()`, `bug_no_grnas_assigned_crashes_run_qc()`? Also, to run the functions, I should first define the function `make_mock_base_data_for_testing_run_qc()` contained within `temp-test-run_qc.R` (and call `load_all(helpers=TRUE)`)?

**[2024-02-14 17:57] Timothy Barry** _(thread reply)_
Great, thanks. Note that the bug related to no gRNA being assigned to any cell is fixed only in the nf-pipeline branch.

**[2024-02-20 18:43] Eugene Katsevich** _(thread reply)_
&gt; Also, we of course already have C++ code for computing a GLM score test, so we might want to discuss use of that in this context.
You are referring to `compute_observed_full_statistic_v2()`, right? I am aware of this function and was planning on using it. So all of the heavy lifting has already been done; it’s just a matter of extracting normal tail probabilities and modifying the high-level and medium-level functions accordingly. This task would probably take you about 10 minutes, me about 30 minutes, and the undergrad about 2 months. It’s not super time-sensitive (as far as I know), which is why I decided to assign it to the undergrad.

**[2024-02-22 14:41] Louis Deutsch**
hi I was going to ask about this in the meeting but maybe it'll be helpful to write here first. Pretty much all of `assign_grnas()` 's unit tests involving thresholding are failing and it seems some of the object internals have changed. What exactly changed? This is on my `api-unit-tests` branch so maybe some of this will change once it's merged

**[2024-02-22 14:48] Louis Deutsch** _(thread reply)_
I haven't tried running those tests in a couple weeks, and I think the only changes to actual sceptre code on this branch came from your `fix bug related to threshold` commit a week ago

**[2024-02-24 15:26] Timothy Barry**
Also, I located the HTTP directory containing the Gasperini gRNA UMI counts, so we should be good to go wrt the Gasperini data for the paper.

**[2024-02-24 17:51] Timothy Barry**
By perturbation present and absent do you mean treatment group and control group, respectively? Also, is the idea that we would be tracking four variables: the number of cells with the perturbation, the number of cells without the perturbation, the number of cells where the gene is expressed, and the number of cells where the gene is not expressed? Of course, given a fixed number of cells, there are only three degrees of freedom here: number of cells, number of cells with nonzero gene expression, and number of cells in the treatment group. So is the idea that we would require all three of these variables to exceed some threshold?

 I think it would interesting to explore this idea, but to me this feels like a post-sceptre3-project issue.

**[2024-02-24 18:00] Timothy Barry**
Finally, I do think that this proposal potentially could greatly simplify construction of negative control pairs, especially in the high-MOI setting. (This is a bit technical; I could explain later.)

**[2024-02-24 19:42] Eugene Katsevich**
&gt; By perturbation present and absent do you mean treatment group and control group, respectively?
Yes.
&gt; Also, is the idea that we would be tracking four variables: the number of cells with the perturbation, the number of cells without the perturbation, the number of cells where the gene is expressed, and the number of cells where the gene is not expressed?
Yes.
&gt; Of course, given a fixed number of cells, there are only three degrees of freedom here: number of cells, number of cells with nonzero gene expression, and number of cells in the treatment group. So is the idea that we would require all three of these variables to exceed some threshold?
One way or another, we would want to check all four, but indeed there are relationships among these quantities than may be exploited.
&gt; I think it would interesting to explore this idea, but to me this feels like a post-sceptre3-project issue.
You’re probably right, but this issue is important enough that I wanted to put it on our radars.
&gt; I will say that several users have expressed concern that powerful CRISPRi perturbations can cause gene expression to go to zero. I was under the impression that, in such settings, we would have insufficient power to detect a change, but maybe this is wrong.
I believe this is wrong, and that the users’ concerns are valid. <https://chat.openai.com/share/e/ed93d11a-bebd-45f9-836c-2f5a511eff30|ChatGPT did a quick simulation> of the scenario in my illustration above, finding that the p-value for that specific contingency table is on the order of 1e-15, and that the power to detect an effect at level 0.05 is 100% under random expression of a gene with 50% expression probability in the control group and 0% expression probability in the treatment group if there are 100 treatment and 100 control cells. Of course, we’d want to plug in more realistic values for all these parameters, but the point is that *it is perfectly possible to have high power when one of the cells in your contingency table has a count of zero.*
&gt; Also, doesn’t our current QC amount to a restriction on the top two blocks of the contingency table?
Yes.
&gt; Finally, I do think that this proposal potentially could greatly simplify construction of negative control pairs, especially in the high-MOI setting. (This is a bit technical; I could explain later.)
That’s great to hear!
  - link: https://chat.openai.com/share/e/ed93d11a-bebd-45f9-836c-2f5a511eff30

**[2024-02-25 13:56] Timothy Barry**
I think I am convinced that the strategy you're proposing is better than the one that we are currently using. There are, however, some lingering questions. For example, what should the default thresholds be for these quantities? Furthermore, incorporating this change into the software would be pretty nontrivial. The API would change, the documentation would need to change, the pairwise QC plot likely would need to change, and most challenging of all, about 4-5 low-level C++ functions would have to be rewritten. So yea, I think we'll need to balance the value-add this change would provide against the disruption it would cause.

**[2024-02-26 18:03] Eugene Katsevich**
Hi folks, I did a <https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc.pdf|deep dive into pairwise QC>. This writeup has two main takeaways:
1. I propose a new way of doing pairwise QC that is conceptually and computationally superior to the one we are currently using.
2. The proposed pairwise QC is surprisingly concordant with our current pairwise QC.
I propose we discuss how to proceed with this during our regular weekly meeting on Thursday. We might need to extend this meeting, considering that we also plan to discuss Tim’s pitch to 10x alongside our usual progress updates. I’m pretty flexible on Thursday, so I would be ok with starting earlier than 3pm or finishing later than 4pm. Please let me know what works for you two.
  - link: https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc.pdf

**[2024-02-26 19:11] Eugene Katsevich**
There may be implications of these observations beyond pairwise QC. The simplicity of all of the expressions involved makes it possible to create an experimental design “calculator” that can tell you how many cells you need given a certain number of gRNAs, an MOI, the expression level of a gene of interest, and how many tests you want to run. Or, for a fixed number of cells, how many gRNAs you can afford to assay while preserving power. The calculator can take the form of a Shiny app. This seems to be a missing piece in the field currently, and our users are frequently asking us for information of this kind. The main limitation is that `sceptre` does not actually use the simple 2x2 chi-square test for association analysis, but I think this sort of calculator would be quite good at getting you in the right ballpark.

**[2024-02-27 11:02] Eugene Katsevich**
Hi folks, at the risk of sending you too many messages within a short time, here are some reflections after sleeping on yesterday’s ideas:

*Pairwise QC:*
• *Strengths:* The proposed approach is promising in that it is more easily interpretable and relies on “ancillary statistics” only.
• *Weaknesses:* I overlooked the fact that for low MOI, we use the NT cells as controls, so the total number of cells relevant to a pair is not the same as the total number of cells in the dataset. Therefore, the computations as well as the proposed pairwise QC plot need to be adjusted to accommodate this. Note that the “cells per gene” quantity now becomes a function of the gRNA, which unfortunately hampers the clean decoupling of the computation of cells per gene and cells per gRNA.
• *Outlook:* It still feels to me that the proposed QC metrics are the “right” ones to use, with the adjustments described in the previous bullet point. The sooner we can implement them, the better. In an ideal world, we would include them in the `sceptre3` paper. However, I acknowledge that our bandwidth is extremely limited, and we may not be able to afford it, at least not in an initial submission.
*Power analysis for experimental design:*
• *Strengths:* Having a set of easily interpretable plots that depend on a few knobs that practitioners can adjust in a Shiny app to inform their experimental design seems like it could be very impactful.
• *Weaknesses:* We are assuming for the sake of pairwise QC that the perturbation has the strongest possible effect, and weeding out those pairs whose data are too sparse to detect even this effect. We probably want to add a few more knobs to the experimental design app, like the effect size of the perturbation and the efficacies of the gRNAs, to make it more fully featured. In this sense, the proposed pairwise QC calculations are at most a useful starting point for an experimental design app.
• *Outlook:* This definitely is a completely separate project from `sceptre3`, for which there is currently no bandwidth in the Katsevich Lab. We may revisit this after `sceptre3` is done. We may want to collaborate with another lab on this. Since Xihong is interested in experimental design for single-cell CRISPR screens, she may be interested in this direction.

**[2024-02-27 11:54] Louis Deutsch**
maybe noob question but for low MOI when we use the same cells in every test as the control group, does that induce correlation in the p-values? If so, can that make the false negative and false positive rates higher variance or something?

**[2024-02-27 13:46] Eugene Katsevich**
&gt; You’re saying that, by fixing as much stuff as possible, we are including more information relevant to understanding power?
Yes. Let’s take an example to clarify. Suppose we have 100 treatment cells and 100 control cells, as you assume. However, we have a lowly expressed gene, which is expressed in only 20 of the 200 total cells. Then, the contingency table you present is infeasible. It would yield high power but is not relevant to the situation at hand because it assumes a more highly expressed gene. If we restrict ourselves to a gene expressed in 20 cells, the most extreme contingency table is [0; 100, 20; 80]. So we have to take into account not just the cells/gRNA but also the cells/gene. This is why my plots have these two key quantities as their x- and y-axes.

**[2024-02-28 13:10] Timothy Barry** _(thread reply)_
OK, sounds good. I'll switch up the assignee on Github.

**[2024-02-28 15:04] Louis Deutsch**
Also here's the code coverage report as it currently stands:

```sceptre Coverage: 49.35%
R/import_functs.R: 0.00%
R/mixture_model_functs.R: 0.00%
R/ondisc_functs.R: 0.00%
R/pair_constructor_functs.R: 0.00%
R/plotting_functions.R: 0.00%
R/qq_plot_helpers.R: 0.00%
R/zzz.R: 0.00%
src/check_functions.cpp: 0.00%
src/mixture_functs.cpp: 0.00%
src/pair_construction_functs.cpp: 0.00%
R/s4_helpers.R: 17.65%
R/getters_and_setters.R: 26.09%
R/s4_analysis_functs_2.R: 54.03%
src/generate_samples_functions.cpp: 56.18%
R/association_testing_functions_v2.R: 56.34%
R/assign_grna_functions.R: 61.48%
R/aux_functions.R: 65.57%
R/check_functions.R: 69.71%
R/pairwise_qc_functs.R: 72.22%
R/medium_level_functs_v2.R: 72.44%
R/preprocessing_functions.R: 73.77%
R/precomputation_functions.R: 75.47%
R/neg_control_functions.R: 80.26%
R/s4_analysis_functs_1.R: 87.90%
src/shared_low_level_functions.cpp: 90.62%
src/sparse_matrix_functs.cpp: 96.05%
src/low_level_full_test.cpp: 98.44%
R/cellwise_qc_functs.R: 100.00%
R/test_functs.R: 100.00%
src/compute_nb_size_param_functions.cpp: 100.00%
src/negative_control_functions.cpp: 100.00%
src/pairwise_qc.cpp: 100.00%```

**[2024-02-28 15:44] Timothy Barry**
Hey Joseph, as I mentioned, I am updating the tests to be compatible with the most recent changes to the package. I am getting the following error:

```── Error ('test-check_functions.R:118:3'): check_import_data_inputs ────────────
Error in `inset(valid_grna_target_data_frame, 1, 1, "bad&amp;id")`: could not find function "inset"```
Do you know what the `inset()` function is? Apparently, it is not found.

**[2024-02-28 20:53] Louis Deutsch**
Hi Tim, can you help me interpret the following?
```> scep@grna_assignments
$grna_group_idxs
$grna_group_idxs$t1
 [1]  1  2  3  4  5  6  7  8  9 10

$grna_group_idxs$t2
 [1] 11 12 13 14 15 16 17 18 19 20

$grna_group_idxs$t3
integer(0)


$indiv_nt_grna_idxs
$indiv_nt_grna_idxs$nt1
 [1]  1  2  3  4  5  6  7  8  9 10

$indiv_nt_grna_idxs$nt2
 [1] 11 12 13 14 15 16 17 18 19 20


$all_nt_idxs
 [1] 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50```
My question: this is from an example where cells `1:10` express t1, `11:20` express t2, cells `31:40`  express nt1, and cells `41:50` express nt2. For the idx in `indiv_nt_grna_idxs$nt1`, are those in reference to `all_nt_idxs`? so if i wanted to get the cells expressing nt1, I'd do `scep@grna_assignments$all_nt_idx[scep@grna_assignments$indiv_nt_grna_idxs$nt1]`? I want to confirm that the idx in the `indiv_nt_grna_idxs` are not positions in the `grna_matrix` but rather for subsetting the `all_nt_idxs` element

(edit) also this seems to only be the behavior of this object when the control group is `"nt_cells"`? When the control group is the complement it doesn't seem like I need to do this double indexing (and the `all_nt_idxs` element isn't there so i guess that one goes without saying)

**[2024-02-28 22:30] Louis Deutsch**
Testing update: after adding in my tests for `run_discovery_analysis()` and `method="mixture"` we're only at ~53% coverage. I'm adding tests now for the plotting functions. I'm just checking the basic structures, like the number of panels, titles, etc are there, to at least confirm that future changes don't mangle them. I'm hoping without too much effort, and without it being too artificial, I can get us up to ~70% coverage

From ChatGPT:
```Many organizations and industry best practices suggest that a test coverage of 70-80% is a healthy benchmark for most projects. This range is often considered a good balance between effort and the confidence level in the codebase.```
(update) this is where we are with the tests for the plotting functions
```sceptre Coverage: 69.51%
R/import_functs.R: 0.00%
R/ondisc_functs.R: 0.00%
R/pair_constructor_functs.R: 0.00%
R/zzz.R: 0.00%
src/check_functions.cpp: 0.00%
src/mixture_functs.cpp: 0.00%
src/pair_construction_functs.cpp: 0.00%
R/s4_helpers.R: 17.65%
R/getters_and_setters.R: 26.09%
R/aux_functions.R: 65.57%
R/assign_grna_functions.R: 66.39%
R/mixture_model_functs.R: 67.24%
R/association_testing_functions_v2.R: 69.72%
R/s4_analysis_functs_2.R: 71.09%
R/pairwise_qc_functs.R: 72.22%
src/generate_samples_functions.cpp: 74.16%
R/precomputation_functions.R: 75.47%
R/check_functions.R: 76.00%
R/preprocessing_functions.R: 79.51%
R/neg_control_functions.R: 80.26%
R/medium_level_functs_v2.R: 82.05%
R/plotting_functions.R: 86.26%
R/s4_analysis_functs_1.R: 91.94%
src/shared_low_level_functions.cpp: 94.85%
src/sparse_matrix_functs.cpp: 96.05%
src/low_level_full_test.cpp: 98.44%
R/cellwise_qc_functs.R: 100.00%
R/qq_plot_helpers.R: 100.00%
R/test_functs.R: 100.00%
src/compute_nb_size_param_functions.cpp: 100.00%
src/negative_control_functions.cpp: 100.00%
src/pairwise_qc.cpp: 100.00%```

**[2024-02-28 22:43] Timothy Barry** _(thread reply)_
Yes, `indiv_nt_grna_idxs$nt1` and `indiv_nt_grna_idxs$nt2` are with reference to `all_nt_idxs`. Yes, you would obtain the cells expressing nt1 via `scep@grna_assignments$all_nt_idx[scep@grna_assignments$indiv_nt_grna_idxs$nt1]` (assuming you have not yet run quality control).

Finally, yes, this double indexing business is for the NT cell control group only; it does not apply to the complement set control group.

**[2024-02-28 22:44] Timothy Barry** _(thread reply)_
Yes, exactly, once you have run QC, you would need to do triple indexing (assuming NT cells control group).

**[2024-02-29 11:35] Louis Deutsch** _(thread reply)_
yes, I'll be fixing those fails in thenext couple hours. It's some parts of the `assign_grnas()` tests that fail due to  some of the internals having changed since i wrote them. And ok great, glad to hear it's been merged! I'll just submit a new PR with the same branch when i finish this

**[2024-02-29 13:43] Louis Deutsch**
Tim, is this a bug? From one of my `assign_grnas()` unit tests that currently is failing. This one uses the maximum method if that matters. Having `NA` as the sublist names doesn't seem intended?

```scep_low_all_1@grna_assignments_raw
$grna_group_idxs
$grna_group_idxs$t1_c1_d1
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

$grna_group_idxs$t1_c2_d1
NULL

$grna_group_idxs$t1_c3_d1
NULL

$grna_group_idxs$t2_c3_d1
NULL


$indiv_nt_grna_idxs
$indiv_nt_grna_idxs$&lt;NA&gt;
NULL

$indiv_nt_grna_idxs$&lt;NA&gt;
NULL

$indiv_nt_grna_idxs$&lt;NA&gt;
NULL```


**[2024-02-29 13:57] Louis Deutsch**
It also looks like `scep_low_all_1@initial_grna_assignment_list` has changed (the failures of the `assign_grnas()` tests I think are all due to a couple of the values of these slots having changed). Before it had one element per gRNA now in this case it just has a single element:
```$g1_t1_c1_d1
 [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24```
is that intended?

**[2024-02-29 14:01] Louis Deutsch** _(thread reply)_
every element of the grna_matrix is 1, and the method is maximum. This is the test named `"assign_grnas method=maximum moi=low grna_matrix all 1"`

**[2024-02-29 18:07] Louis Deutsch**
Hi <@U0239H5UC9E> here is a self-contained example of the 5 tests that are failing for `assign_grnas()` which will hopefully help guide your search! Some might actually be good failures where it was wrong for them to have been passing earlier.

```## making the data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
num_grnas &lt;- 6
num_nt &lt;- 3
nrow_grna_matrix &lt;- num_grnas + num_nt
num_cells &lt;- 12
num_responses &lt;- 13

grna_target_data_frame &lt;- data.frame(
  grna_id = c(paste0("grna_", 1:num_grnas), paste0("nt_", 1:num_nt)),
  grna_target = rep(c("t1", "t2", "non-targeting"), c(num_grnas/2, num_grnas/2, num_nt)),
  chr = "", start = 0, end = 1
)
grna_matrix_all_1 &lt;- matrix(1,  nrow = nrow_grna_matrix, ncol = num_cells) |&gt;
  `rownames&lt;-`(grna_target_data_frame$grna_id)

response_matrix &lt;- rpois(num_responses * num_cells, lambda = 1) |&gt; 
  matrix(ncol = num_cells) |&gt;
  `rownames&lt;-`(paste0("response_", 1:num_responses))

discovery_pairs &lt;- data.frame(
  grna_target = "t1",
  response_id = "response_1"
)

## tests with constant `grna_matrix`: 2 failures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

scep_low_all_1 &lt;- import_data(
  grna_matrix = grna_matrix_all_1,
  response_matrix = response_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low"
) |&gt;
  set_analysis_parameters(
    discovery_pairs = discovery_pairs
  ) |&gt;
  assign_grnas(method = "maximum")

## FAILURE 1
# should there be one element per grna_id? This broke in the recent nf-pipeline merge
expect_equal(
  scep_low_all_1@initial_grna_assignment_list,
  lapply(1:nrow_grna_matrix, function(i) if(i == 1) 1:num_cells else integer(0)) |&gt;
    setNames(grna_target_data_frame$grna_id)
)

## FAILURE 2
# this has some NAs as names for sublists which seems wrong
expect_false(scep_low_all_1@grna_assignments_raw$indiv_nt_grna_idxs |&gt; names() |&gt; <http://is.na|is.na>() |&gt; any())


## tests with a `grna_matrix` where some cells have a clear max but others are constant ~~~~~~~~~~
## there are 3 failures here

grna_matrix_clear_max &lt;- grna_matrix_all_1
for(i in 1:nrow_grna_matrix) grna_matrix_clear_max[i,i] &lt;- 100
grna_matrix_clear_max[,(nrow_grna_matrix + 1):num_cells] &lt;- 0

scep_low_clear_max &lt;- import_data(
  grna_matrix = grna_matrix_clear_max,
  response_matrix = response_matrix,
  grna_target_data_frame = grna_target_data_frame,
  moi = "low"
) |&gt;
  set_analysis_parameters(
    discovery_pairs = discovery_pairs
  ) |&gt;
  assign_grnas(method = "maximum")

## FAILURE 3
# Maybe this one actually shouldn't have worked before, and now is actually correct
# even though the test is failing?
expect_equal(
  scep_low_clear_max@initial_grna_assignment_list,
  as.list(1:nrow(grna_target_data_frame))  |&gt;
    setNames(grna_target_data_frame$grna_id)
)

## FAILURE 4
# all of these cells were considered to have multiple gRNAs before;
# if the `grna_matrix` has a constant value of 1, instead of 0, these appear
# correctly. So having columns of all 0 seems specifically bad
expect_equal(scep_low_clear_max@cells_w_multiple_grnas, (nrow_grna_matrix + 1):num_cells)


## FAILURE 5
# maybe this one also should have actually failed before, and is now correct?
expected_grna_assignments_raw &lt;- list(
  grna_group_idxs = list(t1 = c(1,2,3), t2 = c(4,5,6)),
  indiv_nt_grna_idxs = list(nt_1 = 7, nt_2 = 8, nt_3 = 9)
)
expect_equal(scep_low_clear_max@grna_assignments_raw, expected_grna_assignments_raw)
# either way, kinda weird how cell 10:12 got put right after cell 1
# in `scep_low_clear_max@grna_assignments_raw$grna_group_idxs$t1` 
# instead of at the end, where the ordering would naturally be?```
  - link: http://is.na

**[2024-02-29 22:44] Timothy Barry**
I cleared failures 1-2. For failures 3-5, we should update the test; the sceptre output is correct. Nice job catching the bug for tests 1-2!

As a side note, we probably want to deal with the situation where there clearly are no cells that contain a gRNA as part of cellwise QC, but we have not done this yet.

**[2024-02-29 22:44] Timothy Barry**
(I mean this in the context of the maximum assignment strategy.)

**[2024-03-01 02:03] Louis Deutsch**
I just submitted a PR that is to merge the one extra commit I didn't push earlier that fixes a bug in the plotting function tests. The assign_grna tests are unchanged and I won't modify them until notified

**[2024-03-01 20:18] Timothy Barry**
Hey guys (especially <@U0524GR916C>), I made an update to the package that I hope will simplify the gRNA assignment step. This update should enable us to pass all of Joseph's tests.

First, within the `sceptre_object`, I changed the field `cells_w_multiple_grnas` to `cells_w_zero_or_twoplus_grnas`. As the name change implies, I am now flagging cells that contain zero or 2+ gRNAs (in the low-MOI setting). Such cell will be removed as part of cellwise QC.

Within the context of the maximum assignment method, cells with a total UMI count across gRNAs less than `min_grna_n_umis_threshold` (default value 5) are flagged as containing zero gRNAs. As mentioned above, such cells are removed as part of cellwise QC.  `min_grna_n_umis_threshold` can be passed as an argument to `assign_grnas()`.

Additionally, within the context of the thresholding and mixture methods, cells not assigned any gRNA are flagged as containing zero gRNAs and included in the `cells_w_multiple_grnas` vector. Thus, such cells are now removed as part of cellwise QC. (This should not change the results at all when using the NT cells as the control group.)

Finally, all of this is documented in the `large-scale` branch of the sceptre-book repo. I hope all of this is clear; please let me know if I can clarify some aspect of this!

**[2024-03-05 10:25] Eugene Katsevich**
Here’s a draft of the 10x blog post, also available on GitHub at <https://github.com/Katsevich-Lab/sceptre3-project/blob/main/10x-analysis-guide-draft.docx|sceptre3-project>. It’s currently about 750 words, which is within the range of the word counts in the three blog posts Juan shared with us. It would be great if one or both of you could look this over and provide some feedback. I have a couple questions related to what I wrote, which perhaps Tim can answer:
• I wrote “To use `sceptre` on your 10x data, start by calling `cellranger count` on your FASTQ data to obtain a directory of UMI count output files. Supply the path to the Cell Ranger output to `sceptre`, which can directly import your data into `sceptre_object` format.” Is this a sufficiently accurate description of the junction between Cell Ranger and `sceptre`? Do we need to say something about `.mtx` versus `.h5`?
• I wrote “Besides its computational speed, `sceptre` comes with a light memory footprint. The memory required to analyze a dataset is not much more than the memory taken up by the data itself.” Is this an accurate statement? Have we quantified this?
  - link: https://github.com/Katsevich-Lab/sceptre3-project/blob/main/10x-analysis-guide-draft.docx

**[2024-03-07 14:55] Eugene Katsevich**
I’ve added a <https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc-v2.pdf|new writeup> on pairwise QC, in which I treat the issues of sidedness and MOI more carefully. Here are the two main conclusions:
1. The current QC results in a selection bias, leading to inflated p-values for right-sided and two-sided tests. This may explain the inflation we’ve seen in the Replogle 2022 analysis!
2. There is a fundamental asymmetry between left- and right-tailed tests, with left-tailed tests being more vulnerable to low power. Our QC should accommodate this. Due to this asymmetry, the current QC is too stringent for right-tailed tests.
We should check if indeed the current pairwise QC is the culprit for the inflation in the Replogle analysis; if so, this will (fortunately or unfortunately) force our hand to update our pairwise QC. In the writeup linked above, I propose a sidedness-aware pairwise QC that should address both of the issues listed above.
  - link: https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc-v2.pdf

**[2024-03-07 18:20] Eugene Katsevich**
<https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc-v3.pdf|Here> is the third installation in my series on pairwise QC. It’s shorter and more digestible than the second. In short, it confirms the relationship of the `sceptre` p-value with `n_nonzero_trt` and explains why the effect of this selection bias may be less pronounced on less sparse datasets (like those from `sceptre2`) but more pronounced on more sparse datasets (like Replogle).
  - link: https://github.com/Katsevich-Lab/sceptre3-project/blob/main/pairwise-qc-v3.pdf

**[2024-03-07 18:24] Eugene Katsevich**
<@U0239H5UC9E>: Please take a look. Producing a plot like the attached for Replogle would be very illuminating.

**[2024-03-11 17:03] Eugene Katsevich** _(thread reply)_
Yes, I should have filtered on both. I was implicitly assuming that the filter on `n_nonzero_cntrl` doesn’t do much.

**[2024-03-13 21:31] Eugene Katsevich**
There are at least a couple negative control gRNAs that occur more than once in the top ten. It might be the case that some negative controls “didn’t work” in the sense of having some genomic effect. But removing negative control gRNAs that exhibit miscalibration is a slippery slope I view as a last resort.

**[2024-03-13 21:32] Eugene Katsevich**
Regardless, Tim, stratifying the QQ plot by gRNA may be illuminating.

**[2024-03-13 21:33] Eugene Katsevich**
It might be that the problematic gRNAs are the ones most correlated with the omitted confounder. 

**[2024-03-26 12:42] Eugene Katsevich** _(thread reply)_
This seems sensible. Technically `sceptre` performs two regressions (one for the responses and one for the gRNAs). Do we want to choose a more specific name than `regression_method`, like `response_regression_method`?

**[2024-03-26 12:43] Timothy Barry** _(thread reply)_
That's a good point. Yes, I think `response_regression_method` would be the right name (especially if we want to add alternate gRNA regression methods, e.g. XGBoost or whatever, in the future).

**[2024-03-26 16:11] Eugene Katsevich**
<@U0239H5UC9E>: I have found that `sceptre` is working somewhat more slowly than I remember; the calibration check for the low-MOI example took about 2 minutes and the discovery analysis took about 4.5 minutes. Both were run in parallel. Is this slower than `sceptre` was working before, or am I just misremembering?

**[2024-03-26 16:53] Timothy Barry**
Hm. I just ran a calibration check on the example low-MOI data using eight processors. It took 4 seconds.

**[2024-03-28 20:46] Louis Deutsch** _(thread reply)_
I just did some googling and I’m not seeing any canonical references for it. It’s not really a “method” it’s just a hacky thing I’ve seen people do back in my days of fitting black box predictive models. 

<https://www.researchgate.net/profile/Kedar-Potdar-2/publication/320465713_A_Comparative_Study_of_Categorical_Variable_Encoding_Techniques_for_Neural_Network_Classifiers/links/59e6f9554585151e5465859c/A-Comparative-Study-of-Categorical-Variable-Encoding-Techniques-for-Neural-Network-Classifiers.pdf?origin=publication_detail&amp;_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6InB1YmxpY2F0aW9uIiwicGFnZSI6InB1YmxpY2F0aW9uRG93bmxvYWQiLCJwcmV2aW91c1BhZ2UiOiJwdWJsaWNhdGlvbiJ9fQ|This> is the most official reference I’ve found but it barely describes it (it’s “binary encoding”)

The idea is if a categorical has K levels we number them 0 to K-1, but in binary, then treat those as strings with 0 padding so everything is the same length. Then split those strings into their digits and the encoding is the log2(K) columns formed that way. So we get comparisons between log2(K) groups but it’s completely arbitrary which groups we form 
  - link: https://www.researchgate.net/profile/Kedar-Potdar-2/publication/320465713_A_Comparative_Study_of_Categorical_Variable_Encoding_Techniques_for_Neural_Network_Classifiers/links/59e6f9554585151e5465859c/A-Comparative-Study-of-Categorical-Variable-Encoding-Techniques-for-Neural-Network-Classifiers.pdf?origin=publication_detail&amp;_tp=eyJjb250ZXh0Ijp7ImZpcnN0UGFnZSI6InB1YmxpY2F0aW9uIiwicGFnZSI6InB1YmxpY2F0aW9uRG93bmxvYWQiLCJwcmV2aW91c1BhZ2UiOiJwdWJsaWNhdGlvbiJ9fQ

**[2024-03-30 13:10] Timothy Barry**
Also, Joseph, code for running the rd7 analysis is contained within the "rd7_analysis" subdirectory of the `sceptre3-project-v2` repo (<https://github.com/Katsevich-Lab/sceptre3-project-v2>).
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2

**[2024-03-31 22:11] Timothy Barry**
In the interest of getting the ball rolling with respect to the Bioconductor submission, I have updated the main branch of `ondisc`. I labeled the release v1.2.0. Here is the Github <https://github.com/timothy-barry/ondisc|page> and package <https://timothy-barry.github.io/ondisc/index.html|website>. I'd be quite happy to incorporate any edits before Bioconductor submission!
  - link: https://github.com/timothy-barry/ondisc
  - link: https://timothy-barry.github.io/ondisc/index.html

**[2024-04-01 12:07] Eugene Katsevich**
Hey <@U0239H5UC9E>, I took a look at your recent updates to the ondisc website and sceptre book. On the whole, everything is excellent. This thing is _really_ coming together. What a massive piece of work. My comments are pretty minor:

*`ondisc` website:*
• The changelog entry has “v0.1.0” instead of “v1.2.0".
• The vignette says “See the frequently asked questions page for tips on installing ondisc such that it runs as fast as possible.” Currently, the FAQ in our book has an entry called “How can I improve the speed of sceptre?” Do you intend to update this entry to cover ondisc as well? 
• In the vignette, while discussing `create_odm_from_cellranger()`, you state “As part of importing the data, `create_odm_from_r_matrix()` computes the cell-wise covariates.” Did you intend to refer to `create_odm_from_cellranger()`? 
• In the vignette, you state “`.odm` files are portable. Thus, a user can create an `.odm` file on one computer, move the `.odm` file to another computer, and then open the `.odm` file on the second computer.” Perhaps it’s worth mentioning that R `odm` objects are not portable? So if someone wants to move an `odm` from one machine to another, we recommend moving the backing file and reinitializing the `odm` rather than saving the `odm` file as an R object and moving it to another machine.
*`sceptre` book:*
• Please add funding information in the acknowledgments on the front page. I recommend: “We also thank the funders that made `sceptre` possible; the methods and software development were supported by an award from Analytics at Wharton, NSF grants DMS 2113072 and DMS 2310654, and NIH grants…”
• Please change “10X” to “10x” everywhere.
• In Appendix C there is a typo: “couppled”.
• In Section 1.5: “on-disk” should be “on disk”, i.e. the hyphen is unnecessary.
• Section 8.2 states “The gRNA pod size and pair pod size correspond to _p_ and _r_ in Figure 8.2).” Do I understand correctly that _p_ and _r_ are the numbers of pods rather than the sizes of the pods?
• Section 8.4 states “Indeed, we typically are less concerned about stringent calibration of the p-values far out into the tail of the distribution in the context of a massive-scale trans analysis.” I recommend we justify this statement, as it might not be apparent to our users.
• In the note at the end of Chapter 8, `set_analysis_parameter()` is missing an “s”.

**[2024-04-01 13:22] Timothy Barry** _(thread reply)_
Thanks for the comments! I incorporated all of these edits into the `ondisc`<https://timothy-barry.github.io/ondisc/| website> and sceptre book.

(Note that I am planning on updating the `main` branch of `sceptre-book` at the same time we update the `main` branch of `sceptre`. If you'd like to see the latest version of `sceptre-book`, you can check out the `at-scale` branch. Indeed, the latest version of the FAQ page of `sceptre-book` includes some discussion about optimizing the C++ code within `ondisc`.)
  - link: https://timothy-barry.github.io/ondisc/

**[2024-04-01 14:49] Timothy Barry** _(thread reply)_
Hm. Yes, I pushed. <https://github.com/timothy-barry/sceptre-book/blob/large-scale/faq.qmd|Here> is what `faq.qmd` looks like on the `large-scale` branch on Github.
  - link: https://github.com/timothy-barry/sceptre-book/blob/large-scale/faq.qmd

**[2024-04-02 13:14] Eugene Katsevich**
FYI <@U0524GR916C>: Our low-MOI paper was accepted to Genome Biology! <@U0239H5UC9E>: Should we ask Juan to update the reference to this paper in the analysis guide?

**[2024-04-02 14:56] Eugene Katsevich**
Hey <@U0239H5UC9E>: What is the current barrier to running CRT on an ondisc-backed `sceptre_object`? It’s not fully clear to me why any additional code modification is necessary, given that you’re pulling data out of the response matrix using `load_csr_row()`, which works with both matrices and odms. Does it have to do with the way you interact with the gRNA data? I ask because I am getting to the point in my project with the undergraduate where we need to analyze the Gasperini data. I’m thinking it would be very convenient to analyze the Gasperini data using an ondisc-backed `sceptre_object`, but also that, ideally, we’d use CRT as a reference point.

**[2024-04-02 16:01] Eugene Katsevich** _(thread reply)_
I’m scheduled to meet with my student tomorrow morning to give her an assignment for next week, which ideally would include the analysis of the full Gasperini data. Do you think it is doable for me to update the `sceptre` package to accommodate the CRT for ondisc-backed sceptre objects and to import the Gasperini data into the new `ondisc` format during the remainder of today?

**[2024-04-04 13:33] Eugene Katsevich**
I see that the `sceptre_object` <https://github.com/Katsevich-Lab/sceptre/blob/main/R/s4_class.R#L26|class> still has a field called `run_permutations` rather than `resampling_mechanism`. Is this intentional?
  - link: https://github.com/Katsevich-Lab/sceptre/blob/main/R/s4_class.R#L26

**[2024-04-04 15:43] Eugene Katsevich**
<@U0239H5UC9E>: Here’s <https://github.com/Katsevich-Lab/import-gasperini-2019-v3/tree/main|import-gasperini-2019-v3>. So far, I have simply copied and pasted the `import-gasperini-2019-v2` repo, but I replaced `create_odms_3.R` with <https://github.com/Katsevich-Lab/import-gasperini-2019-v3/blob/main/at-scale/create_sceptre_object_3.R|create_sceptre_object_3.R>.
  - link: https://github.com/Katsevich-Lab/import-gasperini-2019-v3/tree/main
  - link: https://github.com/Katsevich-Lab/import-gasperini-2019-v3/blob/main/at-scale/create_sceptre_object_3.R

**[2024-04-04 18:33] Louis Deutsch**
(edit for posterity: this is wrong, ignore this) Hi Tim, forgot to mention this earlier re: example local nextflow pipeline. I think we can't use the `-r` flag locally so it should be

```nextflow run timothy-barry/sceptre-pipeline main \
 --sceptre_object_fp $sceptre_object_fp \
 --response_odm_fp $response_odm_fp \
 --grna_odm_fp $grna_odm_fp \
 --output_directory $output_directory \
 --grna_assignment_method mixture \
 --pair_pod_size 1000 \
 --grna_pod_size 25 \
 --trial```
in §7.5 of the ebook

**[2024-04-04 20:10] Louis Deutsch** _(thread reply)_
Yes I am. Before you updated `sceptre-pipeline` I was getting the error `Revision option cannot be used when running a local script` and removing the `-r` fixed it. I was able to fully run the trial nextflow pipeline with the few changes I mentioned to you (I also needed to allocate more time to the assign_grnas() step).

Now I'm getting some weird errors that I'm still making sense of. I'll follow up in a few min as I work thru them

**[2024-04-04 20:45] Timothy Barry** _(thread reply)_
OK. Is there anything that we should change within the pipeline or book, do you think? For example, should we request more time for gRNA assignment?

**[2024-04-06 18:59] Eugene Katsevich**
Hey folks, FYI I <https://github.com/Katsevich-Lab/import-gasperini-2019-v3|imported> the full Gasperini data into both ondisc-backed and non-ondisc-backed `sceptre_object` form. These are both located at `~/data/external/gasperini-2019-v3/at-scale/processed`. The README of the import directory linked above has some documentation.
  - link: https://github.com/Katsevich-Lab/import-gasperini-2019-v3

**[2024-04-08 11:29] Timothy Barry**
We have our <https://github.com/Katsevich-Lab/sceptre/discussions/106|very first discussion post>.
  - link: https://github.com/Katsevich-Lab/sceptre/discussions/106

**[2024-04-10 00:21] Timothy Barry**
<https://github.com/Katsevich-Lab/sceptre/issues/108|This> is rather concerning to me. Has anyone encountered this in any form?
  - link: https://github.com/Katsevich-Lab/sceptre/issues/108

**[2024-04-10 22:46] Louis Deutsch** _(thread reply)_
`launch_script.sh` is this one, with the only change being one of the paths
<https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/launch_script.sh>
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/launch_script.sh

**[2024-04-11 13:12] Louis Deutsch**
here's a qqplot for each NT gRNA separately

**[2024-04-12 16:32] Timothy Barry** _(thread reply)_
I do think it would be interested to remove those NT gRNAs (by changing their label from "non-targeting" to something else) and rerun your analysis.

**[2024-04-12 18:25] Louis Deutsch** _(thread reply)_
what's the best way to do this? Will I need to use `import_data_from_cell_ranger` with the modified `grna_target_data_frame`? There's no direct way to initialize a sceptre_object with ondisc response and grna matrices + in-memory `grna_target_data_frame` right?

**[2024-04-13 13:15] Timothy Barry** _(thread reply)_
My recommendation would be to directly update the gRNA target data frame (`sceptre_object@grna_target_data_frame`) after creating the `sceptre_object`.

**[2024-04-13 13:16] Timothy Barry** _(thread reply)_
Also, note that the gRNA target data frame is always in memory, regardless of whether the data are stored on disk or not.

**[2024-04-13 13:30] Timothy Barry** _(thread reply)_
So in fact you could pass a modified `grna_target_data_frame` to `import_data_from_cellranger`, but this would required you to call `import_data_from_cellranger` again (which is kind of slow).

**[2024-04-16 12:58] Eugene Katsevich**
Hey <@U0239H5UC9E>: `sceptre::get_grna_assignments()` appears to return a matrix of class `ngRMatrix`. Is this intentional? I’m currently updating the `sceptreIGVF` package to accommodate our latest push to `sceptre` and this weird class is causing some trouble.

**[2024-04-16 17:05] Louis Deutsch** _(thread reply)_
I'm telling it to use `batch` via
```set_analysis_parameters(
    formula_object = formula(~ batch + log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                               log(response_n_nonzero) + log(response_n_umis))
  ) ```
and then saving the resulting sceptre_object. Is that right? This works when I run it on my laptop for a small number of calibration pairs

**[2024-04-18 16:55] Eugene Katsevich**
<@U0239H5UC9E>: In the process of writing a `sceptre` inference <https://github.com/IGVF-CRISPR/CRISPR-JAMBOREE/blob/sceptre_pipeline/single-cell/inference/sceptre/nextflow/main.nf|Nextflow process> for IGVF, I found that the syntax
```    script:
    """
    inference_sceptre.R \
      $mudata_fp \
      $params.side \
      $params.grna_integration_strategy \
      $params.resampling_approximation \
      $params.control_group \
      $params.resampling_mechanism \
      '${params.formula_object}'
    """```
resolves the issue with special characters in the formula object being incorrectly interpreted by `bash`. I think that using this syntax to allow users to specify formula objects as strings rather than R objects in our `sceptre` pipeline is cleaner. What do you think?
  - link: https://github.com/IGVF-CRISPR/CRISPR-JAMBOREE/blob/sceptre_pipeline/single-cell/inference/sceptre/nextflow/main.nf

**[2024-04-18 18:24] Eugene Katsevich**
&gt; So does this allow for spaces in the formula object?
Yes. I tested it with something like `"~ log(response_n_nonzero) + log(response_n_umis) + prep_batch"`.
&gt; And does the formula object need to be specified as the final parameter passed to the script?
I don’t think so. Just make sure to put `'${params.formula_object}'` instead of `$params.formula_object`, wherever this parameter appears.

**[2024-04-20 11:58] Eugene Katsevich**
<@U0239H5UC9E> I just submitted a <https://github.com/Katsevich-Lab/sceptre/pull/113|pull request> including the vanilla score test functionality as well as the CRT functionality for ondisc-backed `sceptre` objects. Please review!
  - link: https://github.com/Katsevich-Lab/sceptre/pull/113

**[2024-04-25 01:43] Louis Deutsch**
I ran this code and just got this warning which I haven't seen before.
```sceptre_object &lt;- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = path_to_sceptre_object,
  response_odm_file_fp = path_to_gene,
  grna_odm_file_fp = path_to_grna
)  |&gt;
  set_analysis_parameters(
    formula_object = formula(~ log(grna_n_nonzero + 1) + log(grna_n_umis + 1) +
                               log(response_n_nonzero) + log(response_n_umis)),
    resampling_mechanism = "crt", resampling_approximation = "no_approximation"
  ) |&gt;
  assign_grnas() |&gt;
  run_qc()

Warning messages:
1: In max(sceptre_object@n_ok_discovery_pairs, sceptre_object@n_ok_positive_control_pairs) :
  no non-missing arguments to max; returning -Inf
2: In run_qc_pt_2(run_qc_pt_1(sceptre_object, n_nonzero_trt_thresh,  :
  NAs introduced by coercion to integer range```
This seems concerning? This is Gene's fork of `sceptre` on the `score-test` branch so I can use CRT with ondisc-backed sceptre. For my data I ran `import_data_from_cellranger()` with this fork, and then I'm using the `gene.odm`, `grna.odm`, and `sceptre_object.rds` that resulted from that.

**[2024-04-28 11:45] Timothy Barry** _(thread reply)_
Have you played around with the transpose operation? I was looking at the `DelayedArray` manual <https://www.bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf|here>, and I did not find anything about transpose... Also, do you know how recently these features were implemented? Perhaps they came online after we implemented the core `ondisc` functionality.

All that said, yes, I do think we should do comparisons of `ondisc` to alternatives, starting with `HDF5Array`/`DelayedArray` and `LoomR`. Hopefully we will see some sort of advantage in one or more respects.

More broadly, we I think we should try to come up with a concrete plan re benchmarking in the context of this paper. Which packages will we benchmark against? What metrics will we be tracking? Which datasets will we use? Ideally, we would keep benchmarking to a minimum (given our limited bandwidth), but yea, some benchmarking probably is unavoidable.
  - link: https://www.bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf

**[2024-04-28 12:17] Eugene Katsevich** _(thread reply)_
&gt; Have you played around with the transpose operation?
Not yet.
&gt; I was looking at the `DelayedArray` manual <https://www.bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf|here>, and I did not find anything about transpose... Also, do you know how recently these features were implemented? Perhaps they came online after we implemented the core `ondisc` functionality.
I found a reference to transposition on slide 17 <https://www.bioconductor.org/packages/devel/bioc/vignettes/DelayedArray/inst/doc/01-Working_with_large_arrays.pdf|here>, dated 2017 (see screenshot below).
&gt; All that said, yes, I do think we should do comparisons of `ondisc` to alternatives, starting with `HDF5Array`/`DelayedArray` and `LoomR`. Hopefully we will see some sort of advantage in one or more respects.
&gt; 
&gt; More broadly, we I think we should try to come up with a concrete plan re benchmarking in the context of this paper. Which packages will we benchmark against? What metrics will we be tracking? Which datasets will we use? Ideally, we would keep benchmarking to a minimum (given our limited bandwidth), but yea, some benchmarking probably is unavoidable.
I agree.
  - link: https://www.bioconductor.org/packages/devel/bioc/manuals/DelayedArray/man/DelayedArray.pdf
  - link: https://www.bioconductor.org/packages/devel/bioc/vignettes/DelayedArray/inst/doc/01-Working_with_large_arrays.pdf

**[2024-04-29 13:53] Louis Deutsch**
hi all, for the smaller Replogle dataset that I'm going to analyze, where exactly is it?

I found 12 files described as follows:

&gt; This dataset includes data from three Perturb-seq experiments described in Replogle et al. 2022: 
&gt; 1. K562 genome-scale perturb-seq sampled at day 8 post-transduction (K562_gwps)
&gt; 2. K562 essential-scale perturb-seq sampled at day 6 post-transduction (K562_essential)
&gt; 3. RPE1 essential-scale perturb-seq sampled at day 7 post-transduction (rpe1)
&gt; 
&gt; For each dataset, there are four processed Perturb-seq files in AnnData format.
&gt; 1. Raw, single-cell expression data for genes expressed at &gt;0.01 UMI per cell (named $pop_raw_singlecell_01.h5ad)
&gt; 2. Raw, pseudo-bulk expression data for genes expressed at &gt;0.01 UMI per cell (named $pop_raw_bulk_01.h5ad)
&gt; 3. gemgroup Z-normalized single-cell expression data for genes expressed at &gt;0.01 UMI per cell (named $pop_normalized_singlecell_01.h5ad)
&gt; 4. gemgroup Z-normalized pseudo-bulk expression data for genes expressed at &gt;0.01 UMI per cell (named $pop_normalized_bulk_01.h5ad)
is it one of these?

**[2024-04-29 14:07] Eugene Katsevich**
The dataset is K562_essential, and the file you want is probably the raw, single-cell expression data for genes expressed at &gt;0.01 UMI per cell (named $pop_raw_singlecell_01.h5ad).

**[2024-04-29 14:08] Eugene Katsevich**
You should probably add your code to <https://github.com/Katsevich-Lab/import-replogle-2022/tree/main|Katsevich-Lab/import-replogle-2022>; I’ve just given you push access. You can use Tim’s code there as a template.
  - link: https://github.com/Katsevich-Lab/import-replogle-2022/tree/main

**[2024-04-29 14:15] Eugene Katsevich**
Since you’ll be looking at Replogle data import, I’d like you to be on the lookout for *another potential source of miscalibration: incorrect gRNA labeling during import.* Please pay special attention to how we label a gRNA as non-targeting, and triple check the correctness of these labels.

**[2024-05-01 11:38] Eugene Katsevich**
They apparently model variable gRNA efficacy as well. 

**[2024-05-01 11:47] Timothy Barry**
Finally, I appreciate how transparent their power simulation results are.

**[2024-05-01 11:53] Timothy Barry**
I don't see anything about gRNA-to-cell assignment or QC.

**[2024-05-01 12:03] Eugene Katsevich**
The presentation is tomorrow. I can ask about gRNA assignment and QC.

**[2024-05-01 12:11] Timothy Barry** _(thread reply)_
We should keep in mind that our running time on the Gasperini data includes assigning gRNAs to cells using the mixture model.

**[2024-05-01 12:13] Timothy Barry** _(thread reply)_
Do you mean for gRNA assignment? We used the mixture model.

**[2024-05-01 17:19] Timothy Barry**
I have been working with the dataset of one of our clients, Veronica. Veronica is carrying out a cis analysis of a high-MOI, CRISPRi enhancer-targeting screen. Initially, there was miscalibration: 9 negative control pairs were falsely rejected (see plot).

**[2024-05-01 17:31] Timothy Barry**
So, in conclusion, correcting for the confounding effect due to cell cycle seems like a necessary step, at least on some datasets. I think this implies two action items.
1. `sceptre` should provide an option for computing the top k PCs of the response matrix. I added an issue to this effect <https://github.com/Katsevich-Lab/sceptre/issues/125|here>.
2. We should investigate this issue on the Replogle rd7 data.
  - link: https://github.com/Katsevich-Lab/sceptre/issues/125

**[2024-05-01 17:33] Eugene Katsevich**
&gt; 1. `sceptre` should provide an option for computing the top k PCs of the response matrix. I added an issue to to this effect <https://github.com/Katsevich-Lab/sceptre/issues/125|here>.
Would `sceptre` or `ondisc` handle this? It seems like it fits conceptually within the “computation of cellwise covariates” step, which `ondisc` handles.
  - link: https://github.com/Katsevich-Lab/sceptre/issues/125

**[2024-05-01 18:25] Eugene Katsevich** _(thread reply)_
I’m still holding out hope that there was some sort of gRNA mislabeling at the import step, but unobserved confounders definitely are a potential culprit here.

**[2024-05-03 17:11] Eugene Katsevich**
Folks, I wanted to share a quick thought on the confounding issue. The model I always had in mind is that _biological_ variability impacts gene expression but not gRNA presence, because gRNAs assort themselves randomly into cells. _Technical_ variability, on the other hand, can impact our measurements of both gene expression and gRNA presence. This is why, in the past, we have focused exclusively on adjusting for technical sources of variation (e.g., batch, sequencing depth) but not biological sources of variation (e.g, cell cycle).

I wonder if we need to reconsider the statement that “gRNAs assort themselves randomly into cells.”  I could imagine that cells in different biological states (e.g. cell cycle states) may have different probabilities of gRNA infection. This may not cause issues if cells in different biological states have equal chances of receiving any gRNA, conditional on them receiving a gRNA. We would really have an unmeasured confounding problem if cells in different biological states have preferences for different gRNAs. This seems like a less plausible scenario, since all gRNAs come in the same viral packaging.

In summary, it remains unclear to me why biological variability, like cell cycle, would have a confounding effect on single-cell CRISPR screens.

**[2024-05-03 19:01] Eugene Katsevich** _(thread reply)_
Let’s say we have small cells and big cells, and only big cells get gRNAs due to the phenomenon you hypothesize. Then we compare gene expression in big cells with an NT gRNA and big cells with a targeting gRNA. As long as _which_ gRNA a big cell receives is random, there shouldn’t be a confounding issue. We might worry that our conclusion holds for big cells rather than for all cells, but I wouldn’t expect spurious discoveries in this situation. 

**[2024-05-03 19:02] Eugene Katsevich** _(thread reply)_
To clarify, I am thinking about a low MOI setting (like Replogle), where the control group is the NT cells rather than the complement set. 

**[2024-05-04 19:06] Eugene Katsevich**
Just as an FYI, my undergraduate student has completed her <https://github.com/Katsevich-Lab/sceptre-vs-score/blob/main/final-report/final-writeup.pdf|final report> comparing `sceptre` to the vanilla score test on the Gasperini data. The main findings are as follows:
• *Computational:* After precomputation, the score test takes about 0.07 seconds/pair, whereas the CRT-based `sceptre` takes about 2.8 seconds/pair. 
• *Statistical:* The left-sided score test has reasonable calibration and power, while the right-sided score test has nontrivial inflation, especially for moderate effective sample sizes.
  - link: https://github.com/Katsevich-Lab/sceptre-vs-score/blob/main/final-report/final-writeup.pdf

**[2024-05-04 19:17] Eugene Katsevich** _(thread reply)_
This issue slipped through the cracks. <@U0524GR916C>: Did you ever figure this out? If not, do you know if this warning is unique to the `score-test` branch? Which `path_to_sceptre_object`, `path_to_gene`, and `path_to_grna` are you using?

**[2024-05-05 19:28] Timothy Barry**
I appreciate this writeup from Winnie. Perhaps you can tell her that I think she did a good job. :slightly_smiling_face: A few comments:

1. The reported precomputation times seem a bit strange to me. Do I understand correctly that the per-gene precomputation time for `sceptre` was 11.3 seconds? I just used `sceptre` to analyze a single negative control pair on the full, 200,000+ cell Gasperini dataset. The total running time --- which included (1) constructing the negative control pair, (2) running the gene precomputation, (3) running the gRNA precomputation, and (4) running the association test was 0.52 seconds. So for some reason it seems Winnie's code is running about an order of magnitude slower than mine?
2. I noticed that Winnie analyzed 100 negative control pairs. I think this is reasonable, but we should keep in mind that we do not benefit much from sharing of information across gRNAs in this setting. Thus, the CRT runtime (relative to the score test runtime) is a bit conservative.
3. The statistical results are intriguing. I think it would be helpful to plot the negative control results for sceptre as well, i.e. to repeat Figure 1 but using sceptre instead of the score test. If I recall correctly, sceptre makes ~2-3 false discoveries in the right tail out to ~100K pairs.
4. I like the stratification by effective sample size. We see pretty strong evidence that sparsity is culpable (at least in part) for miscalibration of the score test.
5. I took a look at the PR, and briefly, it looks like the score test is implemented as follows: (1) estimate the per-gene size parameter using a Poisson GLM; (2) fit the NB GLM via `MASS`; (3) compute the score test using our C++ score test function. Is this roughly accurate? 
6. Back in the day, we assessed the left-tailed NB GLM Wald test on the Gasperini data out to ~600,000 pairs, and we observed moderate miscalibratation. Do we believe that, had we assessed the left-tailed NB GLM score test out to ~600,000 pairs, we also would see moderate miscalibration?

**[2024-05-05 19:33] Timothy Barry** _(thread reply)_
Oh I see. Veronica’s data are high-MOI, so that's the setting I mostly had in mind. I agree that it's harder to cook up a story to explain why including PCs as covariates would help improve calibration in the low-MOI setting.

**[2024-05-11 11:56] Louis Deutsch**
hi all, I've figured out that warning I was getting. This is all on the `score-test` branch of Gene's fork of `sceptre`.

The warning can be recreated by this code:
```read_ondisc_backed_sceptre_object(
  sceptre_object_fp = path_to_sceptre_object,
  response_odm_file_fp = path_to_gene,
  grna_odm_file_fp = path_to_grna
)  |&gt;
  set_analysis_parameters(
    resampling_approximation = "no_approximation"
  ) |&gt;
  assign_grnas() |&gt;
  run_qc()```
The particular warning is
```Warning messages:
1: In max(sceptre_object@n_ok_discovery_pairs, sceptre_object@n_ok_positive_control_pairs) :
  no non-missing arguments to max; returning -Inf
2: In run_qc_pt_2(run_qc_pt_1(sceptre_object, n_nonzero_trt_thresh,  :
  NAs introduced by coercion to integer range```
The problem is this code in `run_qc_pt2()`:
```if (sceptre_object@resampling_approximation == "no_approximation") {
  mult_fact &lt;- if (sceptre_object@side_code == 0L) 10 else 5
  sceptre_object@B3 &lt;- ceiling(mult_fact * max(
    sceptre_object@n_ok_discovery_pairs,
    sceptre_object@n_ok_positive_control_pairs
  ) / sceptre_object@multiple_testing_alpha) |&gt;
    as.integer()
}```
since I'm doing this without positive control pairs or discovery pairs, `sceptre_object@n_ok_discovery_pairs` and `sceptre_object@n_ok_positive_control_pairs` are both still at their default of `integer(0)` . The `max` then returns `-Inf` with the warning that there are no non-missing arguments, and then the cast to integer gives the second warning.

**[2024-05-12 20:50] Louis Deutsch**
On the main branch of your fork the corresponding bit of `run_qc()` is
``` if (!sceptre_object@fit_parametric_curve) {
    side &lt;- sceptre_object@side_code
    mult_fact &lt;- if (side == 0L) 10 else 5
    sceptre_object@B3 &lt;- ceiling(mult_fact * sceptre_object@n_ok_discovery_pairs/sceptre_object@multiple_testing_alpha) |&gt;
      as.integer()
  }```
On this branch `set_analysis_parameters()` requires `discovery_pairs`, but if I pass a data.frame with 0 rows I can still have `sceptre_object@n_ok_discovery_pairs` be `integer(0)`. Arithmetic with `integer(0)` becomes `numeric(0)` but then with the final cast to int we have `@B3` as `integer(0)`. This doesn't raise any warnings but subsequent parts throw errors.

I'm now running this with the low moi example data and I get an error on `run_calibration_check()`:
```Error: Expecting a single value: [extent=0].```
which goes away if I set `@B3` to something like 10L.

On `katsevich-lab/sceptre` we have
``` if (sceptre_object@resampling_approximation == "no_approximation") {
    mult_fact &lt;- if (sceptre_object@side_code == 0L) 10 else 5
    sceptre_object@B3 &lt;- ceiling(mult_fact * max(
      sceptre_object@n_ok_discovery_pairs,
      sceptre_object@n_ok_positive_control_pairs
    ) / sceptre_object@multiple_testing_alpha) |&gt;
      as.integer()
  }```
so we'll get the previous issue with an `NA` value for `@B3`

I think this might be an easy fix if we just add something to ensure that `B3` never ends up with an empty or missing value?

**[2024-05-15 16:15] Eugene Katsevich**
Hi folks, that’s my bad for not clarifying when we’d be resuming our meetings. I propose that Tim and I meet on Friday, as he proposes, and that Joseph sends us any updates via Slack. We can then resume our regularly scheduled meetings next week.

**[2024-06-08 21:28] Louis Deutsch**
hi all, I'm trying to add the features from `crispat` to our code. I have a data.frame per batch that looks like
```              cell    S_score    G2M_score Batch total_UMI_counts percent_mito
1 AAACCCAAGCTATCTG -0.4246705 -0.364349057     1            11120     8.732014
2 AAACCCAAGGCATGCA -0.1870530  0.122235984     1            22006     9.547396
3 AAACCCAAGTATCCTG  0.1551133 -0.228593455     1            16424     8.085728
4 AAACCCACAATGAAAC -0.4142903 -0.384319798     1            10776    10.282108
5 AAACCCACAGCTACAT  0.1788355  0.009699748     1             9072    10.912698
6 AAACCCAGTATGCGTT -0.2268778  0.484569411     1            14253     8.180734```
Can I match this to the contents of our sceptre object when using ondisc-backing? Or will I need to use HPC3 so I can just do this all in-memory and not use ondisc? I know with `get_response_matrix()` I can subset by gene but I'm not aware of a way to get cell-level information.

Also another QC step they do is, before they create their sceptre object, they filter these covariate data.frames to `total_UMI_counts &gt; 3000 &amp; percent_mito &lt; 20`.

**[2024-07-03 17:50] Louis Deutsch**
hi all, I'm trying to compare our grna_target_data_frame to the Crispat team's but I'm unsure about this part:

They make the grna_target_data_frame in python, and they do so by loading a file they call `K562_essential_library_annotation.csv` in a step described as
```create annotation df from the Supplementary Information file of Replogle et al. and deal with the duplicated gRNAs```
I can't find a file that exactly matches this. There is a file here that has "annotated" in the name, but the column names do not agree with the column names that their file is expected to have <https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_-_commonly_requested_supplemental_files/21632564|https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenoty[…]_et_al_2022_-_commonly_requested_supplemental_files/21632564>

Is there a quick answer to this? Otherwise I'll just wait and see what they say in the morning.
  - link: https://plus.figshare.com/articles/dataset/_Mapping_information-rich_genotype-phenotype_landscapes_with_genome-scale_Perturb-seq_Replogle_et_al_2022_-_commonly_requested_supplemental_files/21632564

**[2024-07-04 09:57] Eugene Katsevich**
List of questions:
• Please walk us through the beginning of the flow of the data through your code (Python and/or R).
• Your gRNA matrix has more rows than ours (yours has 2283) and ours has (2184). Do you know why? Perhaps relatedly, some of your rownames look like “2546_EIF3CL_ENST00000380876.4,ENST00000398943.3,ENST00000398944.3_ENSG00000205609”. Why is this?
• Just to make sure, is your `K562_essential_library_annotation.csv` file the same as the `mmcl.xlsx` from Cell?
• What is the “for duplicated gRNAs (same target sequence) we only keep one (the first) of multiple target gene names” code doing?
• Is there a way that you could provide us with the sceptre objects from successful applications on both datasets, so that we can compare directly.
• Supplementary Figure 3 refers to which dataset?

**[2024-09-14 12:13] Eugene Katsevich**
Hey folks, I’m serving on the thesis committee for a bioengineering PhD student at Penn who is running a CRISPRi screen of enhancers implicated in Alzheimer’s disease. She successfully applied `sceptre` to her data and is happy with the results. Hopefully the preprint will be coming out this year. The one question she had was about how `construct_cis_pairs()` works when `grna_integration_strategy = "union"`, since in that case there are multiple gRNA positions to deal with. Upon inspecting our code, I found that we define distance based on the midpoint between the minimum of the `start` positions and the maximum of the `end` positions. I think this would be worth clarifying in our function documentation and/or e-book. Could one of you please add this clarification?

**[2024-09-14 12:33] Timothy Barry**
I added the following sentence to the "<https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-grna_integration_strategy|gRNA integration strategy>" section of the book: "When `grna_integration_strategy` is set to `"union"`, `construct_cis_pairs()`defines the “position” of a given target as the midpoint between the minimum of the `start` positions and maximum of the `end` positions of the gRNAs within the target." Please let me know if you think this is reasonable.
  - link: https://timothy-barry.github.io/sceptre-book/set-analysis-parameters.html#sec-grna_integration_strategy

**[2024-09-17 17:50] Eugene Katsevich**
Hey <@U0239H5UC9E>, could you please comment on this discussion when you have a chance? Is there something we can change about our documentation (or function itself) to prevent or catch such user mistakes in the future? (I confess to not fully understanding what the mistake was, beyond the distinction between `grna_id` and `grna_target`.)

**[2024-09-28 21:19] Eugene Katsevich**
We’ve gotten three different, but related feature requests:
• In high MOI, define control cells as those receiving only non-targeting perturbations. 
• For high-MOI singleton gRNA analyses, define control cells as the complement set except cells receiving other gRNAs with the same target.
• For low- or high-MOI, define control cells as the complement set except cells receiving other gRNAs targeting the same chromosome.
I was thinking we (or someone…) could implement these options in a unified framework. We already implicitly have the concept of gRNA groups, which are groups of gRNAs we test together (specified by gRNA target, or they could be singletons). Let me call these “treatment gRNA groups.” We could designate another column of our gRNA data frame as specifying “control exclusion gRNA groups.” So each gRNA belongs to one “treatment group” and one “control exclusion group.” The control exclusion group of a gRNA is the set of gRNAs that should not belong to the control group when that gRNA’s treatment group is being tested. So for a given treatment gRNA group, we define the set of control cells as those containing no gRNAs in the control exclusion groups of each gRNA in the treatment group. This accommodates for all of the cases we currently implement and all of the cases we wish to implement. In particular:

• NT control group, accommodating high-MOI: The control exclusion group of a gRNA is its targeting indicator.
• Usual complement control group (singleton gRNA integration strategy): The control exclusion group of a gRNA is its gRNA ID.
• Distinct target complement control group (singleton gRNA integration strategy): The control exclusion group of a gRNA is its gRNA target.
• Distinct target complement control group (union gRNA integration strategy): The control exclusion group of a gRNA is again its gRNA target.
• Chromosome-based control group: The control exclusion group of a gRNA is its chromosome.
We could slightly extend this idea by designating some gRNAs for exclusion from all control groups, like gene-targeting gRNAs in enhancer-targeting screens, or gRNAs targeting near TFs.

This unified framework could accommodate a large variety of control group constructions, including many we cannot even yet envision. It would also provide a route to clean implementation without handling each case separately. Of course, there are some considerations like how this would impact our pairwise QC. One observation is that the main pairwise QC metric we currently lean on is `n_nonzero_trt`, which remains unchanged regardless of the control group. If we did away with `n_nonzero_ctrl`, then we need not compute this quantity for each definition of the control group. (Eventually, we’ll probably want to move away from `n_nonzero_trt`, too…) The other consideration is how this change would impact our actual association tests. In our current implementation, I believe we are recomputing the model fits for each pair when control group is the NT cells. If we switched to always doing this, then I think we wouldn’t have to change much (at the cost of inefficiency when control group is the complement set). Another possibility is to always do model fits on the entire group of cells, regardless of what control group we use.

I’m curious what you guys think of this approach, separately from the question of who has the bandwidth to actually implement it.

**[2024-09-30 11:31] Timothy Barry**
I think this is an interesting suggestion. Here are a few preliminary thoughts and questions.
• I think it would be useful to have a visual to explain this idea to users (or even to the person who will be implementing this idea in the code).
• For "NT control group, accommodating high-MOI," could you explain what you mean by, "The control exclusion group of a gRNA is its targeting indicator"?
• I would be concerned about ditching `n_nonzero_cntrl` but not `n_nonzero_trt` due to selection bias issues.
• The code is indeed complex, but barring a massive rewrite, I am not sure we will be able to make it simpler. It would be very cool and useful if we could use an LLM for this purpose.
• There are a couple relevant optimizations in the code related to control group. First, for the complement set control group, we share a fitted regression model across all targets paired to a given gene. Second, we are careful about generating the synthetic treatment indices for use in the permutation test (or CRT). For instance, when using the NT cells as the control group, we use IWOR sampling to generate the synthetic treatment indices.
With respect to the proposal itself, I think I mostly understand and agree, but I would modify the proposal in the following way. I would add another option `"custom"` for the `control_group` argument. `"custom"` would enable users to manually specify the control group via the `control_exclusion` column of the gRNA target data frame. This would enable us to keep our code backward compatible and our optimizations in place. Additionally, it would keep the API simple for users who wish to use a simple control group (i.e., most users). We would then be able to separately handle the `"custom"` option, likely throwing out most or all optimizations due to the lack of structure.

**[2024-10-01 18:58] Eugene Katsevich**
&gt; • I think it would be useful to have a visual to explain this idea to users (or even to the person who will be implementing this idea in the code).
Agreed.
&gt; • For “NT control group, accommodating high-MOI,” could you explain what you mean by, “The control exclusion group of a gRNA is its targeting indicator”?
The control group in my framework is defined as all cells containing only gRNAs not sharing the “control exclusion gRNA group” with the gRNA (or gRNA group) being tested. If all targeting gRNAs are put in the “targeting” control exclusion group and all non-targeting gRNAs are put in the “non-targeting” control exclusion group, then the control group for each targeting gRNA (group) would end up being all of the non-targeting gRNAs, matching the current behavior.
&gt; • I would be concerned about ditching `n_nonzero_cntrl` but not `n_nonzero_trt` due to selection bias issues.
I already forgot some of the finer insights I had before about the extent of the selection bias introduced by our pairwise QC and how this selection bias depended on filtering on both or just one of our metrics. But I agree with you that ditching `n_nonzero_cntrl` may have unintended statistical consequences.
&gt; • The code is indeed complex, but barring a massive rewrite, I am not sure we will be able to make it simpler. It would be very cool and useful if we could use an LLM for this purpose.
I agree that a massive rewrite may be the only way to make the code more modular and easier to reason about. That’s a tough position for us to be in, because otherwise even apparently small additional features are quite hard to add because of how they interact with the rest of our codebase. But a massive rewrite is probably not in the cards for any human to do in the near future. I’m going to keep hoping that AI gets to the point where you can upload a huge codebase to an LLM and ask it to rewrite the code from scratch in the way that has X, Y, and Z desirable properties. That would be a remarkable feat for an AI to do, and I don’t think it’s there yet.
&gt; • There are a couple relevant optimizations in the code related to control group. First, for the complement set control group, we share a fitted regression model across all targets paired to a given gene. Second, we are careful about generating the synthetic treatment indices for use in the permutation test (or CRT). For instance, when using the NT cells as the control group, we use IWOR sampling to generate the synthetic treatment indices.
&gt; With respect to the proposal itself, I think I mostly understand and agree, but I would modify the proposal in the following way. I would add another option `"custom"` for the `control_group` argument. `"custom"` would enable users to manually specify the control group via the `control_exclusion` column of the gRNA target data frame. This would enable us to keep our code backward compatible and our optimizations in place. Additionally, it would keep the API simple for users who wish to use a simple control group (i.e., most users). We would then be able to separately handle the `"custom"` option, likely throwing out most or all optimizations due to the lack of structure.
Yes, the key challenge is to respect these optimizations while making the code more flexible. Having a more general but less optimized option alongside less general but more optimized options is probably the best short- to medium-term solution, although it would make our codebase more fragmented.

**[2024-10-01 19:01] Eugene Katsevich**
<@U0239H5UC9E>: I’ll let you handle <https://github.com/Katsevich-Lab/sceptre/discussions/148|this> discussion question about the pipeline.
  - link: https://github.com/Katsevich-Lab/sceptre/discussions/148

**[2024-11-02 18:48] Eugene Katsevich**
Hey <@U0239H5UC9E>, could you please look into <https://github.com/Katsevich-Lab/sceptre/issues/146#issuecomment-2451940537|this> user inquiry? Thank you!
  - link: https://github.com/Katsevich-Lab/sceptre/issues/146#issuecomment-2451940537

**[2024-11-05 17:20] Eugene Katsevich** _(thread reply)_
I took care of this. It appears that the problem here was that adding the gRNA covariates when fitting the mixture gRNA assignment model caused perfect multicollinearity. I added fixing this to our <https://github.com/orgs/Katsevich-Lab/projects/7/views/1|backlog>.
  - link: https://github.com/orgs/Katsevich-Lab/projects/7/views/1

**[2024-11-05 17:26] Eugene Katsevich**
<@U0239H5UC9E>: For <https://github.com/Katsevich-Lab/sceptre/issues/157|this issue>, is the problem that running the Nextflow pipeline destroys some aspects of the `sceptre_object`, so running `assign_grnas()` on that messed up `sceptre_object` causes the error?
  - link: https://github.com/Katsevich-Lab/sceptre/issues/157

**[2024-11-05 17:27] Eugene Katsevich** _(thread reply)_
So would the solution be to run `assign_grnas()` on a new `sceptre_object` that is not the output of the Nextflow pipeline? If so, I’d be happy to respond to Nan myself.

**[2024-11-07 22:34] Timothy Barry** _(thread reply)_
I just responded to Nan to clarify this point. But many users (i.e., at least three) would like to run `plot_response_grna_target_pair()` after running the NF pipeline, and we do not currently support this apparently. Indeed, I think this is because we are removing important information from the `sceptre_object` at intermediate steps of the pipeline. We probably should retain this information in the `sceptre_object` so that users simply can load the `sceptre_object` into R and call `plot_response_grna_target_pair()`.

**[2024-11-08 21:27] Eugene Katsevich** _(thread reply)_
I just <https://github.com/Katsevich-Lab/sceptre/discussions/156#discussioncomment-11170027|told this> to `jpulecio`:
&gt; We are looking into this issue, which other users have also encountered. Our current theory for the origin of this error is that the `sceptre_object` returned by the Nextflow pipeline is not meant to be operated on in the same way as regular `sceptre_objects`, because certain fields are removed by the pipeline for memory efficiency purposes. You might have to rerun `assign_grnas()` and `run_qc()` on a newly created `sceptre_object` instead of the one output by the pipeline. Could you please try this?
  - link: https://github.com/Katsevich-Lab/sceptre/discussions/156#discussioncomment-11170027

**[2024-11-09 12:41] Timothy Barry** _(thread reply)_
I think we should just add support for calling `plot_response_grna_target_pair()` after running the NF pipeline. This seems to be our most popular feature request.

**[2024-12-07 21:56] Eugene Katsevich**
Hey folks, I got a ping from Andreas about the question of whether it is possible to do a low-MOI analysis without throwing out cells with more than one gRNA. I’ve decided to try to implement a preliminary solution because I feel bad about keeping him waiting that long. I have some thoughts about how to go about this, which I wanted to run by you.

*Conceptual.* One way to generalize `sceptre` to handle more flexible assignments of cells to treatment and control groups is through the arguments `treatment_group` (which does not yet exist) and `control_group` in `set_analysis_parameters()`.  We can have the following options for each:
• For `treatment_group`, we can have the options `exclusive` (define a cell as being in the treatment group if it contains a gRNA in the targeting group of interest and no gRNAs in any other targeting group) and `inclusive` (define a cell as being in the treatment group if it contains a gRNA in the targeting group of interest and any other gRNAs). We would use `exclusive` as the default for low-MOI and `inclusive` as the default for high-MOI. 
• For `control_group`, we can have the options `complement` (as before), `nt_cells` (define a cell as being in the control group if it contains only NT gRNAs), and later we can have fancier ones like modified versions of the `complement` strategy that excludes gRNAs with the same target. 
In this philosophy, we don’t incorporate gRNA information into cell-wise QC. We just define treatment and control cells appropriately, which may or may not have the side effect of excluding certain cells. If Andreas wants to retain the maximum number of cells in his analysis, he might choose `treatment_group = inclusive` and `control_group = nt_cells`. This allows cells with multiple gRNAs into both treatment and control groups.

*Implementation.* At present, the above strategy seems to be relatively easy to implement. Here is what needs to be done:
• `sceptre_object`: Add a new slot `treatment_group`, and perhaps remove `cells_w_zero_or_twoplus_grnas`.
• `set_analysis_parameters()`: Add the argument `treatment_group` and set appropriate defaults.
• `process_initial_assignment_list()`: Implement the logic above to define `grna_group_idxs` and `indiv_nt_grna_idxs`. 
• `determine_cells_to_retain()`: Remove the filter based on `cells_w_zero_or_twoplus_grnas`.
*Outlook.* I could implement these changes on a new branch and then point Andreas to this branch. This branch would be experimental. I admit that more thought would be required to make sure that our calibration check is doing the right thing for some of the above options. However, if Andreas wants to check calibration properly, I can suggest that he manually relabel one NT gRNA at a time to “undercover” and run discovery analyses for those undercover gRNAs. In the longer term, hopefully we can merge some of this new functionality into the main branch after more exhaustive validation, documentation, testing, etc.

What do you folks think about such a course of action?

**[2024-12-30 15:58] Eugene Katsevich**
Hey <@U0239H5UC9E>, I’ve implemented the above proposal in <https://github.com/Katsevich-Lab/sceptre/pull/163|this pull request>. It was harder than I thought but I believe I’ve overcome the main challenges that presented themselves. Please let me know what you think. Happy to discuss if you have any questions. Hopefully we can merge this into `dev` soon, and then point Andreas to this version of the package, which he’s been awaiting since August.
  - link: https://github.com/Katsevich-Lab/sceptre/pull/163

**[2024-12-30 16:17] Eugene Katsevich** _(thread reply)_
<@U0239H5UC9E>: It appears you deleted `.github` <https://github.com/Katsevich-Lab/sceptre/commit/cc899e7dc55b68c0875eac8d3d752b226817d5a5|in April>, which suspended our `R-CMD-CHECK` CI? Do I understand this correctly?
  - link: https://github.com/Katsevich-Lab/sceptre/commit/cc899e7dc55b68c0875eac8d3d752b226817d5a5

**[2025-01-02 12:13] Timothy Barry**
Sorry, Gene. for the delay; K99 has been all-consuming, but I should have a chance to look into these issues this afternoon.

**[2025-01-02 21:00] Timothy Barry**
Hey Gene, I am reviewing your pull request and have some basic questions. First of all, very nice job with this; the solution you found is quite elegant, and I think our users will appreciate this added flexibility.

As a starting point I want to make sure I understand what is going on with the API. Consider the example below. I am assuming `nt_cells` is the control group.
1. If `treatment_group` is `inclusive` and `remove_cells_w_zero_or_twoplus_grnas` is `FALSE`, then we retain cells 1-3.
2. If `treatment_group` is `inclusive` and `remove_cells_w_zero_or_twoplus_grnas` is `TRUE`, then we retain only cell 3.
3. If `treatment_group` is `exclusive` and `remove_cells_w_zero_or_twoplus_grnas` is `FALSE`, then we retain cells 1 and 3.
4. If `treatment_group` is `exclusive` and `remove_cells_w_zero_or_twoplus_grnas` is `TRUE`, then we retain only cell 3.
Is this accurate?

**[2025-01-03 16:09] Eugene Katsevich**
Action items for Gene:
• Ask Andreas if the “inclusive” option is really necessary, or if the “exclusive” option is sufficient for his purposes. Copy Tim.
• Make `treatment_group` and `remove_cells_w_zero_or_twoplus_grnas` into just one argument, with either two or three options.
• Check that `n_nonzero_trt` and `n_nonzero_cntrl` are still being computed correctly in these fancier cases.
• Disable fancier options when using maximum assignment strategy.

**[2025-01-06 16:56] Eugene Katsevich**
Folks, our vignette and examples each take upwards of ten minutes to run. This makes running devtools::check() annoyingly long, both during development and during CI. In my opinion, one of the main reasons for this is the large number of cells in our example datasets (46,000 in high MOI and 21,000 in low MOI). The reason we need such a large number of cells is that the effect sizes for discovery pairs in real data are fairly low, so our example outputs would not be very interesting if the number of cells were substantially lower. However, if we were to create *synthetic example data*, then we could simulate much larger effect sizes and therefore get away with much smaller numbers of cells (e.g. one or a few thousand), which would greatly accelerate our examples and vignette. What do you think of this idea?

**[2025-01-06 18:30] Timothy Barry**
It doesn't seem that there is any way to avoid this... OK, perhaps we can go head and use simulated data for the examples and vignette.

**[2025-01-06 20:08] Eugene Katsevich**
Should I put those simulated data into `sceptre` or `sceptredata`?

**[2025-01-17 15:41] Eugene Katsevich**
Hey <@U0239H5UC9E>, I have now submitted a <https://github.com/Katsevich-Lab/sceptre/pull/166|pull request> in which I replace the real example datasets with smaller, synthetic ones. I also restored our continuous integration on GitHub and established branch protection rules. Now, `R-CMD-CHECK` runs on pull requests to `main` and `dev`, and must be passed before merging. You’ll note that all the checks have already been passed for this pull request.

I know you must be busy with your K99, but this pull request is pretty short and sweet (especially compared to the others we’re contemplating), so hopefully it will not take much of your time to review and merge it.
  - link: https://github.com/Katsevich-Lab/sceptre/pull/166

**[2025-01-17 18:07] Timothy Barry** _(thread reply)_
Agreed. A lot of mature packages (e.g., <https://github.com/hail-is/hail?tab=readme-ov-file|hail>) have contributor guidelines.
  - link: https://github.com/hail-is/hail?tab=readme-ov-file

**[2025-01-28 17:19] Eugene Katsevich**
Hey <@U0239H5UC9E>, you might have seen Andreas’s email, in which he says that the biggest difference in his analyses would be made if we were able to avoid filtering out cells with multiple targeting gRNAs in the case `control_group = "nt_cells"`. While I have already implemented this functionality, it is problematic in the sense that we cannot easily implement a calibration check for this option. What I am leaning towards now is to do what you had suggested last time we discussed this: simply add one logical argument to `run_qc()` called `remove_cells_w_zero_or_twoplus_grnas`, which if set to `FALSE`, includes cells with multiple gRNAs if at most one of these gRNAs is targeting. This would unfortunately exclude the case Andreas cares about the most, but would lower the complexity of our software and preserve our calibration check functionality regardless of what options our users choose. What do you think?

**[2025-01-31 17:12] Eugene Katsevich**
Hey <@U0239H5UC9E>, I’ve just submitted a <https://github.com/Katsevich-Lab/sceptre/pull/175|pull request> that adds a single new argument `remove_cells_w_zero_or_twoplus_grnas`, as we have discussed. The one issue is that `compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc()` and `compute_n_ok_pairs_ondisc()` and perhaps `compute_n_trt_cells_matrix_ondisc()` from `ondisc` need to be updated as well. Is there a reason why these functions live in `ondisc` rather than in `sceptre`?
  - link: https://github.com/Katsevich-Lab/sceptre/pull/175

**[2025-01-31 17:38] Eugene Katsevich** _(thread reply)_
I’ve converted the pull request to a draft because I think that `remove_cells_w_zero_or_twoplus_grnas` is the wrong argument name. For low-MOI, we _always_ want to remove _at least some_ cells with zero or two-plus gRNAs. Those cells are the ones where we cannot reliably infer a perturbation identity (cells with no gRNAs or gRNAs with multiple targeting gRNAs). The question is just whether we want to remove _all_ cells with zero or two-plus gRNAs or just those where we cannot reliably infer a perturbation identity. Furthermore, `plot_run_qc()` should show the cells removed in either case; currently, it does not show any cells as being removed if `remove_cells_w_zero_or_twoplus_grnas = FALSE`.

**[2025-02-01 13:14] Timothy Barry** _(thread reply)_
There's no particular reason that `compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc()`, `compute_n_ok_pairs_ondisc()`, and `compute_n_trt_cells_matrix_ondisc()` are in `ondisc`, but these functions are different from their in-memory counterparts. For the time being, perhaps we can just throw an error if the user wants to use this specialized feature with an ondisc-backed sceptre object.

Have you had a chance to take a crack at updating the book or documentation?

I'm not so sure. I think `remove_cells_w_zero_or_twoplus_grnas` seems reasonable to me. I prefer this argument name to `keep_more_cells`, which seems kind of vague. Perhaps we can just make it clear in the documentation that `remove_cells_w_zero_or_twoplus_grnas` is relevant only if `control_group` is `"nt_cells"`? :man-shrugging:

**[2025-06-12 14:43] Eugene Katsevich**
Hi folks, I’m starting to poke around `sceptre3` again and have realized that I do not see any code for analyzing the Replogle data in <https://github.com/Katsevich-Lab/sceptre3-project|our repo>. I remember we put a lot of work into (unsuccessfully) analyzing those data. Do you know where all of that code went?
  - link: https://github.com/Katsevich-Lab/sceptre3-project

**[2025-06-12 16:22] Timothy Barry**
Yes. <https://github.com/Katsevich-Lab/sceptre3-project-v2> contains the launch scripts for the Replogle data, and <https://github.com/Katsevich-Lab/import-replogle-2022> contains the import scripts for the Replogle data.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2
  - link: https://github.com/Katsevich-Lab/import-replogle-2022

**[2025-06-23 18:46] Eugene Katsevich**
Hey folks, quick update: I’ve identified the reason why our Replogle analysis was not working. The full report is <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalibration_exploration.pdf|here>. The short version is that some of the NT gRNAs in Replogle apparently are not actually NT, i.e. they have effects on some genes. This renders the data too low-quality for us to use. I therefore propose we try another, higher-quality dataset. Tim has advocated for using the Schnitzler data, which has about 500k cells. We should decide whether these data are big enough.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalibration_exploration.pdf

**[2025-06-27 16:41] Eugene Katsevich**
<@U0239H5UC9E>, I’ve taken a closer look at the Schnitzler paper and they say on p. 3 that “In total, we obtained data for 214,449 cells expressing a single gRNA…” You said before that these data have about 500k cells. Am I missing something?

**[2025-06-29 17:50] Eugene Katsevich**
Hi folks, I’m having second thoughts on the Schnitzler data, due in part to its apparently smaller size and in part because of the issue that the gRNA count information is unavailable. I asked Jesse about the latter issue and he pointed me to the read data on SRA, seemingly suggesting that we can simply rerun Cell Ranger to get whatever count data we want. Which we certainly could do, but that would take time. I could try to prod him further to get Cell Ranger outputs for the gRNA data, but I admit I don’t want to be too annoying to him.

This makes me return to the question: *Would it be so bad for us to simply throw out the bad NT gRNAs and analyze Replogle that way?* To put the question another way: Do the “off-target” effects of some NT gRNAs render this gigantic, rich dataset so flawed that no one should trust any conclusions derived from that dataset? That would be an extreme position, and I would say probably not. Certainly, this dataset should be approached with caution, and there will inevitably always be an asterisk. But can we still carry out a more reliable analysis of that dataset than in the original paper? Very likely. Are there still biological insights that can be gleaned from these data? Also very likely. So is it really justified for us to abandon this dataset entirely?

Given these considerations, and the fact that we have already gone to the effort of importing this dataset, why don’t we just roll with Replogle? This would make our lives a lot easier and, hopefully, we can wrap up this manuscript much sooner than if we set out to analyze another dataset.

**[2025-06-30 04:17] Timothy Barry** _(thread reply)_
Yea, I did see that. When I analyzed and QC'ed the data, I obtained &gt;500,000 cells with a single gRNA. I'm not quite sure where this incongruence is coming from.

**[2025-06-30 04:25] Timothy Barry** _(thread reply)_
I think this would be a reasonable way to proceed, even if it is not ideal. In presenting the results for rd7, we should state that we excluded the NT gRNAs that Replogle et al identified as having off-target effects. Also, we probably should present the "full" results (i.e., the results derived from the entire set of NT gRNAs) in the appendix.

Also, what is our plan for the other two Replogle datasets? We are planning to analyze these as well, right?

**[2025-06-30 09:59] Eugene Katsevich** _(thread reply)_
&gt; In presenting the results for rd7, we should state that we excluded the NT gRNAs that Replogle et al identified as having off-target effects. Also, we probably should present the “full” results (i.e., the results derived from the entire set of NT gRNAs) in the appendix.
Certainly! This is important.

**[2025-07-01 11:23] Eugene Katsevich**
Hi folks, it appears that in our Replogle import code, we <https://github.com/Katsevich-Lab/import-replogle-2022/blob/d4cc83e101231cf94e75321fe0888a4f15beafcc/create_grna_table_4.R#L25|renamed the NT vectors>:
```mutate(vector_id_2 = factor(vector_id,
                            levels = unique(vector_id),
                            labels = paste0(ifelse(non_targeting[1], "non-targeting_", "targeting_"),
                                            "vector_", seq_along(unique(vector_id)))) |> as.character())```
This makes it challenging to know which of the vectors/gRNAs in our sceptre object need to be thrown away, following Replogle’s analysis. <@U0239H5UC9E>: How do you suggest we resolve this issue? The one solution that comes to mind is that we should update the import code to remove this renaming and then reimport the data. One lesson here is that our import code should strive to modify the original data as little as possible.
  - link: https://github.com/Katsevich-Lab/import-replogle-2022/blob/d4cc83e101231cf94e75321fe0888a4f15beafcc/create_grna_table_4.R#L25

**[2025-07-01 12:30] Timothy Barry**
Does Replogle report the individual NT gRNAs that were excluded or the vectors that were excluded? I cannot seem to remember.

**[2025-07-01 12:43] Timothy Barry** _(thread reply)_
From what I can tell looking at the code, it seems that we renamed the vector but not the individual gRNA, right? So, if individual NT gRNAs were excluded, then we would not need to know the vector IDs.

**[2025-07-01 20:41] Eugene Katsevich** _(thread reply)_
Hi <@U0239H5UC9E>: They excluded _vectors_, rather than individual gRNAs. Do you agree that this necessitates re-importing? Is there any other shortcut we can take?  You mentioned we didn’t rename the individual gRNAs but are those individual gRNAs stored somewhere in our sceptre object? If so, I guess we could try to match up the individual gRNAs with the vectors that were excluded? This is somewhat hacky, though.

**[2025-07-01 20:57] Eugene Katsevich**
Hi <@U0524GR916C>: Following <https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1751233855927769|this thread>, we’re going to set aside Schnitzler and focus our attention on Replogle. As I discussed in <https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1750718792512919|this other thread>, the issue with Replogle that was holding us back before was that some of the NT gRNAs were actually not NT. Replogle had in fact identified and removed these non-NT gRNAs from their analysis. Velten did the same, and got good calibration. *So our next task is to redo the Replogle calibration check upon throwing out the non-NT gRNAs.* At a high level, this involves three steps:
1. Identify the gRNAs that Replogle excluded;
2. Create a sceptre object without those problematic gRNAs;
3. Rerun the calibration check on the sceptre object from step 2.
I would like for <@U0239H5UC9E> to weigh in on the best way to carry out step 2. For step 1, I recommend you try to figure out how Velten did this and then try to do something similar. The relevant line in their codebase is <https://github.com/velten-group/crispat_analysis/blob/d7c1ba19cc437663da092e64f2af3fb42b5bcdf7/python/utils_guide_calling.py#L104|this one>. You will need to scan around their codebase and take a look at some of the underlying files to really understand what is going on. I recommend you clone the `crispat_analysis` repo and use Claude Code to help you. Just try not to break the bank on token usage; codebase understanding can be somewhat of an open-ended task and can lead to scanning many files. One way to try to alleviate this is to have Claude first scan the codebase once and then map it out in its `CLAUDE.md`, and then make sure it refers to `CLAUDE.md` to figure out what files are relevant as it scans.
  - link: https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1751233855927769
  - link: https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1750718792512919
  - link: https://github.com/velten-group/crispat_analysis/blob/d7c1ba19cc437663da092e64f2af3fb42b5bcdf7/python/utils_guide_calling.py#L104

**[2025-07-02 05:29] Timothy Barry** _(thread reply)_
Regarding 2, I have an idea! In the gRNA target data frame (i.e., the data frame mapping each gRNA to its target), assign the faulty NT gRNAs to a "dummy target," e.g. `dummy_target`. Then, do not use `dummy_target` anywhere in the analysis, either in the calibration check or discovery analysis. This will enable us to exclude the faulty NT gRNAs without having to update the import code too much and while keeping the `sceptre` and `ondisc` code totally intact.

**[2025-07-02 17:11] Timothy Barry** _(thread reply)_
I don't feel very strongly. Importing the data should take about half an hour. I'd probably just import the data using the main branch of `ondisc`, setting the target of the faulty NT gRNAs to `dummy` and then rolling from there.

**[2025-07-03 01:41] Louis Deutsch**
hi all, I've been using CC to explore crispat, and I've needed to take some time to reacquaint myself with this project, SCEPTRE, Replogle, and what we're doing here.

My plan is to use CC to step through the crispat code until i can see what exactly is in the data structures used for the specific filtering line. In order to do this, I'll need their data but I'm pretty rusty on what the data here is and how it's set up.
• on the cluster we have `data/external/replogle-2022/raw/kd6/K562_essential_other`, which I also have on my laptop. 
• I also need `K562_essential_library_annotation.csv` but I'm not finding this. I see other similar-looking things on the cluster, but I can't find this one. Where is it?

**[2025-07-03 01:53] Louis Deutsch**
And I guess more generally I’d really appreciate some guidance / reminders on what files I need to have in order to run the crispat code, and what they mean. I’m working on figuring this out myself / remembering how this goes, but I figured I’d ask anyway! 

**[2025-07-03 17:04] Eugene Katsevich**
Hi <@U0524GR916C>, here are the steps I recommend you take once Tim finishes including code to download Velten’s objects (please confirm <@U0239H5UC9E> that the steps below make sense):
1. Log onto the HPC3 desktop <https://hpc3-desktop.wharton.upenn.edu/|via browser>.
2. Fire up RStudio via Applications &gt; Programming &gt; RStudio.
3. Install latest versus of `sceptre` and `ondisc` via `devtools::install_github("katsevich-lab/sceptre")` and `devtools::install_github("timothy-barry/ondisc")`.
4. Clone <https://github.com/Katsevich-Lab/import-replogle-2022|import-replogle-2022>.
5. Run the last bit of `download_raw_1.R`, which Tim will have added by then, that download’s Velten’s stuff.
6. Rerun `create_grna_table_4.R`, `create_sceptre_object_kd6_5.R`, and `create_sceptre_object_rd7_5.R`. This will create the new sceptre objects.
7. As a sanity check, make sure that `gene.odm`, `grna.odm`, and `sceptre_object.rds` within `~/data/external/replogle-2022/processed/rd7` were updated recently (i.e. after you completed step 6).
8. Clone or pull <https://github.com/Katsevich-Lab/sceptre3-project-v2|sceptre3-project-v2>. 
9. Check out `sceptre3-project-v2/rd7_analysis/launch_script_remove_nts.sh`, which I have recently created, and see if this script makes sense to you.
10. Go to the terminal, navigate to `sceptre3-project-v2/rd7_analysis`, and submit this script to the cluster via `qsub launch_script_remove_nts.sh`. There will be a log file created in the same directory. Check periodically to make sure things are going ok.
11. Once the job completes, the outputs will be in `$LOCAL_SCEPTRE3_DATA_DIR"/replogle-2022/rd7-remove-nts/`. The most important one will be `plot_run_calibration_check.png`. Check it out by navigating to that directory in the terminal and running `xdg-open plot_run_calibration_check.png`.
  - link: https://hpc3-desktop.wharton.upenn.edu/
  - link: https://github.com/Katsevich-Lab/import-replogle-2022
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2

**[2025-07-07 09:27] Eugene Katsevich**
Hey <@U0239H5UC9E>, what do you make of <https://github.com/Katsevich-Lab/sceptre/pull/183|this> pull request, which claims “This package fails checks and is preventing deployment” and suggests to “update the DESCRIPTION”?
  - link: https://github.com/Katsevich-Lab/sceptre/pull/183

**[2025-07-07 13:10] Timothy Barry** _(thread reply)_
Please see the new <https://katsevich-lab.github.io/sceptre/|website>. I believe that I updated your name throughout. Aside from the book cover, do we need to update your name in any other location?
  - link: https://katsevich-lab.github.io/sceptre/

**[2025-07-07 14:47] Louis Deutsch**
I'm doing the calibration check now but just hit a problem. <@U0239E8QPH6> there was a typo in the submit script that I fixed, where in one place `$project_directory` was used instead of `$output_directory`. But it's also crashing on the first step now. I'm looking at what happened and it seems that there's an issue with parameters not being interpolated. In my `.command.sh` it tried to run
``` set_analysis_parameters.R sceptre_object_fp   response_odm_fp   grna_odm_fp   default   default   default   default   default   default   default   formula_object   discovery_pairs   positive_control_pairs   false```
instead of having the paths actually be there. I don't think this is a problem with the pipeline, but while I go following ChatGPT's suggestions, has anyone had this issue before? ChatGPT doesn't have an immediate simple fix to try so it'll be a bug hunt

**[2025-07-07 18:04] Eugene Katsevich** _(thread reply)_
As a sanity check, could you check our sceptre object versus Velten’s sceptre object (see e.g. `replogle_miscalibration_exploration.Rmd` for how to load the latter) and make sure that we’re using the same gRNAs as they are?

**[2025-07-07 18:13] Louis Deutsch** _(thread reply)_
<@U0239E8QPH6> am i comparing the `grna_target_data_frame`s for each object?

**[2025-07-07 18:13] Louis Deutsch** _(thread reply)_
specifically if `grna_id`s are the same?

**[2025-07-07 18:19] Louis Deutsch**
(starting a new thread) I'm comparing the `grna_target_data_frame`s of our new sceptre object, made with removing some NTs, and the crispat one. They have 152 gRNAs and we have 2666. All 152 of theirs do appear in ours. Here's what it looks like:

**[2025-07-07 18:20] Eugene Katsevich** _(thread reply)_
Could you specifically check their non-targeting gRNAs by doing |&gt; filter(grna_target == “non-targeting”)?

**[2025-07-07 18:37] Louis Deutsch** _(thread reply)_
ok yes, after filtering and sorting, our NT gRNAs are the exact same as theirs. `all(ours$grna_id == theirs$grna_id)` returns TRUE

**[2025-07-08 14:45] Louis Deutsch**
<@U0239E8QPH6> I installed your fork of `ondisc` and that particular branch of `sceptre`, and now I'm trying to run `import-replogle-2022/create_sceptre_object_rd7_remove_nts_add_cell_cycle_5.R` (which is the new script I made to reflect the exact operation we're doing. The only change is that this saves the sceptre_object to `processed/rd7/remove-nts-add-cell-cycle` just to keep things clear).

I'm getting an error from this line:
```sceptre_object &lt;- import_data_from_cellranger(directories = directories,
                                              moi = "low",
                                              grna_target_data_frame = grna_target_data_frame,
                                              use_ondisc = TRUE,
                                              directory_to_write = directory_to_write,
                                              compute_cell_cycle = TRUE)```
the error is
```Error in ondisc::create_odm_from_cellranger(directories_to_load = directories,  : 
  unused argument (compute_cell_cycle = compute_cell_cycle)```
I think the exact issue comes from `sceptre:::import_data_from_cellranger_use_ondisc` which does
```out &lt;- ondisc::create_odm_from_cellranger(directories_to_load = directories, 
        directory_to_write = directory_to_write, write_cellwise_covariates = FALSE, 
        grna_target_data_frame = if (vector_supplied) 
            grna_target_data_frame
        else NULL, compute_cell_cycle = compute_cell_cycle)```
and we get an error. I'm definitely using `ekatsevi/ondisc` and that particular branch of sceptre. What's the fix?

**[2025-07-09 13:31] Timothy Barry**
Book <https://timothy-barry.github.io/sceptre-book/|cover updated>.
  - link: https://timothy-barry.github.io/sceptre-book/

**[2025-07-09 17:01] Louis Deutsch**
I'm running the calibration check with Velten's sceptre object and I'm curious about the best way to do this.

I have an R script that contains this:
```library(sceptre)
library(ondisc)

crispat_sceptre_object <- readRDS("~/data/external/braunger-2024/sceptre_object_RPE1.rds")
crispat_sceptre_object <- crispat_sceptre_object |>
  run_calibration_check( n_calibration_pairs = 1000000, calibration_group_size = 1,
                         parallel = TRUE, n_processors = 20)
saveRDS(crispat_sceptre_object,
        file = "/home/stat/jdeu/data/projects/sceptre3/replogle-2022/velten-calib/crispat_sceptre_object_run_calib_check.rds")```
and then a launch script
```#$ -pe openmp 2
#$ -l m_mem_free=4G

Rscript run_calib_velten.R```
which i submit via `qsub`. It's running, but I only have 1 HPC slot used by this. How can i tell if this is actually running in parallel? Should I use nextflow even though it's a smaller sceptre object?

(also the hard-coded paths are just temporary while I'm figuring this out)

**[2025-07-09 18:00] Eugene Katsevich** _(thread reply)_
Good question. It’s a bit hard, and I’m not sure we have those row names in our object (Tim, let me know if I’m wrong). What I’ve done before to this end is matched on `response_n_nonzero` and `response_n_umis`, which did a pretty good job. See lines 99-195 <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalibration_exploration.Rmd|here>.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalibration_exploration.Rmd

**[2025-07-09 18:39] Louis Deutsch** _(thread reply)_
it's not as generic as yours, but with ChatGPT's help I got a nice way to do this just for this particular case

```our_covariates <- ours |> get_cell_covariates() |> as.data.frame()
crispat_covariates <- theirs |> get_cell_covariates() |> as.data.frame()

our_counts <- our_covariates %>%
  count(response_n_nonzero, response_n_umis, name = "n_matches")

crispat_covariates_matched <- crispat_covariates %>%
  left_join(our_counts, by = c("response_n_nonzero", "response_n_umis")) %>%
  filter(n_matches == 1) %>%
  select(-n_matches) %>%        
  left_join(our_covariates, by = c("response_n_nonzero", "response_n_umis")) ```

**[2025-07-09 21:17] Louis Deutsch**
I'm comparing the `covariate_data_frame`s and have some differences. I used Gene's method to filter Velten's `covariate_data_frame` down to just the NT cells, and then I left join'd our sceptre_object's `covariate_data_frame` to this using `response_n_nonzero` and `response_n_umis` (and only keeping the cells with a unique match). I have 12,172 cells after this. (my code <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalib_compare_remove_nts.Rmd|here>)

G2M and S score are not quite the same, but are very correlated (> 97%), and have a nice symmetric football-shaped scatterplot which looks unconcerning. We also have the same `batch` (modulo a renumbering) and `response_p_mito`.

The one that stands out is`grna_n_nonzero`, and `grna_n_umis` to a much lesser extent. I'm guessing this is because we computed ours from a larger data set than they did? I've included a scatterplot of `grna_n_nonzero` below.

I'm running the calibration check on Velten's sceptre_object and that's about 70% done. I'll update as that gets done, but in the meantime it seems like maybe this difference in covariates could explain the difference?
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/rd7_analysis/replogle_miscalib_compare_remove_nts.Rmd

**[2025-07-10 17:38] Eugene Katsevich**
<@U0239H5UC9E>: Any thoughts on how to proceed? I guess we have to ask ourselves: Why are we doing this calibration check? Probably to justify some big trans analysis we’re going to do next. If we are to follow sceptre’s philosophy, we like having good calibration for as many calibration pairs as we have discovery pairs. If we try a gajillion trans pairs (to showcase our computational performance), then following that philosophy, we should do a gajillion calibration pairs, which we know won’t work. However, as Tim has noted before, technically you need good calibration out to roughly whatever is the BH threshold applied to your discovery pairs. We could run the trans analysis and find out…

**[2025-07-10 20:08] Timothy Barry**
Nice work! Ideally, we would show calibration out to the number of pairs in the trans analysis. However, we are unlikely to achieve that on this dataset. The next best thing would be to show calibration out to the BH rejection threshold. We even could draw the BH rejection threshold onto our calibration check plot. I think we should replicate this analysis within our own codebase/pipeline and then move onto the trans analysis, finding where the BH threshold is.

**[2025-07-10 20:09] Timothy Barry**
I wouldn't be surprised if the BH threshold were around 1e-4 or so.

**[2025-07-10 20:37] Timothy Barry**
To put it differently, attaining calibration out to the Bonferroni threshold is ideal, but attaining calibration out to the BH threshold is probably adequate (assuming we are controlling FDR).

**[2025-07-17 15:30] Eugene Katsevich**
- Tim: Outline figures
- Gene: Construct some kind of a roadmap
- Louis: Get started on playing with CLEANSER and crispat, e.g. on Gasperini

**[2025-07-18 09:36] Eugene Katsevich** _(thread reply)_
In general, testing so many pairs highlights that our practice of matching the number of negative control pairs to the number of discovery pairs can be somewhat conservative, particularly if the number of detectable discovery associations is nontrivial. If we were to rethink the logic of the calibration check from scratch, maybe we would run the discovery analysis _first_, then extract the BH threshold, then run the calibration check _out to that BH threshold_. On the other hand, if the calibration check doesn’t work out well, and you start doing things like adding covariates, you would need to run the discovery analysis and calibration check in a loop, which seems suboptimal. The current practice of (1) iteratively getting your calibration check to work followed by (2) running your discovery analysis once with the right covariates seems more natural. So…:man-shrugging::skin-tone-2:

**[2025-07-18 11:54] Eugene Katsevich**
<@U0524GR916C>: <https://github.com/Gersbachlab-Bioinformatics/CLEANSER|Here>’s the CLEANSER repo.
  - link: https://github.com/Gersbachlab-Bioinformatics/CLEANSER

**[2025-07-18 19:29] Timothy Barry** _(thread reply)_
>  Do we recommend users “run the calibration check by running the pipeline normally and then switch to the massive trans-analysis”?
Yes, I think so. For reasons we've discussed, it's probably not necessary to analyze more than ~10,000,000 negative control pairs, and it should be feasible to do this using the standard R interface (ideally with parallelism).

And I agree with your comments. One thing we may consider doing is drawing a horizontal line on the negative control QQ plot indicating the BH threshold. (Of course, this only would be available after running the discovery analysis.) We automatically could add this line whenever the user calls `plot_run_calibration_check()` after having called `plot_run_discovery_analysis()`. This functionality may help users determine whether their calibration is OK if there is mild departure from uniformity or if there are not enough negative control pairs to match the number of discovery pairs.

**[2025-07-21 18:22] Eugene Katsevich** _(thread reply)_
The most relevant would be Tim’s <https://github.com/Katsevich-Lab/undercover-grna-pipeline|Nextflow pipeline> to run calibration checks on a bunch of methods for a bunch of datasets, where the methods are implemented using a unified interface in <https://github.com/Katsevich-Lab/lowmoi|this R package> (the R package uses `reticulate` to run Python-based methods). The invocation of this pipeline can be found <https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/main/pipeline_launch_scripts/undercover/grp_size_1_0523/grp_size_1.sh|here>, as part of that manuscript’s reproduction repo.
  - link: https://github.com/Katsevich-Lab/undercover-grna-pipeline
  - link: https://github.com/Katsevich-Lab/lowmoi
  - link: https://github.com/Katsevich-Lab/sceptre2-manuscript/blob/main/pipeline_launch_scripts/undercover/grp_size_1_0523/grp_size_1.sh

**[2025-07-24 13:25] Eugene Katsevich**
Hi folks, I’ve scoped out <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/competitors.md|our competitors>, jotted down a rough <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/benchmarking.md|benchmarking outline>, started some sort of <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/todo.md|to do list> that possibly would be better as a GitHub project, and started a <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/software-development.md|list of software development items>, which for now just has some comments on our pipeline from Andreas Gschwind. Let’s discuss during our 2pm meeting today.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/competitors.md
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/benchmarking.md
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/todo.md
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/main/manuscript/software-development.md

**[2025-07-30 10:00] Eugene Katsevich**
Hi <@U0524GR916C>, I see that our `environment.yml` files got much more minimal in <https://github.com/Katsevich-Lab/sceptre3-project-v2/commit/c798fe3500c0003d67694f12f1605bc5de774b52#diff-c3d3cf362d090d4b12740a251306fd0f06b6c3902cb706c0fcfad9a239b7b76f|a recent commit>. That is pretty cool! Could you explain why it’s ok to do this? Is it because the `pip` command will pick up all those dependencies automatically?
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/commit/c798fe3500c0003d67694f12f1605bc5de774b52#diff-c3d3cf362d090d4b12740a251306fd0f06b6c3902cb706c0fcfad9a239b7b76f

**[2025-07-30 13:02] Eugene Katsevich**
<@U0239H5UC9E>: I guess `sceptre v0.3.0` (accompanying the 2024 GB paper) did not have the mixture model gRNA assignment, so we should exclude it from gRNA assignment benchmarking, right? Do you think we should include the current sceptre R package’s mixture model in our benchmarking, in addition to the Nextflow pipeline’s mixture model gRNA assignment functionality?

**[2025-07-30 13:07] Eugene Katsevich**
<@U0524GR916C>: I suspect that you will eventually have to create two Docker images: one for `sceptre v0.3.0` (the previously published version of our software) and one for the current version of the sceptre R package. Because these packages aren’t on CRAN or Bioconductor, there’s not a clean way to use conda to manage these packages. On the other hand, it’s easy to create a Docker image containing these packages, and these images can be read by Nextflow via Singularity, which is on HPC3.

In the context of the gRNA assignment specifically, I suspect we will want to benchmark the sceptre R package’s gRNA assignment performance, as distinct from the sceptre pipeline’s performance. So for the gRNA assignment piece, you’ll already need to make a Docker image with the current sceptre package. Hopefully we can get confirmation from Tim during today’s meeting.

**[2025-07-30 14:29] Eugene Katsevich** _(thread reply)_
Ideally, I guess we would include the thresholding approach as a comparison point for the more sophisticated mixture methods? It would be disappointing if the mixture methods didn’t outperform the baseline of thresholding. From crispat’s benchmarking, this was not super clear. I noticed in general that the crispat folks steer clear of making any definitive conclusions about which methods are better and which methods are worse. This is surprising for a paper whose entire purpose is comparing different methods.

**[2025-07-30 14:44] Eugene Katsevich**
Hi folks, here’s a proposed agenda for today’s meeting at 4pm:
• Update from Louis on guide assignment benchmarking
• Discuss versions of sceptre we will benchmark
• Discuss whether to include thresholding in comparison
• Discuss next steps for Louis
    ◦ Adding datasets
    ◦ Adding methods (e.g. sceptre package)
    ◦ Adding statistical metrics
    ◦ Adding computational metrics
• If time permits, discuss user issues

**[2025-07-30 18:24] Eugene Katsevich**
For Replogle rd7, make the summed-count gRNA matrix available to the methods. This data is available at `~/data/external/replogle-2022/processed/rd7/grna.odm` on HPC.

**[2025-08-12 07:00] Louis Deutsch**
Hi all, in advance of our meeting, here's how I'm currently running sceptre in the nextflow pipeline.
• I'm loading the `response_matrix`, `grna_matrix`, and `grna_target_data_frame` into memory so I can make a `sceptre_object` via `import_data()`. 
• I have a nextflow module that calls an R script, which does this
• version control is via a lock file: when present, those exact versions of R / required packages are used. If not present, it installs the latest versions and makes a lockfile to record what happened. I think this is fine for the R packages themselves. I'm not sure if this is good enough for underlying C++ things. But maybe that also should apply to other packages and we routinely ignore this?

**[2025-08-14 15:37] Eugene Katsevich**
Hey folks, some todos from today’s meeting:
• Tim
    ◦ Figure out if there are a bunch of commits missing. Hopefully not! Also, look into the index out of bounds issue. Did you fix it? There’s a commit from last November that might be relevant.
    ◦ Do more targeted check of NAs in the grna_target_data_frame.
• Louis
    ◦ Figure out Nextflow portions of Boris’s issue
    ◦ Look into issue 162 (specifically, the gRNA assignment part, not the plotting part). Note that the data are linked directly in the GitHub issue.
    ◦ For issue 146: Add a check that covariates being passed into gRNA assignment mixture model are not multicollinear and a test that makes sure that multicollinear inputs yield appropriate error message.
    ◦ Get sceptre into a container and then into the guide assignment benchmarking pipeline.
• Gene
    ◦ Make small bug fix of plotting associated with issue 162

**[2025-08-20 10:58] Timothy Barry**
Hey guys, I fixed our check for NAs within the gRNA target data frame, discovery pairs data frame, and positive control pairs data frame. We are now <https://github.com/Katsevich-Lab/sceptre|passing the command check>!
  - link: https://github.com/Katsevich-Lab/sceptre

**[2025-08-21 13:16] Louis Deutsch**
In advance of our meeting, here's ChatGPT's summary of Boris' issue (emojis and all)

:mag: ChatGPT Report: User [Boris] Issue with `run_discovery_analysis_trans`
:technologist: User context
• Running SCEPTRE pipeline on SLURM
• The user requested *76 processors* on a *single exclusive node*
• Job wall time set to *10 hours*
• `discovery_pairs = "trans"` with a large-scale pooled CRISPR screen dataset (similar to Replogle et al.)
• Pipeline version: `sceptre-pipeline@main`
• SCEPTRE version: `0.10.2`
• ONDISC version: `1.2.0`
• R version: `4.4.2`
:bar_chart: Dataset Characteristics
• *643,752 cells* (325,049 post-QC)
• *36,601 response genes*
• *2,536 targeting gRNAs*
• *130 non-targeting gRNAs*
• ~90 million total discovery pairs in `trans` mode
:exclamation:️ Reported Problems
1. :warning: Subscript out of bounds
```subscript out of bounds (index 530024 >= vector size 530024)```
• Occurs inside `ondisc:::compute_n_ok_pairs_ondisc()`
• Likely a *1-based R to 0-based C++ indexing bug*
• Very likely triggered at the *end of a pod*, e.g., task 48
2. :stopwatch: Task timeout
```Process exceeded running time limit (6h 56m 40s)```
• Individual pod (e.g., pod 48) takes *>6 hours* on (up to) 76 CPUs
• This suggests a *performance bottleneck* or *overly large pods*
3. :warning: Not all pods are running
• 175 discovery pods are created
• Only ~97 appear to be scheduled (`0 of 97` in logs before crash)
• Remaining pods never get launched — likely because the early pod crash halts the process
:hammer_and_wrench: Debug Analysis
:jigsaw: Resource usage & scheduling
• Each pod handles ~514,000 discovery pairs (`90M ÷ 175`)
• At ~6+ hours per pod, *entire pipeline would take 44–50 hours* if run serially
• Nextflow likely stops scheduling further pods once any pod fails
• The user may *misinterpret the scheduling stall* as a config error — but it’s a downstream result of the early task crash
:bomb: Indexing bug in C++ backend


:female-technologist: Recommendations for the SCEPTRE dev team
1. *Fix indexing bug*
2. *Improve pod scheduling behavior*
    ◦ Add an `errorStrategy = 'ignore'` or `retry` to avoid halting all pods on single failure
3. *Add runtime guidance*
    ◦ Warn users about expected time/memory when running `trans` analysis on datasets with >10M pairs
4. *Improve logging*
    ◦ Log max indices before launching C++ code to help detect off-by-one errors preemptively
5. *(Optional)* Break large pods earlier for better runtime distribution

**[2025-08-22 23:16] Louis Deutsch** _(thread reply)_
ok, sounds good. In the meantime, I'm having trouble making sense of the contents of the various directories there. There are 3 Gasperini directories in `data/external` : `gasperini-2019-sra`, `gasperini-2019-v2` and `gasperini-2019-v3`. I'm guessing v3 is the right one? In `data/external/gasperini-2019-v3/at-scale/processed` there are the usual odm files and sceptre_object.rds. For my benchmarking, I will need to have the data available as a h5ad file for crispat and pertpy, and mtx for CLEANSER. Is there a way to convert these odm files to those formats? Otherwise, I'm not sure which parts of the less-processed data are the things I need.

**[2025-08-26 01:37] Louis Deutsch**
Note the relative sizes! (this is Gasperini)
```150M Aug 26 00:56 grna.h5ad
 42M Aug 26 00:53 grna.odm```

**[2025-08-26 10:35] Eugene Katsevich**
Hey <@U0524GR916C>, could you give me an update on where you’re at with gRNA benchmarking? I’ve only partially followed the discussion with Tim. Thank you!

**[2025-08-26 11:43] Louis Deutsch** _(thread reply)_
I needed to have the Gasperini gRNA expression data as a h5ad file for crispat and pertpy (CLEANSER uses mtx), so much of the above is just me converting an odm to h5ad. Tim suggested that the conversion function could go in ondisc hence the discussion. Although it seems that writing a stream to a h5ad is tricky in pure R, so for now I’m doing this conversion in-memory by just loading every row of the odm into a dgRMatrix. 

Right now I’m trying to get the pipeline to run on the HPC for pertpy and crispat. I’m currently troubleshooting some conda things that are happening there that didn’t happen on my laptop 

**[2025-08-26 13:52] Eugene Katsevich** _(thread reply)_
Thanks for the update!
&gt; I needed to have the Gasperini gRNA expression data as a h5ad file for crispat and pertpy (CLEANSER uses mtx), so much of the above is just me converting an odm to h5ad. Tim suggested that the conversion function could go in ondisc hence the discussion. Although it seems that writing a stream to a h5ad is tricky in pure R, so for now I’m doing this conversion in-memory by just loading every row of the odm into a dgRMatrix.
I would say we prioritize speed, so if doing streaming stuff or adding functions to ondisc is going to take longer than a brute force approach to getting data in the right format, then we should choose the brute force approach.
&gt; Right now I’m trying to get the pipeline to run on the HPC for pertpy and crispat. I’m currently troubleshooting some conda things that are happening there that didn’t happen on my laptop
Sounds good! Let me know if you find yourself not making rapid progress on troubleshooting the conda things.

**[2025-08-27 03:11] Louis Deutsch** _(thread reply)_
I've decided to use singularity / apptainer for the python modules in our nextflow pipeline, in addition to the sceptre module. pertpy has a ton of dependencies and it kept crashing on hpc3 while trying to resolve the environment.yml via conda. I switched to a faster version called mamba but that didn't help. I made a lockfile on my laptop, but using a lockfile is a little more complicated. At this point the best answer is basically to build these environments elsewhere and then make them available to the nextflow pipeline. That is getting annoying enough that I think singularity is better, especially since we're using it already. The sceptre singularity image is ~500MB right now, so we'd likely have a few GB of images to store somewhere.

**[2025-08-27 17:04] Louis Deutsch** _(thread reply)_
Yes this image is just for pertpy. I am planning on doing separate images for each module. pertpy does have a ton of dependencies, including things we likely don't need like matplotlib, so it definitely can be trimmed   down. It'll just take more work.

And yes i am able to use this image locally, so I'm pretty sure it's a HPC thing.

**[2025-08-27 17:07] Eugene Katsevich** _(thread reply)_
> pertpy does have a ton of dependencies, including things we likely don’t need like matplotlib, so it definitely can be trimmed  down. It’ll just take more work.
That’s their fault, not ours. Unless having a large image poses insurmountable challenges, I’d say keep it. We are not responsible for pruning the dependencies of another software.

**[2025-08-28 02:21] Louis Deutsch** _(thread reply)_
i think i got it! it's having an issue writing the output but I think i'll be able to figure that out pretty quickly in the morning. At this point it is running pertpy on a subset of the Gasperini data using a singularity image that i made.

**[2025-08-28 02:53] Louis Deutsch** _(thread reply)_
Ok yep it’s fully working now! So today I made a pertpy singularity image, and ran the pertpy part of the nextflow pipeline on a chunk of the gasperini data. Tomorrow I’ll get the images made for the other python packages, and then I’ll get sceptre working on HPC with its singularity image. 

**[2025-08-28 15:10] Timothy Barry**
gRNA assignment metrics and code

**[2025-08-28 15:50] Eugene Katsevich**
I didn’t find any “experimental MOI” reported in Gasperini or Replogle, so we might have to exclude that metric from our list.

**[2025-08-28 16:03] Timothy Barry**
Could we assume the MOI in Replogle to be 1?

**[2025-08-29 00:12] Louis Deutsch**
I've gotten a better sense of what was going wrong with the image yesterday, and I'm posting here for posterity.

By default the singularity image mounts $HOME, which meant `~/.local` was available to the container. This led to the image importing an old version of `anndata` from there, which then didn't have a submodule expected by other packages. This caused the error `ModuleNotFoundError: <http://anndata.io|anndata.io>`. This can be fixed via `autoMounts  = false` (preventing $HOME from mounting) and `runOptions = "--cleanenv --containall"`  (ensuring the image only uses packages actually _in_ it when invoked) in the `singularity` section of the `nextflow.config`.

But then this prevented the image from having access to any of our files at all. The solution to that is to use `bind` to mount just the exact directories that we need. This is what I added to the nextflow.config to get it to work (added at the top level):
```singularity {
  enabled     = true
  autoMounts  = false
  cacheDir    = "${projectDir}/.apptainer-cache"
  runOptions  = """
    --cleanenv --containall \
    --bind ${projectDir}:${projectDir} \
    --bind ${params.dataset_base_dir}:${params.dataset_base_dir}:ro \
    --bind ${params.out_base_dir}:${params.out_base_dir}
  """.trim()
}```
The `:ro` in the second `bind` makes that one read-only. This only applies to modules that use the `container` command to use a `.sif` image, as pertpy does.
  - link: http://anndata.io

**[2025-08-30 02:17] Louis Deutsch**
in the sceptre help page `?assign_grnas`, it never actually says what the default method is. I'd find that very helpful to have there, especially in the description of the argument `method`:
```method	
(optional) a string indicating the method to use to assign the gRNAs to cells, one of "mixture", "thresholding", or "maximum"```

**[2025-08-30 03:19] Louis Deutsch**
As I am focusing on my thesis proposal for maybe the next week, I wrapped up the benchmarking a bit while it is fresh.

Where we are now: the following run on both my laptop and HPC, using a subset of Gasperini:
• sceptre with a singularity image
• pertpy with a singularity image
• cleanser with conda 
• crispat with conda
    ◦ the crispat module does fail on the test data, with an error in the main function `ga_poisson_gaussian()` thrown by `plot_univariate_histogram()`. The actual error is `ValueError: `bins` must be positive, when an integer` . I'm pretty hopeful this is just from the test data being too small, and will go away with real datasets. Regardless, the python environment is set up correctly by nextflow
Beyond that, I don't have any more TODOs for the benchmark pipeline as far as the infrastructure goes, with the caveat it's never been tested on anything big.

Unrelated, how do we want to handle cells that get assigned to no gRNA? That happened with a smaller subset of Gasperini.

**[2025-09-01 11:19] Timothy Barry** _(thread reply)_
Good catch. Would you mind specifying the default in the documentation? Mixture for high-MOI, maximum for low-MOI.

**[2025-09-01 11:21] Timothy Barry** _(thread reply)_
Nice. In high-MOI, cells that do not contain a gRNA are retained for use in the control group. In low-MOI, cells that do not contain a gRNA are filtered out.

**[2025-09-03 18:58] Timothy Barry**
I see we bumped to version `0.10.3`. I just pushed a patch for the `seq.default(start[i], stop[i])` error Louis identified. I got slightly impatient and <https://github.com/Katsevich-Lab/sceptre/releases/tag/v0.10.3|issued the release> on Github. We can update the release documentation tomorrow if we choose. Note that the release includes this latest bug fix.

Also, for tomorrow's meeting, I wanted flag the following two unresolved items:

1. From Veronica: <https://github.com/Katsevich-Lab/sceptre/discussions/180#discussioncomment-13595688>
2. From Seth: <https://github.com/Katsevich-Lab/sceptre/discussions/181>
Both relate to plotting and could be considered feature requests.
  - link: https://github.com/Katsevich-Lab/sceptre/releases/tag/v0.10.3
  - link: https://github.com/Katsevich-Lab/sceptre/discussions/180#discussioncomment-13595688
  - link: https://github.com/Katsevich-Lab/sceptre/discussions/181

**[2025-09-04 14:51] Louis Deutsch**
Does `crispat` have any parallelization at all? I'm trying to see what a "fair" way to run it is. They don't seem to have any places where parallelization can be specified. The only place I/ChatGPT can come up with is MKL / OpenBLAS threads under the hood, so it's suggesting to pin those for a fair comparison. Any thoughts on this?

**[2025-09-04 16:07] Louis Deutsch** _(thread reply)_
Ok I think what was confusing me is that I don’t think there is any parallelization in the particular function we have in `run_crispat.py`. They have others that do parallelize by default with dask, but not `ga_poisson_gauss()`. 

**[2025-09-04 16:09] Louis Deutsch** _(thread reply)_
They do have some parallelization that I missed, but it’s not in the gRNA assignment strategy that I was using I think. I’ll update later as I learn more!

**[2025-09-18 10:56] Louis Deutsch** _(thread reply)_
I can update asynchronously today. I’m not stuck on anything. I looked into adding a `num_grnas` flag for downsampling, but that turned out to be hard, since CLEANSER is a command line tool and crispat expects a path to the data, not a data object already. I decided just to make a single downsampled version of the data instead, for testing. I also have restructured the code some since I hadn’t been thinking about association testing at all when I wrote it before. I’m making good progress on all the pieces we talked about on Monday 

**[2025-09-19 09:28] Eugene Katsevich** _(thread reply)_
Here’s a message from Boris:
&gt; I am reaching out regarding issue #177 on SCEPTRE’s Github. I am running the SCEPTRE pipeline to do a massive-scale trans analysis on the Jurkat data from <https://urldefense.com/v3/__https://www.nature.com/articles/s41588-025-02169-3__;!!IBzWLUs!Vn4YOzunYpiT54p2LKwQG1zqr0XD1m6xlUrQ3rHqdzDlPdBcb7548Nhmv_sWpVu1fIp7EbaCq0mdz1dsAOYEMPLVo-TjpkZ8ESbuEU1W$|Nadig et al. (2024)>, but the pipeline stalls during run_discovery_analysis_trans.
&gt;  
&gt; As advised, I’ve placed the sceptre_object and the code in a shared drive here <https://drive.google.com/drive/folders/1R3QFsekMnPUD8n_ztuY5UKi0oB0OVECV?usp=sharing>
&gt;  
&gt; The code I used to import the data is available here: <https://github.com/boris-vasilev/SCEPTRE_CRISPR_screen/tree/master/src/import_nadig_2024>, which I modeled after your repo <https://github.com/Katsevich-Lab/import-replogle-2022>.
&gt;  
&gt; Would you be able to advise on the execution of the pipeline? From what I understood, you were able to run a similar analysis on the Replogle RPE1 dataset in ~6 hours, so I suspect I may have mis-specified some arguments to the Nextflow pipeline.
  - link: https://urldefense.com/v3/__https://www.nature.com/articles/s41588-025-02169-3__;!!IBzWLUs!Vn4YOzunYpiT54p2LKwQG1zqr0XD1m6xlUrQ3rHqdzDlPdBcb7548Nhmv_sWpVu1fIp7EbaCq0mdz1dsAOYEMPLVo-TjpkZ8ESbuEU1W$
  - link: https://drive.google.com/drive/folders/1R3QFsekMnPUD8n_ztuY5UKi0oB0OVECV?usp=sharing
  - link: https://github.com/boris-vasilev/SCEPTRE_CRISPR_screen/tree/master/src/import_nadig_2024
  - link: https://github.com/Katsevich-Lab/import-replogle-2022

**[2025-09-25 15:05] Eugene Katsevich**
Hey <@U0239H5UC9E>, Louis and I took the liberty of meeting without you today. Louis’s update: He is closing in on running the guide assignment benchmarking pipeline on 4 methods x 2 datasets (subsetted for an initial run).

**[2025-09-25 15:16] Eugene Katsevich**
<@U0524GR916C>: Please let me know when you implement the sceptre guide assignment with dummy response matrix, then when you kick off the first run of the benchmarking pipeline, then when it completes.

**[2025-09-25 16:47] Timothy Barry**
Hi all, I am sorry for missing today's meeting. I was out of office Monday (hiking) and so I lost track of what day it was. It sounds like things are going well. Thank you for the update! Please let me know if there is something I could do to help with the gRNA benchmarking.

**[2025-10-01 23:58] Louis Deutsch**
hey, i'm having trouble getting sceptre to run with the mixture assignment method. I keep getting the errors related to "redundant information".

For example, with the full gasperini dataset, i start with this:
```sceptre_object &lt;- import_data(
  response_matrix = response_matrix,
  grna_matrix = grna_matrix,
  grna_target_data_frame = grna_target_df,
  moi = "low"  # Assume low MOI for now
)```
if i then do
```sceptre_object &lt;- set_analysis_parameters(sceptre_object)```
i get the error `The `formula_object` contains redundant information.`

I don't care about extra covariates for now, so i can get past this by just doing
```sceptre_object &lt;- set_analysis_parameters(sceptre_object, formula = ~1)```
But then when i do
```assign_grnas(sceptre_object, method = "mixture")```
I get the error `The `formula_object` contains redundant information`.

What should i be doing here to get the mixture method to work? I'm just running with maximum for now

**[2025-10-02 09:05] Eugene Katsevich** _(thread reply)_
The mixture model for gRNA assignment has its own formula object (distinct from that used in the association test), specified through `assign_grnas()` rather than through `set_analysis_parameters()`. See <https://timothy-barry.github.io/sceptre-book/assign-grnas.html#sec-mixture_method|the book>.
  - link: https://timothy-barry.github.io/sceptre-book/assign-grnas.html#sec-mixture_method

**[2025-10-02 09:10] Eugene Katsevich** _(thread reply)_
Below are some thoughts on how to diagnose what redundant information is being referred to, you can do the following. First, figure out what formula object is being used by default, you can use
```formula_object = auto_construct_formula_object(
cell_covariates = sceptre_object@covariate_data_frame,
include_grna_covariates = TRUE
)```
Then, you can construct the actual covariate matrix used via
```  covariate_matrix &lt;- convert_covariate_df_to_design_matrix(
    covariate_data_frame = cell_covariate_data_frame,
    formula_object = formula_object
  )```

**[2025-10-02 14:10] Louis Deutsch** _(thread reply)_
I think i have the actual pipeline itself done and ready. The issues I'm troubleshooting now are (1) that sceptre issue i posted about for the mixture method; (2) cleanser's MCMC chains are wandering off to infinity and causing errors; (3) pertpy's caching doesn't seem to be working so it takes more than an hour to run pertpy even on tiny datasets since it's recompiling JAX every time or something like that

**[2025-10-02 15:08] Louis Deutsch** _(thread reply)_
hi Gene, no problem, I figured it was that! Tim and I met anyway. We addressed the sceptre issue above, and figured out some things about cleanser (it's running now with the changes). If cleanser doesn't work with these changes we will ping one of the main devs who Tim knows

**[2025-10-02 15:19] Eugene Katsevich** _(thread reply)_
What was causing the sceptre issue? What about the pertpy issue?

**[2025-10-02 15:20] Louis Deutsch** _(thread reply)_
For sceptre, it was the dummy response matrix. So it’s nice that it’s a quirk of how I’m going this, not sceptre itself. I don’t know about pertpy yet

**[2025-10-02 15:25] Eugene Katsevich** _(thread reply)_
&gt; I don’t know about pertpy yet
Is this the next step then?

**[2025-10-02 15:35] Louis Deutsch** _(thread reply)_
Yeah, pertpy does run, it just has this startup run cost so making changes is slow.

Here’s the exact state of each method:
• crispat is ok, aside from crashing on small test datasets due to some unimportant error with its automatic plotting functions 
• Sceptre seems ok now
• Pertpy runs but isn’t caching correctly
• cleanser is being tested right now. It works on small datasets but is having an MCMC Stan issue on medium+ ones. Not sure if the changes we made will fix that. If I don’t get it working today we will ping Thomas 

**[2025-10-03 13:32] Louis Deutsch** _(thread reply)_
got cleanser working. I needed to adjust some arguments a little but it has now successfully run on replogle-rd7 medium.

**[2025-10-06 14:47] Louis Deutsch**
hi all, i want to make sure i'm collecting the results from sceptre correctly. Right now, i have this:

```assignment_matrix &lt;- get_grna_assignments(sceptre_object)
cells_with_assignments &lt;- assignment_matrix |&gt; apply(2, which)```
so `cells_with_assignments` is a list with one element per cell, and one entry for every grna that got assigned.

On my small replogle-rd7 test data (5000 cells, 500 grnas, low moi, mixture method) i have this result:
``` sapply(cells_with_assignments, length) |&gt; table()

   0    1    2    3
3879 1022   96    3```
so most cells get assigned to 0 or 1 grnas, but some get assigned to 2+. How do we want to handle this? Right now I am just ignoring the cells that got nothing assigned to them for the outputs, but I'mm not sure about the ones with multiple. Should i just consider all of those to be expressed, or do we somehow pick one of the multiple expressed grnas?

**[2025-10-06 15:19] Timothy Barry**
Could you just save the entire outputted gRNA-to-cell assignment matrix for now? We'll need to decide exactly what metric we want to use to assess the quality of the assignments.

**[2025-10-06 18:57] Louis Deutsch**
follow-up: this will still bake the probability threshold into the results. Are we ok with that for now, or would it maybe be better to save the outputs where it's still a sparse matrix but for grna-cell pairs where prob &gt; threshold, we save the probability instead of just the boolean for it exceeding the threshold?

**[2025-10-06 19:37] Eugene Katsevich** _(thread reply)_
For this first run, I think it’ll be simpler to just have boolean assignments. We can then reconsider this question later.

**[2025-10-06 20:07] Louis Deutsch**
another thing: i think we may need to use a GPU to get a realistic run time for pertpy. It is absolutely glacial on CPUs. My understanding (what i've pieced together from chatgpt and the docs) is that pertpy uses JAX + XLA, and we end up compiling tons of slightly different computational graphs for doing gRNA assignment, which is waaaay slower on a CPU. In order to get realistic speeds I need to use JAX with CUDA. I'm currently waiting to get a GPU assigned to me on HPC to test replogle-rd7 small (500 grnas, 5000 cells) to see the difference. As it stands right now, it takes over 2 hours to run pertpy on this dataset with a CPU because of how much time it spends compiling. Which is absolutely ridiculous. I was thinking this was an issue with me not caching stuff correctly, and maybe for reruns of the exact same analysis that is true, but i think in general this is a deeper issue than caching.

But maybe this is getting ahead of ourselves, and this actually is an interesting result, namely that pertpy is bad if all we have is a CPU?

**[2025-10-06 22:38] Louis Deutsch**
i just finished running pertpy on replogle-rd7_small (the 500x5000 one) and even with a GPU it took 1.25 hours. sceptre takes 5 minutes for the same thing. The problem is that JAX/XLA is doing many, many small compilations over and over. I've been trying ChatGPT's different suggestions to address this, but it seems like there isn't really anything we can do short of precomputations / grouping the gRNAs by shape [the number of cells with UMI count > 0 for that gRNA] in advance so that bigger batches get compiled (since the tasks for gRNAs with different shapes are compiled separately). When we haven't done the exact same analysis already, I don't think caching will help. Since we're not helping these methods perform better, I'm going to leave this for now.

I am now kicking off a full CPU run of everything on my medium subsets of gasperini and replogle.

**[2025-10-06 22:41] Louis Deutsch** _(thread reply)_
i'm not really sure what to do about this in order to get a "fair" comparison with pertpy. I am not trying to do parallelization with this GPU. I was just using it so that JAX could use CUDA to speed up this compilation, and it did reduce the time from around 2.5 hours to 1.25, but that's still terrible when the actual algorithm is a tiny fraction of that time

**[2025-10-06 22:42] Louis Deutsch** _(thread reply)_
so maybe pertpy will just look terrible until we do our comparison where we give it more resources?

**[2025-10-07 09:57] Eugene Katsevich**
Hi <@U0524GR916C>, thanks for keeping us updated on your trials and tribulations with pertpy. Generally speaking, we should run a software in the way that the developers recommend, making sure to follow their documentation, but not bending over backwards. If they recommend running on a GPU, then we should run on a GPU. If they have any resources (e.g. documentation or example code) to help with caching, then we should use those resources. If not, then we should run it “out of the box” and the runtime is what it is. At this point, you have invested enough time (probably more…) in figuring out pertpy. I recommend you exclude it from your CPU-based list of methods and get on with running the pipeline with whatever methods remain. Then, I suggest you throw pertpy on a GPU (with parallelization) when we get to the “at-scale” comparison.

**[2025-10-07 10:14] Timothy Barry** _(thread reply)_
I do think that code review would be useful. For instance, in my last meeting with Louis, we noticed that Louis used the "direct capture" rather than "CROP-seq" option in CLEANSER for the Gasperini data. This is quite subtle and could easily be overlooked if one has less familiarity with these data.

**[2025-10-07 10:17] Timothy Barry** _(thread reply)_
This also seems reasonable to me. I think it would be useful at some point to benchmark pertpy on a CPU, but maybe this is lower priority.

**[2025-10-07 10:20] Eugene Katsevich** _(thread reply)_
Later we might or might not decide this is a meaningful comparison, but I’m not against throwing it in there if it doesn’t cost us too much. My stronger feeling is that we shouldn’t bend over backwards to get pertpy working well on a CPU.

**[2025-10-07 10:23] Timothy Barry** _(thread reply)_
&gt; My stronger feeling is that we shouldn’t bend over backwards to get pertpy working well on a CPU.

Well, I think this would be one of our selling points. (If you are working on a laptop or standard desktop, Pertpy will not get you very far. If, on the other hand, you are working on a cluster, the comparison becomes more interesting.)

**[2025-10-07 10:24] Timothy Barry** _(thread reply)_
But maybe the Pertpy developers would object that this is obvious.

**[2025-10-07 10:25] Eugene Katsevich** _(thread reply)_
&gt; If you are working on a laptop or standard desktop, Pertpy will not get you very far. If, on the other hand, you are working on a cluster, the comparison becomes more interesting.
I am not necessarily against delivering this message. I’m just saying that if pertpy takes 2.5 hours out of the box on a CPU, we shouldn’t bend over backwards to try to make it faster.

**[2025-10-07 10:26] Eugene Katsevich** _(thread reply)_
&gt; But maybe the Pertpy developers would object that this is obvious.
If it’s not in their paper (which I don’t think it is), then even if it is obvious to them, it may not be obvious to others.

**[2025-10-07 10:28] Timothy Barry** _(thread reply)_
&gt; That’s helpful context. Was this just a default setting that Louis didn’t touch or did Louis decide proactively which option to choose?

You have to choose one or the other. I agree that, moving forward, it would be good to flag required arguments unrelated to the gene expression matrix, gRNA expression matrix, or covariates.

**[2025-10-07 10:30] Eugene Katsevich** _(thread reply)_
I’m just not sure we need to make such decisions now. If we think there’s a chance we want to include pertpy in the final CPU comparison, let’s just throw it in there. We have lots of other things to do (e.g. we haven’t even started association testing benchmarking), so I want to make sure we do not let these “corner cases” stall our progress.

**[2025-10-07 10:32] Eugene Katsevich** _(thread reply)_
&gt; I agree that, moving forward, it would be good to flag required arguments unrelated to the gene expression matrix, gRNA expression matrix, or covariates.
Yes.

**[2025-10-07 10:34] Timothy Barry** _(thread reply)_
I mean, to me, assuming it's not onerous to do this, it seems worthwhile to throw Pertpy CPU into the mix, at least for the small dataset. I agree we should not be bending over backward to optimize either their CPU or GPU implementation.

**[2025-10-07 11:45] Louis Deutsch** _(thread reply)_
Gene, for context, all I’m thinking of is the python and R scripts in `guide-assignment-pipeline/bin` so it’s not actually that much. I don’t need any nextflow code review 

**[2025-10-07 11:50] Louis Deutsch** _(thread reply)_
And re: pertpy, got it. That all makes sense 

**[2025-10-07 13:00] Louis Deutsch** _(thread reply)_
No, I’m doing that today. Sceptre was so much faster that it finished while I was still watching the qsub log file to make sure things were working. Crispat and cleanser take significantly longer than sceptre, and pertpy even far longer than those. At least with the parameters as I have them now 

**[2025-10-08 01:20] Louis Deutsch**
For both medium datasets (2.5k gRNAs, 50k cells) pertpy is *still* running. It's been around 27 hours. Other people can't have performance this slow, right? Or no one would use this? I asked ChatGPT why other people might have better results than me, and this is what it says:

> • They run on much smaller subsets (the official demos subset to a few thousand cells and a handful of guides).
> • They're on faster GPUs (A100/H100), which cut compile time and execute faster.
> • Their data have fewer distinct shapes (e.g., many guides have similar support sizes), so cache hits are common.
> • Some use simpler assignment methods in pertpy (threshold / max-guide) for screening, and reserve the mixture for a small set.
How valid is the last one? Are people just not trying to use pertpy mixture models even for this medium size? It's all entirely due to the second line here
```pertpy_obj = pt.pp.GuideAssignment()
pertpy_obj.assign_mixture_model(adata, assigned_guides_key="assigned_guide")```
and I don't think there are subtle arguments I'm missing (docs <https://pertpy.readthedocs.io/en/latest/api/preprocessing/pertpy.preprocessing.GuideAssignment.html#pertpy.preprocessing.GuideAssignment.assign_mixture_model|here>). Just wanted to share this so y'all could see it before the meeting!
  - link: https://pertpy.readthedocs.io/en/latest/api/preprocessing/pertpy.preprocessing.GuideAssignment.html#pertpy.preprocessing.GuideAssignment.assign_mixture_model

**[2025-10-08 02:32] Louis Deutsch**
Figured it out. Even on the GPU, where I confirmed JAX was using the GPU, the singularity image was not using the GPU. So these have all been CPU runs. They give a warning that JAX is slow if not using CUDA, so I think this is valid data as far as benchmarking CPU methods, but it seems we shouldn't include pertpy as a CPU method.

**[2025-10-08 15:36] Louis Deutsch**
ok as a preliminary result, for the gasperini_small, which is 500 gRNA and 5000 cells (the .mtx takes up 127K on disc) the peak RAM usage was 1.7 GB. I didn't realize that nextflow tracing only is for each process, and doesn't give me insight into the processes as they run, so i need to add more tracing and rerun the bigger ones

**[2025-10-08 16:51] Eugene Katsevich** _(thread reply)_
&gt; the peak RAM usage was 1.7 GB
For which method? Is that CLEANSER?

**[2025-10-08 17:32] Louis Deutsch**
(doing this as a new post so it's easier to see) This is all the output that i got from the nextflow tracing that i enabled:
```$ cat trace.tsv
task_id process tag     cpus    time    rss     peak_rss        %cpu    %mem    read_bytes      write_bytes
1       CLEANSER_ASSIGN gasperini_small 1       -       3.4 GB  3.7 GB  186.5%  0.7%    221.5 MB        1.2 GB```
So it just traced at the level of the nextflow processes, but not within them. I was imagining a more nuanced time series, which chatgpt tells me I'll need to do a more refined thing to get, involving wrapping my bin scripts in functions that can track usage. I tried a quick version of that just to assess peak memory usage for cleanser, and this is the result:
``` Command being timed: "python /home/mnt/weka/jdeu/code/sceptre3-project-v2/benchmarking/guide-assignment/guide-assignment-pipeline/bin/run_cleanser.py cleanser/grna_matrix.mtx gasperini_small"
        User time (seconds): 646.76
        System time (seconds): 39.60
        Percent of CPU this job got: 186%
        Elapsed (wall clock) time (h:mm:ss or m:ss): 6:08.75
        Average shared text size (kbytes): 0
        Average unshared data size (kbytes): 0
        Average stack size (kbytes): 0
        Average total size (kbytes): 0
        Maximum resident set size (kbytes): 1778192
        Average resident set size (kbytes): 0```
plus some other stuff that doesn't seem important.

Is there a different way that i should be tracing?

**[2025-10-08 17:32] Louis Deutsch** _(thread reply)_
Yes this is all just for cleanser.

**[2025-10-08 18:42] Louis Deutsch**
that means i don't need to rerun the medium datasets. Here are the results for peak rss:

```           process                  tag  peak_rss
0   SCEPTRE_ASSIGN     gasperini_medium    293 MB
1   SCEPTRE_ASSIGN  replogle-rd7_medium  302.5 MB
2   CRISPAT_ASSIGN     gasperini_medium    1.4 GB
3   CRISPAT_ASSIGN  replogle-rd7_medium    1.4 GB
4  CLEANSER_ASSIGN     gasperini_medium   30.9 GB
5  CLEANSER_ASSIGN  replogle-rd7_medium   27.8 GB```

**[2025-10-08 18:43] Louis Deutsch**
this is for 2.5k grnas, 50k cells. (1) nice work sceptre! (2) that's pretty ridiculous for cleanser, right?

Does this seem like enough justification to email that person, Tim?

**[2025-10-08 21:15] Louis Deutsch**
We are using CLEANSER to do gRNA assignment on two datasets, the first one being a subset of the Gasperini data, and the second being a subset of the Replogle rd7 data. Both subsets contain the first 2500 gRNAs and the first 50,000 cells, and the mtx files are under 7M. Using essentially all default arguments, we end up needing around 30GB to run CLEANSER. This is a lot more than we expected, so we are wondering if we are using the method correctly! Our code is <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/guide-assignment-pipeline/bin/run_cleanser.py|here>.

-----

One question: is it ok to link to the code now? Is that like sending a work in progress preprint to someone? [edit: i just noticed it's a public repo so i guess that answers that!] I'm happy to do so if that makes sense. Otherwise, the full python script that runs cleanser is pretty easy to share:

```#!/usr/bin/env python3
import sys
import pandas as pd
import subprocess
import os

input_mtx = sys.argv[1]
dataset_id = sys.argv[2] # used to determine which dataset is present
output_dir = "cleanser_output/"
os.makedirs(output_dir, exist_ok=True)

# Choose flag based on dataset
# Replogle uses --dc (default cell), Gasperini uses --cs (cell-specific)
if "replogle" in dataset_id.lower():
    flag = "--dc"
elif "gasperini" in dataset_id.lower():
    flag = "--cs"
else:
    raise ValueError(f"Unknown dataset '{dataset_id}'. Expected 'replogle' or 'gasperini' in dataset name.")

# Run CLEANSER guide assignment
subprocess.run([
    "cleanser", "-i", input_mtx, "-o", f"{output_dir}/posteriors.csv", flag
], check=True)

# Process CLEANSER output to standardized format
df = pd.read_csv(f"{output_dir}/posteriors.csv", skiprows=1, sep='\t', 
                 names=['grna_id', 'cell_id', 'posterior'])

# Assign each cell to gRNA with highest posterior probability
assignments = df.loc[df.groupby('cell_id')['posterior'].idxmax()]

# Write standardized output
pd.DataFrame({
    'cell_id': assignments['cell_id'],
    'grna_id': assignments['grna_id']
}).to_csv("assignments_cleanser.csv", index=False)```
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/guide-assignment-pipeline/bin/run_cleanser.py

**[2025-10-08 23:20] Louis Deutsch**
Also I will be traveling tomorrow and only available on my phone (not laptop). I will see and respond to messages here. I am running the full gasperini and replogle-rd7 right now for just crispat and sceptre. If they finish overnight I can share those usage stats in the morning!

**[2025-10-14 12:59] Louis Deutsch**
Hi all! I am back online and everything. I just responded to the email with Alex, and I have crispat and sceptre running on the full gasperini and replogle-rd7 datasets now. I will keep you updated as that works. In the meantime, I'm starting to code up what we've talked about in terms of assessing the grna assignments, so we can talk through that later this week maybe

**[2025-10-16 12:08] Louis Deutsch**
here are the basic results from the full gasperini and replogle dataset for crispat and sceptre:

```          process           tag  peak_rss elapsed_time
0  SCEPTRE_ASSIGN  replogle-rd7    1.2 GB       1h 40m
1  SCEPTRE_ASSIGN     gasperini  924.5 MB       2h 14m
2  CRISPAT_ASSIGN  replogle-rd7    1.5 GB        3h 1m
3  CRISPAT_ASSIGN     gasperini    2.5 GB       10h 9m```
a good look for sceptre!

**[2025-10-16 17:11] Timothy Barry**
Hey Louis, to provide a bit more context for the simulation analysis, we will want to model the simulated data on the Gasperini screen. Suppose we have $n$ cells, and consider a given gRNA. Yet $Y_i$ be the UMI count of the gRNA in the $i$th cell. Let $b_i \in \{0,1\}$ be the batch of the $i$th cell. Finally, let $l_i$ be the gRNA sequencing depth (`n_grna_umis`) and $m_i$ the number of gRNAs with nonzero expression (`n_nonzero_grnas`) in the $i$th cell. We model $Y_i$ as a negative binomial model of the covariates:

$$ Y_i \sim NB_\theta(\mu_i); \log(\mu_i) = \beta_0 + \beta_1 b_i + \beta_2 \log(l_i) + \beta_3 \log(m_i).$$

We estimate $\theta$ and $\beta_0, ..., \beta_3$ via NB regression. (As a starting point, we can use `glm.nb` from the `MASS` package in R, although more sophisticated robust estimation procedures also are available.)

Finally, to simulate the data, we introduce a new variable $X_i$ indicating whether the perturbation is present or absent in cell $i$. We can simulate $X_i \sim Bern(0.02)$, where $X_i$ is independent of everything else. We simulate data from the following model. (Note that it is the same as the above model, but we have added $X_i$).

$$ Y_i \sim NB_\theta(\mu_i); \log(\mu_i) = \beta_0 + \gamma X_i + \beta_1 b_i + \beta_2 \log(l_i) + \beta_3 \log(m_i).$$

$\gamma$ is the log-fold change in the gRNA UMI count as a function of whether the cell contains the perturbation $X_i = 1$ or not $X_i = 0$. We will need to select this parameter ourselves. We may consider varying this parameter over a grid, but a reasonable starting point would be $\gamma = log(3),$ i.e., a fold change of 3. Finally, we can apply the methods to analyze the simulated data. ($X_i$ would be a latent variable in the simulations.)

**[2025-10-22 21:09] Louis Deutsch**
On `gasperini_small` (500 gRNAs, 5000 cells, and the mtx file is 127K on disc) my laptop also needs around 3.5GB to run CLEANSER. So this does not seem to be a HPC issue!

**[2025-10-23 16:19] Timothy Barry**
It was good. Louis shared the similarity results and MOI results for crispat and sceptre on the Gasperini and Replogle data. As next steps, we will:
• Finish obtaining the results for cleanser and pertpy on the Replogle/Gasperini data.
• Compute additional metrics, including (1) concordance with the maximum method assignments on the Replogle data, and (2) unimodality of the gRNA UMI counts within cells designated as perturbed/unperturbed.
• Create additional plots to compare assignments across methods, including upset plots and Venn diagrams.
• Benchmark the Nextflow version of sceptre.
• Potentially get started with the simulation study.

**[2025-10-23 16:27] Eugene Katsevich**
Thanks for the update! Where can I find the similarity and MOI results already obtained? 

**[2025-11-05 18:50] Louis Deutsch** _(thread reply)_
i'd like to meet anyway if you have time, to talk about the simulated datasets!

**[2025-11-06 16:43] Timothy Barry** _(thread reply)_
Louis is nearly done with getting Pertpy up-and-running. We had a discussion about generating simulated gRNA count data for benchmarking. We ironed out some kinks and it looks good. One issue is that, on the Gapserini data, batch seems to have a negligible impact on gRNA count, which obviates SCEPTRE's batch adjustment functionality. We therefore may want to consider "artificially" increasing the strength of the batch effect. Louis will apply the methods to analyze the simulated data in the context of the existing pipeline.

**[2025-11-06 16:45] Louis Deutsch** _(thread reply)_
The challenge for pertpy has been that each test run sits in the queue for hours since there aren’t a lot of available resources, and that’s even when I only request a tiny amount (like 1GB and 1hr run time). I also went back to conda so no more singularity for pertpy. That was too much of a complication. I’m at the point where every test run could be the one! So hopefully in the next day or so it’s working!

**[2025-11-10 12:49] Louis Deutsch** _(thread reply)_
I have pertpy running and using the GPU, but now I’m dealing with some secondary issues involving slowdowns with lots of small JIT compilations. I’m hoping to get it fixed this afternoon and then pertpy is in business! 

**[2025-11-10 13:59] Eugene Katsevich** _(thread reply)_
Remember that it’s not our job to optimize pertpy. If these secondary issues occur “naturally” then we should keep them as is, particularly if their impact is not enormous (e.g. if they slow things down by 25% relative to some “ideal” variant). If there’s something where we’re actively doing something suboptimal or running their software in a way that was not intended, especially if it’s 10x-ing the runtime, then yes this probably deserves our attention.

**[2025-11-10 14:19] Louis Deutsch** _(thread reply)_
Yes I hear that! Unfortunately this is bad enough that gasperini small takes well over two hours to finish, so I think I need to figure this out in order to be able to benchmark pertpy at all

**[2025-11-11 01:23] Louis Deutsch**
Here's my update on how pertpy is going. I am trying to get it to work using an interactive GPU node, with a conda environment I made for this test, using gasperini_medium. When I enable full logging, I see a constant stream of output like this:

```Working... ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   1% 4:52:32
WARNING:2025-11-11 01:18:42,672:jax._src.dispatch:198: Finished XLA compilation of jit(convert_element_type) in 0.009664774 sec
2025-11-11 01:18:42 | [WARNING] Finished XLA compilation of jit(convert_element_type) in 0.009664774 sec
WARNING:2025-11-11 01:18:42,675:jax._src.interpreters.pxla:1865: Compiling real with global shapes and types []. Argument mapping: ().
2025-11-11 01:18:42 | [WARNING] Compiling real with global shapes and types []. Argument mapping: ().
WARNING:2025-11-11 01:18:42,680:jax._src.dispatch:198: Finished jaxpr to MLIR module conversion jit(real) in 0.004164219 sec
2025-11-11 01:18:42 | [WARNING] Finished jaxpr to MLIR module conversion jit(real) in 0.004164219 sec
WARNING:2025-11-11 01:18:42,689:jax._src.dispatch:198: Finished XLA compilation of jit(real) in 0.008306026 sec
2025-11-11 01:18:42 | [WARNING] Finished XLA compilation of jit(real) in 0.008306026 sec
WARNING:2025-11-11 01:18:42,695:jax._src.interpreters.pxla:1865: Compiling _reduce_max with global shapes and types [ShapedArray(float32[199,2]), ShapedArray(float32[], weak_type=True)]. Argument mapping: (UnspecifiedValue, UnspecifiedValue).
2025-11-11 01:18:42 | [WARNING] Compiling _reduce_max with global shapes and types [ShapedArray(float32[199,2]), ShapedArray(float32[], weak_type=True)]. Argument mapping: (UnspecifiedValue, UnspecifiedValue).
WARNING:2025-11-11 01:18:42,702:jax._src.dispatch:198: Finished jaxpr to MLIR module conversion jit(_reduce_max) in 0.006636143 sec
2025-11-11 01:18:42 | [WARNING] Finished jaxpr to MLIR module conversion jit(_reduce_max) in 0.006636143 sec
Working... ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━   1% 4:52:32
WARNING:2025-11-11 01:18:42,771:jax._src.dispatch:198: Finished XLA compilation of jit(_reduce_max) in 0.067770720 sec
2025-11-11 01:18:42 | [WARNING] Finished XLA compilation of jit(_reduce_max) in 0.067770720 sec
WARNING:2025-11-11 01:18:42,776:jax._src.interpreters.pxla:1865: Compiling isfinite with global shapes and types [ShapedArray(float32[199])]. Argument mapping: (UnspecifiedValue,).
2025-11-11 01:18:42 | [WARNING] Compiling isfinite with global shapes and types [ShapedArray(float32[199])]. Argument mapping: (UnspecifiedValue,).
WARNING:2025-11-11 01:18:42,782:jax._src.dispatch:198: Finished jaxpr to MLIR module conversion jit(isfinite) in 0.005702019 sec
2025-11-11 01:18:42 | [WARNING] Finished jaxpr to MLIR module conversion jit(isfinite) in 0.005702019 sec
WARNING:2025-11-11 01:18:42,827:jax._src.dispatch:198: Finished XLA compilation of jit(isfinite) in 0.044776440 sec
2025-11-11 01:18:42 | [WARNING] Finished XLA compilation of jit(isfinite) in 0.044776440 sec
WARNING:2025-11-11 01:18:42,833:jax._src.interpreters.pxla:1865: Compiling _reduce_sum with global shapes and types [ShapedArray(float32[199,2])]. Argument mapping: (UnspecifiedValue,).
2025-11-11 01:18:42 | [WARNING] Compiling _reduce_sum with global shapes and types [ShapedArray(float32[199,2])]. Argument mapping: (UnspecifiedValue,).
WARNING:2025-11-11 01:18:42,838:jax._src.dispatch:198: Finished jaxpr to MLIR module conversion jit(_reduce_sum) in 0.005149841 sec
2025-11-11 01:18:42 | [WARNING] Finished jaxpr to MLIR module conversion jit(_reduce_sum) in 0.005149841 sec```
It got through about 7% of gasperini_medium in 20 minutes, so we're looking at a few hours for that dataset, and that's not even the full one. I am pretty stumped at this point. It seems like this might just actually be how it works...? Like maybe the GPU doesn't actually help that much for assignment? Even if it needs to recompile for every unique gRNA, there are only 2500 different gRNAs in this dataset. I don't understand why it's so slow, and if this is an us problem (like a HPC issue) or a them problem

**[2025-11-11 13:15] Louis Deutsch** _(thread reply)_
Ok sounds good! I do have this nagging concern that this is just a HPC issue, and that someone more familiar with GPUs, Wharton’s HPC, and JAX/JIT could immediately fix this, but I shall plug ahead as if this is just how pertpy works. I’ll run it on the full data now 

**[2025-11-11 21:11] Eugene Katsevich** _(thread reply)_
&gt; nagging concern that this is just a HPC issue, and that someone more familiar with GPUs, Wharton’s HPC, and JAX/JIT could immediately fix this
Well, what we could do is send them an email analogous to the one we sent the CLEANSER folks. This will cover our butts. Want to compose a draft?

**[2025-11-11 22:57] Louis Deutsch**
Here's my draft of the body of the email:

> I am trying to use pertpy to do gRNA assignment running on Wharton's high-performance cluster, but it's taking a very long time even on small datasets. I'm using `assign_mixture_model()`, and I'm doing this on a GPU node with JAX installed with CUDA support. But it's taking hours even for a test dataset with 2500 gRNAs and 50,000 cells. I think the issue is that JAX is constantly compiling, and that overhead is massively slowing down the analysis. I would be grateful for any guidance on how to speed this up!
Edit: I saved my steps pretty precisely for the tests i was doing yesterday on an interactive GPU node, to confirm this wasn't a nextflow problem. In case it's helpful to have those, this is every single thing I did to get this problem:

*Confirming there even is a GPU available*
```$ nvcc --version
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2024 NVIDIA Corporation
Built on Wed_Apr_17_19:19:55_PDT_2024
Cuda compilation tools, release 12.5, V12.5.40
Build cuda_12.5.r12.5/compiler.34177558_0
$ nvidia-smi
+-----------------------------------------------------------------------------------------+
| NVIDIA-SMI 555.42.02              Driver Version: 555.42.02      CUDA Version: 12.5     |
|-----------------------------------------+------------------------+----------------------+
| GPU  Name                 Persistence-M | Bus-Id          Disp.A | Volatile Uncorr. ECC |
| Fan  Temp   Perf          Pwr:Usage/Cap |           Memory-Usage | GPU-Util  Compute M. |
|                                         |                        |               MIG M. |
|=========================================+========================+======================|
|   0  NVIDIA L40S                    Off |   00000000:30:00.0 Off |                    0 |
| N/A   20C    P8             20W /  350W |       4MiB /  46068MiB |      0%      Default |
|                                         |                        |                  N/A |
+-----------------------------------------+------------------------+----------------------+

+-----------------------------------------------------------------------------------------+
| Processes:                                                                              |
|  GPU   GI   CI        PID   Type   Process name                              GPU Memory |
|        ID   ID                                                               Usage      |
|=========================================================================================|
|  No running processes found                                                             |
+-----------------------------------------------------------------------------------------+```
*Making my conda environment and setting shell variables*
```conda create -n pert-jax-test python=3.11 -y
conda activate pert-jax-test

python -m pip install -U pip setuptools wheel
pip install -U pertpy scanpy anndata h5py filelock

pip install -U "jax[cuda12]" # above output shows cuda12
unset LD_LIBRARY_PATH # prevents use of system CUDA which lacks cuSPARSE

export JAX_PLATFORM_NAME=gpu
export JAX_LOG_COMPILES=1 # give detailed logs```
*Actually doing gRNA assignment in python*
```import time
import sys

import jax, jax.numpy as jnp
import pertpy as pt
import scanpy as sc

# This is a subset of the (high MOI) Gasperini data, 
# with 2500 gRNAs and 50,000 cells
dataset_name = "grna_matrix.h5ad" 

# Checking that we see the GPU 
print("JAX backend:", jax.default_backend()) # shows `JAX backend: gpu`
print("JAX devices:", jax.devices()) # shows `JAX devices: [CudaDevice(id=0)]`
# quick warmup to force GPU work, just in case this helps
x = jnp.ones((4096, 4096), dtype=jnp.float32)
t0 = time.time()
(x @ x.T).block_until_ready()
print(f"Warmup matmul completed in {time.time() - t0:.2f}s")

# Loading the data and running gRNA assignment
adata = sc.read_h5ad(dataset_name)

ga = pt.pp.GuideAssignment()
print("Running mixture-model assignment...")
t0 = time.time()
ga.assign_mixture_model(
   adata,
   max_assignments_per_cell=100, # I tried this with 40 too
   show_progress=True
)
dt = time.time() - t0
print(f"Finished in {dt:.1f}s")```
Then this shows the many compilation notices, and the end result takes a long time.

**[2025-11-12 14:13] Eugene Katsevich** _(thread reply)_
Yeah! Definitely putting that on my CV, under “Personal endorsements:” Fabian Theis said “thanks from my side for trying pertpy.”

**[2025-11-12 15:32] Louis Deutsch** _(thread reply)_
ok i just posted the issue! I reran everything just now just to be *extra* sure that I'm using the right versions of everything and all that jazz. <https://github.com/scverse/pertpy/issues/881>
  - link: https://github.com/scverse/pertpy/issues/881

**[2025-11-14 16:52] Louis Deutsch**
Hey all, the simulations are going to need more work. Unfortunately the plot I showed yesterday is correct.  Both sceptre and crispat are making _far_ fewer assignments per gRNA than the truth, so the jaccard similarities are tiny because the cardinality of the union (the denominator term) is massive.

The first plot shows the number of assignments we make vs the truth, for each method. On average both sceptre and crispat make around 1000 assignments per gRNA, when actually most simulated gRNAs have over 4000 cells expressing them. Crispat also has that bimodality which is weird.

The second plot shows that, despite the tiny jaccard similarities, we at least do see an increase in similarity when the perturbation effect size increases.

**[2025-11-14 16:54] Eugene Katsevich**
So maybe you inadvertently chose the simulation parameters such that the gRNA assignment problem is very hard!

**[2025-11-14 16:55] Louis Deutsch**
Also does it seem like the simulated data has way too many cells expressing each gRNA? Or should that not matter and we should be able to detect that just fine, given enough signal?

**[2025-11-14 16:56] Louis Deutsch**
We have `Binom(num_cells_in_gasperini, 0.02)` cells expressing each simulated gRNA, so the mean of 4000 is from the ~200k cells in Gasperini x 0.02

**[2025-11-14 17:50] Eugene Katsevich**
&gt; Also does it seem like the simulated data has way too many cells expressing each gRNA? Or should that not matter and we should be able to detect that just fine, given enough signal?
Best to just get these numbers from the data. If the Gasperini data have an average of X cells per gRNA, then use X in your simulation.

**[2025-11-20 13:36] Louis Deutsch**
hi all, one thing i want to discuss today is the simulated data. I've been digging into it way more, and I think the issue is actually statistical.

When I model the actual target gRNA with a Poisson glm, and then also simulate the new fake gRNA from a Poisson, then we are able to perfectly recover the perturbation indicators provided the effect size is large enough (the perturbation probability matters too, but in my experiments the effect size mattered more).

When I instead use a NB regression to model the actual gRNA, I get $\hat\theta \approx 0.015$, so it is very overdispersed. I took a random sample of 15 gRNAs and all of them had \hat\theta in this range. I think this overdispersion leads to a blurring of the classes (perturbed or not) since the variance is absolutely massive, and the mixture model isn't able to separate them out. Even with massive effect sizes (like an effect size of 20) I was not able to get good assignment performance from sceptre on this simulated data with this \theta value.

If instead I use $c \times \hat\theta$ for the NB data I simulate, I am able to get more sensible outputs. The plot below shows the results when i use `c=100` (so I'm simulating with a \theta value of around 1.5). This is still with a NB regression and NB distribution for the new gRNAs, just with a decreased amount of overdispersion. Now we see that when the effect size is large enough (and the perturbation probability too), we get good agreement between the actual and estimated perturbation indicators, but for smaller effect sizes (with this `c`) we still aren't able to recover them. I tried it with `c=10` and that was not enough to get good results.

This data is no longer very realistic though -- only 92% of the entries are 0, vs ~99.5% in the actual Gasperini data.

**[2025-11-24 14:52] Louis Deutsch**
hi all, I would like to just do a slack update for this week if that works, as I'll be traveling off and of throughout the week.

I am continuing to get the association nextflow pipeline going, and to modify the simulation for using rd7. The main challenge here is that we are now simulating low MOI data, so I can't just do iid bernoulli perturbation indicators.

Here is the procedure I'm trying right now. I'll let you know the results as I analyze them.

First, we determine which cell-gRNA combos are expressed.
• for each cell i, sample iid `K_i ~ Pois(0.6) + 1`. This is the number of gRNAs each cell expresses, and is chosen so that the mean of 1.6 matches the empirical MOI we observed, while also ensuring that every cell expresses at least one gRNA.
• now for each cell i where K_ i > 0 (which with the above truncation will be all of them) we sample without replacement the K_i indices that are the expressed gRNAs.
    ◦ optional improvement: sample the K_i expressed gRNAs with probability proportional to the overall activity of that gRNA
This gives us a large binary matrix X with the true perturbation statuses for each cell-gRNA pair.

Next, we need the UMI counts. For this, I am trying the following:
• pick a target gRNA from rd7, call it `y`. Do a Poisson regression of `y ~ batch + log(grna_n_nonzero+1) + log(grna_n_umis+1)` to get covariate effects. Let eta_i be the fitted log means for each cell
•  For cell i and each gRNA g, sample Y_{ig} ~ Pois(lambda_{ig}) where \lambda_{ig} = \exp(\eta_i + \gamma X_{ig}), so \gamma is the effect of a perturbation
    ◦ optional improvements: add an overdispersion to this
So the main change is that now we use a shifted Pois distribution to determine how many gRNAs each cell gets, and I'm using a Poisson regression to estimate the covariates.

**[2025-11-25 10:09] Eugene Katsevich** _(thread reply)_
Thanks for the update! A couple follow-up questions:
1. Before proceeding with rd7, did you verify that the gRNA assignment on those data is easier? Did you see some bimodality? On a similar note, after obtaining the data-generating model, you might want to eyeball what gRNA counts generated from that model look like, to verify that there’s some bimodality we can hang our hat on. 
2. You are including batch in your Poisson regression for the UMI counts. Does that cause any issues because batch has so many levels?

**[2025-11-25 14:51] Louis Deutsch**
Hey, I'm looking at Replogle and I do see some nice bimodal behavior for the UMI counts of the gRNA I'm picking. The first plot here shows that.

I am able to see this bimodality without using SCEPTRE if I condition on the UMI counts being positive. The 0s overwhelm the picture otherwise. However, when I use this same thing of filtering out the zeros, then I also see clear bimodality for the gRNA I used with Gasperini! That's plot 2 here. So it seems that there was some bimodality in Gasperini all along... regardless, I'll keep going ahead with trying to simulate from RD7 instead of Gasperini

**[2025-11-25 15:00] Louis Deutsch** _(thread reply)_
Re point #2, it takes ~12 Fisher scoring iterations to converge, but the model itself seems fine as far as a basic skim of `summary()` indicates. Every level of batch has at least 48 cells with non-zero UMI counts, and some levels have &gt; 1000 cells. I will also try without Batch as a sanity check though!

**[2025-11-25 15:02] Louis Deutsch** _(thread reply)_
at the same time, we did have actual bimodality from the chosen gRNA for Gasperini. Does this change anything?

**[2025-11-25 19:48] Eugene Katsevich** _(thread reply)_
I want to point out that we have two kinds of plots here:
1. Histograms of gRNA counts for cells assigned to have gRNA and not have gRNA.
2. Histogram of gRNA counts for all cells.
The fact that the first plots suggests different distributions for both groups of cells may be misleading: it might just be the algorithm making an arbitrary distinction (analogous to a clustering algorithm obediently giving you the number of clusters you asked for regardless of how many are present in the data). The second type of visualization is more suggestive of bimodal structure.

**[2025-12-02 03:13] Louis Deutsch**
hi all, I'm still making sense of the problem but here is my thinking so far.

For the first plot, I simulated data as we have been doing:
1. fit a NB GLM with `y ~ 1 + log(grna_n_umis + 1) + log(grna_n_nonzero + 1)` (both are significant, and theta is approx 0.004). Let `eta` be the linear predictor here
2. simulate `X = rbinom(num_cells, 1, pert_prob)` with `pert_prob = 0.02`
3. draw new NB data using a log mean of `eta + gamma * X`
The first plot has the marginal distribution of the simulated data and the real gRNA in the upper row, and the lower row has the simulated data conditioned on perturbation status. This is for a large effect size to exaggerate the effect. For the non-perturbed cells it's roughly the right shape, with the mean and variance matching the real gRNA (marginally). But the perturbed cells are way higher variance due to their increased mean, and so overall the marginal distributions are very dissimilar (the range is orders of magnitude different). The perturbed and non-perturbed distributions are also sort of "proportional" so they can't be separated very well.

The next plot repeats this except with a Poisson GLM and sampling new data from a Poisson distribution. Now we are able to actually get two components, although the variability and range are not quite big enough in this case.

I think the overdispersion is the problem, and maybe more generally a lack of fit for the model. For the NB data, we increase the mean of the perturbed cells, and with our huge amount of overdispersion, we get a massive variance which prevents the formation of a component concentrated on larger UMI counts. Even the variance of the unperturbed cells is ~30% greater than the marginal variance of the actual gRNA UMI counts. The last plot gives a sense of the poor fit: this is the actual UMI counts vs the fitted values from the NB GLM [for a random sample of 5000 observations], and they are basically uncorrelated... I'm not sure how bad this actually is, since we know that this model with just two covariates won't explain that much of the data, but it doesn't seem good. Regardless, though, I think the overdispersion is really the killer. If we had less overdispersion then we'd get a more Poisson-like situation where the perturbed component could be concentrated on larger values without being so massively high variance.

So, to summarize so far, I think the NB GLM fits very poorly, leading to a huge overdispersion, which then causes each component of the simulated mixture to be really high variance and therefore effectively inseparable.

**[2025-12-02 09:09] Eugene Katsevich** _(thread reply)_
Thanks for illuminating this vexing issue, <@U0524GR916C>. Is there any value of the dispersion that gets the marginal distributions of the real and simulated gRNA counts to look roughly similar? If so, choosing this value would be a quick-and-dirty solution.

**[2025-12-03 00:43] Louis Deutsch**
I modified the mean a little and now I'm able to get much more realistic shapes! I introduced a new parameter, so now I'm doing this:
```log_mu_pert <- log_mu + pert_effect * X - non_pert_effect * (1-X)```
The extra control over each component is helping. This one does a decent job of matching the range; the sparsity; and the bimodality

This one has `pert_effect = 6` and `non_pert_effect = 4`. I'm also using a much smaller `pert_prob` (the probability a cell expresses this gRNA) since this is Replogle, so with low MOI I figure that I can still do iid Bernoulli indicators as I'm just modeling a single gRNA, but even fewer cells will likely be expressing it. I'm doing `pert_prob = 0.001` here, so we expect ~620 cells to express this gRNA. Does that seem realistic?

If we want to make this seem less hacky, I think I can frame this as the result of some optimization problem since that's intuitively what I was doing anyway.

**[2025-12-03 21:05] Eugene Katsevich** _(thread reply)_
&gt; log_mu_pert &lt;- log_mu + pert_effect * X - non_pert_effect * (1-X)
Can you help me understand what this is doing? What is log_mu_pert? What is log_mu? What relationship does this have to your previous model, `y ~ 1 + log(grna_n_umis + 1) + log(grna_n_nonzero + 1)`?

**[2025-12-03 21:10] Louis Deutsch** _(thread reply)_
Ah yes sorry, I'm using `log_mu` for the linear predictor from the regression, so that's the log mean as modeled by the covariates and isn't affected by the perturbation indicator. Then `log_mu_pert` is the log mean of the simulated gRNAs

**[2025-12-04 17:21] Louis Deutsch**
Hi Gene, we discussed the simulation and agree that it's good enough to move on to other things. I'm going to write up the procedure as it stands, add it to the pipeline, and then move on fully to association. 

Also Lukas responded and said he should be able to get to the pertpy GPU issue in two weeks. 

**[2025-12-04 17:32] Eugene Katsevich** _(thread reply)_
Thanks for the update!
&gt; we discussed the simulation and agree that it’s good enough to move on to other things. I’m going to write up the procedure as it stands, add it to the pipeline, and then move on fully to association.
Sounds good! I think last time we spoke, we agreed you would create a writeup with the simulation results. Is that still on the agenda?
&gt; Also Lukas responded and said he should be able to get to the pertpy GPU issue in two weeks.
:+1::skin-tone-3:

**[2025-12-04 17:33] Louis Deutsch** _(thread reply)_
Yes, now that I've got the simulation somewhere that I’m happy with, my next TODO is that writeup 

**[2025-12-05 12:38] Louis Deutsch**
Lukas just posted an update! It seems he's concerned that MCMC in general won't be fast here. This seems like good news for us!

 <https://github.com/scverse/pertpy/issues/881|https://github.com/scverse/pertpy/issues/881>
  - link: https://github.com/scverse/pertpy/issues/881

**[2025-12-12 01:09] Louis Deutsch**
I just simulated 25 iid copies of the gRNA in <https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1764740580888479|this> post (except with theta=4.5). The Jaccard similarities of Sceptre's assignments with the truth are all above 96%. These are indeed very recoverable!
  - link: https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1764740580888479

**[2025-12-12 01:14] Louis Deutsch**
For the simulations, do we want to consider a variety of parameter settings? Or is one good exemplar like this probably fine (for now)?

**[2025-12-12 08:54] Eugene Katsevich** _(thread reply)_
In an ideal world, we would vary the degree of separation between the two classes over some range. To assess the priority of this task, I would want to see a report of our benchmarking of gRNA assignment methods, including this simulation, to gauge what overall picture is emerging and whether this picture needs augmenting.

**[2025-12-15 14:37] Louis Deutsch**
With GPU resources Pertpy is not able to finish running assignment on the full Replogle before hitting the 4hr limit for GPU resources on HPC and getting killed. Gene, it seems that we could move to running this on that other cluster if we need more than 4 hours, but that seems really annoying to do... what should I do? Run pertpy on a subset with GPU, or switch back to CPU?

**[2025-12-15 14:50] Louis Deutsch**
Also here are the simulation results. First is the histogram of Jaccard similarities between the assignments for each simulated gRNA and the ground truth (the true perturbation indicators). Sceptre and crispat both basically recover the truth, with sceptre slightly better. Pertpy fully recovers some, but leaves many cells unassigned and totally fails on other gRNAs. This is weird as each of the 25 simulated gRNAs are iid. And CLEANSER basically totally fails.

Next, here are the time + space usage results. CLEANSER used all available memory (240GB, the max i can give a single job) but this is the same performance as when i ran it on a subset where it didn't use all available memory, so i think this is a real result and not just a failure due to memory overflow.

Finally, the last plot gives a little more insight into why CLEANSER's jaccard similarities are so low. This plot compares the number of cells assigned to each gRNA. sceptre and crispat are in near total agreement, pertpy is pretty random, and cleanser assigned waaay too many cells to each gRNA

**[2025-12-15 14:54] Louis Deutsch**
On the reduced dataset (100k cells instead of the full 600k for Replogle, with still 25 iid simulated gRNAs) here are the results, just for completeness, and to show that cleanser wasn't failing just due to memory issues:

**[2025-12-15 14:55] Louis Deutsch**
pertpy is also running with a GPU for all of this, so pertpy + gpu is still waaay slower than sceptre and crispat.

**[2025-12-15 15:01] Louis Deutsch**
i have one lingering question with these simulations. I have a simulated dataset that is `(num sims) x (num cells)`. Each row is an iid simulation and should be treated independently. Could there be issues where these methods fail because this isn't how real data behaves, so like if a cell is assigned to one simulated gRNA it won't be assigned to another? Since replogle is low MOI after all.

I don't _think_ this is what's causing pertpy to do badly, since there isn't even an MOI argument, they say it works for both, and each gRNA seems to be modeled independently, but I want to make sure i've thought through this before trumpeting these results!

For reference, this is basically the entirety of the pertpy code:
```import sys
import pandas as pd
import anndata as ad
import pertpy as pt

input_h5ad = sys.argv[1]

# Load data and run pertpy guide assignment
adata = ad.read_h5ad(input_h5ad)
pertpy_obj = pt.pp.GuideAssignment()
pertpy_obj.assign_mixture_model(adata, assigned_guides_key="assigned_guide")```
so there really isn't much wiggle room in that

Edit: ok yes actually there is an argument in there that I overlooked! `max_assignments_per_cell: int = 5` !! I am going to rerun the pertpy simulations with this fixed and post the updates!

Edit 2, some time later: `max_assignments_per_cell` has no effect on this. Assignment is independent for each gRNA in pertpy, and `max_assignments_per_cell` only matters so far as if a cell has more than this many guides assigned to it, in the output it is simply labeled as 'multiple' rather than with each gRNA it actually expresses. So not actually a very interesting parameter.

**[2025-12-15 16:56] Louis Deutsch**
For cleanser with these simulations, would `direct-capture` or `crop-seq` make more sense?

**[2025-12-15 17:57] Timothy Barry**
Do you mean for CLEANSER?

**[2025-12-15 17:57] Louis Deutsch** _(thread reply)_
yeah, i just want to make sure that i actually am doing a fair assessment for these simulations since it's not actually a real gRNA matrix!

**[2025-12-15 17:58] Timothy Barry** _(thread reply)_
Well, our simulated data are modeled on the Repologle data, right? I believe Replogle used direct capture, so maybe we should use that.

**[2025-12-15 18:24] Louis Deutsch**
I'm struggling to get pertpy to work for the simulations... the simulated data has 25 fake gRNAs, each of which is an iid draw from the same process, so i want these to be handled independently for assignment. That means a cell should be able to be assigned to up to 25 gRNAs, otherwise these will compete and mess up the results.

When i run pertpy with the default of `assign_mixture_model(adata, assigned_guides_key="assigned_guide",  max_assignments_per_cell = 5)` then i get the results i posted earlier, where we do a decent-to-great job recovering the ground truth for about 3/4 of the simulated gRNAs.

But now i just finished a run where i used `assign_mixture_model(adata, assigned_guides_key="assigned_guide",  max_assignments_per_cell = 25)` . I thought this might fix things, but instead only ~1.7% of cells get assigned to any gRNAs at all. This seems like a miscalibration issue for the underlying probabilities? I'm rerunning it now so i can save the info needed to compute these, since it only gives me hard assignments, but I'm posting this just in case anyone has any insights in the meantime!

(Edit for posterity: this was a caching issue. I have fixed it)

**[2025-12-15 18:32] Timothy Barry** _(thread reply)_
Does pertpy analyze the gRNAs jointly?

**[2025-12-15 18:45] Louis Deutsch** _(thread reply)_
I'm pretty confident it does each gRNa independently, and this `max_assignments_per_cell` argument is only used at the end. The code that makes me this confident is this (starts at line 246 <https://github.com/scverse/pertpy/blob/main/pertpy/preprocessing/_guide_rna.py|here>). Paraphrasing, it does this:

``` for gene in adata.var_names:
            data = adata[is_nonzero, gene]
            data = np.log2(data)
            assignments = mixture_model.run_model(data)```
i've also used Claude to go through this to double-check my sense
  - link: https://github.com/scverse/pertpy/blob/main/pertpy/preprocessing/_guide_rna.py

**[2025-12-15 18:46] Timothy Barry** _(thread reply)_
OK. In that case, maybe we setting `max_assignments_per_cell` to 25 makes sense to me. Have you taken a look at the docs for this parameter?

**[2025-12-15 18:47] Louis Deutsch** _(thread reply)_
yeah, all it says is "The maximum number of gRNAs that can be assigned to a cell." Which is why i thought setting it to 25 should solve my problem!

**[2025-12-15 18:51] Louis Deutsch** _(thread reply)_
i'll check my new rerun of pertpy in case this is just randomness in the method itself! maybe the MCMC gets weird with this simulated data or something

**[2025-12-15 19:51] Louis Deutsch** _(thread reply)_
ah also i better understand `max_assignments_per_cell` now. It doesn't induce some dependency. All it does is that if a cell exceeds this many gRNAs being assigned to it, that cell is just tagged with "multiple" instead of having the specific gRNAs that it expresses. So that shouldn't actually be changing the behavior like this.

I'm now investigating if it's an issue with caching

**[2025-12-15 23:40] Louis Deutsch** _(thread reply)_
ok i have confirmed that this was indeed a caching issue! There was some sneaky stuff happening. The reason `max_assignments_per_cell = 25` led to that major change was because of these caching issues and not any actual statistical reason

**[2025-12-16 04:54] Louis Deutsch**
ok i have verified everything with the simulations so i am sharing this again. The first plot gives the similarities of each simulated gRNA with the ground truth perturbation status. sceptre is the best, crispat is close behind, pertpy recovers some gRNAs well but fails on others, and cleanser generally is bad.

-------------------------------------------

*Pertpy*
Pertpy does a great job on around half the gRNAs, but it also has 0 cells assigned to 6 of them, and then massively overassigns for a few others. Those are the ones with low similarities. This variability suggests an MCMC issue to me. The perturbed component is modeled as log(y) ~ N(mu, sigma^2), and by default mu ~ N(3, 2^2). The actual mean of the perturbed simulated data is around 1000, so these are pretty far apart. The default fraction perturbed is 15%, when it's only 0.1% for the simulated data. So that's another big difference. By default it's a 50 sample burn-in and 100 samples kept, for a 4 parameter model. Maybe it's not mixing well enough and that's why we get such inconsistent results here?

I reran it with mu ~ N(10, 3^2) and `fraction_positive_expected = 0.01`, hoping that would fit better, but that didn't actually improve things at all. I'm trying with a bigger burn-in and number of samples, just to see if that makes the difference.

something Claude pointed out -- there is a default argument `seed = 0` buried in an internal function, so we can't change it, which causes every gRNA's MCMC to be run with the same seed. Maybe that contributes to it being a coin toss if the parameters settle somewhere good or wander off?


*Cleanser*
It's also weird to me that cleanser is so bad. I believe they use a 2-component NegBin mixture for direct capture (which is what I'm doing). The data exactly are a NegBin mixture, but I guess it doesn't work well. As before, could this be due to MCMC issues / variability? CLEANSER does a warmup of 300 samples and then keeps 1000 samples. That's a lot more, and the model seems comparable in complexity to pertpy's, so I'm surprised that CLEANSER does so bad here. So maybe it's not so simple as just sampling more?

-----------------------------------

I have verified in the code that both pertpy and cleanser are indeed analyzing each gRNA independently, so that is not the concern. I believe these are genuine results at this point.

**[2025-12-16 04:59] Louis Deutsch**
So is this good enough to wrap up the simulations for now? I wanted to be more confident that a failure of a method on the simulations didn't mean that the simulations were terrible, rather than the method itself being suboptimal. Are we at the point where we can actually attribute the failings of pertpy and cleanser on this dataset to those methods, vs "user error" on my part?

**[2025-12-16 05:07] Louis Deutsch**
This is a repeat of the plot in the previous post, with the only change being pertpy. For this run i did this:
```pertpy_obj.assign_mixture_model(
        ...,
        gaussian_mean_prior=(10,3), 
        poisson_rate_prior=1,
        fraction_positive_expected=0.01,
        num_warmup = 200,
        num_samples = 400
)```
so the prior fits better (although that alone was not enough to improve performance earlier), and i increased the warmup and number of samples. This made a big improvement! I think it was indeed just not mixing enough. But now it's slower. I don't have an accurate run-time for this one since i do some persistent caching, but I can rerun it with no cache tomorrow to see how much longer it takes to get this better performance.

**[2025-12-16 10:02] Timothy Barry** _(thread reply)_
OK, nice. So we can just go with `max_assignments_per_cell = 25`?

**[2025-12-16 10:15] Timothy Barry** _(thread reply)_
Nice. Indeed, it's better, but still not as good as crispat or sceptre. I do think it would be good to present the runtimes alongside the results. CLEANER still lags behind. My view is that these simulation results are probably good for now, and that we could move onto other tasks (e.g., association testing).

**[2025-12-16 12:38] Louis Deutsch** _(thread reply)_
Yeah, and if we always want to know every single gRNA that got assigned to a given cell, we can just set this to be massive. This has no effect on the actual assignment process. Purely post-hoc processing 

**[2025-12-16 12:50] Louis Deutsch** _(thread reply)_
Ok sounds good! I’ll get those run times but I’ll otherwise move on. 

Part of the reason I was doing this is I wanted to get extra confirmation that I'm doing the real data analysis with good parameter choices, so that sceptre’s victory is clear. In light of this, do you think I should use these settings for real data? Or just do defaults for now and we can make this improvements to pertpy later?

**[2025-12-16 14:34] Timothy Barry** _(thread reply)_
> In light of this, do you think I should use these settings for real data?
Do you mean the pertpy hyperparameters? Yes, that seems reasonable to me. Ultimately, as Gene has said, it's not our job to optimize their hyperparameters. I think my preference would be to use the default hyperparameters if possible, and if not, to use the parameters that work reasonably well on the simulated data.

**[2025-12-16 20:13] Louis Deutsch**
Pertpy is not able to finish running assignment on replogle in the 4hr time limit for GPUs. I've tried modifying the things in my control (turning off all logging and caching [for future runs] to reduce I/O time, etc) but it's not enough. What should I do?
• run on a subset: this is unsatisfying to me because then our benchmarking results are limited by wharton HPC limits, not pertpy itself
• analyze on that other cluster that the HPC people told us about. This sounds time-consuming so I hope this isn't the answer.
• run on CPU: not ideal because of their statements to run this on a GPU
• My vote: chunk replogle by gRNA and run on pieces. I could treat each chunk as a different dataset, so this would be easy to do with the pipeline (in parallel over chunks, too). Measuring time and space is more complicated now, but still reasonable 
for now, I just used replogle-rd7_medium as a placeholder and I'm doing association.

**[2025-12-23 13:21] Louis Deutsch**
I am delighted to report that the new GPU queue `gpuv100.q` has allowed me to run pertpy on replogle for assignment

**[2026-01-08 13:06] Louis Deutsch**
one thing i'd like to get guidance on today is what exactly the default pertpy differential expression testing approach actually is.

In the tutorial <https://github.com/scverse/pertpy-tutorials/blob/86069936d9d811e23558ddefb462ed78ba3eaaa9/differential_gene_expression.ipynb|here>, they say this:
&gt; Pertpy provides an API to access several types of models for differential expression analysis. The first group of models comprises the T-test and Wilcoxon test as simple statistical tests for comparing expression values between two groups without covariates. The second group includes models of the linear model family that allow modeling complex designs and contrasts. Currently included are <https://academic.oup.com/bioinformatics/article/39/9/btad547/7260507|PyDESeq2>, <https://academic.oup.com/bioinformatics/article/26/1/139/182458|edgeR> as well as a wrapper around statsmodels <https://www.statsmodels.org/|Statsmodels>. which provides access to a wide range of regression models, including ordinary least squares regression, robust linear models and generalized linear models.
I picked edgeR since that's what they do in the tutorial, but I want to confirm if this is what i should be doing!
  - link: https://github.com/scverse/pertpy-tutorials/blob/86069936d9d811e23558ddefb462ed78ba3eaaa9/differential_gene_expression.ipynb
  - link: https://academic.oup.com/bioinformatics/article/39/9/btad547/7260507
  - link: https://academic.oup.com/bioinformatics/article/26/1/139/182458
  - link: https://www.statsmodels.org/

**[2026-01-16 17:20] Louis Deutsch**
Hey all, quick check on the covariates I'm using for running sceptre on a subset of replogle for positive controls:
• i am using the `response_n_nonzero` and `response_n_umis` that sceptre computes from the provided response_matrix subset
• i am computing `grna_n_umis` and `grna_n_nonzero` from the corresponding subset of the original grna_matrix, and passing those in as extra_covariates
• no covariates are used from the grna_matrix that i provide to sceptre, since that's a dummy indicator one
does that seem reasonable for now?

**[2026-01-16 17:23] Louis Deutsch** _(thread reply)_
yeah, except the only complexity is that the gRNA matrix features we use aren't computable from the provided grna_matrix, so they are computed elsewhere and passed in via extra_covariates

**[2026-01-16 17:29] Eugene Katsevich** _(thread reply)_
I thought the plan was to create a sceptre object as usual, pass it through assign_grnas() as usual, then feed that as input into the sceptre portion of the benchmarking pipeline, which just continues by running the discovery analysis.

**[2026-01-16 17:30] Eugene Katsevich** _(thread reply)_
It sounds like what you’re doing instead is extracting the gRNA assignment matrix, and passing that directly into the beginning of sceptre rather than starting with a sceptre object.

**[2026-01-16 17:34] Louis Deutsch** _(thread reply)_
For the positive control comparison, I thought we wanted to have every method run on the same gRNA assignments (edit: and the same data overall)? My approach right now is this:
• determine gRNA assignments in advance using `sceptre-pipeline` on the full dataset. 
• Make a smaller `response_matrix` with all positive control cells, so all these methods can run in-memory, using those pre-computed assignments to know which cells are positive controls for which genes
• Run each method on this subsetted data, using the same gRNA assignments in every case. So for sceptre, the `grna_matrix` argument here is a dummy one recreating the pre-computed gRNA assignments for just the subset of cells currently used. 

**[2026-01-19 04:33] Louis Deutsch**
i've collected the grna assignment results <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/assignment-analysis/assignment-outline.pdf|here>

I think pertpy's long runtime here might be partly due to some GPU settings, so I'm rerunning to confirm. I'll update as soon as they're done
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/assignment-analysis/assignment-outline.pdf

**[2026-01-19 14:31] Louis Deutsch** _(thread reply)_
Cleanser also performs fine (minus the memory and time issues) on real data, so I am investigating its performance on the simulated data 

**[2026-01-20 18:05] Louis Deutsch**
hi all, I've been dotting the i's and crossing the t's for pertpy assignment as I'm doing the writeup, and I want to get some confirmation on what I plan to say about it.

In brief, before I go publicly criticizing it for being glacially slow, I really want to make sure I haven't missed something as silly as `export USE_GPU=true`

Does this seem like enough for now?
• we've made an issue that the authors have looked at
• my code confirms the presence of a GPU
• i have experimented with settings, including caching and memory preallocation, to make sure there isn't a simple cluster setting that makes the difference (for example, trying float32 instead of float64 to speed up XLA since that's known to be a possible problem point)

At this point, I've done all the things I can think of, but it's also my first time ever using a GPU which is why I'm nervous about missing something! I am not spending more time on this now, I just did this last round of double-checking as I make the results writeup.

For context, here is the time and space usage on all medium datasets. Pertpy's not doing great!

**[2026-01-21 02:57] Louis Deutsch**
i just noticed on the <https://github.com/douglasyao/FR-Perturb|main page> for FR-Perturb on github it says this:
> NOTE: This repository is no longer actively maintained. I’m not able to provide updates or bug fixes at this time, but you’re welcome to fork the code and modify it for your own use. Thanks for your understanding!
Does this matter for us?
  - link: https://github.com/douglasyao/FR-Perturb

**[2026-01-21 23:39] Louis Deutsch**
Hi, I've added FR-Perturb to the positive control pipeline, and here are the results on my downsampled low MOI replogle-rd7 association dataset. This dataset has ~50k cells, each of which expresses exactly one of ~80 on-target genes or is NT (with these assignments being computed using `sceptre-pipeline` on the full replogle-rd7 dataset).

The first plot is on-target p-values, the next is on-target log2 FC. Note that the x and y axis scales are free.

They get their p-values by resampling, with a default of B=1000, hence all the p-values of 1/1001. I provided `grna_n_umis` and `grna_n_nonzero` to FR-Perturb as covariates to control for

**[2026-01-22 03:14] Louis Deutsch**
another minor update for assignment: CLEANSER totally fails on the downsampled simulated dataset (basically zero agreement between the ground truth and estimates). I wanted to make sure that this didn't indicate a problem either with the simulated data or how I'm using CLEANSER. I just ran CLEANSER on the full simulated data (except I had to do a loop over guides separate from nextflow, because it does run out of memory if I present it with the whole simulated grna matrix). CLEANSER does indeed recover the ground truth decently well on the bigger dataset, so that previous failure is not a problem on our end.

The issue is that the downsampled simulated dataset has too few perturbed cells for cleanser to get meaningful signal. The other 3 methods could still handle it, but not cleanser. So that's another negative for it. Missing an infrequently perturbed guide entirely is not great. And even with more perturbed cells, as in the plot, some of the simulated guides still totally fail

**[2026-01-22 13:05] Eugene Katsevich**
This may not be our top priority but I think that capturing _scaling_ with some basic problem parameters like #cells or #gRNAs or #pairs for each of the in-memory methods and then with our pipeline would be very striking. Then we could say something like “sceptre package is best among in-memory methods but if you look at all these in-memory methods as a group they don’t scale to large problem sizes but pipeline does and look how much better it is and how big of a problem size we were able to scale it to.” I attach an inspiration from pertpy and a simple mock-up of what our plot may look like.

**[2026-01-27 01:16] Louis Deutsch**
I've made some changes to the simulations. I'll add this to the writeup tomorrow in more detail, but here is the gist. My goal: to get good simulations with 100k cells instead of the 600k in Replogle, and to vary the fraction of perturbed cells while still keeping a reasonable enough UMI distribution.

As a reminder, i pick a highly expressed guide from Replogle, and do a NB regression of it on `grna_n_nonzero` and `grna_n_umis` to estimate how those covariates affect cell UMI counts.

I had been sampling the true perturbation indicators iid from a Bernoulli. But I found it hard to get enough control over the results due to the variability in the NB fitted means, especially when increasing the perturbation rate while downsampling cells.

My new approach: a cell's perturbation probability does depend on those covariates, at least a little, so I now do a logistic regression of the sceptre-pipeline's assignments on  `grna_n_nonzero` and `grna_n_umis` . Then I sample perturbation indicators using the estimated probabilities, with an intercept adjustment to change the overall perturbation rate.

The first plot shows a simulated UMI vector using this method (gamma1 - gamma0 is the mean difference on the log scale for perturbed vs unperturbed), and the second plot is each method's performance. This simulated guide is expressed in around 1.5% of cells.

**[2026-01-27 10:37] Eugene Katsevich** _(thread reply)_
I do want to stress that your main focus should be getting a first “complete” set of results for both assignment and association. There’s still quite a decent chunk of work left on the association front, so I would prioritize that.

**[2026-01-27 13:40] Eugene Katsevich**
Hi folks, there is a new <https://www.biorxiv.org/content/10.64898/2026.01.22.701179v1|preprint> out about guide assignment, which includes a new method and benchmarking study. The method is too new for us to compare against, but we may be interested in both the setup and results of their benchmarking study, which included sceptre, CLEANER, crispat, and a few others methods.
  - link: https://www.biorxiv.org/content/10.64898/2026.01.22.701179v1

**[2026-01-27 13:47] Louis Deutsch** _(thread reply)_
Ok yes I hear that! This was my last todo for assignment until association is basically done, so full steam ahead there now 

**[2026-02-04 23:03] Louis Deutsch**
Hey, i just clicked a link in the sceptre manual and got a 404 error.

&gt; This <http://rdrr.io|rdrr.io> page can’t be found
&gt; No webpage was found for the web address: *<https://rdrr.io/pkg/sceptre/man/run_discovery_analysis.html>*
&gt; HTTP ERROR 404
it's from clicking the `run_discovery_analysis()` link on this page <https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html>. Not urgent at all, just FYI
  - link: http://rdrr.io
  - link: https://rdrr.io/pkg/sceptre/man/run_discovery_analysis.html
  - link: https://timothy-barry.github.io/sceptre-book/run-power-check-and-discovery-analysis.html

**[2026-02-07 22:22] Louis Deutsch**
Hi all, I'm digging into calibration on Replogle and I'd appreciate thoughts on this! First I'll detail the exact data I'm using, then I'll show what I've found.

In brief: with no QC i have clear miscalibration, but even just setting `n_nonzero_trt_thresh = 1` is enough to restore calibration!

The data:
• gene selection: I pick 150 genes at random
• cell selection: I take all cells that: express no targeting guides AND express at least one NT guide AND have at least one non-zero expression for the chosen response genes. I end up with ~13k cells from this.
• I then remove any genes that have all 0 expression in this response matrix subset. I end up with ~100 genes, so my response matrix is essentially 100 x 13k
• I treat each guide as having a unique target, so I end up with 113 targets, for ~113*100 total "discovery pairs"
Here is the sceptre code:
```fmla <- ~ 1 + log(response_n_nonzero_full + 1) +
      log(response_n_umis_full + 1) + log(grna_n_nonzero_full + 1) + 
      log(grna_n_umis_full + 1) # no batch for replogle

scep_nt <- import_data(
  response_matrix = response_matrix_subset,
  grna_matrix = dummy_grna_matrix,
  grna_target_data_frame = nt_grna_target_df,
  moi = "low",
  extra_covariates = nt_cov_df
)

scep_nt <- scep_nt |>
  set_analysis_parameters(
    discovery_pairs = nt_discovery_pairs,
    formula_object = fmla, 
    control_group = "complement"
  ) |>
  assign_grnas(method = "thresholding", threshold = 1)

scep_nt_no_qc <- scep_nt |>
  run_qc(
    n_nonzero_trt_thresh = 0,
    n_nonzero_cntrl_thresh = 0,
    response_n_umis_range = c(0, 1),
    response_n_nonzero_range = c(0, 1),
    p_mito_threshold = 1
  ) |>
  run_discovery_analysis(parallel = TRUE, n_processors = 10)

scep_nt_qc <- scep_nt |>
  run_qc(
    n_nonzero_trt_thresh = 1, # the only difference!
    n_nonzero_cntrl_thresh = 0,
    response_n_umis_range = c(0, 1),
    response_n_nonzero_range = c(0, 1),
    p_mito_threshold = 1
  ) |>
  run_discovery_analysis(parallel = TRUE, n_processors = 10)```
what do you make of this?

I noticed this because initially I had good calibration with Mixscale (with the caveat that it's not really Mixscale at this point, because no "target" has at least 5 DE genes, so no mixscale scores are actually computed or used), but that's because it has this argument:
> `min.pct`	
> only test genes that are detected in a minimum fraction of min.pct cells in either of the two populations. Meant to speed up the function by not testing genes that are very infrequently expressed. Default is 0.1. Same as logfc.threshold, if split.by is set and more than 1 split.by group exists, thiswill be applied to each group and if any of them satisfies this criteria, the gene will be selected.
When running at the default of `min.pct = 0.1`, things looked calibrated, except it was only testing ~30% of genes. But when I did `min.pct = 0`, so that every (nt-target, gene) pair is tested, then there was poor calibration! Although in this case it's poor calibration in the wrong way, as Mixscale computes many p-values of exactly 1.... but regardless, sensible calibrated behavior only emerged from using this filtering argument. I tried to mimic it with sceptre and found what i shared above.

**[2026-02-07 22:35] Louis Deutsch**
*EDIT: this can be ignored. I answered my own question here*

I guess I don't understand why we are having these pairs though. I thought that, by construction, this dataset would have treatment cells for every (target, response) pair.

For example, in `scep_nt_no_qc@discovery_result` I see the one row like this:
```       response_id     grna_target n_nonzero_trt n_nonzero_cntrl pass_qc p_value fold_change se_fold_change log_2_fold_change significant
            <char>          <char>         <int>           <int>  <lgcl>   <num>       <num>          <num>             <num>      <lgcl>
1: ENSG00000170289 non-targeting93             0               5    TRUE   8e-05           0              0              -Inf        TRUE```
but when i check how many cells express this target I have this:
```> sum(dummy_grna_matrix[nt_grna_target_df$grna_id[nt_grna_target_df$grna_target == "non-targeting93"], ])
[1] 605```
so I thought the complement analysis here should have 605 treatment cells versus 13k - 605 cells! Why are there only 5 cells being used for this pair?

------------------------------------

Edit: some more diagnostics. The gRNA assignment is being handled correctly:
```> all.equal(dummy_grna_matrix, sceptre::get_grna_assignments(scep_nt_no_qc) + 0)  # doing +0 to convert from bool to int
[1] TRUE```
so sceptre does consider there to be 605 cells assigned to the target `"non-targeting93"`.

**[2026-02-08 01:14] Louis Deutsch** _(thread reply)_
Well, unsurprisingly, it wasn't this easy. I reran calibration on the full replogle-rd7 using the pipeline (with thresholding instead of mixture as we did before), with those bad NTs removed, and `n_nonzero_trt = 1`, but no luck. It had been run with no qc for n_nonzero_trt before, but that wasn't the issue 

**[2026-02-08 16:34] Louis Deutsch**
Hey, I was just going over Gene's analysis from some months back comparing Velten's crispat object to ours, and I'm realizing the dataset I've been using this whole time actually has the questionable NT guides labeled as targeting. I have 113 NT guides in my data which is what they have. I think the original amount is 130 NT guides. I'm not sure if this means I should change anything I'm doing, but I wanted to flag it!

A consequence of this is that these calibration plots already have the questionable NT guides removed, so there isn't that low-hanging fruit to try

**[2026-02-09 02:13] Louis Deutsch**
FR-Perturb is not calibrated on Replogle. I am using the same dataset that I describe <https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1770549399131619|two posts ago>. I am providing the covariates `'response_n_nonzero_full_log', 'response_n_umis_full_log', 'grna_n_nonzero_full_log', 'grna_n_umis_full_log'` , so all computed on the full data and the same that sceptre gets. I use the skew-t approximation, and I don't specify a NT group so each (target, gene) test uses the cells expressing `target` vs the average of all cells for `gene`.
  - link: https://katsevichlab.slack.com/archives/C05MA6XN5LZ/p1770549399131619

**[2026-02-09 11:13] Eugene Katsevich**
Hi <@U0524GR916C>, I think the conclusion we arrived at back in the day is that it is hard to have good calibration far into the tail of Replogle. My recollection is that Velten simply used fewer pairs, i.e. did not look as far into the tail, and so their calibration appeared better. So I don’t think Velten had any kind of “secret sauce” that lead to awesome calibration. Despite Replogle’s messiness, I don’t want to discard it yet. What I want to check is *the extent to which each method is able to separate positive control pairs from negative control pairs, as measured by ROC curves and areas under those curves.* To make the comparison fair, you should evaluate these metrics based on a common set of pairs analyzed. I did not follow every detail of your messages but it sounded like different methods were doing different QC and looking at different sets of NTs. I would turn off QC for all methods, run all methods on the same Cartesian product of pairs, and then after the fact, if you want, exclude pairs not passing some threshold for some QC metric you think is sensible, but uniformly across methods.

**[2026-02-22 10:27] Timothy Barry** _(thread reply)_
Nice, this result is very clear. This is the high-MOI data, right?

**[2026-02-24 06:30] Louis Deutsch**
I'm working on my combined positive + negative control report, but in the meantime here is one interesting plot. For this plot, each x axis point `x` gives a number of false positives (FPs) on the neg. controls (with the caveat that some of these are actually true positives due to the NT guide issue in Replogle). I then found the p-value threshold for each method such that we got that many FPs, and i counted the number of true positives (TPs) obtain from that threshold (so the threshold varies by method).

This shows that, for these pairs, when a p-value threshold is used such that there is only 1 FP, then mixscale rejects more on-targets than sceptre (~50 to sceptre's ~20), but for every x &gt;= 2 (ie for every accepted number of FPs &gt;= 2) then sceptre has more on-targets rejected. x = 1 is probably somewhat noisy, so in general sceptre has the best density of TPs in its predicted positive sets for these other methods. I think this is pretty good evidence for sceptre having the best FDR!

Thoughts on this?

**[2026-02-28 03:09] Louis Deutsch**
hey, I'd like to discuss the QC that I do. Gene, on Thursday you proposed doing QC based on the total number of cells with positive expression, regardless of trt or ctrl group, so filtering on `n_nonzero_trt + n_nonzero_cntrl &gt;= threshold`. I'm experimenting with that now. Sceptre has the best FDR control regardless, but there are some rough p-values that sneak through with this QC, even when I take what feels like an excessively large threshold. sceptre still has some of those repeat 8e-5 p-values, and mixscale has a bunch exactly equal to 1. The first plot shows calibration with `n_nonzero_trt + n_nonzero_cntrl &gt;= 20`. The second plot shows ROC curves. This QC removes 93060 pairs.

All of those bad p-values are all filtered out if I just do `n_nonzero_trt &gt; 0` . With this filtering we remove 119344 pairs. Based on these, I think the latter QC does a better job at getting rid of poorly performing pairs without unnecessary filtering. But Gene you mentioned that this QC is not typical in the community? Also, if they're not calibrated anyway, maybe it doesn't matter? Sorry for the overlapping plot labels, i just shoehorned that label into existing code.

**[2026-03-01 21:02] Eugene Katsevich**
Hi <@U0524GR916C>, the sceptre p-values of 8e-5 are not inherently problematic; they signal a failure of the skew-normal fit, which triggered the backup option of direct resampling-based p-values, which bottom out at 8e-5. As we see from your plots, this does not impact sceptre calibration. It seems that Mixscale FDR inflation is alleviated by moving from n_nonzero_trt &gt; 0 to thresh = 20, presumably because of the apparently conservativeness that emerges in the bulk (kinda like two wrongs make a right). I am not sure where all those p-values = 1 come from. Before we decide what to do, could you please see what happens when you choose thresholds = 5 and 10? This may help build intuition for what’s going on. <@U0239H5UC9E>: It would be great to hear your thoughts as well.

**[2026-03-01 21:09] Louis Deutsch**
QC `n_nonzero_trt + n_nonzero_cntrl &gt;= threshold` with `threshold = 5`

**[2026-03-01 21:14] Louis Deutsch**
QC `n_nonzero_trt + n_nonzero_cntrl &gt;= threshold` with `threshold = 10`

**[2026-03-02 10:24] Timothy Barry** _(thread reply)_
Thanks for the update. As a result of having missed last week's meeting, I'm not totally following the plot. It looks like we are comparing various pairwise QC thresholds? One option is to filter on `n_nonzero_trt`, while the other is to filter on `n_nonzero_trt + n_nonzero_cntrl >= threshold`? Do I understand correctly that neither of these options are implemented in the sceptre package? And we want to decide which is best? Is the advantage of `n_nonzero_trt + n_nonzero_cntrl >= threshold` is that it mitigates the selection bias issue?

**[2026-03-02 10:41] Eugene Katsevich** _(thread reply)_
&gt; It looks like we are comparing various pairwise QC thresholds? One option is to filter on `n_nonzero_trt`, while the other is to filter on `n_nonzero_trt + n_nonzero_cntrl &gt;= threshold`? And we want to decide which is best?
Essentially correct, though the specific instance of the first option Louis is considering is filtering based on `n_nonzero_trt &gt; 0` . I am not sure how Louis arrived at this choice, but it was what he had been doing prior to last week’s conversation.
&gt; Do I understand correctly that neither of these options are implemented in the sceptre package?
I guess the first option is implemented via `n_nonzero_trt_thresh = 1L, n_nonzero_cntrl_thresh = 0L`, right?
&gt; Is the advantage of `n_nonzero_trt + n_nonzero_cntrl &gt;= threshold` is that it mitigates the selection bias issue?
Yes. Overall, I felt averse to further tying ourselves to a style of QC that we do not believe in. I thought that filtering on the number of cells in which a gene was expressed was more standard. But as you can see from Louis’s results, these different choices have unexpected effects on the behavior of Mixscale, with the filtering on gene expression apparently causing Mixscale p-value to shift towards being less significant (including a big point mass at 1), improving its calibration. On the other hand, `n_nonzero_trt &gt; 0` gives good calibration in the bulk for Mixscale but worse miscalibration in the tail.

**[2026-03-03 02:04] Louis Deutsch**
Here are the results with no filtering. Mixscale has a lot of near exact 1 p-values. I looked into these a little and I think they're genuine, as far as I can tell.

Sceptre still has more positives for each number of false discoveries, but now Mixscale has better estimated FDR.

I think the issue is that mixscale assigns essentially all of these pairs a p-value of 1, which overall helps with these negative pairs, but the actual performance is not as advertised. The final plot shows that when `n_nonzero_trt == 0` mixscale basically guarantees a p-value of 1, so it's not even trying to be calibrated / uniform, but that does help the metrics in this case. Sceptre still has a sensible looking null distribution, and has some power on the positive controls in this setting while mixscale doesn't.

**[2026-03-05 20:10] Louis Deutsch**
hey, in the past i mentioned that i was unsure if it was fair to use the full covariates for sceptre vs if i should use the ones computed from the relevant subset. Mixscale only has access to the subset, so their perturbabtion signature calculation is trying to remove technical factors only using the subset, not the full data. that does make the comparison slightly unfair. I'm not going to worry about this for now, but just noting it!

Edit: mixscale also controls for the total expression of each cell (our `response_n_umis`), which they are computing from just the subset. So i think this is definitely unfair, since both methods use the same covariate but computed from different things!

**[2026-03-05 21:26] Louis Deutsch**
*sceptre with and without covariates*

On the replogle subset with 100 genes, 100 targets, and 50k cells, I got the following:
• with all 4 default covariates it took *31 min* to test all 10k discovery pairs
• with only `response_n_nonzero_full` and `response_n_umis_full` (the "full" denotes that I'm using the ones from the full data, not just computed from the subsample) it took *28 min*
• With only `response_n_umis_full`, which most closely matches Mixscale's regression (aside from full vs subset covariates), it took *27 min*
*Conclusion*: the number of features we use is not the reason for the time discrepancy between sceptre and mixscale

**[2026-03-06 03:21] Louis Deutsch**
*Question: was the problem having all-zero genes in the data?*

In my current test dataset with 100 genes, 100 targets, and 50k cells, ~20% of genes had all zero expression. I had been filtering out all-zero cells but not all-zero genes. That's also why FR-Perturb didn't calculate every pair.

I just reran except now I only took genes that had at least one positive UMI. Same size, just no all-zero genes now.

Unfortunately, mixscale's time is unchanged (~9.5min), and sceptre's decreased by only 2 minutes (32min to 30min)

*Conclusion: the discrepancy was not caused by many genes with ZERO expression* (but maybe lowly expressed vs highly expressed still could matter)

**[2026-03-06 10:15] Eugene Katsevich**
Hi <@U0524GR916C>, the vectorized NB regression idea of `glmGamPoi` is a nice trick I hadn’t heard of before that I can imagine helps a lot. It’s analogous to how a matrix-matrix multiplication is probably a lot faster than the numerically equivalent set of matrix-vector multiplications. So I concur with your assessement that “mixscale might be faster because it’s just faster? no resampling and many model fittings are batched.” If we had infinite time, I could imagine reworking `sceptre` under the hood to exploit the same vectorized NB regression insight. However, we don’t. So I think we’re going to need to come to terms with the fact that `sceptre` is slower than Mixscale and play up the part of the narrative that it’s a lot more memory-intensive as well as the possibility for parallelization via the pipeline. *As a last sanity check of the above conclusion, I do recommend you run sceptre with its usually pairwise QC thresholds and see if the per-pair runtime significantly decreases.* Based on these results, we can determine next steps.

**[2026-03-06 16:42] Timothy Barry** _(thread reply)_
I do not think the other low-MOI methods we were considering would be computationally performant.

**[2026-03-06 19:57] Louis Deutsch**
*Question: does running sceptre with default QC change things?*

I reran sceptre on the same dataset I've been using (100 targets, 100 genes that all have at least one positive UMI count, 50k cells).

QC results: *~3000 pairs pass* out of the 10k discovery pairs

Time: only *6.5min* now. We are testing ~30% of the pairs but it ran in ~20% of the original runtime (originally 30min), so there's a slight gain here.

Stage: for the pairs that pass QC, *only 0.2% (7 out of 3000) of pairs have a stage 3 p-value!* Versus with no QC, 10% of pairs had stage 3 p-values. So we did get a big reduction there, and that probably accounts for the relative speedup beyond just the reduced number of pairs.

*Conclusion: sceptre's default QC reduces the fraction of stage 3 p-values, and provides a speed-up beyond just reducing #pairs* 

**[2026-03-06 21:44] Eugene Katsevich**
I think the challenge we face with applying our pairwise QC (beyond the fact that we don’t like it) is that it doesn’t result in a Cartesian product of pairs, which complicates the computational comparison with other methods. If there was a filtering we could do just at the gRNA and/or gene levels to get a smaller Cartesian product while having an effect similar to that of our pairwise QC (at least in terms of decreasing the number of stage-3 p-values), I think that would be ideal. If not, we can apply Tim’s suggestion of proportionally scaling down the runtimes of the other methods, but this is a bit messy and not fully principled, in my opinion. <@U0239H5UC9E>, what do you think?

**[2026-03-08 22:43] Louis Deutsch**
I just reran where now I'm *picking the 100 genes with the most non-zero UMI counts for the 50k cells that I take*, and I'm back to no QC so it's a cartesian product for everything. Still 100 genes, 100 targets, 50k cells.

Results: *mixscape runs in 10min, sceptre runs in 13min*, so pretty close! Only *0.1% (10 out of 10,000) p-values are stage 3.*

This shows that there is a cartesian product where there are very few stage 3 p-values, and for which the two methods are pretty close!

**[2026-03-09 10:51] Eugene Katsevich**
Indeed, that’s a very encouraging result! So it seems fair to say that sceptre’s computation time generally increases as gene expression decreases. So the result you obtained is a best-case scenario. I think our computational comparison may be more fair if we consider something closer to an average-case scenario, probably of the form “randomly sample genes from among those satisfying a minimum expression threshold.” Hopefully the minimum expression threshold will cut off genes for which sceptre computation is the highest (and very lowly expressed genes are often dropped anyways), while getting a representative selection from the remainder of the distribution.

**[2026-03-09 10:59] Timothy Barry** _(thread reply)_
>  So the result you obtained is a best-case scenario.
Is this the case? Louis seems to have convincingly demonstrated the the difference in computational performance between sceptre and Mixscale was driven primarily by the large proportion of "type 3" pairs (i.e., pairs for which we compute an empirical p-value). In practice, when applying the standard pairwise QC thresholds, I would expect the proportion of type 3 pairs to be very small, e.g. on the order of 0.1% or so. In this sense, Louis has put us into to the more standard problem setting (computationally) with his latest simulation run.

**[2026-03-11 23:49] Louis Deutsch**
Here are the results for *Gasperini*

We basically don't have stage 3 p-values for this! There are 10k discovery pairs and the most any of these datasets has is 6 stage 3 p-values.

So the stage 3 p-value issue seems specific to low MOI, or at least Gasperini is sufficiently dense so that we don't need to resort to them

**[2026-03-12 10:41] Timothy Barry** _(thread reply)_
That's a helpful result.
&gt; So the stage 3 p-value issue seems specific to low MOI.
I think it's probably specific to the Replogle dataset.

**[2026-03-13 10:11] Eugene Katsevich**
Based on this plot, I would take 60% as the threshold and run with it. As I mentioned before, right now the goal is to do a first-pass analysis; this is a higher priority than choosing the threshold exactly right. What do you think, <@U0239H5UC9E>?

**[2026-03-20 06:23] Louis Deutsch**
I've added the plots for all 3 QC methods to this <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/replogle-report-highly-expressed.pdf|pdf>
pass_qc_75 = pairs with the gene in the top 25% of activity
pass_qc_trt = pairs with n_nonzero_trt &gt; 0
pass_qc_def = n_nonzero_trt &gt;= 7 &amp; n_nonzero_cntrl &gt;= 7 ["def" for sceptre default]
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/replogle-report-highly-expressed.pdf

**[2026-03-20 18:04] Louis Deutsch**
<https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/gasperini-report-highly-expressed.pdf|here> are the key plots for gasperini

Same QC definitions:
pass_qc_75 = pairs with the gene in the top 25% of activity [in terms of non-zero entries]
pass_qc_trt = pairs with n_nonzero_trt > 0
pass_qc_def = n_nonzero_trt >= 7 & n_nonzero_cntrl >= 7 ["def" for sceptre default]

I haven't done all the computational benchmarking for gasperini, so I don't know if this top 25% threshold is appropriate, but my takeaway from these is that all 3 QC methods basically produce the same results

I'm thinking that we can then make the decision based on replogle, and gasperini will be fine for whatever we choose
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/gasperini-report-highly-expressed.pdf

**[2026-03-23 15:08] Louis Deutsch**
ok i've just pushed an update to the <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/replogle-report-highly-expressed.pdf|replogle pdf> and now it has the key plots for the following QCs:
• pass_qc_75: top 25% for n_nonzero
• pass_qc_50: top 50%
• pass_qc_30: top 70%
• pass_qc_20: top 80%
Here are the percentages of pairs passing each of these QCs (in that pdf, this appears at the beginning with total counts instead of percents).
```      pass_qc_75 pass_qc_50 pass_qc_30 pass_qc_20
FALSE     60.74%     21.89%      1.85%       0.35%
TRUE      39.26%     78.11%     98.15%      99.65%```
So pass_qc_30 and passs_qc_20 are doing very little filtering, as is reflected in their results. They both look essentially the same.

Also, as a reminder, I'm taking my *existing* pos and neg datasets and filtering down to the pairs passing this QC.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/association/writeups/replogle-report-highly-expressed.pdf

**[2026-03-23 18:06] Timothy Barry**
Great, thanks. I wonder if we could match the mean expression of the Replogle dataset to that of the Gasperini dataset? I.e., if the mean gene expression of Gasperini is 2 UMIs/cell, then we could select the top X% of genes of Replogle such that the average expression is 2 UMIs/cell. Would that be reasonable?

**[2026-03-24 11:35] Eugene Katsevich** _(thread reply)_
But I like the idea of having the same *absolute* (rather than *relative*) QC applied to Gasperini and Replogle, e.g. UMIs/cell &gt;= X. For a fixed X, we’ll retain a greater fraction of genes for Gasperini than for Replogle, which is what we want.

**[2026-03-24 12:10] Eugene Katsevich** _(thread reply)_
Looking at Seurat, the <https://satijalab.org/seurat/reference/findmarkers|gene-level QC thresholds used with DE> are defined in terms of cell numbers expressing feature rather than average UMIs/cell, probably to exclude bursty expression situations where you have one cell that have a UMI count of 100 and then all other cells with a UMI count of 0. So I think we should follow their lead and also use `n_nonzero` rather than `n_umis` type thresholds.
  - link: https://satijalab.org/seurat/reference/findmarkers

**[2026-03-24 13:05] Eugene Katsevich** _(thread reply)_
It is hard to get around the fact that, for low MOI, the relevant set of cells changes for each pair. One way to reason about this is on average over perturbations. For example, let’s say we are considering 1000 perturbations. If we subset to genes expressed in at least 1000 cells, then a gene passing this filter can be expressed in as few as one cell per perturbation.

**[2026-03-24 13:41] Eugene Katsevich**
<@U0239H5UC9E>: Do we anticipate any daylight between the statistical performance of the sceptre package and pipeline on gRNA assignment or association analysis? I am asking because I am wondering whether the pipeline needs to be included in the statistical benchmarking. One place where some daylight might emerge is if for high MOI we use CRT in the package and permutations in the pipeline. Or do we want to just use permutations everywhere so as not to create confusion? But that’s not the default so we’d have to justify it.

**[2026-03-24 18:16] Louis Deutsch** _(thread reply)_
i just added them here <https://github.com/Katsevich-Lab/sceptre3-project-v2/tree/singularity-containers/benchmarking/data>
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/tree/singularity-containers/benchmarking/data

**[2026-03-24 18:35] Louis Deutsch** _(thread reply)_
ok on HPC I put them in `/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/input_data/&lt;dataset name&gt;/gene_summary_stats.csv`, where `&lt;dataset name&gt;` is either `replogle-rd7` or `gasperini`

**[2026-03-25 10:41] Eugene Katsevich** _(thread reply)_
Then we have two choices: either benchmark sceptre pipeline’s statistical performance for association, as distinct from the sceptre package, or use permutations within the sceptre package. Both of these choices have pros and cons. Perhaps the lesser of two evils is to keep the sceptre package at its default (CRT) for high MOI, and then include statistical benchmarking of the permutation-based pipeline in the appendix? Including it in the main text would probably be distracting.

**[2026-03-26 16:48] Eugene Katsevich**
Folks, I think the right version of the QC metric I proposed is
```(number of cells in which the gene is expressed) x MOI / (number of gRNAs) x (gRNAs / target)```
I think this quantity makes sense for low and high MOI. This is the best we can do without bringing the treatment indicators in. If we want to bring the treatment indicators in, then probably we should define the QC metric for a gene as the average of `n_nonzero_trt` for all pairs involving that gene. I think the metric above is trying to get at the latter metric. <@U0239H5UC9E>: What do you think?

**[2026-03-26 17:56] Eugene Katsevich** _(thread reply)_
Ok <@U0524GR916C> I propose for now you use the QC criterion
```(number of cells in which the gene is expressed) x MOI / (number of gRNAs) x (gRNAs / target) &gt;= 7```
to select genes and run with this unless/until there is some big issue.

**[2026-03-27 14:41] Louis Deutsch**
Things are looking good so far! 

For MOI, would this be the “official” MOI, or would I use the average as computed from sceptre’s assignments? 

**[2026-03-27 19:15] Louis Deutsch**
I was just googling about this and saw something about the *MOI of the infected cells only*, vs the *marginal MOI*. I'm using the marginal MOI right now.

In code, I think this would be
```moi_observed &lt;- Matrix::colSums(grna_assignment_matrix) |&gt; mean()

is_infected &lt;- Matrix::colSums(grna_assignment_matrix) &gt; 0
moi_observed_given_infected &lt;- Matrix::colSums(grna_assignment_matrix[, is_infected]) |&gt; mean()```
For Replogle, I have `moi_observed = 1.47` while `moi_observed_given_infected = 1.61`.

For Gasperini, it's `34.46` vs `34.48`.

So not a huge difference, but this isn't a distinction I'd considered before.  Any thoughts? I'm going ahead with the marginal `moi_observed` for now, since I figure the more conservative one is safer, and it doesn't matter much either way, but I'm just curious.

**[2026-03-28 12:24] Timothy Barry** _(thread reply)_
I'd agree with the marginal observed MOI. But quick question: are we not filtering out cells that do not contain a gRNA on the Replogle data? I'm not sure this matters very much, just curious.

**[2026-03-28 14:00] Eugene Katsevich** _(thread reply)_
I wouldn’t pay much attention to this difference but I agree with Tim we should be filtering out cells with no detected gRNAs.

**[2026-03-28 16:41] Timothy Barry** _(thread reply)_
(This is assuming the rest of the results using this QC look reasonable.)

**[2026-03-30 05:10] Louis Deutsch**
For posterity, if anyone is running on the `16xl.q` queue, there can actually be differences in behavior there. The `sceptre-pipeline` by default has a limit of 5min for each `assign_grnas` process. The standard `short.q` wasn't enforcing this, so I didn't realize that this default existed, but `16xl.q` did, so these were all getting killed at 5min even though the exact same code worked on the default queue.  I'm just saying this since it took me some time to figure out why changing queues was giving me errors, and i did not expect different enforcement of things like this on different queues!

**[2026-04-01 13:44] Eugene Katsevich** _(thread reply)_
These awkward choices (and the awkward 98.9k) seem to be manifestations of the fact that varying the number of cells is the wrong choice, at least for Replogle (low MOI). It now seems the right thing to do is to vary the number of genes and/or targets, and to just include whatever cells are involved in those targets. It seems originally that varying the number of cells simultaneously with varying the number of pairs was changing two things at once, but in fact, the number of cells _per pair_ isn’t changing that much, so I am not very concerned by the fact that the _total_ number of cells is changing. What do you think? I am available between now and 2:30pm to discuss via Zoom, if that would be helpful.

**[2026-04-02 05:22] Louis Deutsch**
*Association results:*

All key association results can be found <https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/manuscript/reports/association/association-analysis.pdf|here>, aside from whatever extra datasets we'll be running the pipeline on.

Outline:
1. Replogle
    a. statistical
    b. computational
2. Gasperini
    a. statistical
    b. computational
Gasperini only has 377 on-targets, so that's why the downsampled BH plots look like they do.
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/manuscript/reports/association/association-analysis.pdf

**[2026-04-02 07:36] Louis Deutsch**
*Guide assignment simulation*

I am rethinking the simulation, and I'd appreciate thoughts on this tomorrow if we have time! As a heads up I haven't had time to actually try this, so it's all pending seeing if it actually works. I'll be doing that tomorrow after our meeting.

As a reminder of the methods:
• pertpy's `assign_mixture_model()`: a Poisson-Gaussian mixture model fit via MCMC to the log-transformed counts. 
• CLEANSER: either a Poisson-NB or NB-NB mixture, fit to the raw counts via MCMC, with a scaling factor for total library size. 
• crispat's `ga_poisson_gauss()`: a Poisson-Gaussian mixture fit to the log-transformed counts. MAP estimates are obtained via stochastic variational inference.
• SCEPTRE: a latent variable Poisson GLM, fit via a modified EM algorithm, controlling for some guide covariates and possibly batch
My goal with the simulation: to confirm recovery of ground truth, and to test sensitivity to certain aspects of the data

I propose the following:
1. I fit a NB-NB mixture on a "reasonable" guide vector from Replogle
    a. I do this so it's not too extreme in any direction, so this serves as the default 
    b. no bells and whistles on this one. It's fine that this is what CLEANSER uses, because it won't be the model that CLEANSER uses once the modifications are applied. 
2. I sample datasets from that NB-NB mixture, and modify them in some of the following ways:
    a. sparsity
        i. for a sparsity fraction \pi, i can either zero out \pi fraction of the dataset or sample from the NB-NB conditioned on being positive, to fill in \pi fraction of existing zeros
        ii. how different is this from varying the mixing fraction?
    b. batch 
        i. I'll divide the cells into two even groups, independent of perturbation status
        ii. i will then add a constant to the log mean of one group while subtracting from the other, to sweep out a batch effect
        iii. I won't include the other covariates SCEPTRE uses. When batch_effect is small, sceptre's features don't help, but then as the batch_effect grows, sceptre's use of batch helps out more and more. 
        iv. what about library size? worth trying to do too?
    c. separation
        i. decrease difference in means until they overlap
        ii. increase variance until they overlap
        iii. or some other mixture model quantities
I would do this in the "main effects only" style, where I use the default settings and vary just one parameter up and down at a time, and maybe only for ~3 unique values. I will also use 20k cells max. In the Fishash paper they do 20k cells so I'm comfortable not going above 20k. That is really helpful for getting this done fast because I was doing way more cells before.

To summarize: I'm proposing to shift from trying to make a good model of a gRNA (like using regressions to estimate covariate effects) to using a simpler model, and then putting the effort into varying that simple model and seeing what happens.

I can test this out with SCEPTRE on my laptop so, pending our discussion, I'll do a sceptre-only prototype tomorrow.

**[2026-04-02 10:41] Timothy Barry** _(thread reply)_
If we are going to revisit simulating gRNA data, it may be worthwhile (at some point) to consider using the <https://github.com/williamaeberhard/robnb|robNB package> to fit the gRNA count model in a more robust way.
  - link: https://github.com/williamaeberhard/robnb

**[2026-04-04 20:05] Louis Deutsch** _(thread reply)_
Ok glad to hear it's not a fundamental limit! It is running, it's just causing the pipeline to be hours slower than it actually needs to be if things were running concurrently.

i think this has also happened with the assignment jobs, but it's not as noticeable because they run faster and i don't need as many of them

And every job does eventually get queued, and since I'm on 16xl.q I don't have long times waiting IN the queue (like where `qstat` shows `qw`). That's what's weird about this problem. It is taking a long time for every single job to even appear in `qstat` at all, regardless of whether or not it's running or waiting once it's there. I've been all over the place with claude + chatgpt on this one and i can't pin it down to something i can control.

If there's nothing obvious here, I will just not worry about it and let things take as long as they need. This is making my sceptre-pipeline runs take 5+ hours when it'd be ~1hr if they were concurrent, but it's manageable

**[2026-04-04 20:10] Eugene Katsevich** _(thread reply)_
My first impulse is to have you reach out to the sysadmin folks with some documentation of your issue to try to get it resolved. But recall that these runs are optional improvements to your association benchmarking story. I would prioritize getting all the assignment stuff done before spending another minute on this issue.

**[2026-04-04 20:11] Louis Deutsch** _(thread reply)_
ok will do! And yes I am not spending much time on it, it's mostly background! I just need to get these done before i can do some assignment runs, so that's why them taking 8 hours when it should be 1 is slowing me down

**[2026-04-04 20:17] Eugene Katsevich** _(thread reply)_
If your assignment runs are on hold because of this, I’d recommend killing existing association jobs, submitting assignment jobs, and then returning to association jobs if you finish assignment jobs and have time left over. If the issue recurs with assignment as well, then we’ll have to go to the sysadmins.

**[2026-04-04 20:20] Louis Deutsch** _(thread reply)_
ok sounds good, I won't be spending any more time on association than this. I am making good progress with assignment in the meantime, since there is plenty i can do on my laptop for that

**[2026-04-06 19:55] Louis Deutsch**
*gRNA assignment simulation update*

I have decided not to pursue making my own guide assignment simulation, or at least not for my defense.

• for one thing, that recent Fishash paper did a simulation that included sceptre, crispat, and cleanser, and I don't have much of a value add beyond that. Their simulation is also more biologically appropriate than the mixture models I'd been trying.
• Furthermore, the “simple” models I've tried haven't fit well, so I think the gap between simulation results and reality are large enough that it's not very interesting 


The main type of model I was trying to fit was a 2 component mixture model with a hurdle:
• let Z1, … , Zn be iid Bernoulli(pi), giving our unobserved true perturbation status 
• Y1, …, Yn are the guide UMI counts
• X1, …, Xn are vectors of cell covariates. I only used grna_n_nonzero
• Hurdle: model P(Yi > 0 | Xi, Zi) = inv-logit(a0 + a1 Zi + a2 log(Xi + 1))
• This captures the idea that perturbed cells will be less likely to be zero, as will cells with large grna_n_nonzero 
• Mixture: Sample Yi | Yi > 0, Xi, Zi from a NegBin truncated to 1,2,3,…. with mean log mu_i = b0 + b1 Zi + log(Xi + 1), and each component has its own overdispersion theta0 or theta1 
The problems:
• fitting by maximum likelihood, I was always getting solutions where the non-perturbed class was basically a point mass at 1, and then the other component captured everything else. So the hurdle absorbed the zeros, the NP class absorbed most of the 1s, and then the rest was a single NB. I tried various modifications of this, but nothing simple worked in the time I allotted to this. 
• I also wasn't able to get a good way to smoothly vary the hardness of the separation without making it too biologically unrealistic. I can get nice marginal distributions, but nothing with an interesting but non-trivial class overlap 
• Trying to base this on replogle is hard because of the sparsity
I do think there still could be something useful here, and not redundant with Fishash, but I won't be thinking about this more until post-defense. 

**[2026-04-07 10:08] Timothy Barry** _(thread reply)_
Sounds good. Did you try just fitting a two-component NB mixture model with crispat, as we discussed?

**[2026-04-09 00:24] Louis Deutsch**
I'm following up on the simulation here.

Tim, you asked about how my results were if I used the other assignment methods. I used CLEANSER since (with `--dc`) it fits a NB-NB mixture to the raw non-zero counts, which is very similar to the generative models I was trying.

I have the same problems with CLEANSER that I did with my homebrew simulations. The included plot shows this. The top row is CLEANSER's assignments for my chosen "nice" Replogle guide, and the bottom row is SCEPTRE's assignments for the same guide. The vertical line is at `umi count = 5`, showing that CLEANSER is equivalent to thresholding with a threshold of 5 here.

CLEANSER results:
• One component has a tiny support and absorbs just the smallest values
• Then the other component captures everything else, even though we can see that the "natural" point to draw the line between the two modes in the real data is larger. The "natural" gap there is around `umi count = 15`.
I think this is reasonable enough in terms of a way to get a classification, but it's not great for producing realistic data to test these on. The simulated data would share this same lack-of-fit problem, so our results would look artificially good.

**[2026-04-09 09:02] Louis Deutsch**
hi Tim, here's what I'd like to discuss today.

1. guide simulator [30min]
2. guide QC [10 min]
3. what subsample of full grna matrices to run on [20 min]
-------------------------------------------------------------------------------------------------------------------------------

*1. Guide simulator*

I'm thinking it might work now! Here's how I'm doing it.

Data:
1. Y_1, ..., Y_n are the real counts
2. X_1, ..., X_n are the real `grna_n_nonzero` (I picked this as the most interesting single covariate)
Model:
1. Latent Z_1, ..., Z_n ~ iid Ber(\pi) indicate true perturbation status
2. Hurdle: model logit P(Y_i &gt; 0 | Z_i, X_i) = a0 + a_pert Z_i + log1p(X_i)
    a. Offset for covariate to keep simpler (its coefficient was estimated at 1.1 when i did allow it to vary)
3. Positive part: model Y_i | Y_i &gt; 0, Z_i, X_i as:
    a. non-pert group: Y_i | Y_i &gt; 0, Z_i = 0  ~ Zeta(b0) 
        i. note: as I'm writing this i realize i forgot to include X_i in the non-pert group. I'll do that if we move ahead with this
    b. pert group: Y_i | Y_i &gt; 0, Z_i = 1, X_i ~ NegBin(\log \mu_i = b_pert + log1p(X_i), \theta) 
        i. offset here too
I fit this via maximum likelihood, and then I tweaked the parameters a bit to get the first of the 3 plots below. That's my "very separated" baseline. I then add more overlap from there.

I added overlap by varying 3 parameters jointly:
1. For non-pert: decrease b0 to make the non-pert's tail heavier
2. For pert: decrease b_pert while also decreasing theta to flatten out the dist without changing the range too much 
The 3 plots together show a sequence of datasets made this way. They are ordered so using your arrow keys to go back and forth between them shows the changes really clearly.

Potential problems with this:
• it's still pretty manual. I did base in on an MLE, and it's still pretty close to those values, but it's manual adjustments from there based on my feelings about the results
• a power law non-pert tail might be a bit much? The parameters I'm using all have the first few moments finite. 
*2. Real data analysis: guide QC*

No matter what I do, I will not be running every guide. Do we want to apply a similar guide QC to choose our guides?


*3. Real data analysis: pin down the datasets we want to use*

I am currently subsampling both cells and guides. I do this because:
• cleanser is memory hungry: with 50k cells and 2.5k guides, I hit 50GB
• pertpy is really slow: with 50k cells and 2.5k guides, it takes 10 hours
    ◦ this is on GPU though. It might actually be faster not on GPU
Subsampling allows me to get everything done in a reasonable time frame, and I still get coherent results.

*Main problem with this: the sparsity of replogle.* Taking a &lt;10% sample leaves very few true positives to learn from.

The alternative I am considering, and will be making a decision on today:
*Use many/all cells, use only a few guides*
• Statistical benchmarking
    ◦ A lot of our analyses involve pairwise similarities between methods for what they assign to each guide
    ◦ We will have way fewer datapoints (guides) in those analyses now
    ◦ but each datapoint might be better quality (better model fit) and more realistic (accurate # true positives)
• Computational benchmarking
    ◦ Simpler: only one parameter to vary (#guides)
    ◦ But now we can't estimate how runtime/memory depend on #cells
Thanks!

**[2026-04-11 03:18] Louis Deutsch**
*I figured out why the Fishash paper authors report high memory usage for SCEPTRE*

I'm going through their code to double check some things for my simulations, but along the way I figured this out. Of all the methods they benchmark, SCEPTRE is the only one that requires a gene expression matrix.

In their nextflow pipeline, they have this
```if (!params.skipSceptre) {
  get_gex_matrix_rds(sims_ch)
  run_sceptre_mixture(get_gex_matrix_rds.out)
  out_ch = out_ch.mix(run_sceptre_mixture.out)
}```
which is only triggered for SCEPTRE.

The function `get_gex_matrix_rds()` is how they simulate a gene expression matrix (GEX) for SCEPTRE. The key lines there:
```sce <- splatSimulate(
  batchCells = ncol(guide_counts),
)
colnames(sce) <- colnames(guide_counts)
gex_mat <- counts(sce)```
They then save `gex_mat` as an RDS with no further processing. Well, *the matrix that `counts(sce)` returns is dense!!!*

They are using the default number of genes here, which is 10k. And they keep the number of cells fixed at 20k. So *every single run, SCEPTRE's run script is loading a useless 10k x 20k dense matrix!* (useless because it's just a dummy gene matrix). And no other method uses this.

Below is the plot from their paper with memory use in their simulations, and SCEPTRE hovers around 50GB every time. Reflecting this large constant being added to every run. They also only report memory for simulations, so they wouldn't have noticed that this doesn't happen with real data using a sparse gene matrix.

So, in conclusion, (1) benchmarking is hard, and (2) we should probably tell them since they're unfairly reporting SCEPTRE as being very high memory use!

**[2026-04-16 05:55] Louis Deutsch**
i just found something that seems important re: assignment:

i think sceptre can indeed suffer from using a poisson distribution: for the guides with the heaviest tails, sceptre ends up making a ton of positive predictions. i think this is because getting a poisson to cover all of that means it'll eat up more ambient noise too

the x axis in this plot is `mean(guide_umis[guide_umis &gt;= 10])` for each guide. sceptre breaks down for the highest ones where a poisson fit is the worst!

**[2026-04-16 05:56] Louis Deutsch**
and pertpy can break down with the lightest tails, which we see there

**[2026-04-16 05:57] Louis Deutsch**
the entirety of our association analysis is conditioned on the sceptre assignments. Nothing i can do about that now, but I'm probably gonna have to address that somehow!

**[2026-04-16 05:59] Louis Deutsch**
in general i have not found that the overdispersion is sufficiently explained by the covariates. if it was, then sceptre would be ok. and related, crispat doesn't use any and seems to do really well.

**[2026-04-16 06:19] Louis Deutsch**
for comparison, this plot has the umi counts for a lowly and highly expressed guide (in terms of grna_n_nonzero) for each. gasperini is super light tailed compared to replogle. orders of magnitude different for the highly expressed one here. it makes sense that gasperini wouldn't show this problem

**[2026-04-16 08:36] Louis Deutsch**
hey, another random thought on this: i wonder if the approximation SCEPTRE makes when fitting the mixture breaks down some when the UMI distribution is really heavy tailed. with Gasperini it's ok, but maybe even a few observations that are 500x from the mean could pull a likelihood with a tail that decays as fast as a poisson's does

**[2026-04-16 09:10] Eugene Katsevich**
Hi <@U0524GR916C>, thanks for all these insights. It does seem like sceptre’s assignment model is worth revisiting (after your defense). By the way, <@U0239H5UC9E>, we will be skipping today’s regularly scheduled meeting because you are at a conference and also to accommodate Louis’s thesis defense prep.

**[2026-04-21 15:03] Louis Deutsch**
Also, I have been experimenting with a NegBin instead of Poisson and the results are promising! Each point in this plot is a comparison between `Jaccard(simulated guide ground truth, what this method estimated on that guide)` for the current Poisson GLM SCEPTRE and my new NegBin version. We can see the points that SCEPTRE does poorly on for the largest value of the perturbed group's means (and those are the simulated guides that got the heaviest tails), and the NegBin fixes that!

I did these using the same offset trick, and then fitting a mixture of NB(exp(gamma + offset), disp1) and NB(exp(offset), disp0) via a pure R EM algo.

Something else to consider: I wonder if the approximation used for computing the offset doesn't work as well with Replogle. Gasperini isn't that heavy tailed, but Replogle has values far enough out in the tail that they do distort the fit, and especially of the intercept, which is the most important parameter here since we don't have super strong covariate effects (also as evidenced by crispat doing so well). For my experiments here, I got my offsets by fitting on the data *excluding the top 0.2% of responses in magnitude*. Those are either true perturbations, multiplets, or something else weird, so either way I think it makes sense to do this!

This NegBin model also outperforms crispat in average Jaccard similarities with the truth:
``` pert_mean    avg_scep_jacc   avg_scep_nb_jacc  avg_crispat_jacc
1 121               0.773            0.805            0.778
2 242               0.784            0.833            0.808
3 485               0.810            0.856            0.837
4 970               0.745            0.871            0.855```
Granted this NegBin model is very close to the data generating process, but it's still promising. 

The downside: it's a lot slower to fit.

Note that this is all on simulated data, but since I observed this exact same thing on the real data, I think this is meaningful!

**[2026-04-21 15:15] Louis Deutsch**
To further justify some light filtering of the responses that we use to fit our offset model:

on a simulated guide, the model *WITH* the most extreme 0.2% of observations (so, all the data): beta0_hat = -1.5
on the same guide *WITHOUT* these, beta0_hat = -3.5

those 0.2% of most extreme observations are really distorting the intercept estimate!

This is still just simulated data, but I think the principle applies!

**[2026-04-21 15:38] Louis Deutsch**
Also, the fact that pertpy could struggle on the lowest mean replogle guides suggests that we don't need tails quite as heavy as exp(Pois). So NegBin seems to be heavy enough for Replogle, at least, and then we also have the covariates on top of that

**[2026-04-22 03:04] Louis Deutsch** _(thread reply)_
Unless we hear from Gene, would you be down to meet for a bit anyway to talk about sceptre gRNA assignment? And maybe priorities for getting the benchmarking paper out 

**[2026-04-22 15:11] Louis Deutsch**
Hi Gene, Tim and I just met since I have to proctor an exam during our meeting time tomorrow.

Main takeaways:
• we discussed what our "minimum viable" updates to gRNA assignment could be. I will be updating with some new experiments shortly. I had been fitting  pi NB(exp(gamma + offset), theta1) + (1-pi) NB(exp(offset), theta0) but that is too complex. Tim suggested estimating a single overdispersion theta, and plugging that in, to keep the EM as being over just two parameters. I am coding this up right now and will share shortly
• We discussed whether or not we should do some kind of robustification to the initial GLM fit that gives us our offsets. In my experiments the most extreme UMI counts, which are either perturbed or otherwise not ones we want to use, _do_ exert a meaningful effect on the intercept of this model. But for now we will leave all of this as-is

**[2026-04-22 16:37] Louis Deutsch**
Also, I'd like to check in about our “must haves” for the benchmarking paper!
• what other datasets should we include? Is barnyard the main one?
• How important is repeating the association analysis with different assignments? I’m hoping that this seems less important if we can do a simple change to sceptre that fixes the current issue 
• I am going to do more simulation settings, in particular (1) vary batch effect; (2) vary noise distribution
• Array jobs
• larger problem sizes to really highlight the pipeline’s performance?
• Anything else I should have in mind as I change gears to this?

**[2026-04-22 16:46] Timothy Barry** _(thread reply)_
• What is Barnyard? My feeling is that we have the key datasets already.
• IMO, just use the same set of assignments throughout all of association.
• That sounds like a good idea, especially taking a grid over batch effects.
• That sounds like a very good idea.
• The trans analyses of the full Gasperini and Replogle datasets should be sufficient.
• Beyond adding overdispersion to the gRNA assignment mixture model, we should keep in mind benchmarking the ondisc matrix initialization step. We already have code for creating the ondisc objects, so it's really just a matter of computing the runtime.

**[2026-04-22 16:50] Louis Deutsch** _(thread reply)_
barnyard is a dataset that I learned about from the Fishash paper where there is some amount of ground truth knowledge for guide assignments. So it's real data + ground truth. Somehow. 

**[2026-04-22 16:51] Louis Deutsch** _(thread reply)_
And the Fishash people kindly made their code available so I'm hoping that would make it a lot faster for me to add this dataset to my benchmarking 

**[2026-04-22 17:01] Timothy Barry** _(thread reply)_
Is this mainly for gRNA assignment? I think the cleanser paper used a barnyard dataset as well.

**[2026-04-22 17:14] Timothy Barry** _(thread reply)_
I might personally deprioritize adding new datasets at the moment. A more extensive simulation study would be higher yield than a new dataset IMO.

**[2026-04-22 17:17] Louis Deutsch** _(thread reply)_
ok cool, sounds good! I also am much more curious to see the simulation results, and especially how other methods break down in the presence of strong covariate effects.

and what if some of Replogle's overdispersion comes from marginalizing over a high cardinality batch effect? Granted we're not capturing that either, but who's to say

**[2026-04-22 21:50] Eugene Katsevich**
Hi <@U0524GR916C>, I apologize for missing your earlier message about rescheduling tomorrow’s meeting to today. I’m preparing for a talk on Monday, which has been taking up a lot of my bandwidth recently.

Regarding next steps, I think expanding your simulation and implementing array jobs would both be useful. Another thing I recall from your defense is that your benchmarking plots were missing points, perhaps due to experiments that had not finished running, or experiments you did not run with the pipeline, or experiments you could have run with in-memory methods but didn’t because you thought they would be too slow (e.g. sceptre in-memory on gRNA assignment). As we near manuscript-writing, having a more complete set of benchmarking plots will be useful. I also don’t remember how you chose the degree of parallelization for each dataset and each setting, but it might be good to standardize that as well.

Furthermore, I think it is not too early to start sketching out an outline for the paper, and then systematically filling it in. Such an outline could help us identify and prioritize remaining to-do items. Louis: Perhaps you can start a document somewhere on GitHub called manuscript and start us off with a rough outline. We can then iterate on the outline and shape the content into figures. Note, Louis, that biology-oriented papers are structured around their figures, which have many panels each.

**[2026-04-30 13:04] Eugene Katsevich**
Hi folks, <https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html|here> is a plan for us to get sceptre onto Bioconductor.
  - link: https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html

**[2026-04-30 14:48] Timothy Barry**
What I would do.

• Place a lower bound of about 0.1 on size parameter
• Use one size parameter across conditions
• Do robust estimation of the gRNA UMI count model (maybe)

**[2026-05-06 10:35] Eugene Katsevich**
Hey <@U0239H5UC9E>, I've responded to Laura on <https://github.com/Katsevich-Lab/sceptre/issues/202|Issue #202>. I suspect that Laura's R session is still pointing to an older version of `parallely`. If that ends up being the case, do you think we should revert to the backwards-compatible `parallely` code we had before to avoid other users running into this issue?
  - link: https://github.com/Katsevich-Lab/sceptre/issues/202

**[2026-05-14 14:38] Louis Deutsch**
Here's the possible bug:
<https://github.com/Katsevich-Lab/sceptre/blob/e803b71068d8b55f7a8bc4a8b8c1c6780c1ebe79/src/mixture_functs.cpp#L46>
  - link: https://github.com/Katsevich-Lab/sceptre/blob/e803b71068d8b55f7a8bc4a8b8c1c6780c1ebe79/src/mixture_functs.cpp#L46

**[2026-05-21 16:39] Eugene Katsevich**
Hey <@U0239H5UC9E>, thanks for merging PR 200. I have brought <https://github.com/Katsevich-Lab/sceptre/pull/203|PR 203> up to date with main, and it has passed checks and is waiting for your review (when you have time). Thank you!
  - link: https://github.com/Katsevich-Lab/sceptre/pull/203

**[2026-05-25 01:02] Louis Deutsch**
Lukas updated pertpy <https://github.com/scverse/pertpy/pull/987|here> to address the GPU speed concerns.

Apparently they completely abandoned their current approach (MCMC-based) and now just implement `crispat`. So `pertpy` is now equal to accelerated `crispat`.

I'll have to repeat all my pertpy benchmarking, but I do support this decision on their part from what I've seen.
  - link: https://github.com/scverse/pertpy/pull/987

**[2026-05-25 01:10] Louis Deutsch**
beating crispat is now more important

**[2026-05-25 11:20] Timothy Barry**
I see. Couldn't we still include benchmark comparisons to the original pertpy implementation? As that is what they published? I don't think we necessarily should be expected to accommodate the latest update to every method in our benchmarking analysis.

**[2026-05-25 15:21] Timothy Barry**
That sounds reasonable to me. IMO it would be a waste to throw out our current pertpy results.

**[2026-05-27 19:03] Louis Deutsch**
Yeah, and this kind of simplifies the comparison landscape. 

For laptop methods, it’s basically just sceptre vs crispat, and now for distributed it’s basically nextflow-powered sceptre vs GPU-powered crispat 

**[2026-05-28 12:43] Eugene Katsevich**
Hi folks, <https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html|here>'s an updated summary on where we stand with respect to Bioconductor submission.
  - link: https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html

**[2026-05-28 14:42] Louis Deutsch**
```josep at louis-laptop in ~/katsevich-lab/sceptre/R (glmgampoi-experiment●)
$ grep -r "ondisc:::" .
./import_functs.R:  response_matrix &lt;- ondisc:::create_odm_from_r_matrix_internal(
./import_functs.R:  grna_matrix &lt;- ondisc:::create_odm_from_r_matrix_internal(
./neg_control_functions.R:      n_nonzero_m &lt;- ondisc:::compute_n_trt_cells_matrix_ondisc(
./pairwise_qc_functs.R:      out &lt;- ondisc:::compute_n_ok_pairs_ondisc(
./pairwise_qc_functs.R:      out &lt;- ondisc:::compute_nt_nonzero_matrix_and_n_ok_pairs_ondisc(
./assign_grna_functions.R:      ondisc:::threshold_count_matrix_cpp(```

**[2026-06-03 18:44] Eugene Katsevich**
Hey <@U0239H5UC9E>, <https://github.com/Katsevich-Lab/sceptre/pull/210|PR #210> is ready for your re-review; I removed the deprecated alias.
  - link: https://github.com/Katsevich-Lab/sceptre/pull/210

**[2026-06-04 13:19] Eugene Katsevich**
Hi folks, <https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html|here>'s an updated summary on where we stand with respect to Bioconductor submission.
  - link: https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/singularity-containers/sceptre-dev/sceptre-bioc-submission.html

**[2026-06-05 12:28] Timothy Barry**
Here is my paragraph for our Bioc cover letter.

The primary container for single-cell data in the Bioconductor ecosystem is `SingleCellExperiment`. Instead of using `SingleCellExperiment` to store and manipulate the expression data, we implemented several specialized data structures and algorithms better suited to the specific access patterns required of perturb-seq analysis. After data import (e.g., from a collection of 10x Genomics .mtx files), sceptre computes various cell-wise covariates relevant to perturb-seq analysis, including the gRNA with the greatest UMI count in a given cell and the fraction of UMIs attributable to this gRNA. `sceptre` does so using optimized C++ routines. Next, `sceptre` converts the sparse gene and gRNA expression matrices from column-accessible format into row-accessible format. This transformation is crucial, as all subsequent steps -- including QC, gRNA assignment, and gRNA-to-gene association testing -- require rowwise access to the gene and gRNA expression matrices. `SingleCellExperiment` (and allied packages, such as `scran`) can store these matrices and covariates, but they do not implement the specific, optimized perturb-seq computations of `sceptre`.

Another reason I didn't use  `SingleCellExperiment` is that it has like 50 dependencies and imports, and it seemed totally unnecessary to explode our dependency list just to use their `Matrix`-wrapper storage container, but I'm not sure we need to highlight this. Maybe we could add a sentence like, "Moreover, using `SingleCellExperiment` would substantially increase the dependency footprint of `sceptre`".)

**[2026-06-12 13:23] Eugene Katsevich**
Hi folks, <https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/e565a1c/sceptre-dev/sceptre-bioc-submission.html|here>'s an updated summary on where we stand with respect to Bioconductor submission.
  - link: https://htmlpreview.github.io/?https://raw.githubusercontent.com/Katsevich-Lab/sceptre3-project-v2/e565a1c/sceptre-dev/sceptre-bioc-submission.html

**[2026-06-13 03:45] Louis Deutsch**
hi all, quick update on which method is currently the winner:

*offsets*: two very similar ones I'll be trying more extensively
• trimmed Pois GLM: remove top 0.1% currently
• thresholded Pois GLM: only use obs with `y <= Y_MAX`. Currently using `Y_MAX=100`
On Replogle, the offsets are always fairly correlated, but sometimes are tightly clustered on the line y=x while other times they are shifted.

The NB GLM is a tiny bit better on the simulations, but 3x slower, so I am sticking to Pois. The difference in performance doesn't seem significant.

*mixture*:
• (1-pi) Pois(exp(o_i)) + pi NB(exp(o_i + gamma), theta), fit with EM jointly over (gamma, theta, pi)
*major speedup:*
•  in the EM algo, I force any cells with `y == 0` to be non-perturbed, so the updates of (gamma, theta) grow in the number of positive obs, not the number of obs.
• This provides a 5x speedup, at basically no cost to performance!
With that speedup, this is actually _faster_ than my pure R implementation of current sceptre! And also performs much better.

Tomorrow I am going to benchmark this on more extensive simulations.

**[2026-06-13 03:58] Louis Deutsch**
There's another speedup we feel tantalizingly close to, but I don’t think it actually works:
• if we didn't filter any cells out, as we do with trim or thresh, then we are fitting many GLMs with the exact same X and only varying y
• One way to achieve this is to cap y, so every cell is still present, but there is a max to how outlying any UMI counts can be. This doesn't perform as well as trim or thresh though
• But now this opens up the option of precomputing regression quantities 

My brilliant idea that didn't work well: use OLS to regress log(y+1) on X. Every one of these has the same hat matrix H, so we can get fitted values that are like an offset by doing `H %*% log1p(y)` for the same H every time. Unfortunately this doesn't work well enough statistically. Maybe it could be improved, because this is biased, but also there are so many exact zeros that the GLM is probably worth the compute time. 

If we change cells between each regression, then I'm not sure if there's anything we can precompute that's meaningful. Either way, this can be added later so I'm not worrying about this more for now. 

**[2026-06-13 10:49] Eugene Katsevich** _(thread reply)_
Nice job on finding this speedup! Do I understand correctly that this is not numerically equivalent to the full mixture model fit, because cells with zero counts would have a small but still positive probability of being assigned to the positive class? I am willing to believe that it's close.

**[2026-06-13 17:08] Louis Deutsch**
QUESTION: if we had to pick, is it better to have a mixture that makes false positives, or false negatives?

My instinct is that it’s better to have a smaller sample size of actually perturbed cells, vs a larger sample size with some non-perturbed mixed in.

(Edit: this is relevant to a simulation I'm making with a heavier tailed non-pert distribution, to probe if using a Pois for that is going to cause false positives vs a NB or something. Although i don't think the NB will make much of a difference for NP unless theta_nonpert is tiny, since exp(o_i) < 1)

**[2026-06-13 17:25] Louis Deutsch** _(thread reply)_
oh yes this isn't blocking me. It's just for helping me narrow down the search space for failure modes with the simulations I'm making today

**[2026-06-15 08:57] Louis Deutsch**
hi all, here are the results i want to show today
• i've generalized the simulations, inspired by Cellbender, and the non-pert distribution is a lot more realistic now
• My new proposed one is close to current sceptre, but it is shocking hard to beat the bugged sceptre
• the new method and sceptre are both generally better than crispat here. 
• the new method fits in basically the same time as pure R current sceptre. Both are 2x faster than crispat.

main takeaways:
• the new proposed model is computationally fine, but in most simulations it is slightly worse than sceptre
• the only big difference is for the largest mean for the old simulations. That is a very overdispersed one (theta = 0.1) so that could be it. But I haven't yet found the exact reason why bugged does bad there, and not my new approach.
• i am close to understanding better what the bugged one is doing that makes it so much better, but i'm not quite there
<https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/sceptre-nb/nb-bench_v3.pdf>
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/sceptre-nb/nb-bench_v3.pdf

**[2026-06-15 09:48] Eugene Katsevich** _(thread reply)_
Hi Louis, thank you for sharing this writeup ahead of our meeting. Here are some immediate questions:
1. One of the big changes in your writeup is a new simulation model, because you find our current simulation model is insufficiently faithful to real data. Is it true that this new simulation model has not been applied to gRNA counts before? Did you consider adopting any of the existing gRNA count simulation models instead, such as the FishHash simulation model?
2. Could you help me parse the statement "My new proposed one is close to current sceptre, but it is shocking hard to beat the bugged sceptre." Is there a distinction being made between "current sceptre" and "bugged sceptre"? Isn't the current sceptre bugged? Didn't we find that your new gRNA assignment method was substantially better than current/bugged sceptre (without adding trimming to the latter)?

**[2026-06-15 14:48] Louis Deutsch** _(thread reply)_
1.  I'm now using basically a simplified version of what Fishash does, which is a variant of Cellbender. I'm doing this since I find it more interpretable, but after this meeting I am going to fully use their simulation too. 

**[2026-06-15 14:52] Louis Deutsch** _(thread reply)_
2. Current sceptre == bugged sceptre. Overall, I found one simulation setting where current sceptre did pretty poorly and the new proposed method beats it. For other cases, current sceptre does better but not a ton better. So right now it seems like the new method is more robust in general and avoids the big failure mode, but isn't uniformly better. 

I haven't finished pinning down exactly why current sceptre failed so clearly on that one simulation setting. I am going to figure this out very soon and use that to make more simulations exploring that failure mode.  

**[2026-06-15 14:56] Louis Deutsch** _(thread reply)_
I am still working on this, but I think that the log likelihood ratio of the Pert to NonPert densities is approximately linear in o_i with the bug, and is exponential in o_i when implemented correctly. I suspect that the bug reduces the effect of the offsets, guarding against poor offsets, and makes it closer to thresholding, which works surprisingly well. I included a threshold of 10 in my simulations for reference. 

**[2026-06-17 11:00] Timothy Barry** _(thread reply)_
Got it. Are we going to ignore for the time being that our gRNA assignment house is not quite in order?

**[2026-06-18 02:44] Louis Deutsch**
hi all, I have found a model that is kinda principled and is either equivalent to SCEPTRE or better, for all 3 of my simulations + Replogle + Gasperini! Computationally it's similar too.

*Basic idea: try to reverse engineer something that's more statistics parole than statistics jail, but has nearly identical performance, so we can swap this in with minimal noticeable changes.* 

I discussed the Cellbender model on Monday. We can use a simple version of that to get our additive likelihood: do
```Y_i ~ Pois(exp(o_i)) + Ber(pi) Pois(exp(delta))```
assuming independence, so perturbed is Pois(exp(o_i) + exp(delta)). We are using all Pois; don't have an endogenous component; and there is no cell-level scaling in the pert component. But it is a special case of Cellbender, so there is a decent amount of justification for this form. Also I think this form makes sense. Perturbed cells aren't free from background noise. 

-------

The more controversial part: this model doesn't work well if we fit it normally. It has the same problems that the correctly-implemented offset model has, because e^\delta is huge, and so our decision boundary ends up being way too large.

But if we compute the decision boundary using the log of the mean shift, ie \delta instead of \exp(\delta), then we get extremely similar results to SCEPTRE. I've been going over the math and I think there are some arguments that can be made, but I'm not spending much time on that until I hear if this replacement idea seems plausible. My intuition is that e^\delta is a terrible estimate for most data points, because its dominated by the tail, and really we care about the classifications for y~10. Every method gets y >100 right. It's y ~ 10 that's hard, and using the log mean shift of \delta makes us less sensitive to the tail, and do better on small values.

<https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/sceptre-nb/additive-model.pdf|This pdf> has the metrics plots for the 3 simulations I showed on Monday, plus Replogle. Gasperini is not included because every method is basically identical there. This new method is basically the same or better than SCEPTRE across the board!

I think this method also has the same robustness to the offsets, so we don't seem to need to robustify those with this. That can be another justification for the log. And it does better in the simulation case where SCEPTRE falls apart, and it also has much less of the analogous issue on Replogle where SCEPTRE starts to make way too many positive predictions for the heaviest tailed guides!
  - link: https://github.com/Katsevich-Lab/sceptre3-project-v2/blob/singularity-containers/benchmarking/guide-assignment/sceptre-nb/additive-model.pdf

**[2026-06-18 02:45] Louis Deutsch**
To summarize this ^ post: I am now proposing a more principled model, that performs identically or better (in my 3 simulation settings + Replogle + Gasperini, at least), so we can swap this in and justify this instead of justifying or keeping the bug.

**[2026-06-18 02:53] Louis Deutsch**
My next steps:
1. get the full Fishash simulations running, to further avoid overfitting to the few cases I'm looking at
2. Barnyard?
Pending confirmation that we are interested in this model switch approach.

**[2026-06-18 11:21] Timothy Barry**
Nice! That seems like great progress.

Still, I wonder whether we may consider something even simpler, e.g. the following model:
```G_i \sim Pois(\mu_i) (or alternately NB)
\log(\mu_i) = \beta_0 + \beta_1 X_i + \log(l_i)```
Here, `l_i` is some notion of "library size," either `grna_n_nonzero` or `grna_n_umis`. (It may be necessary to robustify `grna_n_umis` by truncating the top and bottom 1% of cells.)

My hunch is that the initial regression step may be throwing us off. The intuitive advantage of the regression step is that it enables us to "absorb" the intercept and all covariates into the offset terms. But this may be causing problems.

Would it be feasible to assess a simple alternative such as this? I believe the EM algo should have closed form updates for the Poisson distribution and for the NB distribution with fixed/known theta. (Theta probably could be estimated or set outside the EM algo.) We probably could use much of the same C++ code.

**[2026-06-18 11:24] Timothy Barry**
To explain my reasoning a bit more, sceptre is basically an approximation to a GLM mixture model. (The full GLM mixture would be too slow to fit.) At the time I implemented the gRNA assignment module, I thought it would be important to account for all covariates, including things like batch. However, Louis has shown this might not be the case; accounting for a single covariate, i.e. some notion of "gRNA library size," may be adequate.

**[2026-06-18 14:06] Eugene Katsevich**
<@U0239H5UC9E>, does your proposal model the gRNA-negative and gRNA-positive cells using the same dispersion? I recall Louis mentioning at some point that he thinks the two need to be modeled using different dispersions?

**[2026-06-18 14:42] Eugene Katsevich** _(thread reply)_
If we don't feel at peace with any solutions (let's say, after trying for another few days) but want something legit for our submission, we should consider just implementing crispat or something and calling it a day.

**[2026-06-18 17:03] Timothy Barry** _(thread reply)_
I think this is reasonable, although it would be a bit sad. We'd need to figure out what to do with the gRNA assignment benchmarking. Perhaps implementing crispat in optimized C++ code could count as a contribution. (We have a good starting point, and Claude could help.) But one step at a time.

**[2026-06-18 17:11] Timothy Barry** _(thread reply)_
Although part of the utility of the benchmarking analyses is that they've revealed our gRNA assignment house is not quite in order...

**[2026-06-18 20:24] Louis Deutsch** _(thread reply)_
That would turn us vs pertpy into nextflow crispat vs GPU crispat basically 

**[2026-06-18 21:03] Louis Deutsch**
Ok I will experiment with that simplified model, and I am going to think more deeply about the mathematics of the bug. 

Because these additive log-mean models have significantly better precision than crispat in the more realistic simulations. I’m not ready to just use crispat   

**[2026-06-18 21:06] Louis Deutsch**
I'm also exploring a simple extension of using the crispat model but incorporating library size  

**[2026-06-18 21:53] Timothy Barry** _(thread reply)_
This would basically be a Poisson-Gaussian mixture, with offsets for the Poisson and a shifted mean for the Gaussian? I guess that would be equivalent to just subtracting off library size from the count itself for each cell?

**[2026-06-22 18:11] Louis Deutsch**
Hi all, I’m putting together something collecting my attempts and insights so far. I will share tomorrow, and hopefully then we'll be ready to pick a direction for guide assignment in sceptre for now 

**[2026-06-24 09:46] Louis Deutsch**
hi all, i'm still working on putting everything together, but i have a lot more insight into the offsets and covariates. I'll have all of this ready for us to make a decision in our meeting.

Here's one interesting thing: because of Replogle's heavy tails, there is a meaningful amount of circularity in using `grna_n_umis` for guide g without subtracting out Y_{ig} for each cell. In this plot, both models are regular Poisson GLMs with `y ~ log1p(grna_n_umis) + log1p(grna_n_nonzero)` . The _only_ difference is this:
```df_circular <- scep@covariate_data_frame |>
  transmute(
    y = grna_odm[guide,],
    log1p_grna_n_umis = log1p(grna_n_umis ),
    log1p_grna_n_nonzero = log1p(grna_n_nonzero)
  )
df_not_circular <- scep@covariate_data_frame |>
  transmute(
    y = grna_odm[guide,],
    log1p_grna_n_umis = log1p(grna_n_umis - y), # <-- remove y_i's influence
    log1p_grna_n_nonzero = log1p(grna_n_nonzero - (y > 0)) # <-- remove I(y_i>0)'s influence
  )```
I also added one that both fixes the circularity and also only uses cells with `y_i <= 1000`.

Main takeaways
• the circularity makes it look like the offsets are doing a good job in terms of describing the mean for large y_i, but it's really just echoing `y_i`. 
• additionally thresholding to `y <= 1000`  removes essentially all of the apparent signal in this covariate
Other conclusions so far:
• These covariates generally don't have much signal
• The GLMs are dominated by the tail (on Replogle, I am seeing >85% of the deviance, and 98.5% of the total counts, coming from the largest 0.1% of obs). It seems like the tail is heavy enough that it overpowers the sparsity. 
• grna_n_nonzero explains the most, so if we had to pick one, that's it. But it doesn't work as an offset on Replogle, and can have a negative coefficient, so its relationship must be estimated
• I tried a bunch more models to no avail.

**[2026-06-24 09:48] Louis Deutsch**
So i think part of both crispat and bugged sceptre's success is that the covariates don't help much. Crispat ignores them, and the bug guards against them.

**[2026-06-24 10:38] Timothy Barry** _(thread reply)_
Apologies for my lack of understanding, but what do you mean by "circular"? For a given gRNA, letting Y_i denote its cellwise UMI count and l_i denote the (log-transformed) gRNA library size, Y_i and l_i are highly correlated?

&gt; grna_n_nonzero explains the most, so if we had to pick one, that's it. But it doesn't work as an offset on Replogle
Why is this?

**[2026-06-24 10:42] Timothy Barry** _(thread reply)_
Interesting. And "covariate" includes `n_grna_umis`, `n_grna_nonzero`, batch, etc?

**[2026-06-24 10:51] Louis Deutsch** _(thread reply)_
i checked using the 2 grna ones and the 2 gene ones (grna and response total and num_nonzero)

I think batch could make a big difference, because it would affect smallish counts like y ~ 10. That's really the only hard regime in this problem (y_i &gt;&gt; 20 is probably perturbed, y_i &lt;&lt; 20 is probably non-perturbed). But for now, I'm happy to fall back on the "low effort" solution of just doing everything conditioned on batch, so I'm not considering it much otherwise.

**[2026-06-24 11:23] Louis Deutsch** _(thread reply)_
I just mean that cell j's grna_n_nonzero includes Y_{jg} for each g. I think that's the reason for the part that looks like a straight line in that plot. The guide counts are heavy tailed enough that some cells will have grna_n_nonzero_j ~ Y_{jg}.

**[2026-06-24 11:36] Timothy Barry** _(thread reply)_
Ah, OK, I see... What if we were to do a "leave one out" strategy and drop the given gRNA when computing its `grna_n_nonzero`? Although I don't want this to be too much of a diversion...

**[2026-06-25 07:31] Louis Deutsch**
Hi all, here are my bullet points in advance of today's meeting.

GOAL: make a decision regarding guide assignment in SCEPTRE

I have not been able to come up with a fully principled method that doesn't have weird failure modes or hard tuning parameters, in the time I have had to try. For my dissertation, I propose doing a pure R implementation of crispat, and coming back to this later. Claude Code can do a pure R crispat easily.

-----------

Here are the things that make this hard:

• the covariates are very subtle:
    ◦ right now, there is a circularity in `grna_n_umis`, since it includes the current cell's `y_i`. (This also applies to `grna_n_nonzero` but the effect is smaller). The first plot shows the effects on offsets for Replogle. The circularity leads to much smaller offsets for small `y_i` cells. The effects of this are pretty delicate and I haven't been able to pin down a clear conclusion so far. Changes in the offset coefficients or linear predictors don't predictably translate to performance differences. I'm not even sure how bad this is. Generally the circular model fits the training data better, but not always. The coefficients behave pretty differently if the response covariates are present, but the offsets aren't as different.
    ◦ On Replogle, if we only pick one covariate, `grna_n_nonzero` is the best in terms of explanatory power. But it also can behave oddly: using it as the only covariate, if we restrict to modeling cells with `y > 0`, then the coefficient is generally estimated to be significant and *negative*. So increasing grna_n_nonzero is associated with a _decrease_ in y. I suspect this is due to low MOI competition, or something like that, although I haven't tried this on Gasperini to check. I am very hesitant to use it as an offset with a fixed coef of 1 because of this. This doesn't seem like it can be used as "library size" scaling the mean up and down.
• the offsets are also subtle:
    ◦ I think there is a narrow range where they can help. For tiny y_i, we should assign to NP no matter o_i; and for large y_i, we should assign to P no matter o_i. It's rare that our decision isn't essentially entirely determined by y_i alone.
    ◦ The model's fit is mostly determined by the largest y, and the offsets for small y are all over the place. A good model probably doesn't want to be very sensitive to them. Crispat doesn't use them at all, and the additive form of bugged sceptre reduces their impact. 
    ◦ Batch could really matter, but I think that there aren't that many cells right on the boundary where library size-type coefficient information can really help
    ◦ The second plot helps form my opinion that the Replogle offsets aren't doing a good job characterizing the non-pert mean. This shows the offsets (linear predictors) versus log y, for 4 Replogle guides: the first column uses all y, the second uses thresholding to remove the tail. The color indicates if the circularity is removed. The numbers in the upper left corner are the spearman correlations between the offsets and log y, for the points above the dashed horizontal line. When both the circularity and upper tail are removed, then there is very little relationship between the offsets and y. 
More thoughts on the SCEPTRE bug
• I think this additive form with a small gamma is actually really reasonable
    ◦ we don't care about modeling the pert tail. We just want a good decision boundary. Bugged sceptre accidentally focuses on this region
    ◦ I got very comparable performance to SCEPTRE by using (1-pi) Pois(exp(o_i)) + pi Pois(exp(o_i) + K), with K=10 provided as a hyperparameter
        ▪︎ If we could estimate K this could work, but it seems pretty equivalent to the problem of picking a guide threshold, which is the whole thing we're trying to avoid with our mixture models
Thoughts on CRISPAT
• In some simulations I can beat crispat, but it's pretty reliable across the board
• ignoring the 0s is a good call, and a Pois-Gaussian mixture on the log scale is very sensible
• it's a nice, simple model and my simulations haven't uncovered any surprise failure modes (unlike with SCEPTRE)

**[2026-06-25 08:53] Louis Deutsch**
Also, here are some more method results.

These are the same simulated datasets from before. 3 different settings, all of the form
```Pois(small * L_i) + Ber(endog prob) Pois(larger * L_i) + Ber(pert prob) NB(mu * L_i, theta)```
• first plot: "NPlowvar" is a very short non-pert tail. Mostly exog noise [the initial `Pert(small)` that every cell receives]. 10 is a very large value. "NPhighvar" has a heavier endog part, so non-pert can attain 20 or so. The plot doesn't show it, but I increased the overdisp as i shrank mu. This is the "normal" one, and these all look pretty realistic to real data, marginally. 
• Second plot: fix mu_pert, vary theta over 10, 1, .1. Much more overdispersed than real data looks. 
• Third plot: `endog prob = 0` so just 2 components. The endog Pois is nearly always &lt;= 5. The largest mean NB dist is quite overdispersed, more than real data. Overdisp decreases as mu_pert decreases. Largest mean is also most overdispersed, and its perts go as low as 0, so more overlap. That setting is harder for precision than the others. 
• I wasn't sure if it makes sense to put the L_i factors in the non-pert part, but it doesn't seem "wrong", and I figured I'd err on the side of overdispering vs under.
Key: the name of each method is of the form `&lt;offset description&gt;_&lt;mixture description&gt;`.
For offsets:
• `glmpois` = poisson GLM fit on everything. What current SCEPTRE does
• `glmpoisgrnafix` = poisson GLM using only grna features (`grna_n_nonzero` and `grna_n_umis`) and a circularity adjustment:
```grna_n_nonzero - (y &gt; 0)
grna_n_umis - y```
• `threshglmpois1000` = poisson GLM fit to only the cells with `y &lt;= 1000`. This is my representative of a robust regression
For mixtures:
• `poisbug` is exactly what SCEPTRE does now, with the bug doing `exp(o_i) + gamma` instead of `exp(o_i + gamma)`
• `poisthresh10` is `(1-pi) Pois(exp(o_i)) + pi Pois(exp(o_i) + 10)`
• `pois0nb` is `(1-pi) Pois(exp(o_i)) + pi NB(exp(o_i + gamma), theta)`, optimized jointly over (pi, gamma, theta), and forcing all cells with y=0 to be non-pert
main takeaways
• Sim 1, the "most realistic" one
    ◦ crispat has noticeably the lowest precision, and therefore jaccard
    ◦ Bugged sceptre does really well, and was very hard to beat in my experiments.
• Sim 2, the "extremely overdispersed" one
    ◦ a pretty clear precision-recall tradeoff between methods
• Sim 3, not the most realistic data, but interesting
    ◦ SCEPTRE has a big and unique failure in precision due to making way too many positive predictions. A unique failure model could be a sign of the fragility of an unprincipled method 
    ◦ crispat is the best on this one


**[2026-06-25 09:52] Eugene Katsevich**
Hi folks, I have completed my portion of the preparation for Bioconductor submission (see summary <https://raw.githack.com/Katsevich-Lab/sceptre3-project-v2/95ba3cef519be07851b4d8c3f674d91d863c6d5b/sceptre-dev/sceptre-bioc-submission.html|here>). In particular, <@U0239H5UC9E>, I have filed three PRs for you to review, two in `sceptre` and one in `sceptre-book`. These are summarized in the document linked above. Once you merge these PRs, we will be ready for submission. Please let me know if any of these PRs are not modular enough; I can try to break them up for your convenience.
  - link: https://raw.githack.com/Katsevich-Lab/sceptre3-project-v2/95ba3cef519be07851b4d8c3f674d91d863c6d5b/sceptre-dev/sceptre-bioc-submission.html

**[2026-06-25 09:57] Eugene Katsevich** _(thread reply)_
Separately, <@U0239H5UC9E>, I filed a <https://github.com/Katsevich-Lab/sceptre/pull/218|tiny PR> that fixes the CI trigger so R-CMD-check runs on every pull request (including forks and PRs stacked onto feature branches, which were previously skipped) while limiting push-triggered runs to the `main` and `dev` branches.
  - link: https://github.com/Katsevich-Lab/sceptre/pull/218

**[2026-06-25 12:29] Eugene Katsevich**
Thoughts on *how our mixture model bug interacts with Bioconductor submission:*

Our current high-MOI mixture implementation differs from the EM algorithm that we intended to implement. At the same time, this implementation has been around for a while, is the behavior used in prior sceptre analyses, and has performed reasonably well in our benchmarks. So I do not think we should simply remove it or rush to replace it with an unvalidated new implementation.

I think the honest path is to document the situation clearly: the current `method = "mixture"` option should be described as preserving historical sceptre behavior, rather than as a clean implementation of the intended EM algorithm. We may also want to add a clearer alias such as `method = "legacy_mixture"` while continuing to support `method = "mixture"` for backward compatibility.

In parallel, we should work on adding a more principled Poisson-Gaussian mixture option, following the crispat approach. But I do not think that needs to block Bioconductor submission, provided the submitted package documents the current behavior accurately. Once the new implementation is validated, we can decide whether to make it the recommended/default high-MOI assignment method.

Ideally, we would submit to Bioconductor with the current behavior honestly documented, then add and validate the Poisson-Gaussian option during July/August, with the goal of having the intended high-MOI story settled before the first Bioconductor release of sceptre. Since the detailed fall Bioconductor schedule is not posted yet, I would treat mid-to-late September as the practical target for landing any major default/API changes.

**[2026-06-25 13:21] Timothy Barry** _(thread reply)_
So, the idea is to submit the package soon and then update the mixture functionality before it is released on Bioconductor?

**[2026-06-25 16:13] Timothy Barry** _(thread reply)_
Yes, but the documentation does not have a lot to say about the mixture method. (We might need to scratch "principled" in "a principled mixture model.")