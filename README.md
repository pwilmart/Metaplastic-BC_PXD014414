# Metaplastic-BC_PXD014414
Re-analysis of metaplastic breast cancer 27-sample TMT data from PXD014414


---

### Delta-mass histograms are considered first:

The PAW pipeline uses wider parent ion tolerances in Comet (typically 1.25 Da) for two reasons: (1) accurate mass does not distinguish correct and incorrect matches unless incorrect masses are allowed to have inaccurate masses, and (2) Comet scores are converted to PeptideProphet-like discriminant scores. The second reason needs a reliable DeltaCN score, which requires a large number of theoretical peptides be scored in the search. That is guaranteed in a wider tolerance search.

![deltamass 2+](images/mass_2.png)

These are the deltamass histograms for the 2+ peptide top hits. The top panel shows the wide tolerance search full range (+/- 1.25 Da). Matches to the forward protein sequences are in blue. Matches to the sequence-reversed decoys are in red. We have (as expected) a major peak at 0.0 Da. We also have very small peaks at 0.984 Da (deamidation) and 1.003 Da (an isotopic peak mis-trigger). The lower two panels are expanded scale plots of these two regions. The dotted lines define adjustable windows to capture the peak(s) in each region. The windows will be called the 0-Da (left) And 1-Da (right) regions.

The incorrect scores (in red) are uniformly distributed across the +/- 1.25 Da delta-mass space, as are the incorrect matches to the target sequences (there are also flat distributions in the blue histograms). We will also make score histograms for the delta-mass values that were not inside either the 0-Da region or the 1-Da region. The instruments may not be able to measure an accurate mass for all of the very low abundance peptides.

---

![deltamass 3+](images/mass_3.png)

These are the deltamass histograms for the 3+ peptide top hits. The 0-Da peak is slightly wider and about half as tall as the 2+ peptides. The peaks in the 1-Da region are relatively larger. 3+ peptide tend to be longer so they are more likely to contain an asparagine residue (the source of the deamidations). 3+ peptide will also tend to be heavier with a less prominent monoisotopic peak, so mis-triggering on the first isotopic peak is more likely.

---

![deltamass 4+](images/mass_4.png)

These are the deltamass histograms for the 4+ peptide top hits. We have far fewer 4+ peptides with trypsin compared to 2+ and 3+ peptides. The 0-Da peak gets wider with increasing charge state, and the relative chance of mis-triggers gets greater.

---

### Conditional score histograms are checked next:

The PAW pipeline takes a divide and conquer approach to controlling PSM FDR. We have three charge states (2+, 3+, and 4+) in three delta-mass regions (0-Da, 1-Da, and outside windows) for 9 sets of conditional score histograms. Matches are further divided into unmodified peptides and each modification grouped by mass shift. We only specified oxidized methionine in the Comet searches, so we have one +16 mass shift. There are some samples where semi-trpytic peptides are common. Semi-trpytic searches are supported by separating matches by number of tryptic termini (NTT). We required fully-tryptic specificity, so all of the NTT=1 histograms will be empty.

![scores 0-Da 2+ nomods](images/score_0_2_nomod.png)
![scores 0-Da 2+ oxidation](images/score_0_2_ox.png)

These are the scores for 2+ matches that had delta masses inside the 0-Da window. The unmodified peptides are in the top figure and the oxidized methionine peptides are in the lower figure. The target matches are in blue and the decoy matches are in red. The numerical tables below the histograms show the PSM FDR for the selected histogram at the location of the dotted line (the highlighted row). The dotted line in the right histogram corresponds to the 1% FDR cutoff. The majority of the correct target matches will be retained (to the right of the dotted line).

> The thresholds for a 1% FDR are set automatically, but they can be changed interactively.

---

![scores 0-Da 3+ nomods](images/score_0_3_nomod.png)
![scores 0-Da 3+ oxidation](images/score_0_3_ox.png)

These are the similar score histograms for the 3+ peptides inside the 0-Da window.

---

![scores 0-Da 4+ nomods](images/score_0_4_nomod.png)
![scores 0-Da 4+ oxidation](images/score_0_4_ox.png)

Finally, here are the 4+ peptide scores for the 0-Da window. We have a lot fewer 4+ peptides.

---

![scores 1-Da 2+ nomods](images/score_1_2_nomod.png)
![scores 1-Da 2+ oxidation](images/score_1_2_ox.png)

We have relatively few correct target matches for peptides in the 1-Da regions. We simplify things a little by combining deamidated peptides (note that we did not specify deamidation as a variable modification) and C13 mis-triggers into a single score histogram. We have a difference in the relative magnitudes of the target matches (in blue) and the decoy matches (in red) compared to the 0-Da regions. We will need to set slightly higher score cutoffs to maintain the 1% FDR. Only the 2+ peptides are shown for the 1-Da regions. The 3+ and 4+ peptides are similar.

---

![scores outside 3+ nomods](images/score_out_3_nomod.png)
![scores outside 3+ oxidation](images/score_out_3_ox.png)

Last, but not least, are the matches that did not have defined delta masses (outside of the 0-Da or the 1-Da windows). It was not obvious in the delta-mass histograms above that there was any excess of target matches compared to decoy matches outside of the 0-Da or 1-Da regions. However, we see that there are clearly high-scoring target matches that can be recovered. The sensitivity to recover these peptides is not great, but every little bit helps. Only the 3+ peptides are shown for the "outside" matches. The 2+ and 4+ peptides are similar.

---
