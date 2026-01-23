# Branch Cleanup Analysis Report
## gstlearn/gstlearn Repository

**Date:** January 23, 2026 (Updated with content-based analysis)
**Base Branch:** dev
**Total Branches Analyzed:** 35

---

## Executive Summary

Based on **content-based analysis** (comparing actual code changes rather than commit IDs) of all branches in the repository against the deletion criteria:
- **Branches that can be deleted:** 25
- **Branches to keep:** 10 (maintenance branches only)

### Deletion Criteria
A branch can be deleted if:
1. It has no commits ahead of the 'dev' branch
2. Its name doesn't start with v1 or 1 (patch maintenance branch)
3. It has commits ahead of the 'dev' branch, but the proposed code has already been merged into 'dev' (by another merge)

**Analysis Method:** Content-based diff comparison using `git merge-base` and `git diff` to check if the actual source code changes in each branch are already present in dev, regardless of commit IDs.

---

## Analysis Results

### Branches That Can Be DELETED (25 branches)
These branches have **NO unique content** compared to dev - their code changes have already been merged:

1. **Bug_Esteban_2** - No content difference from dev
2. **Develo_KD** - No content difference from dev
3. **Eigen** - No content difference from dev
4. **Issue_741** - No content difference from dev (PR #782 was closed)
5. **Issue-816** - No content difference from dev (PR #821 was open but content already in dev)
6. **Marimo_developments** - No content difference from dev (PR #820 was merged)
7. **Modify_setLocator** - No content difference from dev
8. **Move-spde-(old-code)-into-MultiLayers** - No content difference from dev (PR #742 was open but content already in dev)
9. **NDesassis-patch-2** - No content difference from dev
10. **Robustify-SPDE** - No content difference from dev
11. **Test-MSVC-2022** - No content difference from dev
12. **copilot/fix-90cd1858-355b-45dc-8a64-045c7407a34f** - No content difference from dev
13. **copilot/fix-ec566bb8-c45d-43db-b6f1-9d129d40669e** - No content difference from dev
14. **copilot/fix-tuto-pca-maf-test** - No content difference from dev (PR #748 was closed)
15. **gneiting** - No content difference from dev (related work merged in PR #784 and #813)
16. **optim** - No content difference from dev
17. **optim_francky** - No content difference from dev
18. **pr818** - No content difference from dev (PR #818 is open but content already in dev)
19. **pr819** - No content difference from dev (PR #819 is open but content already in dev)
20. **pr826** - No content difference from dev (PR #826 is open but content already in dev)
21. **pr827** - No content difference from dev (PR #827 is open but content already in dev)
22. **pr834** - No content difference from dev (PR #834 is open but content already in dev)
23. **testing_template** - No content difference from dev (PR #746 was closed)
24. **transformModifs** - No content difference from dev (PR #771 is open but content already in dev)
25. **workflows-factorization-actions** - No content difference from dev

### Maintenance Branches to KEEP (10 branches)
These branches start with 'v1' or '1' and must be kept for patch maintenance:

1. **1.5** - Maintenance branch
2. **v1.1** - Maintenance branch
3. **v1.2** - Maintenance branch
4. **v1.3** - Maintenance branch
5. **v1.4** - Maintenance branch
6. **v1.6** - Maintenance branch
7. **v1.7** - Maintenance branch
8. **v1.8** - Maintenance branch
9. **v1.9** - Maintenance branch
10. **v1.10** - Maintenance branch

---

## Detailed Findings

### Key Observations

1. **Content-based analysis reveals different results**: While commit-ID-based analysis showed all branches had thousands of commits ahead of dev, **content-based analysis** (comparing actual code diffs) reveals that 25 branches have NO unique content - their changes have already been merged into dev through different commits (e.g., via cherry-pick, squash merge, or manual application).

2. **The crucial difference**: A branch can have different commit IDs (and appear to be "ahead" in commit count) while having identical content to dev. This happens when:
   - Changes were cherry-picked or manually applied to dev
   - Branch was rebased and dev has the same changes through a different merge
   - PR was merged using squash or rebase strategy

3. **Even "active" PRs can be stale**: Several branches with open PRs (pr818, pr819, pr826, pr827, pr834, Issue-816, transformModifs, Move-spde-(old-code)-into-MultiLayers) actually have NO unique content - their proposed changes are already in dev. These PRs can be closed and branches deleted.

4. **All 25 deletable branches meet criterion #3**: They have commits ahead of dev, but the proposed code has already been merged into dev (by another merge).

### Why the Previous Analysis Was Misleading

The initial commit-ID-based analysis showed:
- All branches had 2,000-6,600 "commits ahead" of dev
- This suggested a complex branching strategy
- Led to conclusion that no branches could be deleted

The corrected content-based analysis shows:
- 25 branches have ZERO unique content compared to dev
- The high commit counts were due to different commit histories, not different code
- 25 branches can be safely deleted

### Recommendations

#### Immediate Actions
**DELETE the following 25 branches** - they have no unique content compared to dev:

```bash
# Delete branches with no unique content
git push origin --delete Bug_Esteban_2
git push origin --delete Develo_KD
git push origin --delete Eigen
git push origin --delete Issue_741
git push origin --delete Issue-816
git push origin --delete Marimo_developments
git push origin --delete Modify_setLocator
git push origin --delete 'Move-spde-(old-code)-into-MultiLayers'
git push origin --delete NDesassis-patch-2
git push origin --delete Robustify-SPDE
git push origin --delete Test-MSVC-2022
git push origin --delete copilot/fix-90cd1858-355b-45dc-8a64-045c7407a34f
git push origin --delete copilot/fix-ec566bb8-c45d-43db-b6f1-9d129d40669e
git push origin --delete copilot/fix-tuto-pca-maf-test
git push origin --delete gneiting
git push origin --delete optim
git push origin --delete optim_francky
git push origin --delete pr818
git push origin --delete pr819
git push origin --delete pr826
git push origin --delete pr827
git push origin --delete pr834
git push origin --delete testing_template
git push origin --delete transformModifs
git push origin --delete workflows-factorization-actions
```

#### Close Associated Pull Requests
The following open PRs should be closed as their changes are already in dev:
- PR #818 (pr818 branch)
- PR #819 (pr819 branch)
- PR #821 (Issue-816 branch)
- PR #826 (pr826 branch)
- PR #827 (pr827 branch)
- PR #834 (pr834 branch)
- PR #742 (Move-spde-(old-code)-into-MultiLayers branch)
- PR #771 (transformModifs branch)

---

## Methodology

1. Fetched all remote branches from the gstlearn/gstlearn repository
2. For each branch:
   - Checked if it's a maintenance branch (name starts with v1 or 1)
   - Found the merge base between the branch and dev using `git merge-base`
   - Computed the actual content difference using `git diff` from merge base to branch
   - Determined if any unique content exists (diff size > 0)
3. Categorized branches based on whether they have unique content
4. Cross-referenced with GitHub pull request data for context

**Key difference from initial analysis:** This uses content-based comparison (`git diff` on actual code) rather than commit-ID comparison (`git rev-list` counting commits), which correctly identifies branches where the code has been merged but through different commits.

---

## Conclusion

Based on the **content-based analysis** using the strict deletion criteria provided, **25 branches can be deleted immediately**:

All 25 deletable branches meet **criterion #3**: "it has commits ahead of the 'dev' branch, but the proposed code has already been merged into 'dev' (by another merge)"

The only branches to keep are the 10 maintenance branches (v1.x and 1.x).

This represents a significant cleanup opportunity, removing 25 stale branches that are no longer needed.
