# Branch Cleanup Analysis Report
## gstlearn/gstlearn Repository

**Date:** January 23, 2026
**Base Branch:** dev
**Total Branches Analyzed:** 35

---

## Executive Summary

Based on the analysis of all branches in the repository against the deletion criteria:
- **Branches that can be deleted:** 0
- **Branches to keep:** 35

### Deletion Criteria
A branch can be deleted if:
1. It has no commits ahead of the 'dev' branch
2. Its name doesn't start with v1 or 1 (patch maintenance branch)  
3. It has commits ahead of the 'dev' branch, but the proposed code has already been merged into 'dev' (by another merge)

---

## Analysis Results

### Maintenance Branches (KEEP - 10 branches)
These branches start with 'v1' or '1' and should be kept for patch maintenance:

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

### Branches with Active Open PRs (KEEP - 6 branches)
These branches have open pull requests and active development:

11. **pr818** - PR #818 (open) - "Fix memory leaks" - 6,569 commits ahead
12. **pr819** - PR #819 (open) - "Fix memory leak in FFT" - 6,566 commits ahead
13. **pr826** - PR #826 (open) - "Make FFT thread-safe" - 6,606 commits ahead
14. **pr827** - PR #827 (open) - "Allow using stencils to build PrecisionOp" - 6,609 commits ahead
15. **pr834** - PR #834 (open) - "Fix crash at application exit with VS2017" - 6,611 commits ahead
16. **Issue-816** - PR #821 (open) - 6,582 commits ahead

### Branches with Open PRs but Different Names (KEEP - 2 branches)

17. **Marimo_developments** - PR #820 (merged) - 6,565 commits ahead
    - Note: PR was merged, but branch still has commits ahead
18. **Move-spde-(old-code)-into-MultiLayers** - PR #742 (open) - 6,294 commits ahead

### Branches with Closed PRs (REVIEW - 2 branches)
These branches have closed PRs that were not merged. Need review to determine if work is still needed:

19. **Issue_741** - PR #782 (closed without merge) - 6,436 commits ahead
    - Last commit: 2025-11-25
20. **testing_template** - PR #746 (closed without merge) - 6,351 commits ahead
    - Last commit: 2025-11-06

### Branches with Unique Commits and No Associated PR (KEEP - 15 branches)
These branches have unique commits not in dev and no associated PR found:

21. **Bug_Esteban_2** - 5,515 commits ahead (last: 2025-06-24)
22. **Develo_KD** - 2,261 commits ahead (last: 2023-04-04)
23. **Eigen** - 2,332 commits ahead (last: 2023-12-20)
24. **Modify_setLocator** - 4,339 commits ahead (last: 2025-05-22)
25. **NDesassis-patch-2** - 6,157 commits ahead (last: 2025-09-30)
26. **Robustify-SPDE** - 4,418 commits ahead (last: 2024-12-19)
27. **Test-MSVC-2022** - 6,501 commits ahead (last: 2025-12-12)
28. **copilot/fix-90cd1858-355b-45dc-8a64-045c7407a34f** - 6,157 commits ahead (last: 2025-09-30)
29. **copilot/fix-ec566bb8-c45d-43db-b6f1-9d129d40669e** - 6,093 commits ahead (last: 2025-09-29)
30. **copilot/fix-tuto-pca-maf-test** - 6,314 commits ahead (last: 2025-10-30)
    - Note: Associated with PR #748 which was closed
31. **gneiting** - 6,413 commits ahead (last: 2026-01-06)
    - Note: Related work merged in PR #784 and #813
32. **optim** - 5,418 commits ahead (last: 2025-06-05)
33. **optim_francky** - 6,388 commits ahead (last: 2025-11-20)
34. **transformModifs** - PR #771 (open) - 6,445 commits ahead (last: 2025-11-26)
35. **workflows-factorization-actions** - 6,083 commits ahead (last: 2025-09-29)

---

## Detailed Findings

### Key Observations

1. **All branches have commits ahead of dev**: Every non-maintenance branch shows a large number of commits ahead of dev (ranging from 2,261 to 6,611 commits). This suggests:
   - The repository uses a complex branching strategy where feature branches diverge significantly
   - These branches may have rebased histories or different merge bases
   - The 'dev' branch itself may have been rebased or reset at some point

2. **No branches meet the deletion criteria**: None of the analyzed branches:
   - Have zero commits ahead of dev
   - Are maintenance branches (v1.x or 1.x) that should be deleted
   - Have all their changes already merged into dev

3. **Active development branches**: Several branches (pr818, pr819, pr826, pr827, pr834, Issue-816, transformModifs) have active open PRs and should definitely be kept.

4. **Stale branches to review**: 
   - **Issue_741** and **testing_template** have closed PRs and may be candidates for deletion after review
   - **copilot/fix-tuto-pca-maf-test** also had a closed PR (#748)

### Recommendations

#### Immediate Actions
**NONE** - No branches meet the strict deletion criteria at this time.

#### Recommended Reviews
1. **Issue_741** - Review if work from closed PR #782 is still needed
2. **testing_template** - Review if work from closed PR #746 is still needed  
3. **copilot/fix-tuto-pca-maf-test** - Review if work from closed PR #748 is still needed
4. **Marimo_developments** - PR #820 was merged but branch still has commits ahead - investigate why

#### Further Investigation Needed
For the branches with no associated PRs and old last commit dates:
- **Develo_KD** (last commit: 2023-04-04) - Over 2 years old
- **Eigen** (last commit: 2023-12-20) - Over 2 years old
- **Robustify-SPDE** (last commit: 2024-12-19) - Over 1 year old

These may be abandoned feature branches that could potentially be deleted after confirming the work is no longer needed.

---

## Methodology

1. Fetched all remote branches from the gstlearn/gstlearn repository
2. For each branch:
   - Checked if it's a maintenance branch (name starts with v1 or 1)
   - Counted commits ahead and behind the dev branch
   - Checked if branch is an ancestor of dev (fully merged)
   - Searched for associated pull requests
   - Recorded last commit date
3. Categorized branches based on the deletion criteria
4. Cross-referenced with GitHub pull request data

---

## Conclusion

Based on the strict deletion criteria provided, **no branches can be deleted at this time**. All branches either:
- Are maintenance branches (v1.x, 1.x)
- Have active open pull requests
- Have unique commits not yet merged into dev
- Require further review to determine their status

However, there are **4 branches** that warrant review for potential deletion:
1. Issue_741 (closed PR)
2. testing_template (closed PR)
3. copilot/fix-tuto-pca-maf-test (closed PR)
4. Develo_KD (very old, no PR)

It is recommended to consult with the repository maintainers about these branches before making any deletion decisions.
