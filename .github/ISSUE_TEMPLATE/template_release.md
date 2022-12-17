---
name: Release template
about: Template for release related issues.
title: ''
labels: 'release'
assignees: ''

---

### Release v[x.x.x]

Release issue to document testing and bug fixing. Check out the GitHub wiki workflow for the [drafting of a new release](https://github.com/ExploreASL/ExploreASL/wiki/Drafting-a-new-Release).


### General information
Test withing the branch **n/a**.
Don't forget to pull a new version of the following repositiories before testing: `Testing`, `FlavorDatabase`, `TestDataSet`!!

### Pre-test before creating release branch

- [ ] UnitTest: **n/a**
- [ ] TestDatasets: **n/a**
- [ ] Flavors: **n/a**
- [ ] Fix and close last issues: **n/a**
- [ ] Create release branch once all the tasks above are done: **n/a**

### Testing

Only start testing once the release branch is created. Test all again after all major fixes.

#### Run UnitTesting - (xASL_test_Unittesting) - all should be passed

- [ ] Linux: **n/a**
- [ ] Windows: **n/a**
- [ ] MacOS: **n/a**

#### Run TestDatasets - (xASL_test_TestExploreASL) - post the final table here in this issue

- [ ] Linux: **n/a**
- [ ] Windows: **n/a**
- [ ] MacOS: **n/a**

#### Run Flavors test - (xASL_test_Flavors([],0,0) - Run only the import part including BIDS2Legacy and stop before processing. Paste the printed errors.

- [ ] Linux: **n/a**
- [ ] Windows: **n/a**
- [ ] MacOS: **n/a**

### Administration

Do these tasks in the release branch while testing:

- [ ] Document `CHANGES.md` using release notes of all issues: **n/a**
- [ ] Update `VERSION` file: **n/a**
- [ ] Merge and close all branches+issues for this version: **n/a**
- [ ] Test-run docu crawler: **n/a**

### Release

Start the release process only once all the tasks above are done:
- [ ] Update reference table in TestDataSet: **n/a**
- [ ] Merge release branch to main: **n/a**
- [ ] Clean project and milestones: **n/a**
- [ ] Rebase develop with latest commits: **n/a**
- [ ] Update version to BETA in develop: **n/a**
- [ ] Crawl new documentation: **n/a**
- [ ] Update online documentation and delete BETA: **n/a**
- [ ] Do the release, copy release text, tag, and Zenodo: **n/a**
