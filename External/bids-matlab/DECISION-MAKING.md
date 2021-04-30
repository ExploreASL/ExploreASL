# Decision-making rules

## Introduction

These rules have been taken and adapted from those of the
[BIDS specification](https://github.com/bids-standard/bids-specification/blob/master/DECISION-MAKING.md).

### Maintainers Group

| Name                                                | Time commitment | Role / Scope |
| --------------------------------------------------- | --------------- | ------------ |
| Rémi Gau ([@Remi-Gau](https://github.com/Remi-Gau)) | 5h/week         | Admin        |
| Guillaume Flandin ([ ?? ](https://github.com/ ??))  | ?? h/week       | Admin        |
| Robert Oostenvold ([ ?? ](https://github.com/ ??))  | ?? h/week       | Admin        |

Maintainers may declare a limited scope of responsibility. Such a scope can
range from maintaining a modality supported in the specification to nurturing a
welcoming community. One or more scopes can be chosen by the maintainer and
agreed upon by the Maintainers Group. A maintainer is primarily responsible for
issues within their chosen scope(s), although contributions elsewhere are
welcome, as well.

### Contributors Group

Contributors are listed in the [README](./README.md#contributors) using the
all-contributors bot.

## GitHub Workflow

For the day-to-day work on the BIDS-MATLAB, we currently abide by the following
rules with the intention to:

- Strive for consensus.
- Promote open discussions.
- Minimize the administrative burden.
- Provide a path for when consensus cannot be made.
- Grow the community.
- Maximize the [bus factor](https://en.wikipedia.org/wiki/Bus_factor) of the
  project.

The rules outlined below are inspired by the
[lazy consensus system used in the Apache Foundation](https://www.apache.org/foundation/voting.html)
and heavily depend on the
[GitHub Pull Request Review system](https://help.github.com/articles/about-pull-requests/).


## Rules

1. Every modification of the specification (including a correction of a typo,
   adding a new Contributor, an extension adding support for a new data type, or
   others) or proposal to release a new version needs to be done via a Pull
   Request (PR) to the Repository.
1. Anyone can open a PR (this action is not limited to Contributors).
1. A PR is eligible to be merged if and only if these conditions are met:
   1. The last commit is at least 5 working days old to allow the community to
      evaluate it.
   1. The PR features at least two
      [Reviews that Approve](https://help.github.com/articles/about-pull-request-reviews/#about-pull-request-reviews)
      the PR from Contributors of which neither is the author of the PR. The
      reviews need to be made after the last commit in the PR (equivalent to
      [Stale review dismissal](https://help.github.com/articles/enabling-required-reviews-for-pull-requests/)
      option on GitHub).
   1. Does not feature any
      [Reviews that Request changes](https://help.github.com/articles/about-required-reviews-for-pull-requests/).
   1. Does not feature "WIP" in the title (Work in Progress).
   1. Passes all automated tests and checks if the PR is aimed at the `main` branch. This means for example that some checks regarding styling or code quality are allowed will not prevent a merge if the PR is aimed at the `dev` branch see below for more details).
   1. Is not proposing a new release or has been approved by at least one
      Maintainer (that is, PRs proposing new releases need to be approved by at
      least one Maintainer).
1. A Maintainer can merge any PR - even if it's not eligible to merge according
   to Rule 4.
1. Any Contributor can Review a PR and Request changes. If a Contributor
   Requests changes they need to provide an explanation what changes should be
   added and justification of their importance. Reviews requesting changes can
   also be used to request more time to review a PR.
1. A Contributor that Requested changes can Dismiss their own review or Approve
   changes added by the Contributor who opened the PR.
1. If the author of a PR and Contributor who provided Review that Requests
   changes cannot find a solution that would lead to the Contributor dismissing
   their review or accepting the changes the Review can be Dismissed with a vote
   or by a Maintainer. Rules governing voting:
   1. A Vote can be triggered by any Contributor, but only after 5 working days
      from the time a Review Requesting changes has been raised and in case a
      Vote has been triggered previously no sooner than 15 working days since
      its conclusion.
   1. Only Contributors can vote, each contributor gets one vote.
   1. A Vote ends after 5 working days or when all Contributors have voted
      (whichever comes first).
   1. A Vote freezes the PR - no new commits or Reviews Requesting changes can
      be added to it while a vote is ongoing. If a commit is accidentally made
      during that period it should be reverted.
   1. The quorum for a Vote is 30% of all Contributors.
   1. The outcome of the Vote is decided based on a simple majority.

### Stable VS latest versions, releases and fixes

Version number follow semantic versioning. 

<!-- add link -->

The `main` branch holds the stable version of the toolbox.

The `dev` branch is where the latest version can be fetched.

Version bumps and new releases are triggered:
- by hotfixes of bug 
- by a merge of the develop branch in the main branch.

A diagram version of the decision-making flow we are aiming for is shown below. ([source](https://blog.axosoft.com/gitflow/))

![git_flow](commenting_images/gitflow_diagram.png "gitflow_diagram")

#### Conditions for merge into `dev`

We can accumulate a certain level of ["technical debt"](https://en.wikipedia.org/wiki/Technical_debt) on the development branch and therefore, we are more flexible as to what can be merged in there.

Conditions:

- All unit tests of code not changed in the PR must still pass.

<!-- that feels bit tautological -->

#### Conditions for merge into `main`

Eventually though this technical debt must be paid back before a new release and merging into the main branch.

Conditions: 
- All unit and integration tests must pass.
- All checks for code style and quality must pass.

## Comments

1. Researchers preparing academic manuscripts describing work that has been
   merged into this repository are strongly encouraged to invite all Maintainers
   as co-authors as a form of appreciation for their work.

<!-- 1. PRs MUST be merged using the "Create a merge commit" option in GitHub (by using
   the "merge pull request" option). This is necessary for our automatic
   changelog generator to do its work reliably. See the [GitHub help page](https://help.github.com/en/articles/about-merge-methods-on-github)
   for information on merge methods. See the changelog generator implementation
   in our [circleci configuration file](./.circleci/config.yml). -->
