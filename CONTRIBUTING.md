# How to contribute

Thank you for considering contributing to the development of FATES. There are a few guidelines that we need contributions to follow.

## Read and Understand the FATES license:

https://github.com/NGEET/fates/blob/master/LICENSE.txt

## Please read and follow the code of conduct:

https://github.com/NGEET/fates/blob/master/CODE_OF_CONDUCT.md

## Getting Started

Those who wish to contribute code to FATES must have those changes integrated through the developer repository NGEET/fates.  Changes that make it to public releases must go through this repository first, as well.  Please refer to the [developer section](hhttps://fates-users-guide.readthedocs.io/en/latest/developer/developer-guide.html) of the [User's Guide](https://fates-users-guide.readthedocs.io/en/latest/index.html) for more details.  Here are some basic first steps:

* All developers should create a fork of the NGEET/fates repository into their personal space on github
* Follow the [developer work-flow](https://fates-users-guide.readthedocs.io/en/latest/developer/FATES-Development-Workflow.html) 
* Each set of changes should have its own [feature branch](https://fates-users-guide.readthedocs.io/en/latest/developer/Feature-Branch-Naming-Convention.html) that encapsulates your desired changes.
* The work-flow will lead you eventually to submit a Pull-Request to NGEET/fates:master, please follow the template in the Pull Request and communicate as best you can if you are unsure how to fill out the text
* It is best to create an issue to describe the work you are undertaking prior to starting.  This helps the community sync with your efforts, prevents duplication of efforts, and science is not done in a vaccuum!
* Expect peers to interact, help, discuss and eventually approve your submission (pull-request)

## Use of AI Agents

AI agents should not be used to submit pull requests, issues, discussion items or other features available on the developer platforms on which the FATES model is officially hosted (i.e. Github), except on a case-by-case basis.  Approval to do so should be sought from a core member of the FATES development team.  Unapproved submissions will be closed.  Repeated unapproved submissions may culminate in the blocking of the submitting user.

Contributors that use generative AI to help develop code are asked to describe its usage when submitting a pull request. Specifically, we request information such as the platform use and the extent of its use in producing the pull request (for example, specific commands versus chunks of code written with generative AI).

## Joining Discussions and Meetings

In addition to the github discussions, we hold a roughly biweekly call, which covers both scientific and technical issues related to FATES. We use a google group to organize, schedule, and discuss these calls.  Emails to the list are moderated and we try to be pretty ruthless about preventing anything other than this topic from appearing on it.  To join, apply to the group here: https://groups.google.com/forum/#!forum/fates_model

## Things to Remember

* Make commits in logical units (i.e. group changes)
* Changes that are submitted should be limited to 1 single feature (i.e. don't submit changes to the radiation code and the nutrient cycle simultaneous, pick one thing)
* Check for unnecessary whitespace with `git diff --check` before committing
* We have no standard protocol for commit messages, but try to make them meaningful, concise and succinct.
* You will most likely have to test (see workflow above), see the [testing protocols](https://fates-users-guide.readthedocs.io/en/latest/developer/Testing-Protocols.html) for more details.

## Coding Practices and Style

Please refer to the FATES [style guide](https://fates-users-guide.readthedocs.io/en/latest/developer/style.html)

## Trivial Changes

If changes are trivial, it's possible testing will not be required. Conversations via the Pull Request will address if tests are not needed

## Documentation

Yes please!  If you are creating new code, fixing existing code, anything.  Please add comments in the code itself.  Please also follow the style guide for comments.  Please update the tech note, and submit a PR to [the fates-docs repo](https://github.com/NGEET/fates-docs/) alongside your FATES PR.  Also, please create and/or modify existing wiki documentation.  You may be asked to add documtation prior to having a pull-request approved.

## Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [GitHub pull request documentation](https://help.github.com/articles/creating-a-pull-request/)
* [FATES Wiki](https://github.com/NGEET/fates/wiki)
