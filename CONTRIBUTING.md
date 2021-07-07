# How to contribute
We are glad you are reading this.  Thanks for taking the time to contribute!

## Resources
* LIS website: https://lis.gsfc.nasa.gov/
* [How to read LISF documentation](https://github.com/NASA-LIS/LISF/tree/master/docs)
* [LDT Users' Guide](https://github.com/NASA-LIS/LISF/blob/master/docs/LDT_users_guide/LDT_usersguide.adoc)
* [LIS Users' Guide](https://github.com/NASA-LIS/LISF/blob/master/docs/LIS_users_guide/LIS_usersguide.adoc)
* [LVT Users' Guide](https://github.com/NASA-LIS/LISF/blob/master/docs/LVT_users_guide/LVT_usersguide.adoc)
* [Working with GitHub](https://github.com/NASA-LIS/LISF/blob/master/docs/working_with_github/working_with_github.adoc) for LISF
* [How to create a testcase](https://github.com/NASA-LIS/LISF/blob/master/docs/howto_create_lis_testcases/howto_create_lis_testcases.adoc) for LISF

### I am having a problem running LDT/LIS/LVT and would like to ask a question
First, please see the Users' Guides listed in the Resources above to ensure that your config files and other settings are correct.

Next, please check all log files and error messages in the output of your run.  These messages typically inform the user what should be fixed.

If you are still having a problem with your run, please first go to the ["Discussions"](https://github.com/NASA-LIS/LISF/discussions) tab.  On that page,
"Search" in the box in the upper left for previous Discussions, which may already have an answer to your question.

If none of these resources answer your particular question, please open a "New Discussion" using the green box in the upper left.

**NOTE:** Please don't open an "Issue" to ask your question.  You will get faster results by first searching the Guides and Discussions, and then by opening
a "New Discussion", if necessary.

The LIS team will answer your questions on the "Discussion" page, and help to fix your problem.

### After consultation with the LIS team, it has been determined there is an issue that should be fixed
Only after the LIS team has commented on your "Discussion", if both you and the LIS team agree that there is indeed an "Issue" in the code, please go to the
["Issues"](https://github.com/NASA-LIS/LISF/issues) tab and open a "New Issue".

### I would like to contribute either a new feature or a bug fix to the LISF repository
Great!  We welcome contributions (new features, enhancements, bug fixes, etc.).

Please read the [Working with GitHub](https://github.com/NASA-LIS/LISF/blob/master/docs/working_with_github/working_with_github.adoc) page for best practices
on how to submit a "Pull Request" for review by the LIS team before it is merged into the repository.

Generally, the branch for the Pull Request is determined by the below table.  However, the LIS team will help determine the appropriate branch for the
Pull Request within the conversation about the "Issue".  If you are not sure which branch to use, please ask on the "Issue" conversation.

| Branch | Description |
| ------ | ------------- |
| master | New features/enhancements **and** bug fixes for code that is _not_ in the latest release |
| support | Bug fixes for code that _is_ in the latest release |

Please see https://chris.beams.io/posts/git-commit/ for guidance on writing a good commit message.

**NOTE:** Be sure to end your commit message with a line as "Resolves: #????" - where ???? is the "Issue" number associated with this Pull Request.  This line
will automatically close the "Issue" once the code is merged into the repository.

If you are contributing a new feature, we will need a testcase so we can reproduce the result before the pull request is merged.  Please read the
[How to create a testcase](https://github.com/NASA-LIS/LISF/blob/master/docs/howto_create_lis_testcases/howto_create_lis_testcases.adoc) document.
