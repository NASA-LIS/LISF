# Contribution Guidelines

Welcome to the LISF repository! Thanks for taking the time to contribute! Please read this entire document *before* contributing.

## Table of Contents

* [Support](#support)
* [Issues](#issues)
* [Pull Requests](#pull-requests)
    * [Branches](#branches)
    * [Testcases](#testcases)
    * [Documentation](#documentation)
* [Resources](#resources)

## Support

First, please review the Users' Guides listed in the [Resources](#resources) below to ensure that your config files and other settings are correct.

Next, please check **all** log files and standard error outputs for error messages. These messages typically inform the user what should be fixed and they may not be present in *every* log file.

If you are still unable to resolve the issue, please visit our [Discussions](https://github.com/NASA-LIS/LISF/discussions) forum and use the search feature to check for a solution. It is possible that another user has encountered the same issue and we have already provided a solution.

If your particular question has not been addressed yet, please open a "New Discussion" describing your problem and attach your configuration file (e.g., `lis.config`) and a log file (e.g., `lislog.0000`). Attach a log file containing an error message or, if no error message exists, attach the longest log file to start.

**NOTE:** Please don't open an "Issue" to ask your question. You will get faster results by first searching the Users' Guides and Discussions, and then by opening a "New Discussion", if necessary.

## Issues

Please **do not** open a new Issue until the LIS team has commented on your Discussion post and confirmed that there is indeed an Issue in the code. To open an Issue either:

* Navigate to the relevant Discussion thread and click "Create issue from discussion" on the right side of the screen.
* Navigate to the [Issues](https://github.com/NASA-LIS/LISF/issues) tab and click the "New Issue" button. Select the relevant Issue template and provide the information described in the text box. Link the Issue to the original Discussion thread by including `#123` in the Issue description (where `123` is the Discussion number).

If you plan to contribute a fix yourself, please [assign yourself to the Issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/assigning-issues-and-pull-requests-to-other-github-users).

## Pull Requests

We welcome Pull Requests to LISF that add new features, enhance existing capabilities, and fix bugs! **BUT**, we ask that you first open an Issue (and, before that, a Discussion) so that we are aware of your intentions and have an opportunity to provide our feedback.

Before beginning to work on your contribution, please review our [Working with GitHub](https://github.com/NASA-LIS/LISF/blob/master/docs/working_with_github/working_with_github.adoc) document for more information about our GitHub practices and a step-by-step guide to submitting a Pull Request.

**NOTE:** Please [link the relevant Issue to your Pull Request](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue) so that it is automatically closed once the code is merged into the repository.

### Branches

Generally, the target branch for a given Pull Request can be determined by the below table. However, the LIS team will help determine the appropriate branch for the Pull Request within the Issue or Discussion conversation. If you are not sure which branch to use, please ask! Again, please review our [Working with GitHub](https://github.com/NASA-LIS/LISF/blob/master/docs/working_with_github/working_with_github.adoc) document before beginning your work to ensure you are building off of the correct branch.

| Branch  | Description                                                                              |
| ------- | -----------------------------------------------------------------------------------------|
| master  | New features/enhancements **and** bug fixes for code that is _not_ in the latest release |
| support | Bug fixes for code that _is_ in the latest release                                       |

### Commit Messages

Please use descriptive commit messages. See [this post](https://chris.beams.io/posts/git-commit/) for guidance on writing a good commit message.

### Testcases

If you are contributing a new feature, we will need a testcase so we can reproduce the result before the pull request is merged. Please read the [How to Create a Testcase](https://github.com/NASA-LIS/LISF/blob/master/docs/howto_create_lis_testcases/howto_create_lis_testcases.adoc) document.

**IMPORTANT:** Please **do not** commit your testcase into the Pull Request/repository.  We will make arrangements to obtain your input/output files and configuration files.

### Documentation

If you are adding or modifying config entries, please update the appropriate config documentation file (e.g., `lis/configs/lis.config.adoc`) to reflect these changes. Our Users' Guides are written in the [AsciiDoc markup language](https://docs.asciidoctor.org/asciidoc/latest/). If you require assistance with documentation, please ask in the Issue thread.

## Resources

* [LIS Website](https://lis.gsfc.nasa.gov/)
* [LISF Documentation](https://github.com/NASA-LIS/LISF/tree/master/docs)
* [Working with GitHub](https://github.com/NASA-LIS/LISF/blob/master/docs/working_with_github/working_with_github.adoc)
* [How to Create a Testcase](https://github.com/NASA-LIS/LISF/blob/master/docs/howto_create_lis_testcases/howto_create_lis_testcases.adoc)

## License

By contributing to the LISF repository, unless you explicitly state otherwise, you agree that your contributions will be licensed under the [LISF License](https://github.com/NASA-LIS/LISF/blob/master/LICENSE.txt), which is under Apache License 2.0.
