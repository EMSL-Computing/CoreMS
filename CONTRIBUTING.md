# Contributing to CoreMS

Thank you for considering contributing to CoreMS! We appreciate your interest in helping us improve our project. This document outlines the guidelines and steps for contributing to CoreMS.

## Table of Contents

- [Getting Started](#getting-started)
- [Versioning](#versioning)
- [Merge Request Checklist](#merge-request-checklist)
- [Code Style](#code-style)
- [Issue Reporting](#issue-reporting)
- [License](#license)

## Getting Started

To get started with contributing to CoreMS, please follow these steps:

1. Create an issue proposing a fix or expanded functionality and make sure it's substantially different from an existing one.
2. Fork the CoreMS repository. If you are part of the development team, you can forgo a fork and instead make a branch.
3. Install the necessary dependencies. Refer to the [README](./README.md) for detailed installation instructions.
4. Make your changes or additions.
5. Test your changes thoroughly.
6. Commit your changes and push them to your forked repository. Reference your original issue in your commits (i.e. closes #23)
7. Submit a merge request to the main CoreMS repository and select an appropriate reviewer for the changes. Note the merge request checklist below that will be checked before each merge into the master branch. See the merge request checklist

## Versioning

We strive to use semantic versioning. To bump a new version and regenerate documentation, use one of the following make commands (according to version number)  `make major`, `make minor`, or `make patch`.  This should accompany each PiPy release.

## Merge Request Checklist

Before merging *into the master branch*, each of these will be checked by a reviewer.

1. CI/CD pipeline must pass.
2. Each merge request must be accompanied by an appropriate bump in version number, following the major.minor.patch format (semantic versioning). 
    - Major: Incremented when making incompatible API changes.
    - Minor: Incremented when adding new features in a backwards-compatible manner.
    - Patch: Incremented for backwards-compatible bug fixes.
3. Unit tests must be added or updated to cover the changes made.
4. Documentation must be updated and rerendered to reflect any new features or changes.
5. Any relevant issues or pull requests should be referenced in the merge request (i.e. closes #23).

## Code Style

CoreMS follows the [NumPy documentation style guide](https://numpydoc.readthedocs.io/en/latest/format.html). Please ensure that your code adheres to this style to maintain consistency throughout the project.  

Documentation is rendered using the [pdoc package](https://github.com/mitmproxy/pdoc/tree/main).

## Issue Reporting

If you encounter any issues or bugs while using CoreMS, please report them by opening an issue in the issue tracker. Please provide as much detail as possible, including steps to reproduce the issue and any relevant error messages.

## License

By contributing to CoreMS, you agree that your contributions will be licensed as described in the [LICENSE](./LICENSE) file.

We appreciate your contributions and look forward to working with you to improve CoreMS!
