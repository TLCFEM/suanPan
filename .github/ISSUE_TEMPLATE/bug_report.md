---
name: Bug Report
about: Create a report to help us improve
title: "[BUG]"
labels: 'bug'
assignees: ''
body:
  - type: input
    attribute:
      label: Version
      description: Indicate version tags or 'dev' if using the head of 'dev' branch.
      validation:
        required: true
  - type: textarea
    attribute:
      label: Describe the bug
      placeholder: |
        A clear and concise description of what the bug is.
        You can also upload a MWE to reproduce the behaviour.
      validation:
        required: true
  - type: textarea
    attribute:
      label: Additional context
      placeholder: |
        Add any other context about the problem here.
        If there are any references, please also provide them.
---
