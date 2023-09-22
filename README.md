# GENEC
Welcome to the Geneva stellar evolution code.

To copy the code to your computer, run:

git clone git@github.com:GESEG/GENEC.git

git tutorials can be found here:
https://www.atlassian.com/git/tutorials/

The master branch, which you download with the command above will stay the same until the next release.
The "develop" branch is the main development branch.
 
Before making any changes, please create your own branch off the development branch using this command:

git checkout -b feature/my_feature develop

[meaning that you create a new branch called "feature/my_feature" starting from the "develop" branch]
[my_feature should be replaced by the topic of the feature you are working on]

When you want to save changes to the git repo, use this command (after running git status, add, commit):

git push origin feature/my_feature
[this pushes the changes you made to your (local) branch to the remote (origin) repository]

