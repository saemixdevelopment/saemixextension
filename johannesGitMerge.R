# To change stuff in master and carry them into multiOutcome
# work in master branch, then

git checkout multiOutcome
git merge master


# Creating a tag

git tag -d v3.0
git tag v3.0 c3f31db
git push --tags

git tag -d v3.1
git tag v3.1 c65318c
git push --tags
