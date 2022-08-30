mkdir /tmp/R_fixes
git checkout multiOutcome
cp R/* /tmp/R_fixes
git checkout master
cp /tmp/R_fixes/*.R R/
git status
git add R/*.R
git commit -m "Johannes saved my gitday"
git push
