#!/bin/bash
## Convenience script to add, commit and push updates to GitHub 
## Author: Thomas Girke
## Last update: Feb 4, 2021

## (1) Makes sure you are in correct branch (currently there is only one)
git checkout main  

## (2) Build site and copy rendered pages to docs 
rm -rf public/* # clean out public
rm -rf docs/* # clean out docs (assumes dir exists)
Rscript -e "blogdown::build_site()" # Build site with blogdown via command-line call
cp -r public/* docs/ # Copy build/rendered pages to docs

## (3) Commit and push changes
git add -A :/
git commit -am "some edits"
git push -u origin main
echo "Committed/pushed changes to GitHub"
