#!/bin/bash
## Convenience script to add, commit and push updates to GitHub 
## Author: Thomas Girke
## Last update: Feb 4, 2021

## (1) Makes sure you are in correct branch
git checkout main  

## (2) Build site and copy rendered pages to docs 
rm -rf public/*
rm -rf docs/*
Rscript -e "blogdown::build_site()"
cp -r public/* docs/

## (3) Commit and push changes
git add -A :/
git commit -am "some edits"
git push -u origin main
echo "Committed/pushed changes to gh-pages branch on GitHub"
