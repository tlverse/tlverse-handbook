#!/bin/bash

# test PAT and git branch
set -e
[ -z "${GITHUB_PAT}" ] && exit 0
[ "${TRAVIS_BRANCH}" != "master" ] && exit 0

# configure credentials
git config --global user.email "nh@nimahejazi.org"
git config --global user.name "Nima Hejazi"
git config --global http.postBuffer 100000000

# clone the gh-pages branch and copy updated book contents to it
git clone -b gh-pages \
  https://${GITHUB_PAT}@github.com/${TRAVIS_REPO_SLUG}.git \
  book-output
cd book-output
cp -r ../_book/* ./

# stage, commit, push copied files to branch gh-pages
if [ "${TRAVIS_PULL_REQUEST}" = "false" ]
then
  COMMIT_MESSAGE="Update book via ${TRAVIS_COMMIT}."
else
  COMMIT_MESSAGE="Update book via PR #${TRAVIS_PULL_REQUEST} ($TRAVIS_COMMIT)."
fi
git add --all *
git commit -m "${COMMIT_MESSAGE}"
git push origin gh-pages
