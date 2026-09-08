# Pushing this repository

Prepared and committed locally; it has no remote. Two routes.

## With the GitHub CLI

    sudo apt install gh          # not currently installed on zeus
    gh auth login
    cd ~/beacon-repo
    gh repo create beacon-seti --private  --source=. --remote=origin --push
    #                          ^^^^^^^^^ or --public

## Without it

Create an empty repository on github.com (do not add a README, licence or
.gitignore -- this repo has its own), then:

    cd ~/beacon-repo
    git remote add origin git@github.com:<user>/beacon-seti.git
    git push -u origin main

## Before making it public

- The paper says the repository becomes public **on submission**. If that is
  still the plan, create it private and flip it later.
- No data is committed and no credentials appear in tracked files (checked).
- Add a LICENSE. For code accompanying a paper, MIT or BSD-3 is usual; for the
  result JSONs, CC-BY-4.0 is the common choice.
- Once public, assume it is cached and forked immediately; treat it as
  irreversible.
