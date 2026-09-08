# Pushing this repository

gh 2.100.0 is installed on zeus. The repository is committed on branch `main`
with no remote. Authentication is the one step that needs you.

## 1. Authenticate

    gh auth login

Choose GitHub.com, then HTTPS, then "Login with a web browser". gh prints a
one-time code; open the URL it gives on any machine and enter it. Nothing is
stored anywhere but zeus.

## 2. Create and push

    cd ~/beacon-repo
    gh repo create beacon-seti --private --source=. --remote=origin --push

Use `--public` instead of `--private` only when you intend it to stay public.

## Before making it public

- The paper says the repository becomes public **on submission**. If that is
  still the plan, create it private and flip it later:
  `gh repo edit --visibility public --accept-visibility-change-consequences`
- No data is committed; no credentials appear in any tracked file (checked).
- Add a LICENSE. MIT or BSD-3 for the code, CC-BY-4.0 for the result JSONs is
  the usual pairing for a paper repository.
- Replace the data-availability line in the paper with the real URL.
- Once public, assume it is cached and forked at once. Treat it as irreversible.
