You are helping write a git commit message for recent changes in this project (DCPower).

$ARGUMENTS

## Step 1 — Get the diff

Run:
```bash
git diff --cached
git diff
git status
```

Collect all staged and unstaged changes.

## Step 2 — Route to Codex to draft the message

Pass the diff to Codex to draft a commit message:

```bash
git diff HEAD | codex exec \
  --ephemeral \
  --sandbox read-only \
  -o /tmp/dcpower_commit.txt \
  "You are writing a git commit message for a biostatistics R project called DCPower — a simulation-based power analysis framework for circadian rhythm differential expression studies.

   Write a commit message following these rules:
   - First line: imperative mood, under 72 chars, no period
   - Blank line
   - Body: bullet points grouped by theme (methodology, figures, paper, bugfix)
   - Focus on WHY, not just what files changed
   - No mention of AI tools or assistants
   
   The diff is provided via stdin."
```

Read `/tmp/dcpower_commit.txt` and display the drafted message.

## Step 3 — Review

Show the drafted commit message and ask:
> "Use this message? (yes / modify / cancel)"

If the user says modify, ask what to change and revise accordingly.

## Step 4 — Commit

Once approved, run:
```bash
git add -u
git commit -m "<approved message>"
```

Do not push unless the user explicitly asks.
