# Reference papers

This directory used to bundle the reference PDFs cited by `data/library/*.yaml`.
They've been removed from the repository (most are publisher-copyrighted,
e.g. ACM/Elsevier journal articles) and are now `.gitignore`d — see
`data/papers/*.pdf` in the top-level `.gitignore`.

Each catalog entry's `doi:` field links out to the paper instead. Drop a
PDF here locally (matching filename referenced historically, e.g.
`wellrng.pdf`) if you want one for personal reference; it won't be tracked
by git and won't ship in the wheel or CLI installs.
