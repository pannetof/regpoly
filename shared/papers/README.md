# Reference PDFs for the published-generators library

This directory holds PDFs of the academic papers cited by entries in
`docs/library/*.yaml` via their `reference.pdf` field.  Files here are
served by the web app under `/papers/<filename>`.

## Redistribution policy

**Commit a PDF here only when the publisher's license explicitly
permits redistribution.**  In practice:

- ✅ arXiv preprints.
- ✅ Final author-page PDFs where the author owns the copyright and
  has granted distribution rights (or the work is in the public
  domain / under a permissive Creative Commons license).
- ✅ ACM Author-Izer links that resolve to publisher-hosted PDFs are
  fine to **link to**; the file is not copied into this directory.
- ❌ Publisher-typeset PDFs (Math. Comp., ACM TOMS, …) that remain
  under exclusive publisher copyright.  Leave `reference.pdf` empty
  in the YAML and rely on `reference.doi` instead.

When in doubt, check the publisher's green open-access policy on
[SHERPA Romeo](https://v2.sherpa.ac.uk/romeo/), drop the `pdf:`
field, and let the DOI carry the reader to the official copy.
