# RSDS migration handoff — 2026-10-26

The research domain has moved from ODSP source commit d58fb759416fc458f69ef61ea79b7b1d4164f2c5.

- Original source-domain **390 assets / 4,620,492 bytes** have identical raw Git blob hashes in RSDS.
- The initial import identified six root-level N2 Serengeti receipts that were copied and hash-verified in the runner but accidentally omitted from the initial git stage. The import workflow was corrected to stage those six files. The corrected [import run](https://github.com/zuizui0223/RSDS/actions/runs/40260326289) completed successfully.
- The final full 390-file receipt is at [RSDS_TRANSFER_VERIFICATION_20261026.json](RSDS_TRANSFER_VERIFICATION_20261026.json). The offline checker must pass at the RSDS branch head before this PR may merge.
- The negative/undefined external original GPS solar-radiation response is **preserved**. There is no new post-transfer causal or globally universal solar avoidance claim.
- ODSP remains a historical methodological and provenance record; RSDS becomes the primary research repository. Archived source workflows are not runnable workflows; active RSDS provenance CI is separate.
- The one historical Movebank API response-page tar archive is not one of the 390 transferred text/analysis assets. Its parsed public source catalogues, raw registry manifest, provider source identities and DOIs are included; the remaining binary is documented in the frozen manifest.
