# PathoSV
PathoSV: A transcriptome-aware tool for assessing Structural Variant pathogenicity.


PathoSV estimates Structural Variant (SV) pathogenicity using a 'truncated ratio' that quantifies predicted gene transcription disruption in a specific tissue relative to its TPM. Input SV coordinates and tissue; PathoSV outputs the ratio (flagging > 0.25 as potentially pathogenic), SV population frequency, affected gene annotations (OMIM, GO), and offers deepseek interpretation for relating SV impact to a specified disease.
