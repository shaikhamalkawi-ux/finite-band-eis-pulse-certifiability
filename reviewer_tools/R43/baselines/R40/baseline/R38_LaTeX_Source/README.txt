Project 52 — Batteries | IEEE OJVT | R38 EP1 editorial-production source

Science: R38 full-domain conditional external-cell verification (locked).
This is an editorial-only production pass of the approved R38 manuscript.

Main manuscript: main.tex
Supporting Information: supplementary.tex (includes S1–S13 through the appendix TeX files)
Cover letter: cover_letter.tex
Graphics: figures/ (eight original PNG images preserved)

Compile from this folder with a current TeXLive installation containing IEEEtran:
    latexmk -pdf -halt-on-error -interaction=nonstopmode main.tex
    latexmk -pdf -halt-on-error -interaction=nonstopmode supplementary.tex
    latexmk -pdf -halt-on-error -interaction=nonstopmode cover_letter.tex

Do not add the raw RWTH or Bayreuth third-party datasets to this source archive.
The separate Computational_Supplement.zip supplies saved verifiers and inputs;
its external source archives are not redistributed.
Original scientific R38 source and change comparisons are retained internally.
No journal-portal submission or IEEE online PDF checker execution is claimed.
