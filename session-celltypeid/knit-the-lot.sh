#!/bin/sh
cd $HOME/git/excelerate-scRNAseq/session-celltypeid
 R -e "rmarkdown::render('celltypeid.Rmd')"
