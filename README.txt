There are four folders. /Code produces all the plots used in the LaTeX. It contains a library, that is a modification of the StatOrdPattHxC library, called "myLibrary". All data can be found under /data. The bib and xlsx folder is for the bibliometric analysis. The TXT folder contains data for ElNinoHolcene.R. Powerlaw.R uses computer generated data. Weather.R uses Weather.csv under CSV. The plots produced by the code can be found under /Figures. The figures in the report responds to the following plots

Figure 2.1: powerlaw/histogram.pdf - produced by line 20-25 in powerlaw.R

Figure 2.2: ElNino/ArticleEntropyPlot.jpg - produced by ElNinoHolcene.R

Figure 2.3: ElNino/Entropy.pdf - produced by ElNinoHolcene.R

Figure 4.1: is all plots in powerlaw except histogram.pdf - produced by powerlaw.R 

Figure 4.2: Weather/tiesTable.pdf - produced by line 37-73 in Weather.R

Figure 4.3: Weather/entropyTable.pdf - produced by line 78-197 in Weather.R

Figure 4.4: Weather/noiseStochasticTheoretical.pdf  - produced by line 248-325 in Weather.R

Figure 4.5: Weather/random_vs_theoreticalSplit.pdf  - produced by line 248-325 in Weather.R

Figure 4.6: Weather/constantWithWhiteNoiseStochasticTheoretical.pdf - produced by line 329-368 in Weather.R

Figure 4.7: Weather/random_vs_theoreticalSplitWhiteNoise.pdf - produced by line 329-368 in Weather.R

Figure 4.8: Weather/confidenceIntervalPlot.pdf  - produced by line 201-243 in Weather.R

Figure 4.9: Weather/pValuesTheoretical,10=iterations,Sorted.pdf - produced by line 406-441 in Weather.R

Figure 4.10: Weather/pValuesTheoretical.pdf - produced by line 373-401 in Weather.R


Text/articles contains two jabRefs files with referenced. Reference.bib is a raw list of articles pulled for bibliometric analysis, where the ones in data/bib, where the used ones. In the data/bib/BibliomatrixAnalysis.bib each article has a comment saying "Article *number*", where *numbeer* corresponds to its evaluation in the xlsx/reference.xlsx. RefrenceFinal.bib is all references used in the report itself. The figures folders are copied into Text, so they can be included in the report.