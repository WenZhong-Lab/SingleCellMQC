# v0.7.0 (25)

-   Fixed bugs in `Read10XH5Data()`.
-   Fixed `RunLQ_MAD()`.
-   Optimized sample, cell, feature and batch method functions.
-   Optimized `RunReport()` and improved HTML report content.
-   Add `RunBenchmarkDoublet` and `PlotBenchmarkDoublet` for benchmarking of doublet detection methods.
-   Optimized `FindCommonPCTOutlier` function.
-   Delete some functions.

# v0.6.1

-   The default auto-generated sample identifier in `Read10XData()` has been updated from Sample_X to SampleX\_ when no sample name is provided.

# v0.6.0

-   Fixed scCDC correction in `RunDbt_ADT()`.
-   Fixed `PlotADTCutoff()` failure when processing BPCells matrices.
-   Fixed bugs in `FindContaminationFeature()`.
-   Fixed bugs in `RunCorrection()`.
-   Fixed issues with factor variables affecting certain functions.
-   Added `RunVarExplainedPerFeature()` for per-feature variance explained analysis.
-   Optimized `RunReport()` and improved HTML report content. Restructured batch QC content into three levels: sample, cell, and feature.
-   Enhanced overall function efficiency and speed.
-   Updated documentation for clarity and completeness.

# v0.4.0

-   Fixed `RunDbt_DoubletFinder()` compatibility with DoubletFinder (≥ v2.0.6).\
-   Fixed scCDC correction in `RunDbt_ADT()` .\
-   Fixed `PlotADTCutoff()` failure when processing BPCells matrices.\
-   Fixed label misalignment in `FindSampleMetricsWarning()` output plots.
