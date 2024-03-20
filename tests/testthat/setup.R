# Setting up the options for a mock MultiStateQTLExperiment.

set.seed(42)
nQTL <- 100
nStates <- 10

sumstats <- mockSummaryStats(nStates=nStates, nQTL=nQTL, names=TRUE)
qtle <- QTLExperiment(
    assay=list(
        betas=sumstats$betas,
        errors=sumstats$errors,
        pvalues=sumstats$pvalues,
        lfsrs=sumstats$pvalues))


sumstats_noNames <- mockSummaryStats(nStates=nStates, nQTL=nQTL, names=FALSE)
state_ids <- colnames(sumstats$betas)
feature_ids <- gsub("\\|.*", "", row.names(sumstats$betas))
variant_ids <- gsub(".*\\|", "", row.names(sumstats$betas))

## Simulated object with LFSRS -------------------------------------------------

sim <- qtleSimulate(
    nstates=10, nfeatures=100, ntests=1000,
    global=0.2, multi=0.4, unique=0.2, k=2)
sim <- callSignificance(sim, mode="simple", assay="lfsrs",
                        thresh=0.0001, second.thresh=0.0002)
sim_sig <- getSignificant(sim)
sim_top <- getTopHits(sim_sig, assay="lfsrs", mode="state")


## GTEX data -------------------------------------------------------------------

input_path <- system.file("extdata", package="multistateQTL")
state <- c("lung", "thyroid", "spleen", "blood")

input <- data.frame(
    state=state,
    path=paste0(input_path, "/GTEx_tx_", state, ".tsv"))

gtex <- sumstats2qtle(
    input,
    feature_id="molecular_trait_id",
    variant_id="rsid",
    betas="beta",
    errors="se",
    pvalues="pvalue",
    verbose=TRUE)

