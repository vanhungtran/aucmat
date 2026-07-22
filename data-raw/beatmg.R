# ==============================================================================
# data-raw/beatmg.R
#
# Builds the package datasets `beatmg_baseline` and `beatmg_ontreatment` from
# the raw BeatMG (NeuroNEXT NN102) trial exports in data-raw/:
#   - MG_NPX.csv: Olink Explore NPX export (Inflammation + Inflammation_II
#     panels), long format, 149 sample-timepoint records x 738 assays
#   - Static_Demographics_BeatMG.csv: baseline demographics and randomized
#     arm for 52 trial participants
#
# BeatMG was a randomized, placebo-controlled trial of rituximab (anti-CD20,
# B-cell depleting) for myasthenia gravis. The cleaned datasets provide a
# real-data worked example for aucmat(): does the randomized arm show up in
# the plasma proteome, and if so, when and via which proteins?
#
# Re-run this script (`source("data-raw/beatmg.R")`) any time the raw
# exports change; it regenerates data/beatmg_baseline.rda and
# data/beatmg_ontreatment.rda.
# ==============================================================================

npx <- read.csv("data-raw/MG_NPX.csv", stringsAsFactors = FALSE)
demo <- read.csv("data-raw/Static_Demographics_BeatMG.csv", stringsAsFactors = FALSE)

# ---- QC filtering ----
# Sample_Type == "SAMPLE": drop Olink inter-plate CONTROL samples
# QC_Warning == "PASS": drop sample-level QC failures
# Assay_Warning != "EXCLUDED": drop assay-level exclusions (WARN-flagged
#   assays are kept, per Olink's own convention, but noted in documentation)
npx_qc <- subset(npx, Sample_Type == "SAMPLE" & QC_Warning == "PASS" &
                        Assay_Warning != "EXCLUDED")

# ---- Outcome: randomized arm ----
demo_y <- demo[, c("SubjectID", "trt_desc_subjectstatus")]
names(demo_y) <- c("SampleID", "arm")

# ---- Build a wide protein matrix (rows = samples) for one timepoint ----
# One row per SampleID x Assay is expected; aggregate() guards against
# accidental duplicate NPX readings for the same sample/assay.
build_matrix <- function(npx_qc, demo_y, tp) {
  d <- subset(npx_qc, Timepoint == tp)
  d <- aggregate(NPX ~ SampleID + Assay, data = d, FUN = mean)
  wide <- reshape(d, idvar = "SampleID", timevar = "Assay", direction = "wide")
  names(wide) <- sub("^NPX\\.", "", names(wide))
  merged <- merge(wide, demo_y, by = "SampleID")
  merged <- merged[!is.na(merged$arm) & merged$arm %in% c("Rituximab", "Placebo"), ]
  merged$arm <- factor(merged$arm, levels = c("Placebo", "Rituximab"))
  rownames(merged) <- NULL
  merged
}

# Timepoint 1 = baseline (pre-treatment), n = 48 after QC + arm matching.
# Timepoint 3 = latest well-powered on-treatment collection, n = 44.
# Timepoint 4 has only 5 samples after QC -- too sparse for AUC inference,
# excluded from the shipped datasets.
beatmg_baseline    <- build_matrix(npx_qc, demo_y, tp = 1)
beatmg_ontreatment <- build_matrix(npx_qc, demo_y, tp = 3)

save(beatmg_baseline, file = "data/beatmg_baseline.rda", compress = "xz")
save(beatmg_ontreatment, file = "data/beatmg_ontreatment.rda", compress = "xz")
