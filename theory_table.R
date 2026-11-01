# Combine the three theorem-matched experiments into main Table 2.
# Run after the individual study scripts have rebuilt their result CSVs.
out <- Sys.getenv("DFT_RESULTS", "results")
training <- read.csv(file.path(out, "training_summary.csv"))
power <- read.csv(file.path(out, "local_power_results.csv"))
delay <- read.csv(file.path(out, "delay_results.csv"))
fmt <- function(x) ifelse(is.na(x), "---", sprintf("%.2f", x))
row <- function(label, m, values) {
    paste0(paste(c(label, format(m, scientific=FALSE), fmt(values)), collapse=" & "), " \\\\")
}
headrow <- function(title) paste0("\\multicolumn{6}{l}{\\itshape ",title,"} \\\\")
lines <- c("\\begin{tabular}{@{}lr rrrr@{}}\\toprule",
   "Design & $m$ & \\multicolumn{2}{c}{RSMS} & \\multicolumn{2}{c}{SSMS} \\\\",
   "& & Panel & Ref. & Panel & Ref. \\\\ \\midrule",
   headrow("A. Training contamination: alarms by the horizon (\\%)"))
paired <- function(method, what, dat, id, col) {
    # Each lookup must identify exactly one design cell.
    index <- which(dat$method == method & dat[[id]] == what)
    stopifnot(length(index) == 1L)
    dat[[col]][index]
}
for (profile in c("flat","front-cycle")) {
    conditions <- if (profile == "flat") c("clean","subcritical","boundary","diffuse") else c("clean","boundary")
    labels <- c(clean="clean",subcritical="sparse, $m^{1/4}$",boundary="sparse, $\\sqrt m$",diffuse="diffuse")
    for (condition in conditions) {
        d <- training[training$profile==profile & training$m==1600 & training$design==condition,]
        values <- unlist(lapply(c("RSMS","SSMS"), function(method)
            vapply(c("actual","reference"),function(z) paired(method,z,d,"measure","estimate"),numeric(1))))
        lines <- c(lines,row(paste(if(profile=="flat") "Flat," else "Front-cycle,",labels[condition]),1600,100*values))
    }
}
lines <- c(lines,"\\midrule",headrow("B. Loading change: first alarm after the break (\\%)"))
for (i in seq_len(5)) {
    design <- c("local_weak","local_strong",rep("fixed",3))[i]
    size <- c(720,720,80,240,720)[i]
    label <- c("Local, $h=0.7$","Local, $h=1.6$",rep("Fixed, $\\delta=0.13$",3))[i]
    d <- power[power$m==size & power$scenario==design,]
    values <- unlist(lapply(c("RSMS","SSMS"),function(method) {
        a <- paired(method,"fitted_panel",d,"model","post_break")
        b <- if(design=="fixed") NA_real_ else paired(method,"shifted_reference",d,"model","post_break")
        c(a,b)
    }))
    lines <- c(lines,row(label,size,100*values))
}
lines <- c(lines,"\\midrule",headrow("C. Detection delay: median on the predicted time scale"))
for(onset in c("Jump","Linear")) for(size in c(128,8192)) {
    d <- delay[delay$onset==onset,]
    values <- unlist(lapply(c("RSMS","SSMS"),function(method)
        vapply(c(as.character(size),"Limit"),function(z) paired(method,z,d,"grid","median"),numeric(1))))
    lines <- c(lines,row(onset,size,values))
}
writeLines(c(lines,"\\bottomrule\\end{tabular}"),file.path(out,"theory_table.tex"))
