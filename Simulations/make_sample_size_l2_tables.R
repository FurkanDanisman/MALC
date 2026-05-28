## LaTeX summary tables for sample-size simulation L2 errors.
##
## Four tables are generated:
##   1. Group 1 distributions, median (IQR)
##   2. Group 1 distributions, mean (SD)
##   3. Group 2 distributions, median (IQR)
##   4. Group 2 distributions, mean (SD)
##
## The smallest median and smallest mean within each distribution/sample-size
## row are bolded separately.

results_dir <- Sys.getenv(
  "DENSOLOG_RESULTS_DIR",
  unset = "~/Desktop/RA_Documents/RA-DOCUMENTS/PDF/Method"
)
results_dir <- path.expand(results_dir)

output_file <- file.path("Simulations", "sample_size_l2_summary_tables.tex")

sample_cols <- c("X1e2", "X1e3", "X1e4", "X1e5", "X1e6")
sample_labels <- c("$10^2$", "$10^3$", "$10^4$", "$10^5$", "$10^6$")

methods <- list(
  list(label = "MALC", file_suffix = "method1", repo_dir = "DensOLog"),
  list(label = "BK2002", file_suffix = "method2", repo_dir = "BK2002"),
  list(label = "BinnedNP", file_suffix = "method3", repo_dir = "binnednp"),
  list(label = "KernSmooth", file_suffix = "method4", repo_dir = "KernSmooth")
)

table1_distributions <- list(
  list(label = "normal", prefix = "normal"),
  list(label = "Beta", prefix = "beta"),
  list(label = "Gamma", prefix = "gamma"),
  list(label = "logistic", prefix = "logistic"),
  list(label = "Student's $t$", prefix = "t")
)

table2_distributions <- list(
  list(label = "Laplace", prefix = "laplace"),
  list(label = "chi-square", prefix = "chisqr"),
  list(label = "log-normal", prefix = "lnorm"),
  list(label = "Weibull", prefix = "weibull"),
  list(label = "Pareto", prefix = "pareto")
)

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("_", "\\\\_", x)
  x
}

format_num <- function(x, digits = 4) {
  formatC(x, format = "f", digits = digits)
}

bold_if <- function(x, flag) {
  if (isTRUE(flag)) {
    paste0("\\textbf{", x, "}")
  } else {
    x
  }
}

find_result_file <- function(prefix, method) {
  flat_file <- file.path(
    results_dir,
    paste0(prefix, "_sample_sim_data_", method$file_suffix, ".csv")
  )
  if (file.exists(flat_file)) {
    return(flat_file)
  }

  if (!is.na(method$repo_dir)) {
    repo_file <- file.path(
      "Simulations", "Sim_Results", method$repo_dir,
      paste0(prefix, "_sample_sim_data_", method$file_suffix, ".csv")
    )
    if (file.exists(repo_file)) {
      return(repo_file)
    }
  }

  stop(
    "Could not find sample-size result file for prefix='", prefix,
    "' and method='", method$label, "'. Looked in ", flat_file,
    call. = FALSE
  )
}

read_distribution_results <- function(dist) {
  out <- vector("list", length(methods))
  names(out) <- vapply(methods, `[[`, character(1), "label")

  for (m in seq_along(methods)) {
    f <- find_result_file(dist$prefix, methods[[m]])
    dat <- read.csv(f, check.names = TRUE)

    missing_cols <- setdiff(sample_cols, names(dat))
    if (length(missing_cols) > 0L) {
      stop(
        "Missing columns in ", f, ": ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    out[[m]] <- dat[sample_cols]
  }

  out
}

summarize_one <- function(x) {
  c(
    median = stats::median(x, na.rm = TRUE),
    q1 = stats::quantile(x, 0.25, na.rm = TRUE, names = FALSE),
    q3 = stats::quantile(x, 0.75, na.rm = TRUE, names = FALSE),
    iqr = stats::IQR(x, na.rm = TRUE),
    mean = mean(x, na.rm = TRUE),
    sd = stats::sd(x, na.rm = TRUE)
  )
}

make_median_cell <- function(s, bold = FALSE) {
  paste0(bold_if(format_num(s[["median"]]), bold), " (", format_num(s[["iqr"]]), ")")
}

make_mean_cell <- function(s, bold = FALSE) {
  paste0(bold_if(format_num(s[["mean"]]), bold), " (", format_num(s[["sd"]]), ")")
}

make_table_rows <- function(distributions, stat = c("median", "mean")) {
  stat <- match.arg(stat)
  rows <- character()

  for (dist_index in seq_along(distributions)) {
    dist <- distributions[[dist_index]]
    dat <- read_distribution_results(dist)

    for (j in seq_along(sample_cols)) {
      stats_by_method <- lapply(dat, function(d) summarize_one(d[[sample_cols[j]]]))
      medians <- vapply(stats_by_method, function(s) s[["median"]], numeric(1))
      means <- vapply(stats_by_method, function(s) s[["mean"]], numeric(1))

      best_median <- which(medians == min(medians, na.rm = TRUE))
      best_mean <- which(means == min(means, na.rm = TRUE))

      if (stat == "median") {
        cells <- vapply(
          seq_along(stats_by_method),
          function(i) make_median_cell(
            stats_by_method[[i]],
            bold = i %in% best_median
          ),
          character(1)
        )
      } else {
        cells <- vapply(
          seq_along(stats_by_method),
          function(i) make_mean_cell(
            stats_by_method[[i]],
            bold = i %in% best_mean
          ),
          character(1)
        )
      }

      dist_label <- if (j == 1L) dist$label else ""
      row <- paste(
        c(dist_label, sample_labels[j], cells),
        collapse = " & "
      )
      if (j == length(sample_cols)) {
        row <- paste0(row, " \\\\[0.6ex]")
      } else {
        row <- paste0(row, " \\\\")
      }

      rows <- c(rows, row)
    }

    if (dist_index < length(distributions)) {
      rows <- c(rows, "\\midrule")
    }
  }

  rows
}

print_latex_summary_table <- function(distributions, stat = c("median", "mean"),
                                      caption, label) {
  stat <- match.arg(stat)
  header <- paste(
    c("Distribution", "$n$", vapply(methods, `[[`, character(1), "label")),
    collapse = " & "
  )

  rows <- make_table_rows(distributions, stat = stat)

  c(
    "\\begin{table}[p]",
    "\\centering",
    "\\scriptsize",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\renewcommand{\\arraystretch}{1.15}",
    "\\resizebox{\\textwidth}{!}{%",
    "\\begin{tabular}{llcccc}",
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabular}%",
    "}",
    paste0("\\caption{", caption, " Bold entries mark the smallest value across methods for that distribution and sample size.}"),
    paste0("\\label{", label, "}"),
    "\\end{table}",
    ""
  )
}

latex_out <- c(
  "% Generated by Simulations/make_sample_size_l2_tables.R",
  "% Requires booktabs; uses \\resizebox from graphicx.",
  print_latex_summary_table(
    distributions = table1_distributions,
    stat = "median",
    caption = "Median (IQR) sample-size simulation L2 errors for the normal, beta, gamma, logistic, and Student's $t$ distributions.",
    label = "tab:l2_sample_size_group1_median"
  ),
  print_latex_summary_table(
    distributions = table1_distributions,
    stat = "mean",
    caption = "Mean (SD) sample-size simulation L2 errors for the normal, beta, gamma, logistic, and Student's $t$ distributions.",
    label = "tab:l2_sample_size_group1_mean"
  ),
  print_latex_summary_table(
    distributions = table2_distributions,
    stat = "median",
    caption = "Median (IQR) sample-size simulation L2 errors for the Laplace, chi-square, log-normal, Weibull, and Pareto distributions.",
    label = "tab:l2_sample_size_group2_median"
  ),
  print_latex_summary_table(
    distributions = table2_distributions,
    stat = "mean",
    caption = "Mean (SD) sample-size simulation L2 errors for the Laplace, chi-square, log-normal, Weibull, and Pareto distributions.",
    label = "tab:l2_sample_size_group2_mean"
  )
)

writeLines(latex_out, output_file)
message("Wrote ", normalizePath(output_file, mustWork = FALSE))
