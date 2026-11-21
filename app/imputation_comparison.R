compare_imputation <- function(data, mouse_id = NULL, mouse_classes = NULL) {
  numeric <- vapply(data, is.numeric, logical(1))
  categorical <- vapply(data, function(x) is.character(x) || is.logical(x), logical(1))
  mice <- !is.null(mouse_id)
  if (mice) classes <- split(mouse_classes$Souris, mouse_classes$Groupe)

  score <- function(scenario, method, actual, predicted, categorical = FALSE) {
    filled <- !is.na(predicted)
    error <- if (categorical) {
      sum(actual[filled] != predicted[filled])
    } else {
      sum(abs(actual[filled] - predicted[filled]))
    }
    data.frame(scenario, method, tested = length(actual), filled = sum(filled), error)
  }

  runs <- withr::with_seed(42, lapply(seq_len(5L), function(run) {
    if (mice) {
      held_out <- unlist(lapply(classes, sample, size = 2L), use.names = FALSE)
      test_rows <- which(mouse_id %in% held_out)
    } else {
      test_size <- min(nrow(data) - 1L, ceiling(nrow(data) / 5))
      test_rows <- sample.int(nrow(data), test_size)
    }
    reference_rows <- setdiff(seq_len(nrow(data)), test_rows)
    if (length(reference_rows) > 3000L) reference_rows <- sample(reference_rows, 3000L)
    reference <- data[reference_rows, , drop = FALSE]
    validation <- data[test_rows, , drop = FALSE]
    results <- list()

    if (any(numeric)) {
      train <- reference[numeric]
      test <- validation[numeric]
      observed <- !is.na(test)
      observed[, colSums(!is.na(train)) == 0L] <- FALSE
      positions <- which(observed)
      isolated <- matrix(FALSE, nrow(test), ncol(test))
      isolated[positions[sample.int(length(positions), min(length(positions), 500L))]] <- TRUE
      masks <- list("Mesures numériques isolées" = isolated)

      if (mice) {
        by_mouse <- matrix(FALSE, nrow(test), ncol(test))
        for (mouse in held_out) {
          rows <- mouse_id[test_rows] == mouse
          available <- which(colSums(observed[rows, , drop = FALSE]) == sum(rows))
          if (length(available)) by_mouse[rows, available[sample.int(length(available), 1L)]] <- TRUE
        }
        masks[["Toutes les mesures d’une protéine par souris"]] <- by_mouse
      }

      medians <- vapply(train, median, numeric(1), na.rm = TRUE)
      means <- vapply(train, mean, numeric(1), na.rm = TRUE)
      for (scenario in names(masks)) {
        mask <- masks[[scenario]]
        if (!any(mask)) next
        rows <- rowSums(mask) > 0L
        active_mask <- mask[rows, , drop = FALSE]
        actual <- test[mask]
        hidden <- test[rows, , drop = FALSE]
        hidden[active_mask] <- NA
        columns <- which(active_mask, arr.ind = TRUE)[, "col"]
        predictions <- list("Médiane" = medians[columns], "Moyenne" = means[columns],
                            "KNN (5 voisins)" = impute_knn(hidden, train, 5L)[active_mask],
                            "Garder les cellules vides" = rep(NA_real_, length(actual)))
        results <- c(results, lapply(names(predictions), function(method) {
          score(scenario, method, actual, predictions[[method]])
        }))
      }
    }

    if (any(categorical)) {
      train <- as.data.frame(lapply(reference[categorical], as.character))
      test <- as.data.frame(lapply(validation[categorical], as.character))
      modes <- vapply(train, most_frequent, character(1))
      eligible <- !is.na(test)
      eligible[, is.na(modes)] <- FALSE
      positions <- which(eligible, arr.ind = TRUE)
      if (nrow(positions)) {
        masked <- positions[sample.int(nrow(positions), min(nrow(positions), 500L)), , drop = FALSE]
        actual <- as.matrix(test)[masked]
        scenario <- "Catégories isolées"
        results <- c(results, list(
          score(scenario, "Catégorie la plus fréquente", actual, modes[masked[, "col"]], categorical = TRUE),
          score(scenario, "Garder les cellules vides", actual, rep(NA_character_, length(actual)), categorical = TRUE)
        ))
      }
    }
    do.call(rbind, results)
  }))

  scores <- do.call(rbind, runs)
  if (is.null(scores)) return(NULL)
  totals <- aggregate(cbind(tested, filled, error) ~ scenario + method, data = scores, FUN = sum)
  test_order <- match(paste(totals$scenario, totals$method), paste(scores$scenario, scores$method))
  totals <- totals[order(test_order), ]
  errors <- totals$error / totals$filled
  categorical <- totals$scenario == "Catégories isolées"
  labels <- ifelse(categorical, sprintf("%.1f %%", 100 * errors), sprintf("%.3f", errors))
  labels[totals$filled == 0L] <- NA_character_

  data.frame(
    Scénario = totals$scenario, Méthode = totals$method,
    Masquées = as.integer(totals$tested), Estimées = as.integer(totals$filled),
    `Restées vides` = as.integer(totals$tested - totals$filled), Erreur = labels,
    check.names = FALSE
  )
}
