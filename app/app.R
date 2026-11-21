library(shiny)
library(readxl)
library(ggplot2)

source("exploration.R", local = TRUE, encoding = "UTF-8")
source("preparation.R", local = TRUE, encoding = "UTF-8")
source("imputation_comparison.R", local = TRUE, encoding = "UTF-8")

missing_tokens <- c("", "NA", "N/A", "?")
numeric_methods <- c("Médiane" = "median", "Moyenne" = "mean",
                     "Estimation par voisins (KNN)" = "knn", "Garder les cellules vides" = "keep")

normalize_missing <- function(column) {
  if (is.character(column)) {
    column <- trimws(column)
    column[toupper(column) %in% missing_tokens] <- NA_character_
  }
  if (is.numeric(column)) column[!is.finite(column)] <- NA
  column
}

most_frequent <- function(column) {
  observed <- column[!is.na(column)]
  if (!length(observed)) return(column[NA_integer_])
  categories <- unique(observed)
  counts <- tabulate(match(observed, categories))
  categories[which.max(counts)]
}

impute_column <- function(column, numeric_method, categorical_method) {
  observed <- column[!is.na(column)]
  if (!length(observed) || !anyNA(column)) return(column)

  if (is.numeric(column)) {
    if (numeric_method == "median") column[is.na(column)] <- median(observed)
    if (numeric_method == "mean") column[is.na(column)] <- mean(observed)
  } else if ((is.character(column) || is.logical(column)) &&
             categorical_method == "mode") {
    column[is.na(column)] <- most_frequent(observed)
  }
  column
}

impute_knn <- function(data, reference, neighbors) {
  if (!ncol(data)) return(data)
  donors <- as.matrix(reference)
  queries <- as.matrix(data)
  scales <- apply(donors, 2, sd, na.rm = TRUE)
  usable <- is.finite(scales) & scales > 0
  scaled_donors <- sweep(donors[, usable, drop = FALSE], 2, scales[usable], "/")
  scaled_queries <- sweep(queries[, usable, drop = FALSE], 2, scales[usable], "/")
  medians <- apply(donors, 2, median, na.rm = TRUE)

  for (row in which(rowSums(is.na(queries)) > 0)) {
    distances <- rep(Inf, nrow(donors))
    observed <- !is.na(scaled_queries[row, ])
    if (any(observed)) {
      differences <- sweep(scaled_donors[, observed, drop = FALSE], 2,
                           scaled_queries[row, observed], "-")
      distances <- rowMeans(differences^2, na.rm = TRUE)
    }
    nearest <- order(distances[is.finite(distances)])
    candidates <- which(is.finite(distances))[nearest]
    for (column in which(is.na(queries[row, ]))) {
      available <- candidates[!is.na(donors[candidates, column])]
      data[row, column] <- if (length(available)) {
        mean(donors[head(available, neighbors), column])
      } else {
        medians[column]
      }
    }
  }
  data
}

count_outliers <- function(column) {
  if (!is.numeric(column)) return(NA_integer_)
  observed <- column[!is.na(column)]
  if (length(observed) < 4) return(NA_integer_)
  quartiles <- quantile(observed, c(.25, .75))
  margin <- 1.5 * diff(quartiles)
  sum(observed < quartiles[1] - margin | observed > quartiles[2] + margin)
}

panel_card <- function(title, ...) {
  tags$section(class = "result-section", h2(title), ...)
}

metric <- function(label, value) {
  div(class = "metric", tags$span(label), tags$strong(format(value, big.mark = " ", scientific = FALSE)))
}

chart_theme <- function() {
  theme_minimal(base_size = 13) +
    theme(
      text = element_text(colour = "#20292c"),
      axis.text = element_text(colour = "#48565b"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "#e6eaeb"),
      strip.text = element_text(face = "bold", colour = "#285c63"),
      plot.margin = margin(12, 16, 12, 12)
    )
}

ui <- navbarPage(
  "Exploration de données",
  collapsible = TRUE,
  header = tags$head(tags$link(rel = "stylesheet", href = "dashboard.css")),
  tabPanel(
    "Import et diagnostic",
    h1(class = "page-heading", "Diagnostic des données"),
    sidebarLayout(
      sidebarPanel(
        h2("Jeu de données"),
        radioButtons("source", "Source", c("Mice Protein Expression" = "mice",
                                           "Mon fichier" = "upload")),
        conditionalPanel(
          "input.source === 'upload'",
          fileInput("file", "Fichier", accept = c(".csv", ".tsv", ".xls", ".xlsx")),
          selectInput("delimiter", "Format du texte",
                      c("Virgule, décimales avec point" = "csv",
                        "Point-virgule, décimales avec virgule" = "csv2",
                        "Tabulation, décimales avec point" = "tsv")),
          textAreaInput("description", "Description du jeu (facultative)", value = "", rows = 3,
                        placeholder = "Origine des données, signification d’une ligne…")
        ),
        uiOutput("diagnostic_columns"),
        width = 3
      ),
      mainPanel(
        uiOutput("dataset_context"),
        uiOutput("dimensions"),
        tabsetPanel(
          type = "tabs",
          tabPanel(
            "Aperçu",
            div(class = "diagnostic-view",
                p(class = "diagnostic-help", "Les six premières lignes permettent de vérifier les noms des colonnes et la lecture des valeurs."),
                div(class = "table-scroll", tableOutput("preview")))
          ),
          tabPanel(
            "Qualité des variables",
            div(class = "diagnostic-view",
                p(class = "diagnostic-help",
                  "Les colonnes les plus incomplètes sont affichées en premier. Pour les nombres, la règle de 1,5 × IQR signale les valeurs atypiques."),
                div(class = "table-scroll", tableOutput("profile")),
                conditionalPanel(
                  "input.source === 'mice'",
                  tags$details(
                    class = "diagnostic-detail",
                    tags$summary("Mice : vérifier les mesures manquantes par souris"),
                    p(class = "diagnostic-help",
                      "Chaque souris a 15 mesures prévues par protéine. Si les 15 manquent, aucune mesure de cette protéine n’est disponible pour la souris."),
                    selectInput("protein", "Protéine à vérifier", NULL),
                    textOutput("mouse_explanation"),
                    div(class = "table-scroll", tableOutput("mouse_missing"))
                  )
                ))
          ),
          tabPanel(
            "Distribution de la cible",
            div(class = "diagnostic-view",
                p(class = "diagnostic-help", "Comparez les effectifs pour repérer les classes peu représentées. Chaque barre compte des lignes du tableau."),
                textOutput("class_summary"),
                plotOutput("class_plot", height = "290px"),
                div(class = "table-scroll", tableOutput("class_counts")))
          )
        ),
        width = 9
      )
    )
  ),
  tabPanel(
    "Nettoyage et imputation",
    h1(class = "page-heading", "Remplir les cellules vides"),
    conditionalPanel(
      "input.source === 'upload' && !output.dataset_available",
      p(class = "diagnostic-help", "Charge un fichier lisible dans Import et diagnostic pour commencer.")
    ),
    conditionalPanel(
    "output.dataset_available",
    sidebarLayout(
      sidebarPanel(
        h2("Méthode de remplissage"),
        selectInput("numeric_method", "Pour les nombres", numeric_methods),
        conditionalPanel(
          "output.has_categorical_gaps",
          selectInput("categorical_method", "Pour les catégories",
                      c("Catégorie la plus fréquente" = "mode", "Garder les cellules vides" = "keep"),
                      selected = "keep")
        ),
        conditionalPanel(
          "input.numeric_method === 'knn'",
          numericInput("neighbors", "Nombre de voisins (KNN)", value = 5, min = 1, step = 1)
        ),
        tags$details(
          class = "imputation-options",
          tags$summary("Protéger des colonnes"),
          div(
            selectizeInput("protected", "Colonnes à laisser telles quelles", NULL, multiple = TRUE),
            helpText("Les colonnes protégées restent dans le CSV avec leurs cellules vides, quelle que soit la méthode choisie.")
          )
        ),
        width = 3
      ),
      mainPanel(
        div(class = "current-dataset", textOutput("imputation_source")),
        uiOutput("imputation_summary"),
        uiOutput("method_help"),
        tabsetPanel(
          tabPanel(
            "Avant et après",
            div(class = "diagnostic-view",
                selectInput("inspect", "Colonne à comparer", NULL),
                div(class = "diagnostic-help", textOutput("variable_summary")),
                plotOutput("comparison", height = "320px"))
          ),
          tabPanel(
            "Bilan par variable",
            div(class = "diagnostic-view",
                p(class = "diagnostic-help",
                  "Le tableau indique combien de valeurs manquaient au départ, combien ont été remplacées et combien restent."),
                div(class = "table-scroll", tableOutput("imputation")))
          ),
          tabPanel(
            "Comparer les méthodes",
            div(
              class = "diagnostic-view",
              p(class = "diagnostic-help",
                "On masque des valeurs connues, puis on compare les estimations aux vraies valeurs. Plus l’erreur est petite, mieux c’est."),
              actionButton("compare_imputation", "Lancer la comparaison"),
              tags$details(
                class = "diagnostic-detail",
                tags$summary("Protocole du test"),
                p(class = "diagnostic-help",
                  "Cinq tirages. Sur Mice, les mesures d’une souris restent ensemble et deux souris de chaque groupe servent au test. Pour un fichier importé, environ 20 % des lignes servent au test et jusqu’à 3 000 autres lignes servent de référence."),
                p(class = "diagnostic-help",
                  "Chaque tirage masque jusqu’à 500 valeurs isolées par type de variable. Les méthodes utilisent les mêmes valeurs masquées dans chaque scénario. L’erreur indique l’écart absolu moyen pour les nombres et le pourcentage de catégories mal estimées pour les catégories.")
              ),
              div(class = "diagnostic-view",
                  div(class = "table-scroll", tableOutput("imputation_test_scores")))
            )
          )
        ),
        div(class = "imputation-export",
            downloadButton("download", "Télécharger le tableau (CSV)"),
            p(class = "diagnostic-help",
              "Pour tester un modèle, appliquez la méthode à partir du lot d’entraînement de chaque partition.")),
        width = 9
      )
    ))
  ),
  tabPanel(
    "Exploration",
    h1(class = "page-heading", "Explorer les données"),
    conditionalPanel(
      "input.source === 'upload' && !output.dataset_available",
      p(class = "diagnostic-help", "Charge un fichier lisible dans Import et diagnostic pour commencer.")
    ),
    conditionalPanel("output.dataset_available", exploration_ui("eda"))
  ),
  tabPanel(
    "Préparation",
    h1(class = "page-heading", "Préparer les variables"),
    conditionalPanel(
      "input.source === 'upload' && !output.dataset_available",
      p(class = "diagnostic-help", "Charge un fichier lisible dans Import et diagnostic pour commencer.")
    ),
    conditionalPanel("output.dataset_available", preparation_ui("prepare"))
  ),
  tabPanel(
    "Étude Mice",
    h1(class = "page-heading", "Mice Protein Expression"),
    div(class = "study-overview",
          uiOutput("mice_summary"),
          p("Les huit groupes combinent le génotype, le protocole comportemental et le traitement."),
          p(class = "diagnostic-help",
            tags$strong("Comment lire c-CS-m ? "),
            "c indique une souris témoin (t pour trisomique), CS signifie contexte puis choc (SC pour l’ordre inverse) et m indique la mémantine (s pour la solution saline)."),
          p(a("Description du jeu de données (UCI)", href = "https://doi.org/10.24432/C50S3Z"))),
      sidebarLayout(
        sidebarPanel(
          h2("Comparer les groupes"),
          selectInput("study_protein", "Protéine à comparer", NULL),
          p(class = "diagnostic-help", "Choisis une protéine pour comparer son niveau entre les huit groupes."),
          width = 3
        ),
        mainPanel(
          panel_card("Niveau de la protéine par groupe",
                     p(class = "diagnostic-help",
                       "Chaque point représente une souris. Sa valeur est la médiane de ses mesures disponibles pour cette protéine."),
                     plotOutput("study_protein_plot", height = "340px")),
          panel_card("Effectifs des groupes",
                     p(class = "diagnostic-help", "Les effectifs sont comptés en souris. Toutes les mesures d’une souris doivent rester dans le même lot lors d’une future évaluation."),
                     div(class = "table-scroll", tableOutput("mice_groups"))),
          width = 9
        )
      )
  )
)

mice_path <- file.path("..", "data", "Data_Cortex_Nuclear.xls")
mice_data <- as.data.frame(read_xls(mice_path))
proteins <- names(mice_data)[vapply(mice_data, is.numeric, logical(1))]
mouse_id <- sub("_.*", "", mice_data$MouseID)
mouse_classes <- unique(data.frame(Souris = mouse_id, Groupe = mice_data$class))

server <- function(input, output, session) {
  raw_data <- reactive({
    if (input$source == "mice") return(mice_data)
    req(input$file)
    path <- input$file$datapath
    format <- tolower(tools::file_ext(input$file$name))
    validate(need(format %in% c("csv", "tsv", "xls", "xlsx"), "Formats acceptés : CSV, TSV, XLS et XLSX."))

    imported <- tryCatch(
      switch(format,
             xls = read_xls(path, na = missing_tokens),
             xlsx = read_xlsx(path, na = missing_tokens),
             {
               delimiter <- if (format == "tsv") "tsv" else input$delimiter
               separator <- switch(delimiter, csv = ",", csv2 = ";", tsv = "\t")
               read.table(path, header = TRUE, sep = separator,
                          dec = if (delimiter == "csv2") "," else ".",
                          quote = "\"", comment.char = "", row.names = NULL,
                          check.names = FALSE, stringsAsFactors = FALSE,
                          na.strings = missing_tokens, strip.white = TRUE)
             }),
      error = function(error) error
    )
    validate(need(!inherits(imported, "error"), "Lecture impossible. Vérifiez le format et le séparateur."))
    validate(need(nrow(imported) > 0 && ncol(imported) > 0, "Le fichier doit contenir des données."))
    imported <- as.data.frame(imported)
    names(imported) <- make.unique(ifelse(nzchar(names(imported)), names(imported), "Variable"))
    imported
  })

  output$diagnostic_columns <- renderUI({
    if (identical(input$source, "upload") && is.null(input$file)) {
      return(p(class = "diagnostic-help", "Choisis un fichier pour afficher ses colonnes."))
    }
    variables <- names(raw_data())
    tagList(
      selectInput("target", "Variable cible", choices = c("Aucune" = "", variables),
                  selected = if (identical(input$source, "mice")) "class" else ""),
      tags$details(
        tags$summary("Corriger le type d’une colonne"),
        selectizeInput("categorical", "Traiter comme des catégories",
                       choices = variables, multiple = TRUE),
        helpText("Par exemple, des codes 1, 2 et 3 peuvent désigner des groupes.")
      )
    )
  })

  dataset_available <- reactive({
    tryCatch({
      raw_data()
      TRUE
    }, shiny.silent.error = function(error) FALSE)
  })
  output$dataset_available <- dataset_available
  outputOptions(output, "dataset_available", suspendWhenHidden = FALSE)

  observeEvent(raw_data(), {
    data <- raw_data()
    variables <- names(data)
    is_case <- input$source == "mice"
    updateSelectizeInput(session, "protected", choices = variables,
                         selected = if (is_case) c("MouseID", "Genotype", "Behavior", "Treatment") else character())
    incomplete <- variables[vapply(data, function(column) anyNA(normalize_missing(column)), logical(1))]
    updateSelectInput(session, "inspect", choices = variables,
                      selected = if (is_case) "BCL2_N" else head(c(incomplete, variables), 1))
  })

  updateSelectInput(session, "protein", choices = proteins, selected = "BCL2_N")
  updateSelectInput(session, "study_protein", choices = proteins, selected = "DYRK1A_N")

  clean_data <- reactive({
    data <- raw_data()
    data[] <- lapply(data, normalize_missing)
    categorical <- intersect(names(data), input$categorical)
    data[categorical] <- lapply(data[categorical], as.character)
    data
  })

  imputation_columns <- reactive({
    setdiff(names(clean_data()), c(input$target, input$protected))
  })

  has_categorical_gaps <- reactive({
    if (identical(input$source, "upload") && is.null(input$file)) return(FALSE)
    columns <- clean_data()[imputation_columns()]
    any(vapply(columns, function(column) {
      (is.character(column) || is.logical(column)) && anyNA(column)
    }, logical(1)))
  })

  output$has_categorical_gaps <- reactive(has_categorical_gaps())
  outputOptions(output, "has_categorical_gaps", suspendWhenHidden = FALSE)

  imputed_data <- reactive({
    data <- clean_data()
    predictors <- imputation_columns()
    data[predictors] <- lapply(data[predictors], impute_column,
                               numeric_method = input$numeric_method,
                               categorical_method = input$categorical_method)
    if (identical(input$numeric_method, "knn")) {
      neighbors <- input$neighbors
      validate(need(length(neighbors) == 1L && is.finite(neighbors) &&
                      neighbors >= 1 && neighbors == floor(neighbors),
                    "Choisissez un nombre entier de voisins supérieur ou égal à 1."))
      numeric <- predictors[vapply(data[predictors], is.numeric, logical(1))]
      data[numeric] <- impute_knn(data[numeric], data[numeric], neighbors)
    }
    data
  })

  exploration_server("eda", clean_data, raw_data, imputed_data)
  cleaning_method <- reactive({
    req(input$numeric_method)
    method <- names(numeric_methods)[match(input$numeric_method, numeric_methods)]
    if (identical(input$numeric_method, "knn")) method <- sprintf("KNN (%s voisins)", input$neighbors)
    if (identical(input$categorical_method, "mode") && has_categorical_gaps()) {
      method <- paste(method, "mode pour les catégories", sep = " · ")
    }
    method
  })
  preparation_server("prepare", imputed_data, imputation_columns, cleaning_method)

  output$dataset_context <- renderUI({
    raw_data()
    if (identical(input$source, "mice")) {
      return(div(
        class = "dataset-context",
        h2("Mice Protein Expression"),
        p("Ce jeu de données mesure l’expression de 77 protéines dans le cortex cérébral de 72 souris, dont 34 trisomiques. Chaque ligne correspond à l’une des 15 mesures d’une souris."),
        p("La colonne class combine le génotype, le protocole comportemental et le traitement. MouseID permet de retrouver les mesures d’une même souris."),
        a("Description et dictionnaire des variables (UCI)", href = "https://doi.org/10.24432/C50S3Z")
      ))
    }
    description <- trimws(input$description)
    div(class = "dataset-context", h2(input$file$name),
        if (isTRUE(nzchar(description))) p(class = "user-description", description))
  })

  output$dimensions <- renderUI({
    data <- clean_data()
    div(class = "metric-grid", metric("Observations", nrow(data)), metric("Variables", ncol(data)),
        metric("Valeurs manquantes", sum(is.na(data))))
  })
  output$preview <- renderTable(head(clean_data()), na = "NA")

  output$profile <- renderTable({
    data <- clean_data()
    numeric <- vapply(data, is.numeric, logical(1))
    qualitative <- vapply(data, function(x) is.character(x) || is.logical(x), logical(1))
    missing <- colSums(is.na(data))
    profile <- data.frame(
      Variable = names(data),
      Type = ifelse(numeric, "Quantitative", ifelse(qualitative, "Qualitative", "Autre")),
      Modalités = ifelse(qualitative, vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1)), NA),
      Manquantes = missing,
      `Manquantes (%)` = round(100 * missing / nrow(data), 2),
      Atypiques = vapply(data, count_outliers, integer(1)),
      check.names = FALSE
    )
    profile[order(-missing), ]
  }, na = "—")

  class_distribution <- reactive({
    data <- clean_data()
    validate(need(input$target %in% names(data), "Sélectionnez une variable cible."))
    values <- data[[input$target]]
    observed <- as.character(values[!is.na(values)])
    validate(
      need(length(observed) > 0L, "La cible ne contient aucune classe renseignée."),
      need(length(unique(observed)) <= 50L, "Choisissez une cible avec au plus 50 catégories.")
    )
    counts <- sort(table(observed), decreasing = TRUE)
    list(counts = counts, missing = sum(is.na(values)))
  })

  output$class_plot <- renderPlot({
    counts <- class_distribution()$counts
    bars <- data.frame(Classe = factor(names(counts), levels = rev(names(counts))),
                       Effectif = as.integer(counts))
    ggplot(bars, aes(Classe, Effectif)) +
      geom_col(fill = "#347780", width = .66) + coord_flip() + chart_theme() +
      labs(x = NULL, y = "Nombre de lignes")
  })

  output$class_counts <- renderTable({
    counts <- class_distribution()$counts
    data.frame(Classe = names(counts), Effectif = as.integer(counts),
               `% des cibles renseignées` = round(100 * as.integer(counts) / sum(counts), 1),
               check.names = FALSE)
  })

  output$class_summary <- renderText({
    distribution <- class_distribution()
    counts <- distribution$counts
    missing <- sprintf(ngettext(distribution$missing, "%d ligne sans cible.", "%d lignes sans cible."),
                       distribution$missing)
    if (length(counts) < 2L) return(paste("Une seule classe renseignée.", missing))
    sprintf("La classe la plus représentée compte %.1f fois l’effectif de la moins représentée. %s",
            max(counts) / min(counts), missing)
  })

  output$imputation_summary <- renderUI({
    before <- sum(is.na(clean_data()))
    after <- sum(is.na(imputed_data()))
    div(class = "metric-grid", metric("Manquantes au départ", before), metric("Valeurs remplacées", before - after),
        metric("Manquantes restantes", after))
  })

  output$imputation_source <- renderText({
    raw_data()
    if (identical(input$source, "mice")) "Mice Protein Expression" else input$file$name
  })

  output$method_help <- renderUI({
    req(input$numeric_method, input$categorical_method)
    explanation <- switch(
      input$numeric_method,
      median = "La médiane remplace chaque valeur manquante par la valeur centrale de sa colonne. Elle est peu sensible aux valeurs extrêmes.",
      mean = "La moyenne remplace chaque valeur manquante par la moyenne des valeurs observées dans sa colonne.",
      knn = "KNN estime chaque valeur manquante à partir de lignes similaires sur les autres variables quantitatives. La médiane sert de solution de repli si aucun voisin ne dispose d’une valeur observée.",
      keep = "Les valeurs manquantes des variables quantitatives sont conservées."
    )
    div(class = "dataset-context",
        p(explanation),
        if (identical(input$categorical_method, "mode") && has_categorical_gaps()) {
          p("Pour les variables qualitatives, le mode remplit les valeurs manquantes avec la catégorie la plus fréquente.")
        },
        if (identical(input$target %in% names(clean_data()), TRUE)) {
          p(tags$strong(input$target), " est la cible. Ses valeurs restent telles quelles.")
        })
  })

  output$variable_summary <- renderText({
    original <- clean_data()
    req(input$inspect %in% names(original))
    before <- sum(is.na(original[[input$inspect]]))
    after <- sum(is.na(imputed_data()[[input$inspect]]))
    if (before == 0L) return("Cette colonne ne contient aucune cellule vide.")
    if (!input$inspect %in% imputation_columns()) {
      return(sprintf("Cette colonne est protégée. Ses %d cellules vides restent telles quelles.", before))
    }
    if (before == nrow(original)) {
      return("Toute la colonne est vide. Il faut des valeurs observées pour estimer les valeurs manquantes.")
    }
    sprintf("%d valeurs remplies sur %d cellules vides au départ.", before - after, before)
  })

  output$imputation <- renderTable({
    before <- colSums(is.na(clean_data()))
    after <- colSums(is.na(imputed_data()))
    report <- data.frame(Colonne = names(before), `Vides au départ` = before,
                         Remplies = before - after, `Encore vides` = after,
                         check.names = FALSE)
    report[order(-before), ]
  })

  imputation_test <- eventReactive(input$compare_imputation, {
    data <- clean_data()[imputation_columns()]
    validate(need(nrow(data) >= 2L, "Le test nécessite au moins deux lignes."))
    is_mice <- identical(input$source, "mice")
    scores <- withProgress(message = "Comparaison des méthodes", {
      compare_imputation(data, if (is_mice) mouse_id, if (is_mice) mouse_classes)
    })
    list(data = data, source = input$source, scores = scores)
  })

  output$imputation_test_scores <- renderTable({
    test <- imputation_test()
    validate(
      need(identical(test$source, input$source) &&
             identical(test$data, clean_data()[imputation_columns()]),
           "Relance la comparaison pour appliquer la nouvelle sélection."),
      need(!is.null(test$scores), "Il faut des valeurs connues dans les colonnes à traiter pour lancer ce test.")
    )
    test$scores
  }, na = "—")

  output$comparison <- renderPlot({
    before <- clean_data()
    req(input$inspect %in% names(before))
    original <- before[[input$inspect]]
    updated <- imputed_data()[[input$inspect]]
    validate(need(any(!is.na(updated)), "Cette variable ne contient aucune valeur observée."))
    plot_data <- data.frame(
      Valeur = if (is.numeric(original)) c(original, updated) else as.character(c(original, updated)),
      Étape = factor(rep(c("Avant", "Après"), each = length(original)), levels = c("Avant", "Après"))
    )
    plot_data <- plot_data[!is.na(plot_data$Valeur), ]
    chart <- ggplot(plot_data, aes(Valeur, fill = Étape)) + facet_wrap(~ Étape) +
      scale_fill_manual(values = c(Avant = "#87989d", Après = "#347780"), guide = "none") +
      chart_theme() + labs(x = input$inspect, y = "Effectif")
    if (is.numeric(original)) {
      return(chart + geom_histogram(bins = 25))
    }
    validate(need(length(unique(plot_data$Valeur)) <= 30,
                  "Le graphique qualitatif est limité à 30 modalités."))
    chart + geom_bar() + coord_flip()
  })

  output$download <- downloadHandler(
    filename = function() "donnees_nettoyees.csv",
    content = function(path) write.csv(imputed_data(), path, row.names = FALSE, na = "")
  )

  output$mice_summary <- renderUI({
    div(class = "metric-grid", metric("Mesures", nrow(mice_data)), metric("Souris", length(unique(mouse_id))),
        metric("Protéines", length(proteins)), metric("Classes", length(unique(mice_data$class))))
  })

  output$mice_groups <- renderTable({
    counts <- table(mouse_classes$Groupe)
    data.frame(Groupe = names(counts), Souris = as.integer(counts))
  })

  output$study_protein_plot <- renderPlot({
    req(input$study_protein %in% proteins)
    values <- mice_data[[input$study_protein]]
    observed <- !is.na(values)
    validate(need(any(observed), "Aucune mesure disponible pour cette protéine."))
    by_mouse <- aggregate(
      list(Niveau = values[observed]),
      by = list(Souris = mouse_id[observed], Groupe = mice_data$class[observed]),
      FUN = median
    )
    by_mouse$Groupe <- factor(by_mouse$Groupe, levels = sort(unique(mouse_classes$Groupe)))
    ggplot(by_mouse, aes(Groupe, Niveau)) +
      geom_boxplot(fill = "#edf4f4", colour = "#285c63", width = .55, outlier.shape = NA) +
      geom_point(position = position_jitter(width = .13, height = 0, seed = 42),
                 colour = "#347780", alpha = .85, size = 2) +
      scale_x_discrete(drop = FALSE) +
      chart_theme() + labs(x = "Groupe expérimental", y = input$study_protein)
  })

  mouse_missing_report <- reactive({
    req(input$protein %in% proteins)
    report <- aggregate(
      cbind(Manquantes = as.integer(is.na(mice_data[[input$protein]])), Mesures = 1L),
      by = list(Souris = mouse_id, Classe = mice_data$class), FUN = sum
    )
    report[order(-report$Manquantes, report$Souris), ]
  })

  output$mouse_explanation <- renderText({
    report <- mouse_missing_report()
    sprintf("%s : %d souris sur %d n’ont aucune mesure pour cette protéine. %d autres ont des mesures partielles.",
            input$protein, sum(report$Manquantes == report$Mesures), nrow(report),
            sum(report$Manquantes > 0 & report$Manquantes < report$Mesures))
  })

  output$mouse_missing <- renderTable(mouse_missing_report())
}

shinyApp(ui, server)
