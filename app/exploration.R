exploration_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      radioButtons(ns("mode"), "Que veux-tu examiner ?",
                   c("Une colonne" = "one", "Deux colonnes" = "two")),
      selectizeInput(ns("first"), "Colonne à examiner", NULL),
      conditionalPanel(
        "input.mode === 'two'", ns = ns,
        selectizeInput(ns("second"), "Colonne à comparer", NULL)
      ),
      conditionalPanel(
        "input.mode === 'one' && output.first_is_numeric", ns = ns,
        radioButtons(ns("numeric_chart"), "Graphique",
                     c("Histogramme" = "histogram", "Boîte à moustaches" = "boxplot"))
      ),
      width = 3
    ),
    mainPanel(
      uiOutput(ns("interpretation")),
      plotOutput(ns("chart"), height = "400px"),
      div(class = "table-scroll", tableOutput(ns("counts"))),
      width = 9
    )
  )
}

exploration_server <- function(id, dataset, raw_dataset) {
  moduleServer(id, function(input, output, session) {
    observeEvent(raw_dataset(), {
      columns <- names(raw_dataset())
      first <- if ("DYRK1A_N" %in% columns) "DYRK1A_N" else columns[1]
      second <- if ("ITSN1_N" %in% columns) "ITSN1_N" else head(setdiff(columns, first), 1)
      updateSelectizeInput(session, "first", choices = columns, selected = first, server = TRUE)
      updateSelectizeInput(session, "second", choices = columns, selected = second, server = TRUE)
    })

    selection <- reactive({
      data <- dataset()
      validate(need(input$first %in% names(data), "Choisis une première colonne."))
      columns <- input$first
      if (identical(input$mode, "two")) {
        validate(
          need(ncol(data) > 1L, "Il faut au moins deux colonnes pour les comparer."),
          need(input$second %in% names(data), "Choisis une deuxième colonne."),
          need(!identical(input$first, input$second), "Choisis deux colonnes différentes.")
        )
        columns <- c(columns, input$second)
      }
      selected <- data[columns]
      numeric <- vapply(selected, is.numeric, logical(1))
      selected[!numeric] <- lapply(selected[!numeric], as.character)
      selected
    })

    output$first_is_numeric <- reactive({
      data <- dataset()
      isTRUE(input$first %in% names(data)) && is.numeric(data[[input$first]])
    })
    outputOptions(output, "first_is_numeric", suspendWhenHidden = FALSE)

    complete <- reactive({
      selected <- selection()
      selected[complete.cases(selected), , drop = FALSE]
    })

    counts <- reactive({
      selected <- complete()
      req(ncol(selected) == 1L, !is.numeric(selected[[1]]))
      sort(table(selected[[1]]), decreasing = TRUE)
    })

    output$interpretation <- renderUI({
      selected <- selection()
      data <- complete()
      n <- nrow(data)
      omitted <- nrow(selected) - n
      if (ncol(data) == 1L) {
        name <- names(data)[1]
        if (is.numeric(data[[1]])) {
          explanation <- if (identical(input$numeric_chart, "boxplot")) {
            "Le trait marque la médiane. La boîte contient les 50 % de valeurs du milieu. Les points au-delà des moustaches signalent des valeurs atypiques à examiner."
          } else {
            "Chaque barre compte les valeurs comprises dans un intervalle. Une barre haute indique que beaucoup de lignes ont des valeurs dans cet intervalle."
          }
        } else {
          categories <- length(counts())
          explanation <- paste(
            "Chaque barre compte les lignes qui appartiennent à une catégorie.",
            if (categories > 15L) "Le graphique montre les 15 catégories les plus fréquentes.",
            if (categories > 50L) "Le tableau en affiche 50."
          )
        }
        return(div(class = "dataset-context",
                   h2(name),
                   p(explanation),
                   p(class = "diagnostic-help", sprintf("%d valeurs observées · %d manquantes", n, omitted))))
      }

      explanation <- pair_interpretation(data)
      div(class = "dataset-context",
          h2(paste(names(data), collapse = " et ")),
          p(class = "diagnostic-help",
            sprintf("%d lignes avec les deux valeurs · %d lignes écartées faute de mesure", n, omitted)),
          explanation)
    })

    output$chart <- renderPlot({
      data <- complete()
      validate(need(nrow(data) > 0L, "Aucune ligne avec les valeurs nécessaires."))
      if (ncol(data) == 2L) return(pair_plot(data))

      values <- data[[1]]
      if (is.numeric(values)) {
        chart <- ggplot(data.frame(Valeur = values), aes(x = Valeur)) + chart_theme()
        if (identical(input$numeric_chart, "boxplot")) {
          return(chart + geom_boxplot(fill = "#edf4f4", colour = "#285c63", width = .38) +
                   scale_y_continuous(breaks = NULL) + labs(x = names(data)[1], y = NULL))
        }
        return(chart + geom_histogram(bins = 25, fill = "#347780", colour = "white") +
                 labs(x = names(data)[1], y = "Nombre de lignes"))
      }

      largest <- head(counts(), 15L)
      bars <- data.frame(Catégorie = factor(names(largest), levels = rev(names(largest))),
                         Effectif = as.integer(largest))
      ggplot(bars, aes(Catégorie, Effectif)) +
        geom_col(fill = "#347780", width = .66) +
        coord_flip() + chart_theme() + labs(x = NULL, y = "Nombre de lignes")
    })

    output$counts <- renderTable({
      frequencies <- counts()
      shown <- head(frequencies, 50L)
      data.frame(Catégorie = names(shown), Effectif = as.integer(shown),
                 `% des lignes renseignées` = round(100 * as.integer(shown) / sum(frequencies), 1),
                 check.names = FALSE)
    })
  })
}

pair_interpretation <- function(data) {
  if (!nrow(data)) return(p("Aucune ligne ne renseigne les deux colonnes."))
  numeric <- vapply(data, is.numeric, logical(1))
  distinct <- vapply(data, function(x) length(unique(x)), integer(1))
  association <- "Il faut au moins deux valeurs distinctes dans chaque colonne pour calculer une association."

  if (all(numeric)) {
    if (all(distinct > 1L)) {
      rho <- cor(data[[1]], data[[2]], method = "spearman")
      association <- sprintf("Spearman (ρ) = %.2f", rho)
    }
    return(tagList(
      p(sprintf("Un point = une ligne. Par exemple, la première ligne complète est placée à %s sur l’axe horizontal et à %s sur l’axe vertical.",
                format(signif(data[[1]][1], 4)), format(signif(data[[2]][1], 4)))),
      p(tags$strong(association)),
      p(class = "diagnostic-help",
        "Spearman compare le rang des valeurs. Vers +1, les rangs montent ensemble. Vers −1, l’un monte quand l’autre descend. Vers 0, il y a peu de tendance dans un même sens.")
    ))
  }

  if (any(distinct[!numeric] > 20L)) {
    return(p("Le graphique et le calcul sont limités à 20 catégories par colonne. Choisis une colonne avec moins de catégories."))
  }

  if (any(numeric)) {
    values <- data[[which(numeric)]]
    groups <- data[[which(!numeric)]]
    if (all(distinct > 1L)) {
      fitted <- ave(values, groups, FUN = mean)
      eta_squared <- sum((fitted - mean(values))^2) / sum((values - mean(values))^2)
      association <- sprintf("Eta carré (η²) = %.2f. Environ %.0f %% de la variation observée se situe entre les groupes.",
                             eta_squared, 100 * eta_squared)
    }
    return(tagList(
      p("Une boîte résume les valeurs numériques d’une catégorie. Compare les médianes et regarde si les boîtes se recouvrent."),
      p(tags$strong(association))
    ))
  }

  if (all(distinct > 1L)) {
    observed <- table(data[[1]], data[[2]])
    expected <- outer(rowSums(observed), colSums(observed)) / sum(observed)
    chi_squared <- sum((observed - expected)^2 / expected)
    cramer_v <- sqrt(chi_squared / (sum(observed) * (min(dim(observed)) - 1)))
    association <- sprintf("V de Cramér = %.2f", cramer_v)
  }
  tagList(
    p("Chaque barre représente une catégorie de la première colonne. Les couleurs montrent les proportions des catégories de la deuxième colonne."),
    p(tags$strong(association)),
    p(class = "diagnostic-help", "Un V proche de 0 indique une faible association. Plus il se rapproche de 1, plus l’association observée est forte.")
  )
}

pair_plot <- function(data) {
  numeric <- vapply(data, is.numeric, logical(1))
  if (all(numeric)) {
    points <- data.frame(X = data[[1]], Y = data[[2]])
    return(ggplot(points, aes(X, Y)) +
             geom_point(colour = "#347780", alpha = .55, size = 2) +
             geom_point(data = head(points, 1L), colour = "#a95324", size = 3) +
             chart_theme() + labs(x = names(data)[1], y = names(data)[2],
                                  caption = "Le point orange correspond à l’exemple au-dessus."))
  }

  categories <- vapply(data[!numeric], function(x) length(unique(x)), integer(1))
  validate(need(all(categories <= 20L), "Plus de 20 catégories dans une des colonnes."))
  if (any(numeric)) {
    groups <- data.frame(Groupe = data[[which(!numeric)]], Valeur = data[[which(numeric)]])
    return(ggplot(groups, aes(Groupe, Valeur)) +
             geom_boxplot(fill = "#edf4f4", colour = "#285c63") +
             coord_flip() + chart_theme() +
             labs(x = names(data)[!numeric], y = names(data)[numeric]))
  }

  groups <- data.frame(Groupe = data[[1]], Catégorie = data[[2]])
  ggplot(groups, aes(Groupe, fill = Catégorie)) +
    geom_bar(position = "fill", width = .66) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), " %")) +
    coord_flip() + chart_theme() +
    labs(x = names(data)[1], y = "Part dans la catégorie", fill = names(data)[2])
}
