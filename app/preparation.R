preparation_ui <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      h2("Transformations"),
      checkboxInput(ns("standardize"), "Standardiser les nombres"),
      conditionalPanel(
        "output.can_encode", ns = ns,
        checkboxInput(ns("encode"), "Encoder les catégories en 0/1")
      ),
      conditionalPanel(
        "!output.can_encode", ns = ns,
        p(class = "diagnostic-help", "Aucune catégorie éligible à l’encodage.")
      ),
      width = 3
    ),
    mainPanel(
      div(class = "diagnostic-help", textOutput(ns("cleaning_method"))),
      uiOutput(ns("summary")),
      panel_card(
        "Aperçu du tableau transformé",
        div(class = "table-scroll preparation-preview", tableOutput(ns("preview")))
      ),
      div(class = "imputation-export", downloadButton(ns("download"), "Télécharger le tableau préparé (CSV)")),
      width = 9
    )
  )
}

preparation_server <- function(id, imputed, eligible, cleaning_method) {
  moduleServer(id, function(input, output, session) {
    available <- reactive({
      encoding_candidates(imputed()[eligible()])
    })

    settings <- reactive({
      reference <- imputed()
      find_transformations(reference, eligible(), isTRUE(input$standardize),
                           if (isTRUE(input$encode)) available()$categories else list())
    })

    can_encode <- reactive(length(available()$categories) > 0L)
    output$can_encode <- can_encode
    outputOptions(output, "can_encode", suspendWhenHidden = FALSE)
    output$cleaning_method <- renderText(paste("Nettoyage utilisé", cleaning_method(), sep = " · "))

    output$summary <- renderUI({
      transforms <- settings()
      too_many <- available()$too_many
      tagList(
        div(class = "metric-grid",
            metric("Cellules vides après nettoyage", sum(is.na(imputed()))),
            metric("Colonnes standardisées", length(transforms$scales)),
            if (can_encode()) metric("Colonnes encodées", length(transforms$categories)),
            if (can_encode()) metric("Indicateurs créés", sum(lengths(transforms$categories)))),
        if (length(too_many) && (isTRUE(input$encode) || !can_encode())) {
          p(class = "analysis-note",
            paste("Plus de 30 catégories", paste(too_many, collapse = ", "), sep = " · "))
        }
      )
    })

    output$preview <- renderTable(apply_transformations(head(imputed()), settings()), na = "NA")
    output$download <- downloadHandler(
      filename = function() "donnees_preparees.csv",
      content = function(path) {
        write.csv(apply_transformations(imputed(), settings()), path, row.names = FALSE, na = "")
      }
    )
  })
}

encoding_candidates <- function(columns) {
  qualitative <- columns[vapply(columns, function(x) is.character(x) || is.logical(x), logical(1))]
  observed <- lapply(qualitative, function(x) sort(unique(as.character(x[!is.na(x)]))))
  list(categories = observed[lengths(observed) > 0L & lengths(observed) <= 30L],
       too_many = names(observed)[lengths(observed) > 30L])
}

find_transformations <- function(reference, eligible, standardize, categories) {
  columns <- reference[eligible]
  scales <- list()

  if (standardize) {
    numeric <- columns[vapply(columns, function(x) is.numeric(x) && any(!is.na(x)), logical(1))]
    scales <- lapply(numeric, function(x) {
      deviation <- sd(x, na.rm = TRUE)
      c(center = mean(x, na.rm = TRUE), scale = if (is.finite(deviation) && deviation > 0) deviation else 1)
    })
  }

  list(scales = scales, categories = categories)
}

apply_transformations <- function(data, transforms) {
  for (column in names(transforms$scales)) {
    parameters <- transforms$scales[[column]]
    data[[column]] <- (data[[column]] - parameters[["center"]]) / parameters[["scale"]]
  }

  taken <- names(data)
  for (column in names(transforms$categories)) {
    levels <- transforms$categories[[column]]
    labels <- tail(make.unique(c(taken, paste(column, levels, sep = "__"))), length(levels))
    values <- as.character(data[[column]])
    for (index in seq_along(levels)) {
      data[[labels[index]]] <- as.integer(values == levels[index])
    }
    taken <- c(taken, labels)
    data[[column]] <- NULL
  }
  data
}
