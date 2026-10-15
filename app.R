library(shiny)
library(readxl)

col_names <- c("No", "PAM_seq",
               "F1_AA", "F1_avail", "F1_best", "F1_fold",
               "F2_AA", "F2_avail", "F2_best", "F2_fold",
               "F3_AA", "F3_avail", "F3_best", "F3_fold")

sense_data <- read_excel("yu2023_supplS3.xlsx", sheet = "Sense strand",
                         skip = 4, col_names = col_names)
sense_data <- sense_data[!is.na(sense_data$No), ]

antisense_data <- read_excel("yu2023_supplS3.xlsx", sheet = "Anti-sense strand",
                             skip = 4, col_names = col_names)
antisense_data <- antisense_data[!is.na(antisense_data$No), ]

# Reverse complement of a DNA sequence.
revcomp <- function(seq) {
  chartr("ACGTN", "TGCAN", paste(rev(strsplit(seq, "")[[1]]), collapse = ""))
}

# Sequence as displayed: reverse-complemented (coding orientation) for the anti-sense strand.
display_seq <- function(seq, strand) if (strand == "reverse") revcomp(seq) else seq

# Format a 6-char PAM sequence as HTML with the PAM (positions 2-4) in red bold.
# frame: 1 => split [123][456], 2 => [1][234][56], 3 => [12][345][6]
# strand = "reverse": display the reverse complement (coding orientation), so the
#   PAM becomes CCN at positions 3-5; frames split the same way on the displayed sequence.
# fixed_width = TRUE adds min-width:9ch so >> aligns in Q3 options (not wanted in inline text)
# underline: positions (in the displayed sequence) to underline
# complement = TRUE: display the complementary strand (read 3'->5'), aligned base-by-base
pam_html <- function(seq, frame, strand = "forward", fixed_width = FALSE,
                     underline = integer(0), complement = FALSE) {
  pam_pos <- 2:4
  if (strand == "reverse") pam_pos <- 3:5
  seq <- display_seq(seq, strand)
  if (complement) seq <- chartr("ACGTN", "TGCAN", seq)
  chars <- strsplit(seq, "")[[1]]
  spaces_after <- switch(as.character(frame),
    "1" = 3L,
    "2" = c(1L, 4L),
    "3" = c(2L, 5L)
  )
  inner <- ""
  for (i in seq_along(chars)) {
    ch <- chars[i]
    style <- c(if (i %in% pam_pos) "color:red;font-weight:bold",
               if (i %in% underline) "text-decoration:underline")
    if (length(style) > 0) {
      inner <- paste0(inner, '<span style="', paste(style, collapse = ";"), '">', ch, "</span>")
    } else {
      inner <- paste0(inner, ch)
    }
    if (i %in% spaces_after) inner <- paste0(inner, " ")
  }
  width_css <- if (fixed_width) ";display:inline-block;min-width:9ch" else ""
  paste0('<span style="font-family:\'Courier New\',Courier,monospace', width_css, '">', inner, "</span>")
}

# Replace positions 3 and 4 of seq with the two characters of best_pam.
modify_seq <- function(seq, best_pam) {
  chars <- strsplit(seq, "")[[1]]
  bp    <- strsplit(best_pam, "")[[1]]
  chars[3] <- bp[1]
  chars[4] <- bp[2]
  paste(chars, collapse = "")
}

# Both strands of a 6-char sequence, aligned, with 5'/3' ends labelled.
duplex_html <- function(seq, frame, strand, underline = integer(0)) {
  end <- function(x) paste0('<span class="strand-end">', x, "</span>")
  paste0('<div class="duplex">',
         end("5&prime; "), pam_html(seq, frame, strand, underline = underline), end(" 3&prime;"),
         "<br>",
         end("3&prime; "), pam_html(seq, frame, strand, underline = underline, complement = TRUE),
         end(" 5&prime;"),
         "</div>")
}

is_na_val <- function(x) is.na(x) || trimws(as.character(x)) == "n.a."

# ---- UI ---------------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      body {
        max-width: 760px;
        margin: 40px auto;
        padding: 0 20px 120px 20px;
        font-family: 'Helvetica Neue', Helvetica, Arial, sans-serif;
        font-size: 16px;
        line-height: 1.5;
        color: #222;
        background: #fff;
      }
      h3.q-title {
        font-size: 1.45em;
        font-weight: 600;
        margin-top: 30px;
        margin-bottom: 14px;
      }
      p.subtitle {
        color: #888;
        font-size: 0.82em;
        margin-top: -10px;
        margin-bottom: 20px;
      }
      .radio { margin: 9px 0; }
      .radio label {
        font-size: 1.05em;
        cursor: pointer;
        display: flex;
        align-items: center;
        gap: 8px;
        white-space: nowrap;
      }
      .form-group { margin-bottom: 0; }
      .selectize-control { margin-top: 0; }
      .selectize-input, .selectize-input input {
        font-family: 'Courier New', Courier, monospace !important;
        font-size: 1.05em !important;
      }
      .selectize-dropdown-content .option {
        font-family: 'Courier New', Courier, monospace;
        font-size: 1.05em;
        padding: 5px 10px;
      }
      hr.q-divider {
        border: none;
        border-top: 1px solid #eee;
        margin: 18px 0 0 0;
      }
      .result-box {
        margin-top: 20px;
        padding: 18px 22px;
        background: #f4f8ff;
        border-radius: 8px;
        border-left: 4px solid #3498db;
        font-size: 1.05em;
        line-height: 1.9;
      }
      .coedit {
        display: flex;
        flex-wrap: wrap;
        align-items: flex-end;
        gap: 10px 26px;
        margin: 10px 0 14px 0;
      }
      .coedit-label {
        font-size: 0.85em;
        font-style: italic;
        color: #666;
      }
      .coedit-arrow {
        font-size: 1.4em;
        color: #999;
        padding-bottom: 10px;
      }
      .duplex {
        font-family: 'Courier New', Courier, monospace;
        line-height: 1.45;
        white-space: nowrap;
      }
      .strand-end { color: #aaa; }
      .result-box.not-possible {
        background: #fff5f5;
        border-left-color: #e74c3c;
        color: #7b2020;
      }
      #reset_btn {
        position: fixed;
        bottom: 24px;
        right: 24px;
        padding: 10px 22px;
        font-size: 1em;
        font-weight: 700;
        letter-spacing: 0.06em;
        background: #c0392b;
        color: white;
        border: none;
        border-radius: 6px;
        cursor: pointer;
        box-shadow: 0 2px 10px rgba(0,0,0,0.2);
        z-index: 999;
      }
      #reset_btn:hover { background: #a93226; }
      #footer {
        position: fixed;
        bottom: 14px;
        left: 50%;
        transform: translateX(-50%);
        text-align: center;
        font-size: 0.74em;
        color: #bbb;
        line-height: 1.55;
        white-space: nowrap;
      }
      #footer a { color: #bbb; text-decoration: none; }
      #footer a:hover { text-decoration: underline; }
    "))
  ),

  # Q1 -----------------------------------------------------------------
  h3("On which strand is the NGG PAM?", class = "q-title"),
  radioButtons("strand", label = NULL,
    choiceNames  = list(
      "Sense strand (strand with the codons)",
      "Anti-sense strand (strand opposite)"
    ),
    choiceValues = list("forward", "reverse"),
    selected     = character(0)
  ),

  # Q2, Q3, result rendered server-side
  uiOutput("q2_ui"),
  uiOutput("q3_ui"),
  uiOutput("result_ui"),

  # RESET (page reload resets all state cleanly)
  tags$button("RESET", id = "reset_btn",
    onclick = "window.location.reload()"),

  # Footer
  tags$div(id = "footer",
    tags$div(HTML(
      'part of prime-editing in zebrafish <a href="https://dx.doi.org/10.17504/protocols.io.36wgqxx3xlk5/v1" target="_blank">protocols.io</a>'
    )),
    tags$div(HTML('<a href="https://github.com/francoiskroll/PAMcoedit" target="_blank">source code</a>')),
    tags$div(HTML(
      'data from <a href="https://pubmed.ncbi.nlm.nih.gov/37119812/" target="_blank">Yu et al., 2023. <em>Cell.</em></a>'
    ))
  )
)

# ---- Server -----------------------------------------------------------

server <- function(input, output, session) {

  # Q2: appears after strand chosen
  output$q2_ui <- renderUI({
    req(input$strand)
    dat     <- if (input$strand == "forward") sense_data else antisense_data
    pam_seqs <- dat$PAM_seq

    tagList(
      tags$hr(class = "q-divider"),
      h3("What is the PAM sequence context?", class = "q-title"),
      p("Are you editing the PAM itself? Put here the sequence after edit.", class = "subtitle"),
      selectizeInput("pam_seq", label = NULL,
        choices  = pam_seqs,
        selected = NULL,
        options  = list(
          placeholder = "Click to choose a PAM sequence...",
          onInitialize = I('function() { this.setValue(""); }'),
          render = I('{
            option: function(item, escape) {
              if (!item.value) {
                return "<div style=\'color:#ccc;padding:5px 10px\'>Click to choose a PAM sequence...</div>";
              }
              var s = item.value;
              return "<div style=\'font-family:Courier New,Courier,monospace;padding:5px 10px\'>" +
                     s[0] +
                     "<span style=\'color:red;font-weight:bold\'>" + s.slice(1,4) + "</span>" +
                     s.slice(4) + "</div>";
            },
            item: function(item, escape) {
              if (!item.value) return "<div></div>";
              var s = item.value;
              return "<div style=\'font-family:Courier New,Courier,monospace;font-size:1.05em\'>" +
                     s[0] +
                     "<span style=\'color:red;font-weight:bold\'>" + s.slice(1,4) + "</span>" +
                     s.slice(4) + "</div>";
            }
          }')
        )
      )
    )
  })

  # Q3: appears after PAM sequence chosen
  output$q3_ui <- renderUI({
    req(input$strand, input$pam_seq, nchar(input$pam_seq) == 6)

    seq <- input$pam_seq
    dat <- if (input$strand == "forward") sense_data else antisense_data
    row <- dat[dat$PAM_seq == seq, ]
    req(nrow(row) == 1)

    aa1 <- row$F1_AA
    aa2 <- row$F2_AA
    aa3 <- row$F3_AA

    subtitle <- if (input$strand == "reverse") {
      p("Below is the reverse complement sequence of the PAM sequence context above.",
        class = "subtitle")
    } else NULL

    make_label <- function(html_seq, aa) {
      HTML(paste0(html_seq,
                  '<span style="color:#777"> &gt;&gt; </span>',
                  '<span style="font-family:\'Courier New\',Courier,monospace">', aa, '</span>'))
    }

    tagList(
      tags$hr(class = "q-divider"),
      h3("What is the reading frame?", class = "q-title"),
      subtitle,
      radioButtons("frame", label = NULL,
        choiceNames  = list(
          make_label(pam_html(seq, 1, input$strand, fixed_width = TRUE), aa1),
          make_label(pam_html(seq, 2, input$strand, fixed_width = TRUE), aa2),
          make_label(pam_html(seq, 3, input$strand, fixed_width = TRUE), aa3)
        ),
        choiceValues = list("1", "2", "3"),
        selected     = character(0)
      )
    )
  })

  # Result: appears after frame chosen
  output$result_ui <- renderUI({
    req(input$strand, input$pam_seq, input$frame, nchar(input$pam_seq) == 6)

    seq   <- input$pam_seq
    frame <- as.integer(input$frame)
    dat   <- if (input$strand == "forward") sense_data else antisense_data
    row   <- dat[dat$PAM_seq == seq, ]
    req(nrow(row) == 1)

    best <- row[[paste0("F", frame, "_best")]]
    fold <- row[[paste0("F", frame, "_fold")]]

    if (is_na_val(best)) {
      tagList(
        tags$hr(class = "q-divider"),
        div(class = "result-box not-possible",
          "It is not possible to inactivate the PAM without changing the amino acid."
        )
      )
    } else {
      mod     <- modify_seq(seq, best)
      changed <- which(strsplit(display_seq(seq, input$strand), "")[[1]] !=
                       strsplit(display_seq(mod, input$strand), "")[[1]])
      fold_fmt <- sprintf("%.2f", as.numeric(fold))

      tagList(
        tags$hr(class = "q-divider"),
        div(class = "result-box",
          HTML(paste0(
            "You can inactivate the PAM without changing the amino acid by modifying&nbsp;it:",
            '<div class="coedit">',
              '<div><div class="coedit-label">from</div>',
                duplex_html(seq, frame, input$strand, changed), "</div>",
              '<div class="coedit-arrow">&rarr;</div>',
              '<div><div class="coedit-label">to</div>',
                duplex_html(mod, frame, input$strand, changed), "</div>",
            "</div>",
            "It is expected to be <strong>", fold_fmt, "&times;</strong>",
            " better (measured in HEK293T cells)."
          ))
        )
      )
    }
  })
}

shinyApp(ui, server)
