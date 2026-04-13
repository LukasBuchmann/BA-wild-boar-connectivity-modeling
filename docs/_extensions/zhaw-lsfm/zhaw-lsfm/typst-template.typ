
// This is an example typst template (based on the default template that ships
// with Quarto). It defines a typst function named 'article' which provides
// various customization options. This function is called from the 
// 'typst-show.typ' file (which maps Pandoc metadata function arguments)
//
// If you are creating or packaging a custom typst template you will likely
// want to replace this file and 'typst-show.typ' entirely. You can find 
// documentation on creating typst templates and some examples here: 
//   - https://typst.app/docs/tutorial/making-a-template/
//   - https://github.com/typst/templates


#let article(
  title: none,
  subtitle: none,
  authors: none,
  date: none,
  institut: none,
  keywords: none,
  confidential: false,
  thesis-type: none,
  degree-type: none,
  study-year: none,
  submission-date: none,
  study-direction: none,
  supervisors: none,
  cover-image: none,
  abstract: none,
  abstract-title: none,
  cols: 1,
  margin: (x: 1.25in, y: 2cm),
  paper: "us-letter",
  lang: "en",
  region: "US",
  font: (),
  fontsize: 11pt,
  title-size: 1.5em,
  subtitle-size: 1.25em,
  heading-family: "libertinus serif",
  heading-weight: "bold",
  heading-style: "normal",
  heading-color: black,
  heading-line-height: 0.65em,
  sectionnumbering: none,
  number-depth: none,
  pagenumbering: "1",
  toc: false,
  toc_title: none,
  toc_depth: none,
  toc_indent: 1.5em,
  doc,
) = {
  set page(
    paper: paper,
    margin: margin,
    numbering: none,
  )
  set par(justify: true)
  set text(lang: lang,
           region: region,
           size: fontsize,
           ..if font != () { (font: font) } else { (:) })
  // State to track if we're in appendix mode
  let appendix-mode = state("appendix-mode", false)

  set heading(
    numbering: if sectionnumbering != none {
      (..nums) => if nums.pos().len() <= number-depth {
        let is-appendix = appendix-mode.get()
        if is-appendix {
          numbering("A.1.", ..nums)
        } else {
          numbering(sectionnumbering, ..nums)
        }
      }
    } else { none }
  )

  show heading: it => {
    it
    v(0.5em)
  }
  // Pagebreaks are handled by Lua filter to avoid container conflicts
  
  // Title page
  if title != none or authors != none {
    page(numbering: none)[
      #align(center)[
        #if lang == "en" [
          ZURICH UNIVERSITY OF APPLIED SCIENCES \
          DEPARTEMENT LIFE SCIENCES AND FACILITY MANAGEMENT
        ] else [
          ZÜRCHER HOCHSCHULE FÜR ANGEWANDTE WISSENSCHAFTEN \
          DEPARTEMENT LIFE SCIENCES UND FACILITY MANAGEMENT
        ]
        #if institut != none [
          \ #upper(institut)
        ]
      ]
      #v(1fr)
      #if title != none {
        align(center)[#block(inset: 2em)[
          #set par(leading: heading-line-height)
          #text(weight: "bold", size: title-size)[#title]
          #if subtitle != none {
            parbreak()
            text(weight: "bold", size: subtitle-size)[#subtitle]
          }
          #if cover-image != none {
            v(0.5cm)
            image(cover-image.src, width: eval(cover-image.max-width), height: eval(cover-image.max-height), fit: "contain")
            v(0.1cm)
          }
          #if confidential {
            parbreak()
            if lang == "en" [
              #text(weight: "bold", size: subtitle-size)[confidential]
            ] else [
              #text(weight: "bold", size: subtitle-size)[vertraulich]
            ]
          }
          #if thesis-type != none {
            v(0.2cm)
            text(weight: "bold", size: subtitle-size)[#thesis-type]
          }
        ]]
      }

      #v(0.5cm)

      #if authors != none {
        align(center)[
          #if lang == "en" [
            by \
          ] else [
            von \
          ]
          #for author in authors [
            #author.name \
          ]
          #if degree-type != none and study-year != none [
            \ #degree-type #study-year
          ]
          #if submission-date != none [
            \ #if lang == "en" [
              Submission date #submission-date
            ] else [
              Abgabedatum #submission-date
            ]
          ]
          #if study-direction != none [
            \ Studienrichtung #study-direction
          ]
        ]
      }

      #v(1fr)
      
      #if supervisors != none {
        align(left)[
          #if lang == "en" [
            *Supervisors:* \
          ] else [
            *Betreuer / Betreuerinnen:* \
          ]
          #for supervisor in supervisors [
            #if supervisor.title != none [#supervisor.title] #supervisor.name \
            #supervisor.affiliation \
            \
          ]
        ]
      }
    ]
  }

  // Imprint page (second page)
  if title != none or authors != none {
    page(numbering: none)[
      #v(1fr)
      
      #align(left)[
        #if lang == "en" [
          *Imprint*
        ] else [
          *Impressum*
        ]
      ]
      
      #v(2em)
      
      // Citation section
      #if authors != none and title != none {
        align(left)[
          #if lang == "en" [
            *Recommended Citation:*
          ] else [
            *Zitiervorschlag:*
          ]
        ]
        
        align(left)[
          #text(size: 10pt)[
            #for author in authors [
              #author.name
              #if author != authors.last() [, ]
            ]
            #if date != none [(#date). ] else if submission-date != none [(#{
              let date-str = repr(submission-date).trim("[").trim("]")
              let year = date-str.split("-").at(0)
              year
            }). ] else [(n.d.). ]
            #emph[#title#if subtitle != none [: #subtitle]]. 
            #if institut != none {
              if lang == "en" [
                Zurich University of Applied Sciences, Department Life Sciences and Facility Management, #institut.
              ] else [
                Zürcher Hochschule für Angewandte Wissenschaften, Departement Life Sciences und Facility Management, #institut.
              ]
            }
          ]
        ]
        
        v(2em)
      }
      
      // Keywords
      #if keywords != none {
        align(left)[
          #if lang == "en" [
            *Keywords:* #keywords
          ] else [
            *Schlagworte:* #keywords
          ]
        ]
        v(2em)
      }
      
      // Institute
      #if institut != none {
        align(left)[
          #if lang == "en" [
            #institut \
            Department Life Sciences and Facility Management \
            Zurich University of Applied Sciences
          ] else [
            #institut \
            Departement Life Sciences und Facility Management \
            Zürcher Hochschule für Angewandte Wissenschaften
          ]
        ]
      }
    ]
  }

  // Start page numbering for main content with header
  set page(
    numbering: pagenumbering,
    header: [
      #text(size: 0.8em)[
        #grid(
          columns: (1fr, 1fr, 1fr),
          align: (left, center, right),
          gutter: 1em,
          [ZHAW LSFM],
          [#if thesis-type != none [#thesis-type]],
          [#if authors != none [
            #for author in authors [
              #author.name
              #if author != authors.last() [, ]
            ]
          ]]
        )
      ]
      #v(0.3em)
      #line(length: 100%, stroke: 0.5pt)
    ],
    footer: [
      #text(size: 0.8em)[
        #align(center)[
          #context counter(page).display()
        ]
      ]
    ]
  )
  counter(page).update(1)

  if abstract != none {
    block(inset: 2em)[
    #text(weight: "semibold")[#abstract-title] #h(1em) #abstract
    ]
  }

  if toc {
    let title = if toc_title == none {
      auto
    } else {
      toc_title
    }
    block(above: 0em, below: 2em)[
    #outline(
      title: toc_title,
      depth: toc_depth,
      indent: toc_indent
    );
    ]
  }

  if cols == 1 {
    doc
  } else {
    columns(cols, doc)
  }
}

#set table(
  inset: 6pt,
  stroke: none
)
