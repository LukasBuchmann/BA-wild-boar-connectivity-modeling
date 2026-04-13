
// Typst custom formats typically consist of a 'typst-template.typ' (which is
// the source code for a typst template) and a 'typst-show.typ' which calls the
// template's function (forwarding Pandoc metadata values as required)
//
// This is an example 'typst-show.typ' file (based on the default template  
// that ships with Quarto). It calls the typst function named 'article' which 
// is defined in the 'typst-template.typ' file. 
//
// If you are creating or packaging a custom typst template you will likely
// want to replace this file and 'typst-template.typ' entirely. You can find
// documentation on creating typst templates here and some examples here:
//   - https://typst.app/docs/tutorial/making-a-template/
//   - https://github.com/typst/templates
#show: doc => article(
$if(title)$
  title: [$title$],
$endif$
$if(subtitle)$
  subtitle: [$subtitle$],
$endif$
$if(by-author)$
  authors: (
$for(by-author)$
$if(it.name.literal)$
    ( name: [$it.name.literal$],
      affiliation: [$for(it.affiliations)$$it.name$$sep$, $endfor$],
      email: [$it.email$] ),
$endif$
$endfor$
    ),
$endif$
$if(date)$
  date: [$date$],
$endif$
$if(institut)$
  institut: [$institut$],
$endif$
$if(keywords)$
  keywords: [$keywords$],
$endif$
$if(confidential)$
  confidential: $confidential$,
$endif$
$if(thesis-type)$
  thesis-type: [$thesis-type$],
$endif$
$if(degree-type)$
  degree-type: [$degree-type$],
$endif$
$if(study-year)$
  study-year: [$study-year$],
$endif$
$if(submission-date)$
  submission-date: [$submission-date$],
$endif$
$if(study-direction)$
  study-direction: [$study-direction$],
$endif$
$if(supervisors)$
  supervisors: (
$for(supervisors)$
    (
      title: [$supervisors.title$],
      name: [$supervisors.name$],
      affiliation: [$supervisors.affiliation$],
    ),
$endfor$
  ),
$endif$
$if(zhaw-cover-image)$
  cover-image: (
$if(zhaw-cover-image.src)$
    src: "$zhaw-cover-image.src$",
$else$
    src: "$zhaw-cover-image$",
$endif$
$if(zhaw-cover-image.max-width)$
    max-width: "$zhaw-cover-image.max-width$",
$else$
    max-width: "100%",
$endif$
$if(zhaw-cover-image.max-height)$
    max-height: "$zhaw-cover-image.max-height$",
$else$
    max-height: "5cm",
$endif$
  ),
$endif$
$if(lang)$
  lang: "$lang$",
$endif$
$if(region)$
  region: "$region$",
$endif$
$if(abstract)$
  abstract: [$abstract$],
  abstract-title: "$labels.abstract$",
$endif$
$if(margin)$
  margin: ($for(margin/pairs)$$margin.key$: $margin.value$,$endfor$),
$endif$
$if(papersize)$
  paper: "$papersize$",
$endif$
$if(mainfont)$
  font: ("$mainfont$",),
$elseif(brand.typography.base.family)$
  font: $brand.typography.base.family$,
$endif$
$if(fontsize)$
  fontsize: $fontsize$,
$elseif(brand.typography.base.size)$
  fontsize: $brand.typography.base.size$,
$endif$
$if(title)$
$if(brand.typography.headings.family)$
  heading-family: $brand.typography.headings.family$,
$endif$
$if(brand.typography.headings.weight)$
  heading-weight: $brand.typography.headings.weight$,
$endif$
$if(brand.typography.headings.style)$
  heading-style: "$brand.typography.headings.style$",
$endif$
$if(brand.typography.headings.color)$
  heading-color: $brand.typography.headings.color$,
$endif$
$if(brand.typography.headings.line-height)$
  heading-line-height: $brand.typography.headings.line-height$,
$endif$
$endif$
  sectionnumbering: "$section-numbering$",
  number-depth: $number-depth$,
  pagenumbering: $if(page-numbering)$"$page-numbering$"$else$none$endif$,
$if(toc)$
  toc: $toc$,
$endif$
$if(toc-title)$
  toc_title: [$toc-title$],
$endif$
$if(toc-indent)$
  toc_indent: $toc-indent$,
$endif$
  toc_depth: $toc-depth$,
  cols: $if(columns)$$columns$$else$1$endif$,
  doc,
)
