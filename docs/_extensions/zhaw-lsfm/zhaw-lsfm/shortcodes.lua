-- Shortcodes for generating figure and table lists in Typst
return {
  ['list-of-figures'] = function(args, kwargs, meta)
    -- Generate typst code for figure list
    local typst_code = '#outline(title: none, target: figure.where(kind: "quarto-float-fig"))'
    return pandoc.RawBlock('typst', typst_code)
  end,
  
  ['list-of-tables'] = function(args, kwargs, meta)
    -- Generate typst code for table list  
    local typst_code = '#outline(title: none, target: figure.where(kind: "quarto-float-tbl"))'
    return pandoc.RawBlock('typst', typst_code)
  end,
  
  ['table-of-contents'] = function(args, kwargs, meta)
    -- Generate typst code for table of contents with depth awareness
    local toc_depth = 2 -- default depth to match ZHAW requirements
    
    -- Check if toc-depth is specified in metadata
    if meta["toc-depth"] then
      toc_depth = tonumber(pandoc.utils.stringify(meta["toc-depth"])) or 2
    end
    
    local typst_code = '#outline(title: none, depth: ' .. toc_depth .. ')'
    return pandoc.RawBlock('typst', typst_code)
  end,
  
  ['references'] = function(args, kwargs, meta)
    -- Generate references div for bibliography
    return pandoc.Div({}, pandoc.Attr("refs"))
  end,

  ['zhaw-title-block'] = function(args, kwargs, meta)
    if not quarto.doc.is_format("html") then return {} end

    local lang = meta.lang and pandoc.utils.stringify(meta.lang) or "en"

    local L = {}
    if lang == "de" then
      L.university          = "Zürcher Hochschule für Angewandte Wissenschaften"
      L.department          = "Departement Life Sciences und Facility Management"
      L.author_by           = "von"
      L.date_prefix         = "Abgabedatum"
      L.direction_prefix    = "Studienrichtung"
      L.supervisors_label   = "Betreuer / Betreuerinnen:"
      L.confidential_label  = "Vertraulich"
      L.imprint_title       = "Impressum"
      L.citation_label      = "Zitiervorschlag:"
      L.keywords_label      = "Schlagworte:"
    else
      L.university          = "Zurich University of Applied Sciences"
      L.department          = "Departement Life Sciences and Facility Management"
      L.author_by           = "by"
      L.date_prefix         = "Submission date"
      L.direction_prefix    = "Study direction"
      L.supervisors_label   = "Supervisors:"
      L.confidential_label  = "Confidential"
      L.imprint_title       = "Imprint"
      L.citation_label      = "Recommended Citation:"
      L.keywords_label      = "Keywords:"
    end

    local s = pandoc.utils.stringify
    local lines = {}
    local function add(x) table.insert(lines, x) end

    add('<header id="title-block-header">')

    if meta.institut then
      add('<div class="zhaw-header">' .. L.university .. '<br>' .. L.department .. '<br>' .. s(meta.institut) .. '</div>')
    end
    if meta.title    then add('<h1 class="title">' .. s(meta.title) .. '</h1>') end
    if meta.subtitle then add('<p class="subtitle">' .. s(meta.subtitle) .. '</p>') end

    if meta['zhaw-cover-image'] then
      add('<div class="cover-image"><img src="' .. s(meta['zhaw-cover-image'].src) .. '" alt="Cover"></div>')
    end

    if meta.confidential then add('<p class="confidential">' .. L.confidential_label .. '</p>') end
    if meta['thesis-type'] then add('<p class="thesis-type">' .. s(meta['thesis-type']) .. '</p>') end

    if meta.author then
      add('<div class="author-section"><p class="von">' .. L.author_by .. '</p>')
      local authors = meta.author.tag == 'MetaList' and meta.author or {meta.author}
      for _, a in ipairs(authors) do add('<p class="author">' .. s(a) .. '</p>') end
      add('</div>')
    end

    if meta['degree-type'] then
      local yr = meta['study-year'] and (' ' .. s(meta['study-year'])) or ''
      add('<p class="degree-type">' .. s(meta['degree-type']) .. yr .. '</p>')
    end
    if meta['submission-date'] then
      add('<p class="submission-date">' .. L.date_prefix .. ' ' .. s(meta['submission-date']) .. '</p>')
    end
    if meta['study-direction'] then
      add('<p class="study-direction">' .. L.direction_prefix .. ' ' .. s(meta['study-direction']) .. '</p>')
    end

    if meta.supervisors then
      add('<div class="supervisors"><p>' .. L.supervisors_label .. '</p><ul>')
      for _, sup in ipairs(meta.supervisors) do
        local title = sup.title and s(sup.title) or ''
        local name  = sup.name  and s(sup.name)  or ''
        local aff   = sup.affiliation and s(sup.affiliation) or ''
        add('<li>' .. title .. ' ' .. name .. '<br>' .. aff .. '</li>')
      end
      add('</ul></div>')
    end

    add('</header>')

    -- Imprint
    add('<section id="imprint" class="imprint-section"><hr>')
    add('<h3 class="imprint-title">' .. L.imprint_title .. '</h3>')

    if meta.author then
      local authors = meta.author.tag == 'MetaList' and meta.author or {meta.author}
      local author_str = {}
      for _, a in ipairs(authors) do table.insert(author_str, s(a)) end
      local year = ''
      if meta.date then
        year = s(meta.date):match("^(%d%d%d%d)") or s(meta.date)
      elseif meta['submission-date'] then
        year = s(meta['submission-date']):match("^(%d%d%d%d)") or ''
      end
      local full_title = s(meta.title) .. (meta.subtitle and (': ' .. s(meta.subtitle)) or '')
      local inst_str = meta.institut and (' ' .. L.university .. ', ' .. L.department .. ', ' .. s(meta.institut) .. '.') or ''
      add('<div class="imprint-citation"><p class="imprint-label">' .. L.citation_label .. '</p>')
      add('<p class="citation-text">' .. table.concat(author_str, ', ') .. ' (' .. year .. '). <em>' .. full_title .. '</em>.' .. inst_str .. '</p></div>')
    end

    if meta.keywords then
      add('<div class="imprint-keywords"><p><strong>' .. L.keywords_label .. '</strong> ' .. s(meta.keywords) .. '</p></div>')
    end

    add('<hr></section>')

    return pandoc.RawBlock('html', table.concat(lines, '\n'))
  end
}