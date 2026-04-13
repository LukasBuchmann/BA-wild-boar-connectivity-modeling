-- Switch heading numbering to alphabetic (A, B, C, ...) for appendix chapters
-- in Quarto book projects. Uses the file_metadata API to detect appendix chapters
-- rather than matching heading text, so this works regardless of language.

local function is_typst_book()
  local file_state = quarto.doc.file_metadata()
  return quarto.doc.is_format("typst") and
         file_state ~= nil and
         file_state.file ~= nil
end

local appendix_started = false

local header_filter = {
  Header = function(el)
    if not is_typst_book() then
      return nil
    end

    if el.level ~= 1 then
      return nil
    end

    local file_state = quarto.doc.file_metadata()
    local bookItemType = file_state.file.bookItemType

    if bookItemType == "appendix" and not appendix_started then
      appendix_started = true
      -- Emit state change and reset heading counter before the appendix divider
      local state_change = pandoc.RawBlock('typst',
        '#state("appendix-mode", false).update(true)\n#counter(heading).update(0)')
      return {state_change, el}
    end

    return nil
  end
}

return quarto.utils.combineFilters({
  quarto.utils.file_metadata_filter(),
  header_filter
})
