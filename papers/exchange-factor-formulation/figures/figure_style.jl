# Shared figure style for the paper. Figures are drawn at a fixed scale: TEXT_WIDTH canvas pixels correspond to the
# full text width of the manuscript, and every figure is included at (canvas width / TEXT_WIDTH) × \linewidth.
# Equal sizes in this file then give equal printed sizes in every figure, whatever its shape or number of panels.

using LaTeXStrings                                 # L"..." strings for LaTeX-typeset labels

const TEXT_WIDTH   = 1200                          # canvas width of a figure included at the full text width (\linewidth)
const SINGLE_WIDTH = 720                           # canvas width of a single-panel figure, included at 0.6\linewidth
const LABEL_SIZE   = 24                            # axis labels, titles and colour bar labels
const TICK_SIZE    = 20                            # tick labels and legend entries
const MARKER_SIZE  = 24                            # data markers in Makie figures, drawn as thin font glyphs ('+', '×')

if isdefined(Main, :Plots)                         # figures drawn with Plots
    Plots.default(fontfamily = "Computer Modern",  # LaTeX font for every text element
                  framestyle = :box,               # full frame around every axis, as Makie draws it
                  guidefontsize = LABEL_SIZE,      # axis labels
                  titlefontsize = LABEL_SIZE,      # titles
                  tickfontsize = TICK_SIZE,        # tick labels
                  legendfontsize = TICK_SIZE,      # legend entries
                  colorbar_titlefontsize = LABEL_SIZE,   # colour bar label
                  colorbar_tickfontsize = TICK_SIZE)     # colour bar ticks
end

# Theme for figures drawn with Makie. A Makie script applies it explicitly with set_theme!(paper_makie_theme())
# after `using CairoMakie`, so it does not depend on the order in which packages and this file were loaded.
paper_makie_theme() = merge(Theme(fontsize = LABEL_SIZE,   # labels and titles
                                  Axis = (xticklabelsize = TICK_SIZE, yticklabelsize = TICK_SIZE),   # tick labels
                                  Colorbar = (ticklabelsize = TICK_SIZE,),                          # colour bar ticks
                                  Legend = (labelsize = TICK_SIZE,)),                               # legend entries
                            theme_latexfonts())    # LaTeX font for every text element