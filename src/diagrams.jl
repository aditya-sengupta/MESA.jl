using Unitful: K, g, cm, ustrip
using UnitfulAstro: Lsun, Rsun
using Plots: plot, plot!, scatter, scatter!
using Chain
using DataFrames, InMemoryDatasets

"""
Generates a color-magnitude diagram from a Gaia catalog.

Parameters:
    df::DataFrame - the stellar catalog
    band::String - the band whose magnitude we want to use

Returns:
    plot::Plot - the color-magnitude diagram, plotting color (BP - RP) against abs magnitude
"""
function color_magnitude_diagram(df, band="G", color1=:bp_rp, color2=nothing; kwargs...)
    if color1 == :bp_rp && isnothing(color2)
        color = df[!, :bp_rp]
        magname = Symbol(lowercase(band)*"mag")
        if !(string(magname) in names(df))
            magname = Symbol("mag_"*lowercase(band))
        end
        xl = "BP - RP (mag)"
    else
        # SDSS
        color = df[!, color1] .- df[!, color2]
        if "petromag_"*lowercase(band) in names(df)
            magname = Symbol("petromag_"*lowercase(band))
        elseif "mag_"*lowercase(band) in names(df)
            magname = Symbol("mag_"*lowercase(band))
        else
            magname = Symbol(lowercase(band))
        end
        xl = "$(string(color1)[end]) - $(string(color2[end])) (mag)"
    end
    diagram(color, df[!, magname]; xlabel=xl, ylabel="Magnitude ($band)", yflip=true, kwargs...)
end

dfval(df::Dataset, c::Symbol) = getproperty(df, c).val
dfval(df::DataFrame, c::Symbol) = getproperty(df, c)

function diagram(df, xcol::Union{String,Symbol}, ycol::Union{String,Symbol}; kwargs...)
    diagram(
        dfval(df, xcol),
        dfval(df, ycol),
        xlabel=get(kwargs, :xlabel, uppercasefirst(string(xcol))),
        ylabel=get(kwargs, :ylabel, uppercasefirst(string(ycol)));
        kwargs...
    )
end

function diagram(xseries::AbstractArray, yseries::AbstractArray; addplot=nothing, subset=nothing, savepath=nothing, points=nothing, method=scatter!, kwargs...)
    if isnothing(addplot)
        plot()
    end
    if isnothing(subset)
        subset = 1:length(xseries)
    elseif subset isa Int
        subset = sample(1:length(xseries), subset)
    end
    shape = :none
    if occursin("scatter", string(method))
        shape = :star4
    end
    p = method(xseries[subset], yseries[subset], shape=shape, label=nothing; kwargs...)
    if !isnothing(savepath)
        Plots.savefig(p, savepath)
    end
    if !isnothing(points)
        for (descr, n) in points
            scatter!([xseries[n]], [yseries[n]], label=descr, msw=0)
        end
    end
    p
end

"""
Generates a Hertzsprung-Russell diagram from a MESA run.
"""
function hertzsprung_russell_diagram(df; nticks=5, kwargs...)
    if hasproperty(df, :log_Teff)
        # MESA
        temperatures = 10.0 .^ df.log_Teff * K
        luminosities = 10.0 .^ df.log_L * Lsun
    else
        # Gaia
        temperatures = df.teff_val * K
        luminosities = df.lum_val * Lsun
    end

    xtr = tickrange(temperatures, nticks)
    ytr = tickrange(luminosities, nticks)
    diagram(
        temperatures, luminosities,
        xscale=:log10, yscale=:log10, 
        xflip=true,
        label=nothing, xlabel="Effective temperature", ylabel="Luminosity",
        xticks=(xtr, xtr), yticks=(ytr, ytr); kwargs...
    )
end

extract_field(d, u, l) = l ? (10 .^ d * u) : (d * u)
scale(l) = l ? :log10 : :identity

# don't use this, just modify the df first
function diagram(df, xfield::NamedTuple, yfield::NamedTuple; kwargs...)
    xr = extract_field(getproperty(df, xfield.name), xfield.unit, xfield.islog)
    yr = extract_field(getproperty(df, yfield.name), yfield.unit, yfield.islog)
    diagram(xr, yr; xscale=scale(xfield.islog), yscale=scale(yfield.islog), kwargs...)
end

density_temperature_diagram(df; kwargs...) = diagram(
    df,
    (name = "logRho", unit=g/cm^3, islog=true),
    (name = "logT", unit=K, islog=true),
    xlabel="Density", ylabel="Temperature";
    kwargs...
)

function abundances_diagram(df)
    radius = 10 .^ df.logR
    p = Plots.plot()
    plot!(xlabel = "Radius ($Rsun)", ylabel="Abundances", yscale=:log10, legend=:bottomright)
    plot!(radius, df.h1, label="Hydrogen")
    plot!(radius, df.he4, label="Helium")
    plot!(radius, df.c12, label="Carbon-12")
    plot!(radius, df.n14, label="Nitrogen-14")
    p
end

export diagram, color_magnitude_diagram, hertzsprung_russell_diagram, mesa_profile_diagram, density_temperature_diagram, abundances_diagram