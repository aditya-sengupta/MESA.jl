module MESA
    using FileWatching
    using Base.Threads: @spawn
    using Plots
    using Chain
    using DataFrames
    Plots.default(titlefont = ("computer modern"), legendfont = ("computer modern"), tickfont = ("computer modern"), guidefont = ("computer modern"))

    """
    https://stackoverflow.com/questions/17479782/julia-request-user-input-from-script
    """
    function input(prompt)
        print(prompt)
        return chomp(readline())
    end

    function vcat_with_columns(df1, df2)
        columns = intersect(names(df1), names(df2))
        vcat(DataFrames.select(df1, columns), DataFrames.select(df2, columns))
    end

    """
    Reads a MESA history file and outputs it as a DataFrame.
    """
    function read_mesa(fpath::String, sink=DataFrame)
        return CSV.read(fpath, sink, skipto=7, delim=" ", ignorerepeated=true, header=6)
    end

    """
    Reads a set of MESA history files and concatenates them into a DataFrame.
    """
    function read_mesa(fpaths::Vector{String})
        return reduce(vcat, read_mesa.(fpaths))
    end

    """
    Move into the working directory of a MESA run, start the run, and move back out.
    """
    function run_mesa(fpath, runfile="rn", make=false)
        println("Starting a Julia-spawned MESA run with executable $(abspath(runfile)) on path $fpath")
        if make
            cd(() -> run(`./mk` & `./$runfile`), fpath)
        else
            cd(() -> run(`./$runfile`), fpath)
        end
    end

    """
    Watch the history.data file of an ongoing MESA run and plot new points.
    """
    function run_and_watch_mesa(fpath, runfile="rn", make=false; p=nothing, kwargs...)
        mesa_task = Task(() -> run_mesa(fpath, runfile, make))
        @spawn schedule(mesa_task)
        history = abspath(joinpath(fpath, "LOGS", "history.data"))
        while !isfile(history)
            sleep(0.1)
        end
        io = open(history)
        if isnothing(p)
            p = plot()
        end
        headers = []
        while !("model_number" in headers)
            headers = @chain io readline split(r"\s+") filter(!isempty, _)
        end
        inds = map(a -> findfirst(x -> x == a, headers),  ["log_Teff", "log_L"])

        first_point_done = false
        lastp = [0.0, 0.0]
        nlines = countlines(history) # could be inefficient but fine for now - can replace with filesize if need be
        while !istaskdone(mesa_task)
            w = watch_file(history, 5)
            n = 0
            new_nlines = 0
            try
                new_nlines = countlines(history)
            catch SystemError
                break
            end
            if new_nlines < nlines
                # we've switched to a new inlist; have to close and reopen history.data
                println("new inlist, reopening history.data")
                close(io)
                io = open(history)
                headers = []
                while !("model_number" in headers)
                    headers = @chain io readline split(r"\s+") filter(!isempty, _)
                end
            end
            nlines = new_nlines
            if !isempty(peek(io, String))
                while !eof(io)
                    data = @chain readline(io) begin
                        split(r"\s+")
                        filter(!isempty, _)
                        parse.(Float64, _)
                    end
                    if length(data) > 0
                        n = data[1]
                        thisp = data[inds]
                        if first_point_done && (n != 0)
                            plot!([lastp[1], thisp[1]], [lastp[2], thisp[2]], legend=nothing, xflip=true, xlabel="log effective temperature (K)", ylabel="log luminosity (Lsun)"; kwargs...)
                        else
                            first_point_done = true
                        end
                        lastp = thisp
                        display(p)
                    end
                end
            end
        end
        close(io)
        p
    end

    export read_mesa, run_mesa, run_and_watch_mesa

end # module MESA
