using Plots
default(linewidth=2, grid=true, fontfamily="Computer Modern")
# function gen_irf(irf::DataFrame)
#     fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(16, 12))
#     ax[1,1].plot(irf_df.u, label=:u, alpha=0.6)
#     ax[1,1].plot(irf_df.v, label=:v, alpha=0.6)
#     ax[1,1].plot(irf_df.θ, label=:θ, alpha=0.6)
#     ax[1,1].plot(irf_df.e, label=:e, alpha=0.6)
#     ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[1,1].legend()

#     ax[1,2].plot(irf_df.C, label=:C, alpha=0.6)
#     ax[1,2].plot(irf_df.Y_c, label=:Y_c, alpha=0.6)
#     ax[1,2].plot(irf_df.Y, label=:Y, alpha=0.6)
#     ax[1,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[1,2].legend()

#     ax[2,1].plot(irf_df.N_e, label=:N_e, alpha=0.6)
#     ax[2,1].plot(irf_df.N, label=:N, alpha=0.6)
#     ax[2,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[2,1].legend()

#     ax[2,2].plot(irf_df.z, label=:z, alpha=0.6)
#     #ax[2,2].plot(irf_df.w, label=:w, alpha=0.6)
#     #ax[2,2].plot(irf_df.w_R,label=:w_R, alpha=0.6)
#     ax[2,2].plot(irf_df.Y, label=:Y, alpha=0.6)
#     ax[2,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[2,2].legend()
#     display(fig)
# end

"""
Comparison of irfs to highlight model transmission mechanism

"""
# function gen_irf_comp(irf_bas::DataFrame, irf_alt::DataFrame, labels)
#     fig, ax = plt.subplots(ncols=2, nrows=3, figsize=(16, 16))

#     ax[1,1].plot(irf_bas.u, alpha=0.6, label=labels[1])
#     ax[1,1].plot(irf_alt.u, alpha=0.6, label=labels[2])
#     ax[1,1].set_title("u")
#     ax[1,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[1,1].legend()

#     ax[1,2].plot(irf_bas.v, alpha=0.6, label=labels[1])
#     ax[1,2].plot(irf_alt.v, alpha=0.6, label=labels[2])
#     ax[1,2].set_title("v")
#     ax[1,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[1,2].legend()

#     ax[2,1].plot(irf_bas.N, alpha=0.6, label=labels[1])
#     ax[2,1].plot(irf_alt.N, alpha=0.6, label=labels[2])
#     ax[2,1].set_title("N")
#     ax[2,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[2,1].legend()

#     ax[2,2].plot(irf_bas.N_e, alpha=0.6, label=labels[1])
#     ax[2,2].plot(irf_alt.N_e, alpha=0.6, label=labels[2])
#     ax[2,2].set_title("N_e")
#     ax[2,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[2,2].legend()

#     ax[3,1].plot(irf_bas.C, alpha=0.6, label=labels[1])
#     ax[3,1].plot(irf_alt.C, alpha=0.6, label=labels[2])
#     ax[3,1].set_title("C")
#     ax[3,1].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[3,1].legend()

#     ax[3,2].plot(irf_bas.K, alpha=0.6, label=labels[1])
#     ax[3,2].plot(irf_alt.K, alpha=0.6, label=labels[2])
#     ax[3,2].set_title("K")
#     ax[3,2].yaxis.set_major_formatter(matplotlib.ticker.PercentFormatter())
#     ax[3,2].legend()

#     display(fig)
# end

function gen_irf(irf_df::DataFrame)
    # Set up a 3x2 layout with specified size
    p = Plots.plot(layout=(3,2), size=(800, 600), 
            legend=true, alpha=0.6)
    
    # Top left plot: labor market variables
    plot!(p[1], irf_df.u, label="u", subplot=1)
    plot!(p[1], irf_df.v, label="v")
    plot!(p[1], irf_df.θ, label="θ")
    plot!(p[1], irf_df.e, label="e")
    # Format y-axis as percentage
    yticks!(p[1], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Top middle plot: product line variables
    plot!(p[2], irf_df.N_e, label="N_e")
    plot!(p[2], irf_df.N, label="N")
    yticks!(p[2], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Top right plot: vacancy value variables
    plot!(p[3], irf_df.Q, label="Q")
    plot!(p[3], irf_df.K, label="K")
    yticks!(p[3], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom left plot: consumption and data-consistent counterpart
    plot!(p[4], irf_df.C, label="C")
    plot!(p[4], irf_df.C_R, label="C_R")
    yticks!(p[4], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom middle plot: output and data-consistent counterpart
    plot!(p[5], irf_df.Y, label="Y")
    plot!(p[5], irf_df.Y_R, label="Y_R")
    yticks!(p[5], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom right plot
    plot!(p[6], irf_df.z, label="z")
    plot!(p[6], irf_df.labor_prod, label="labor productivity")
    plot!(p[6], irf_df.δ, label="δ")
    plot!(p[6], irf_df.Y, label="Y")
    yticks!(p[6], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Display the plot
    display(p)
    png("clipboard")
end

function gen_irf_comp(irf_bas::DataFrame, irf_alt::DataFrame, labels)
    # Create a 3x2 subplot layout
    p = Plots.plot(
        layout=(3,2), 
        size=(800,800), 
        legend=:topright,
        fmt=:png
    )

    # Plot for u
    plot!(p[1,1], 
        [irf_bas.u irf_alt.u],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="u",
        ylabel="%"
    )

    # Plot for v
    plot!(p[1,2], 
        [irf_bas.v irf_alt.v],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="v",
        ylabel="%"
    )

    # Plot for N
    plot!(p[2,1], 
        [irf_bas.N irf_alt.N],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="N",
        ylabel="%"
    )

    # Plot for N_e
    plot!(p[2,2], 
        [irf_bas.N_e irf_alt.N_e],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="N_e",
        ylabel="%"
    )

    # Plot for C
    plot!(p[3,1], 
        [irf_bas.C irf_alt.C],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="C",
        ylabel="%"
    )

    # Plot for K
    plot!(p[3,2], 
        [irf_bas.K irf_alt.K],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title="K",
        ylabel="%"
    )

    # Display the plot
    display(p)
    png("clipboard")
end