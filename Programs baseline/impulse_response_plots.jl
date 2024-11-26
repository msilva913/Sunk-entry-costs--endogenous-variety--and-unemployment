using Plots
using LaTeXStrings

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
    p = Plots.plot(layout=(2, 4), size=(1200, 500), 
            legend=:right, alpha=0.6)
    
    # Top left plot: labor market variables
    plot!(p[1], irf_df.u, label=L"u", subplot=1, legend=:right)
    plot!(p[1], irf_df.v, label=L"v")
    plot!(p[1], irf_df.θ, label=L"θ")
    plot!(p[1], irf_df.e, label=L"e")
    # Format y-axis as percentage
    yticks!(p[1], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Top middle plot: product line variables
    plot!(p[2], irf_df.N_e, label=L"N_e")
    plot!(p[2], irf_df.N, label=L"N")
    yticks!(p[2], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Top right plot: vacancy value variables
    plot!(p[3], irf_df.Q, label=L"Q")
    plot!(p[3], irf_df.K, label=L"K")
    yticks!(p[3], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom left plot: consumption and data-consistent counterpart
    plot!(p[4], irf_df.C, label=L"C")
    plot!(p[4], irf_df.C_R, label=L"C_R")
    yticks!(p[4], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom middle plot: output and data-consistent counterpart
    plot!(p[5], irf_df.Y, label=L"Y")
    plot!(p[5], irf_df.Y_R, label=L"Y_R")
    yticks!(p[5], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    plot!(p[6], irf_df.w, label=L"w")
    plot!(p[6], irf_df.w_R, label=L"w_R")
    plot!(p[6], irf_df.ls, label="labor share")
    yticks!(p[6], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    plot!(p[7], irf_df.d_f, label=L"d_f")
    plot!(p[7], irf_df.Y_c, label=L"Y_c")
    yticks!(p[7], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom right plot
    plot!(p[8], irf_df.z, label=L"z")
    plot!(p[8], irf_df.labor_prod, label="labor productivity")
    plot!(p[8], irf_df.δ, label=L"δ")
    #plot!(p[8], irf_df.Y, label=L"Y")
    plot!(p[8], irf_df.w, label=L"w")
    yticks!(p[8], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Display the plot
    display(p)
    png("clipboard")
end

function gen_irf_comp(irf_bas::DataFrame, irf_alt::DataFrame, labels)
    # Create a 3x2 subplot layout
    p = Plots.plot(
        layout=(2, 4), 
        size=(1200, 500), 
        legend=:topright,
        fmt=:png
    )

    # Plot for u
    plot!(p[1], 
        [irf_bas.u irf_alt.u],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"u",
        ylabel="%"
    )

    # Plot for v
    plot!(p[2], 
        [irf_bas.v irf_alt.v],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"v",
        ylabel="%"
    )

    # Plot for e
    plot!(p[3], 
        [irf_bas.e irf_alt.e],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"e",
        ylabel="%"
    )

    # Plot for N_e
    plot!(p[4], 
        [irf_bas.N_e irf_alt.N_e],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"N_e",
        ylabel="%"
    )

    # Plot for C_R
    plot!(p[5], 
        [irf_bas.C_R irf_alt.C_R],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"C_R",
        ylabel="%"   
    )

    # Plot for Y_R
    plot!(p[6], 
        [irf_bas.Y_R irf_alt.Y_R],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"Y_R",
        ylabel="%"
    )

    plot!(p[7], 
        [irf_bas.d_f irf_alt.d_f],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"d_f",
        ylabel="%"
    )

    plot!(p[8], 
    [irf_bas.ls irf_alt.ls],
    label=[labels[1] labels[2]],
    alpha=0.6,
    title=L"ls",
    ylabel="%"
)
    
    # Display the plot
    display(p)
    png("clipboard")
end