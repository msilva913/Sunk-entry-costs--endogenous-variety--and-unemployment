using Plots
import Plots:default
default(linewidth=2, grid=true, fontfamily="Computer Modern")

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
    plot!(p[6], irf_df.w_int, label=L"w_{int}")
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
    plot!(p[8], irf_df.w_R, label=L"w_R")
    yticks!(p[8], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Display the plot
    display(p)
    png("clipboard")
end

function gen_irf_comp(irf_bas::DataFrame, irf_alt::DataFrame, labels)
    # Create a 3x2 subplot layout
    p = Plots.plot(
        layout=(3, 3), 
        size=(1100, 700), 
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

      # Plot for p
      plot!(p[5], 
      [irf_bas.p irf_alt.p],
      label=[labels[1] labels[2]],
      alpha=0.6,
      title=L"p",
      ylabel="%"
  )

    # Plot for C_R
    plot!(p[6], 
        [irf_bas.C_R irf_alt.C_R],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"C_R",
        ylabel="%"   
    )

    # Plot for Y_R
    plot!(p[7], 
        [irf_bas.Y_R irf_alt.Y_R],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"Y_R",
        ylabel="%"
    )

    plot!(p[8], 
        [irf_bas.d_f irf_alt.d_f],
        label=[labels[1] labels[2]],
        alpha=0.6,
        title=L"d_f",
        ylabel="%"
    )

    plot!(p[9], 
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