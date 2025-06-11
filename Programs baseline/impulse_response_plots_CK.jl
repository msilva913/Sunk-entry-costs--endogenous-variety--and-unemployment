using Plots
import Plots:default
default(linewidth=2, grid=true, fontfamily="Computer Modern")

function gen_irf(irf_df::DataFrame)
    # Set up a 2x2 layout with specified size
    p = Plots.plot(layout=(2, 2), size=(1000, 500), 
            legend=:right, alpha=0.6)
    
    # Top left plot: labor market variables
    plot!(p[1], irf_df.u, label=L"u", subplot=1, legend=:right)
    plot!(p[1], irf_df.δ, label=L"δ")
    # Format y-axis as percentage
    yticks!(p[1], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
    # Top right plot: vacancy value variables
     plot!(p[2], irf_df.v, label=L"v")
    plot!(p[2], irf_df.θ, label=L"θ")
    #plot!(p[2], irf_df.e, label=L"e")
    yticks!(p[2], :auto, fmt=x->string(round(x*100,digits=1),"%"))
   

    # Bottom left plot: consumption and data-consistent counterpart
    plot!(p[3], irf_df.Q, label=L"Q")
    plot!(p[3], irf_df.K, label=L"K")
    plot!(p[3], irf_df.e, label=L"e")
    yticks!(p[3], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Bottom right plot
    plot!(p[4], irf_df.z, label=L"z")
    plot!(p[4], irf_df.w, label=L"w")
    yticks!(p[4], :auto, fmt=x->string(round(x*100,digits=1),"%"))
    
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