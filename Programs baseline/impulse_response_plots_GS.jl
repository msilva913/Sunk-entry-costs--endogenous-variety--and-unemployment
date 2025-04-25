using Plots
gr()
import Plots:default
default(linewidth=2, grid=true, fontfamily="Computer Modern")

function gen_irf(irf_df::DataFrame)
    # Set up a 3x2 layout with specified size
    p = Plots.plot(layout=(1, 2), size=(1200, 500), 
            legend=:right, alpha=0.6)
    
    # Left plot: labor market variables
    plot!(p[1], irf_df.u, label=L"u", subplot=1, legend=:right)
    plot!(p[1], irf_df.z, label=L"z")
    plot!(p[1], irf_df.θ, label=L"θ")
    # Format y-axis as percentage
    yticks!(p[1], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    
    # Top right plot: vacancy value variables
    plot!(p[3], irf_df.v, label=L"v")
    plot!(p[3], irf_df.e, label=L"e")
    yticks!(p[3], :auto, fmt=x->string(round(x*100,digits=1),"%"))

    # Display the plot
    display(p)
    png("clipboard")
end

