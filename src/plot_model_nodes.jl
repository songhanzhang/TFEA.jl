function plot_model_nodes(Nodes; node_size = 2)

    scatter!(Nodes[:,2], Nodes[:,3], label = "",
             marker = :circle, markerstrokecolor = :dodgerblue, markercolor = :dodgerblue, markersize = node_size)

end
