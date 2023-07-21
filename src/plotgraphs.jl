
using Plots
using LaTeXStrings
using DelimitedFiles


#fid=open("solution.dat","r");


mydata = readdlm("solution.dat", '\t', Float64, '\n')

#display(mydata)

#plot(nodecoords[:,1], nodecoords[:,2], line=(:black,0.9,5,:dot))
gr()
#plotlyjs()
#display(plot(mydata[:,2]*-0.01,   mydata[:,1]*0.01, line=(:black,0.9,3,:solid), lw=2));
p = plot(mydata[:,3]*-0.01,  mydata[:,1]*0.01, linecolor=:blue, linewidth=2, linestyle=:dash, label=L"\widetilde{v}", xlim=(0,3), legendfont=font(18));

display(p)
savefig(p, "myfig.png")
savefig(p, "myfig.pdf")