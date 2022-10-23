using PyPlot
fig, ax = subplots()
ax.plot(2:10, sin.(2:10))
gcf()