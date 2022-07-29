heat.plot.fn = function (mat, ylabs, xlabs, main.title = "N sdreps", shape.size = 3, pch = 22, palette = "heat")
{
    origpar <- par(no.readonly = TRUE)
    par(mar = c(8, 5, 1,5), oma = c(1, 1, 3, 1), xpd = NA)
    z <- as.matrix(mat)

    z.unique = sort(unique(c(z)))
    print(length(z.unique))
    ncolors = length(z.unique)
    colors.unique = hcl.colors(ncolors, palette = palette)
    #print(heat.colors(ncolors))
    z.colors = matrix(colors.unique[match(z, z.unique)], NROW(z))
    legend.z = sort(unique(quantile(c(z), probs = seq(0,1,length = 11), type = 1)))
    legend.colors = colors.unique[match(legend.z, z.unique)]
    xvals = seq(0,1,length = NCOL(z))
    yvals = seq(0,1,length = NROW(z))
    plot(range(xvals), range(yvals), xlim = range(xvals), ylim = range(yvals), 
      xlab = "", ylab = "", type = "n", axes = FALSE)
    axis(1, at = xvals, lab = rep("", length(xlabs)))
    axis(2, at = yvals, lab = rep("", length(ylabs)))#, cex.axis = 0.75, las = 1)
    text(-0.17, yvals, labels = ylabs, srt = 0, cex = 0.75, adj = 0)
    text(xvals, -0.08, labels = xlabs, srt = 90, cex = 0.75, adj = c(1,0))
    box()
    for (j in 1:length(yvals)) {
      #vals = as.integer(z[j,])+1
      points(xvals, rep(yvals[j], length(xvals)),
        cex = shape.size, col = z.colors[j,], bg = z.colors[j,], pch = pch)
    }
    legend(x = 1.05, y = 0.5, xpd = TRUE, legend = legend.z,
        pch = rep(pch, length(legend.z)), pt.cex = shape.size,
        col = legend.colors, pt.bg = legend.colors, y.intersp = 2, x.intersp = 1, bty = "n")
    title(main.title, outer = TRUE, line = 1)
    par(origpar)
}
use.df.ems = df.ems[1:20,]
ynames = paste0(
  c("MR","SR")[match(use.df.ems$SR_model, c(2,3))],
  "_", c("ME","MF")[match(use.df.ems$M_est,c(TRUE,FALSE))],
  "_", use.df.ems$re_config)
xnames = paste0(
  c("HM","MSY")[match(df.oms$Fhist, c("H-MSY","MSY"))],
  "_", c("RsigL","RsigH")[match(df.oms$R_sig,c(0.5,1.5))],
  "_", c("Nsig0", "NsigL", "NsigH")[match(df.oms$NAA_sig, c(NA,0.25,0.5))], 
  "_", c("OEL", "OEH")[match(df.oms$obs_error, c("L","H"))])

png(file.path(here::here(), "Project_0", "results", "naa_om_n_sdrep.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(n_ems_sdrep, xlabs = xnames, ylabs = ynames, palette = "viridis")
dev.off()
#source("c:/work/SSRTWG/SSRTWG/Project_0/code/bubble_plot.R")

png(file.path(here::here(), "Project_0", "results", "naa_om_n_fit.png"), 
  width = 10*144, height = 10*144, res = 144, pointsize = 12)
heat.plot.fn(100 - n_ems_null, xlabs = xnames, ylabs = ynames, main.title = "N fits", palette = "viridis")
dev.off()
