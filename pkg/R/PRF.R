# Compute PRFs for all respondents simultaneously.
# Also, compute a functional data object with all the PRFs.
# See: Ramsay, Hooker, & Graves (2009). Functional data analysis with R and Matlab.
PRF <- function(x, h=.09, N.FPts=101) #x = an object from 'PerFit' class
{
  matrix      <- x$Matrix
  Diffs       <- 1 - colMeans(matrix)
  focal.pts   <- seq(0, 1, length.out=N.FPts)
  GaussKernel <- function(x) {dnorm(x, mean=0, sd=1)}
  KernArg     <- expand.grid(focal.pts, Diffs)
  KernArg     <- matrix(KernArg[,1] - KernArg[,2], nrow=length(focal.pts), byrow=F) / h
  weights     <- GaussKernel(KernArg) / rowSums(GaussKernel(KernArg))
  PRFest      <- weights %*% t(matrix)
  # 
  # Specify a B-spline basis system (Chapter 3).
  # 'basis.bspline' below is a functional basis object of class 'basisfd'.
  # It is based on B-splines (piecewise polinomials, all of the same degree/order [order = degree + 1]).
  # Here we focus of degree three / order four splines (i.e., cubic polinomial segments),
  #   with one knot per break point.
  # This allows any two consecutive splines (piecewise polinomials), sp1 and sp2, with common break point BP,
  #   verifying sp1(BP) = sp2(BP), sp1'(BP) = sp2'(BP), and sp1''(BP) = sp2''(BP).
  # At 0 and 1 (extremes of the x-range), four (= order) knots are used.
  basis.bspline <- create.bspline.basis(rangeval = c(0, 1), norder = 4, nbasis = (4 + 9))
  # 
  # Specify coefficients c for the B-spline basis system computed above and then create functional data objects.
  # Based on smoothing using regression analysis (Section 4.3 in Ramsay et al., 2009).
  x.values     <- focal.pts
  basis.values <- eval.basis(evalarg=x.values, basisobj=basis.bspline)
  y.values     <- PRFest
  basis.coefs  <- solve(crossprod(basis.values), crossprod(basis.values, y.values))
  # Observe that 'basis.values %*% basis.coefs' is the B-spline approximation of PRFest.
  fd.obj       <- fd(basis.coefs, basis.bspline, list("Item difficulty", "Subject", "Probability correct answer"))
  list(PRFdiffs=focal.pts, PRFest=PRFest, FDO=fd.obj)
}

PRF.VarBands <- function (x, h=.09, N.FPts=101, alpha=.05) #x = an object from 'PerFit' class
{
  matrix    <- x$Matrix
  focal.pts <- seq(0, 1, length.out=N.FPts)
  N         <- dim(matrix)[1]; I <- dim(matrix)[2]
  PRFscores <- PRF(x, h, N.FPts)$PRFest
  # Jackknife estimate of the SE:
  PRF.SEarray <- array(NA, c(length(focal.pts), I, N))
  x.jack      <- x
  for (it in 1:I)
  {
    x.jack$Matrix       <- matrix[, -it]
    PRF.SEarray[, it, ] <- PRF(x.jack, h, N.FPts)$PRFest
  }
  PRF.SE <- apply(PRF.SEarray, 3, function(mat)
    {
    sqrt( ((I-1)/I) * rowSums((mat - rowMeans(mat))^2) )
    })
  crit.val         <- qnorm(1-alpha, mean=0, sd=1)
  PRF.VarBandsLow  <- PRFscores-crit.val*PRF.SE
  PRF.VarBandsHigh <- PRFscores+crit.val*PRF.SE
  # 
  # Specify a B-spline basis system.
  basis.bspline <- create.bspline.basis(rangeval = c(0,1), norder = 4, nbasis = (4 + 9))
  # Specify coefficients c for the B-spline basis system computed above and then create functional data objects.
  x.values      <- focal.pts
  basis.values  <- eval.basis(evalarg=x.values, basisobj=basis.bspline)
  # 
  y.values        <- PRF.VarBandsLow
  basis.coefs.Low <- solve(crossprod(basis.values), crossprod(basis.values, y.values))
  fd.obj.Low      <- fd(basis.coefs.Low, basis.bspline, list("Item difficulty", "Subject", "Probability correct answer"))
  # 
  y.values         <- PRF.VarBandsHigh
  basis.coefs.High <- solve(crossprod(basis.values), crossprod(basis.values, y.values))
  fd.obj.High      <- fd(basis.coefs.High, basis.bspline, list("Item difficulty", "Subject", "Probability correct answer"))
  #
  list(PRF.VarBandsLow=PRF.VarBandsLow, PRF.VarBandsHigh=PRF.VarBandsHigh, 
       FDO.VarBandsLow=fd.obj.Low, FDO.VarBandsHigh=fd.obj.High)
}

PRFplot <- function (x, respID, h=.09, N.FPts=101, 
                     VarBands=FALSE, VarBands.area=FALSE, alpha=.05,
                     Xlabel=NA, Xcex=1.5, Ylabel=NA, Ycex=1.5, title=NA, Tcex=1.5)
{
  N <- dim(x$Matrix)[1]; I <- dim(x$Matrix)[2]
  # Sanity check - Class PerFit:
  Sanity.cls(x)
  # Sanity check - Dichotomous data only:
  dico.PFS <- c("Cstar", "C.Sato", "U3", "ZU3", "G", "Gnormed", "D.KB", "r.pbis", "NCI", "Ht", "A.KB", "E.KB", "lz", "lzstar")
  poly.PFS <- c("Gpoly", "Gnormed.poly", "U3poly", "lzpoly")
  if (!(x$PFStatistic %in% dico.PFS))
  {
    stop('The person response function is only computed for dichotomous data. Aborted.')
  }
  # 
  res1 <- PRF(x, h, N.FPts)
  res2 <- PRF.VarBands(x, h, N.FPts, alpha)
  # 
  basis.bspline    <- create.bspline.basis(rangeval = c(0,1), norder = 4, nbasis = (4 + 9))
  x.values         <- seq(0,1,length.out=N.FPts)
  basis.values     <- eval.basis(evalarg=x.values, basisobj=basis.bspline)
  PRF.VarBandsLow  <- basis.values %*% res2$FDO.VarBandsLow$coefs
  PRF.VarBandsHigh <- basis.values %*% res2$FDO.VarBandsHigh$coefs
  for (i in 1:length(respID))
  {
    readline(prompt=paste0("Respondent ", respID[i], ": Press ENTER."))
    par(mar=c(4,4,2,1)+.1, las=1)
    plot(1, type="n", axes=FALSE, ann=FALSE, frame.plot=TRUE, xlim=c(0,1), ylim=c(0,1))
    tmpx <- if (is.na(Xlabel)) {"Item difficulty"} else {Xlabel}
    axis(1,at=seq(0, 1, by=.2)); mtext(side=1, text=tmpx, line=2.5, col="black", cex=Xcex, font=1)
    tmpy <- if (is.na(Ylabel)) {"Probability correct answer"} else {Ylabel}
    axis(2,at=seq(0, 1, by=.2)); mtext(side=2, text=tmpy, line=2.8, col="black", cex=Ycex, font=1, las=3)
    if (VarBands.area)
    {
      polygon(c(x.values,rev(x.values)),
              c(PRF.VarBandsHigh[,respID[i]], rev(PRF.VarBandsLow[,respID[i]])),col = "lightpink1",border=NA)
      par(new=TRUE); plot(res2$FDO.VarBandsLow[respID[i]],ann=FALSE, xlim=c(0,1),ylim=c(0,1), lty=2,lwd=1.5, href=FALSE, axes=FALSE)
      par(new=TRUE); plot(res2$FDO.VarBandsHigh[respID[i]],ann=FALSE,xlim=c(0,1),ylim=c(0,1), lty=2,lwd=1.5, href=FALSE, axes=FALSE)
      
    }
    if (!VarBands.area & VarBands)
    {
      par(new=TRUE); plot(res2$FDO.VarBandsLow[respID[i]],ann=FALSE, xlim=c(0,1),ylim=c(0,1), lty=2,lwd=1.5, href=FALSE, axes=FALSE)
      par(new=TRUE); plot(res2$FDO.VarBandsHigh[respID[i]],ann=FALSE,xlim=c(0,1),ylim=c(0,1), lty=2,lwd=1.5, href=FALSE, axes=FALSE)
    }
    par(new=TRUE); plot(res1$FDO[respID[i]],lwd=2,axes=FALSE,ann=FALSE,frame.plot=T,xlim=c(0,1),ylim=c(0,1), href=FALSE)
    tmp <- if (is.na(title)) {paste("PRF (respID # ",respID[i]
                                    ,")",sep="")} else {title}
    mtext(side=3, text=tmp, line=.5, col="black", cex=Tcex, font=2)
  }
  list(PRF.FDO = res1$FDO, VarBandsLow.FDO = res2$FDO.VarBandsLow,  VarBandsHigh.FDO = res2$FDO.VarBandsHigh)
}

# # PRF vs PRF.fda:
# id <- 55
# p1 <- PRF(x.new,id)
# p2 <- PRF.fda(x.new)
# # plot(p2[[3]])
# p1[[2]] - p2[[2]][,id]
# 
# # PRF.VarBands vs PRF.VarBands.fda
# p1vb <- PRF.VarBands(x.new,id)
# p2vb <- PRF.VarBands.fda(x.new)
# # plot(p2vb[[3]])
# # plot(p2vb[[4]])
# p1vb[[1]] - p2vb[[1]][,id]
# p1vb[[2]] - p2vb[[2]][,id]
# 
# # PRFplot vs PRFplot.fda
# PRFplot2(x.new,id, VarBands.area=TRUE)
# PRFplot.fda(x.new,id, VarBands.area=TRUE)
# # In the same window:
# PRFplot2(x,id)
# par(new=TRUE); PRFplot.fda(x,id)
# #
