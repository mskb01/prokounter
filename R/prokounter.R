.lsume <- function( x, w=rep(1,length(x)) ){
  log(sum(w*exp(x)))
}

#' Get recovered abundance dependent statistics for a given genus in each sample in the data.
#'
#' @param mat a taxa x samples matrix.
#' @param genus a vector indicating the genus category.
#' @param g the genus for which the calculations are to be done.
#' @param div.too if
#' @return A dataframe with components
#'         \itemize{
#'         \item{lr - log abundances of the genus in each sample}
#'         \item{lm - log number of taxa for the genus in each sample}
#'         \item{samples - the sample name to which the row-wise entries correspond to}
#'         \item{g - the genus name given}
#'         \item{chao1 - the chao1 estimator as calculated from the package Vegan }
#'         \item{chao1.se - the standard error fot the chao1 estimate }
#'         \item{ace - the ACE estimator as calculated from the package Vegan }
#'         \item{ace.se - the standard error for the ACE estimate}
#'         }
.getRecovAbPat <- function(mat, genus, g,adt=TRUE){
  require(vegan)

  fi <- which( genus==g )
  sm <- colSums(mat[fi,,drop=FALSE])
  divs <- suppressWarnings(vegan:::estimateR( t(mat[fi,,drop=FALSE]) )) #warnings generated here  are with respect to estimation of asymptotic diversities.
  data.frame(
    "lr"=log(colSums(mat[fi,,drop=FALSE])),
    "lm"=log(colSums(mat[fi,,drop=FALSE]>0)),
    "samples"=colnames(mat),
    "g"=rep(g, ncol(mat)),
    "chao1"=divs["S.chao1",],
    "chao1.se"=divs["se.chao1",],
    "lrot"=log(colSums(mat[fi,,drop=FALSE]))-log(colSums(mat)),
    "ace"=divs["S.ACE",],
    "ace.se"=divs["se.ACE",]
  )
}

#' Obtains the prokounter within-genus taxa accumulation trend.
#'
#' The result objects generated with this function is used for differential richness inference
#' sample-wide (through the function \code{\link{getDRforSampleGroups}}),
#' genus-wide (\code{\link{getDRforGenera}}), and for collections of genera (e.g., phyla, \code{\link{getDRforTaxaCollection}}).
#'
#' @param mat a taxa x samples count matrix obtained from a 16S survey.
#'            We recommend not to post-process or filter 16S data
#'            so as to capture the accumulation rates of error-like low abundant, sparse taxa.

#' @param genus a character vector indicating the genus category for each taxa (i.e., row).
#'              Total length must equal the number of rows in the data matrix \code{mat}. More generally, this can be any annotation category within which one can assume most taxa discoveries are false.
#' @param fit.proc Trend deriving procedure. "ss" when the trend is estimated as smoothing spline (default). "lo" when a loess smoother is used.
#' @param toig Genera categories to ignore. By default this includes unknown or undefined categories \code{c("g__", "unknown", "NA", "", "undef")}.
#' @param fg.fdr False discovery rate threshold for assessing the significance of \code{fG} effects.
#'               By default set to \code{.1}. This parameter determines the significant genus effects incorporated in the technical trend while
#'               assessing differntial richness in functions \code{\link{getDRforSampleGroups}}, \code{\link{getDRforGenera}}, and \code{\link{getDRforTaxaCollection}}.
#' @param plt By default, set to FALSE. TRUE if recovered abundance trend \code{fR} term in manuscript's notation needs to be plotted.
#'            If using \code{fit.proc="ss"}, the random choice of knots may affect the shape of \code{fR}, but the net predicted mean of the linear predictor will be similar.
#' @param kn.rand FALSE by default. TRUE if ties in the choice of spline knots are broken at random.
#' @param ... other parameters to control the fits of smoothing spline \code{gss::ssanova} or the loess (loess) fitting procedures.
#' @return A list with three components.
#'        \itemize{
#'        \item{fit - the fitted trends object.}
#'        \item{
#'        df - a data frame with multiple components
#'        \itemize{
#'         \item{lr - log recovered abundances of the genus across samples}
#'         \item{lm - log number of taxa for the genus across samples}
#'         \item{samples - the sample name to which the row-wise entries correspond to}
#'         \item{g - the genus name given}
#'         \item{fit.proc - fitting procedure used}
#'         \item{chao1 - the chao1 estimator as calculated from the package Vegan }
#'         \item{chao1.se - the standard error fot the chao1 estimate }
#'         \item{ace - the ACE estimator as calculated from the package Vegan }
#'         \item{ace.se - the standard error for the ACE estimate}
#'         \item{ltv - the marginal recovered abundance dependent component (\code{fR} term in the main manuscript) with an additive constant}
#'         \item{ltv.se - the standard error of ltv }
#'         \item{gtv - the genus effect term. \code{fG} in the manuscript}
#'         \item{gtv.se - the standard error of \code{fG}}
#'         \item{gltv - estimates \code{fR + fG} with an additive constant  }
#'         \item{gltv.se - standard error of gltv }
#'         \item{gltv.wonly.gsigs - \code{fR + fG} with an additive constant, but only the statistically significant \code{fG} are incorporated in the sum.  }
#'         }
#'        }
#'        \item{ug - the set of genera over which the calculations were done.  }
#'        }
#' @export
getProkounterTrends <- function(mat, genus,
                       fit.proc="ss",
                       toig=c("g__", "unknown", "NA", "", "undef") ,
                       fg.fdr=.1,
                       plt=FALSE, kn.rand=FALSE,...
) {
  #prokounter has been changed to getProkounterTrends here.

  require(matrixStats)
  require(gss)
  require(RColorBrewer)
  require(MASS)

  stopifnot( nrow(mat) == length(genus) )

  # get genuses
  genus <- as.character(genus)
  ug <- unique( genus ); names(ug ) <- ug;
  g2ft <- ug[ which( !(ug %in% toig) ) ]


  #get recov ab
  df <- do.call( rbind, lapply( g2ft, function(g) .getRecovAbPat( mat, genus, g  ) ) )
  df <- df[ is.finite( df$lr ) & is.finite(df$lm),  ]
  df$g <- factor( df$g )

  if(fit.proc == "ss"){
    qlr <- c(.1,.5, .75, .999)
    q <- quantile( df$lr, prob=qlr, names=TRUE )
    which.min.btr <- function(y) {
      min.indx <- seq_along(y)[y == min(y, na.rm=TRUE)]
      if (length(min.indx) > 1)
        sample(min.indx, size = 1)
      else min.indx
    }
    if(kn.rand){
      id.basis <-  sapply( q, function(x) which.min.btr( ( x-df$lr )^2 ) )
    } else {
      id.basis <-  sapply( q, function(x) which.min( ( x-df$lr )^2 ) )
    }
    fit <- ssanova( lm ~ lr*g, data = df, id.basis = id.basis, ...  )
    ltv <- predict( fit, newdata=df, include=list( "1", "lr" ), se.fit=TRUE  )
    gltv <- predict( fit, newdata=df, include=list( "1", "lr", "g" ), se.fit=TRUE  )
    gtv <- predict( fit, newdata=df, include=list( "g" ), se.fit=TRUE  )

    i <- sapply( ug, function(gxx) match( gxx, df$g ) )
    ugg <- ug[ !is.na(i) ]
    i <- i[!is.na(i)]
    stopifnot( ugg == df$g[i] )
    ps <- 2*pnorm( -abs(gtv$fit[i]/gtv$se[i]) )
    sigs <- which( p.adjust( ps, "BH" ) < fg.fdr )
    gtv.sigs <- rep( 0, nrow(df) )
    for( gx in ugg[sigs] ){
      gtv.sigs[ which( df$g == gx ) ] <- 1
    }
    gltv.wonly.gsigs <- ltv$fit + gtv.sigs*gtv$fit

    df <- cbind( df, "ltv"=ltv$fit, "ltv.se"=ltv$se, "gtv"=gtv$fit, "gtv.se"=gtv$se, "gtv.sigs"=gtv.sigs,
                 "gltv"=gltv$fit, "gltv.se"=gltv$se, "gltv.wonly.gsigs"=gltv.wonly.gsigs
                 )
  } else if( fit.proc == 'lo' ){ #here the genus-wise components are undefined.
    fit <- loess( lm ~ lr, data = df, ...  )
    ltv <- predict( fit, newdata=df, se=TRUE  )
    df <- cbind( df, "ltv"=ltv$fit, "ltv.se"=ltv$se, "gtv"=NA, "gtv.se"=NA,"gtv.sigs"=NA,
                 "gltv"=NA, "gltv.se"=NA,"gltv.wonly.gsigs"=NA
                 )
  } else {
    stop("unkown fitting procedure.")
  }

  RES <- list()
  RES$fit <- fit
  RES$df <- df
  RES$ug <- ug
  RES$fit.proc <- fit.proc

  if(plt){
    .plot.fr( pkobj )
    }

  RES

}

#' Plots the fitted \code{fR} trend from the manuscript.
#' @param pkobj the output object of the getProkounterTrends() function.
.plot.fr <- function(pkobj){
  require(gss)
  df <- pkobj$df
  fit <- pkobj$fit
  g2ft <-pkobj$ug
  fit.proc <- pkobj$fit.proc
  if(fit.proc == "ss"){
    lrvec <- seq(min(df$lr), max(df$lr), .1)
    newd = data.frame( "lr"= rep(lrvec,
                                 each=length(g2ft) ),
                       "g"= rep( g2ft,
                                 each=length( lrvec)
                       )
    )
    ng <- length( df$g )
    ltv <- predict( fit, newdata=newd[,"lr",drop=FALSE], include=c( "1", "lr"), se.fit=TRUE);
    plot( df$lm ~ df$lr, col = 'gray', cex=.5, xlab = "Log Recovered Abundance", ylab = "Log Number of Taxa" )
    abline( h=0 )
    o <- order( newd$lr)
    lines( ltv$fit[o] ~ newd$lr[o], col = 'red', cex=.5 )
    lines( (ltv$fit[o]+1.96*ltv$se[o]) ~ newd$lr[o], col = 'red', cex=.5, lty=2 )
    lines( (ltv$fit[o]-1.96*ltv$se[o]) ~ newd$lr[o], col = 'red', cex=.5, lty=2 )
  } else if( fit.proc=="lo" ){
    lrvec <- seq(min(df$lr), max(df$lr), .1)
    newd = data.frame( "lr"= lrvec )
    ltv <- predict( fit, newdata=newd,se=TRUE ); abline( h=0 )
    plot( df$lm ~ df$lr, col = 'gray', cex=.5, xlab = "Log Recovered Abunadnce", ylab = "Log Number of Taxa" )
    o <- order( newd$lr )
    lines( ltv$fit[o] ~ newd$lr[o], col = 'red', cex=.5 )
    lines( (ltv$fit[o]+1.96*ltv$se[o]) ~ newd$lr[o], col = 'red', cex=.5, lty=2 )
    lines( (ltv$fit[o]-1.96*ltv$se[o]) ~ newd$lr[o], col = 'red', cex=.5, lty=2 )
  }
}


.getAccumControlForSamples <- function(pkobj, sample.names=unique( pkobj$df$samples ), fg.inc=TRUE, fg.sigs.only=TRUE ){

  if(fg.inc){
    if(!all(is.na(pkobj$df$gltv))){
      if(fg.sigs.only)
        netltv <- aggregate( pkobj$df$gltv.wonly.gsigs, by=list(pkobj$df$samples), .lsume )
      else
        netltv <- aggregate( pkobj$df$gltv, by=list(pkobj$df$samples), .lsume )
    } else { # loess trend
      netltv <- aggregate( pkobj$df$ltv, by=list(pkobj$df$samples), .lsume)
    }
  } else {
    netltv <- aggregate( pkobj$df$ltv, by=list(pkobj$df$samples), .lsume)
  }

  indx <- match( sample.names, netltv$Group  )
  y <- netltv$x[ indx ]
  names(y) <- sample.names
  y

}

.getAccumControlForGenera <- function(pkobj, genera.names=unique( pkobj$df$g ),  sample.names=unique( pkobj$df$samples ), ...  ){

  names( genera.names ) <- genera.names
  res <- lapply( genera.names, function(gi){
    so <- subset( pkobj$df, g==gi )
    .getAccumControlForSamples( pkobj = list( "df"=so ), sample.names=sample.names, ... )
  } )
  names(res) <- genera.names
  res

}


.getdrglm <- function( y, y.se=rep(1,length(y)), des, glm.fam, rem.obs.w.nainfzro.ses=FALSE, ... ){
  # Forces infinite variance to maximum non-zero variance and zero variances to minimum non-zero variance.
  # If this behavior is not suitable, subset the input variables according to preference.

  if( rem.obs.w.nainfzro.ses ){
    indx2keep <- !( (!is.finite(y.se)) | is.na(y.se) | (y.se==0) )
    y <- y[indx2keep]
    y.se <- y.se[indx2keep]
    des <- des[indx2keep,,drop=FALSE]
  } else {
    sv <- y.se[ !( !is.finite(y.se) | is.na(y.se) | (y.se==0) ) ]
    y.se[(!is.finite(y.se)) | is.na(y.se)] <- max(sv)
    y.se[y.se==0] <- min(sv)
  }

  y.w <- (1/(y.se^2) )/ ( sum( 1/(y.se^2) ) )
  if( glm.fam=="NB" ){
    fit <- glm.nb( round(y) ~ -1 + des, weights = y.w, ...  )
  } else if( glm.fam=="poisson") {
    fit <- glm(round(y) ~ -1 + des, weights=y.w, ...   )
  }

  tab <- summary( fit )$coef

  list( "fit"=fit,
        "tab"=tab
  )
}



.chao1.fn <- function(mat){
  require(vegan)
  divs <- vegan:::estimateR( t(mat) ) #warnings generated here usually are with respect to estimation of asymptotic diversities.
  list( y=divs["S.chao1",],
        y.se=divs["se.chao1",]
        )
}

.ace.fn <- function(mat){
  require(vegan)
  divs <- vegan:::estimateR( t(mat) ) #warnings generated here usually are with respect to estimation of asymptotic diversities.
  list( y=divs["S.ACE",],
        y.se=divs["se.ACE",]
  )
}

.Chao1.getDRforSampleGroups <- function( mat, des, pkobj, glm.fam="NB",
                                        depth.inc =FALSE,
                                        fg.inc=TRUE,
                                        fg.sigs.only=TRUE, rem.obs.w.nainfzro.ses=FALSE, ... ){
  # Developmental, uses Chao1 as the response variable instead of the observed richness as response variable in sample-wide richness inference.
  #
  # Code can be adapted for other asymptotic estiamtors.
  #
  # Forces any infinite Chao1 variances to maximum non-zero variance and zero variances to minimum non-zero variance.
  # If this behavior is not suitable, call .getdrglm directly with appropriate inputs.

  if(!all(colnames(mat) == rownames(des) ) ){
    stop( "Sample names i.e. column names of input matrix, mat, and rownames of design matrix, des, do not match" )
  }

  chao1 <- .chao1.fn( mat )
  netltv <- .getAccumControlForSamples( pkobj = pkobj, sample.names = colnames(mat), fg.inc = fg.inc, fg.sigs.only = fg.sigs.only )
  des <- cbind( des, netltv[ rownames(des)  ] )
  if(depth.inc)
    des <- cbind(  des, colSums( mat ) )

  .getdrglm( y=chao1$y, y.se=chao1$y.se, des = des, glm.fam = glm.fam, rem.obs.w.nainfzro.ses = rem.obs.w.nainfzro.ses, ...   )

}


.Chao1.getDRforGenera <- function( mat, des, pkobj, glm.fam="NB",
                                        depth.inc =FALSE,
                                        fg.inc=TRUE,
                                        fg.sigs.only=TRUE, rem.obs.w.nainfzro.ses=FALSE, ... ){
  # Developmental. Uses Chao1 as the response variable instead of the observed richness as response variable in sample-wide richness inference.
  #
  # Code can be adapted for other asymptotic estiamtors.
  #
  # Forces infinite variance to maximum non-zero variance and zero variances to minimum non-zero variance.
  # If this behavior is not suitable, call .getdrglm directly with appropriate inputs

  if( colnames(mat) != rownames(des) ){
    stop( "Sample names i.e. column names of input matrix, mat, and rownames of design matrix, des, do not match" )
  }

  gs <- unique( as.character(pkobj$df$g) ); names(gs) <- gs

  fits <- lapply( gs, function(gi){
    so <- subset( pkobj$df, g==gi )
    mmat <- mat[which( genus == gi ),,drop=FALSE]
    tryCatch({
      Chao1.getDRforSampleGroups( mat = mmat, pkobj = list( "df"=so ), des = des, glm.fam = glm.fam, resp.type = resp.type[1],
                            depth.inc = depth.inc,
                            fg.inc = fg.inc,
                            fg.sigs.only = fg.sigs.only, rem.obs.w.nainfzro.ses = rem.obs.w.nainfzro.ses, ...  )
    },error=function(e){
      print(e)
      NA
    })
  } )
  names(fits) <- gs

  #summarizing fits
  .summarizeMultSGFits( sgfits=fits, des=des )

}

#getSampleWideRichness renamed
#' Sample-wide differential richness inference for 16S surveys
#'
#' Obtains differential richness for sample groups and a given design matrix.
#'
#' @param pkobj the output object from \code{\link{getProkounterTrends}} function.
#' @param resp.type Three possible options. By default, \code{"lm"} when the association test is to be done with respect to observed number of taxa (sampled richness).
#'                  \code{"chao1"} or \code{"ace"} when the association test is to be done with the respective asymptotic richness estimators.
#'                  In 16S surveys, we often find these options to be equivalent.
#' @param des the design matrix.
#' @param depth.inc TRUE if logged sample depth (a sample's total count) needs be incorporated as a predictor. Default FALSE as we find prokounter's trend \code{ft} to capture sampling effort well.
#' @param fg.inc TRUE by default. TRUE if prokounter genus specific effects (the fg term in the prokounter trend model) need to be incorporated in deriving the net technical effect.
#'               This is available only when the prokounter trend was derived with \code{fit.proc="ss"} option from the function \code{\link{getProkounterTrends}}.
#' @param fg.sigs.only TRUE if only significant genus effects (fg term in the prokounter trend model) needs to be incorporated.
#'                     The significance of the effects is calculated using the Benjamini-Hoschberg procedure at a default false discovery rate threshold of .1.
#'                     This FDR threshold can be controlled when running \code{\link{getProkounterTrends}} through the option \code{fg.fdr}.
#' @param glm.fam the family function for the generalized linear model. \code{"poisson"} when estimation is performed by \code{glm}.
#'                Default \code{"NB"} for negative binomial, when estimation is performed with \code{MASS::glm.nb}.
#' @param ... other parameters to pass to \code{glm} or \code{glm.nb}
#' @return A list with two entries.
#'\itemize{
#'   \item{ \code{fit} - holds the GLM fit  }
#'   \item{ \code{tab} - the association statistics (effect size, standard error, test statistic, p-value)
#'         corresponding to each column in the design matrix, row-wise for each sample-group.}
#' }
#' @examples
#' #Obtain counts matrix and sample group information
#' library(metagenomeSeq)
#' library(prokounter)
#' data(mouseData)
#' taxaBySamplesCountMatrix <- MRcounts(mouseData, norm=FALSE)
#' des <- model.matrix( ~diet,  pData(mouseData) ) #the rownames of the design matrix must contain corresponding sample names.
#' genus <- as.character(fData(mouseData)$genus)
#' pkt <- prokounter::getProkounterTrends( mat = taxaBySamplesCountMatrix,genus=genus)
#' dr.sw <- getDRforSampleGroups( pkobj=pkt, des = des )
#' dr.sw$tab
#' @export
#' @seealso \code{\link{getDRforTaxaCollection}}, \code{\link{getDRforSampleGroups}}
getDRforSampleGroups <- function( pkobj, des, resp.type=c("lm","chao1","ace"), glm.fam="NB",
                                  depth.inc =FALSE,
                                  fg.inc=TRUE,
                                  fg.sigs.only=TRUE,...  ){
  require(MASS)
  netlr <- aggregate( pkobj$df$lr, by=list(pkobj$df$samples), .lsume)
  resp.type <- resp.type[1]
  if(resp.type=="lm"){
    netlm <- aggregate( pkobj$df$lm, by=list(pkobj$df$samples), .lsume)
  } else if( resp.type=="chao1" ){
    netlm <- aggregate( pkobj$df$chao1, by=list(pkobj$df$samples), function(x) log(sum( x, na.rm=TRUE )) )
  } else if( resp.type=="ace" ){
    netlm <- aggregate( pkobj$df$ace, by=list(pkobj$df$samples), function(x) log(sum( x, na.rm=TRUE )) )
  } else {
    stop("unknown response type")
  }

  if(fg.inc){
    if(!all(is.na(pkobj$df$gltv))){
      if(fg.sigs.only)
        netltv <- aggregate( pkobj$df$gltv.wonly.gsigs, by=list(pkobj$df$samples), .lsume )
      else
        netltv <- aggregate( pkobj$df$gltv, by=list(pkobj$df$samples), .lsume )
    } else { # loess trend
      netltv <- aggregate( pkobj$df$ltv, by=list(pkobj$df$samples), .lsume)
    }
  } else {
    netltv <- aggregate( pkobj$df$ltv, by=list(pkobj$df$samples), .lsume)
  }

  netdf <- merge( merge( netlr, netlm, by="Group.1" ), netltv, by="Group.1" )
  colnames( netdf ) <- c( "sample", "netlr", "netlm", "netltv"  )
  rownames( netdf ) <- netdf$sample


  netdf <- netdf[ rownames(des), ,drop=FALSE ]
  if( depth.inc )
    des <- cbind( des, depth=netdf$netlr, "netltv"=netdf$netltv )
  else
    des <- cbind( des, "netltv"=netdf$netltv )

  if( glm.fam=="NB" ){
    fit <- glm.nb( round(exp(netlm)) ~ -1 + des, data = netdf, ...  )
  } else if( glm.fam=="poisson") {
    fit <- glm(round( exp(netlm)) ~ -1 + des, data = netdf, family=glm.fam, ...   )
  }

  tab <- summary( fit )$coef

  list( "fit"=fit,
        "tab"=tab
        )

}


.summarizeMultSGFits <- function( sgfits, des, covnames=paste0( "des",colnames(des) ), gs=names(sgfits) ){
  gs <- names(sgfits);
  torem <- which(sapply(seq_along(sgfits), function(i){ is.na(sgfits[i]) } ))
  if(length(torem)>0){
    sgfits <- sgfits[ -torem  ]
    print( paste( "No valid fits found for ", length(torem), "features" ) )
  }

  summaries <-  lapply( seq_along(sgfits), function(i){
    gn <- gs[i]
    x <- sgfits[[i]]$tab
    data.frame( x, "g"=gn, "rn"=rownames(x)  )
  })
  full <- do.call( rbind, summaries )
  kx <- unique( full$rn ); names(kx) <- kx
  tab <- lapply( kx, function(u) {
    h <- subset( full, rn==u )
    rownames(h ) <- h$g
    h[,1:4,drop=FALSE]
  } )
  names(tab) <- kx

  list( "fits"=sgfits,
        "tab"=tab
  )
}

#getGenusSpecificRichness renamed
#' Genus-specific differential richness inference for 16S surveys
#'
#' Obtains differential richness for each genus in the study and a given design matrix.
#'
#' @param pkobj the output object from prokounter::getProkounterTrends function.
#' @param des the design matrix.
#' @param resp.type Three possible options. By default, \code{"lm"} when the association test is to be done with respect to observed number of taxa (sampled richness).
#'                  \code{"chao1"} when the association test is to be done with an asymptotic richness estimator Chao1.
#'                  \code{"ace"} when the association test is to be done with an asymptotic richness estimator ACE.
#'                  In 16S surveys, we have found these options to be equivalent because sampled richness is heavily correlated with asymptotic richness estimates.
#' @param glm.fam the glm family. Default "poisson" when estimation is performed by \code{glm}.
#'                "NB" for negative binomial, when eestimation is performed with \code{MASS::glm.nb}.
#'                #' @param depth.inc TRUE if logged sample depth (a sample's total count) needs be incorporated as a predictor. Default FALSE as we find prokounter's trend \code{ft} to capture sampling effort well.
#' @param depth.inc TRUE if logged sample depth (a sample's total count) needs be incorporated as a predictor. Default FALSE as we find prokounter's trend \code{ft} to capture sampling effort well.
#' @param ... other parameters to pass to \code{glm} or \code{glm.nb}
#' @return A list with two entries.
#' \itemize{
#'   \item{ fits - holds the glm fits for each genus.  }
#'   \item{ tab - a list where each entry holds the association statistics (effect size, standard error, test statistic, p-value)
#'         corresponding to each column in the design matrix, row-wise for each genus.}
#' }
#'
#' @examples
#' #Obtain counts matrix and sample group information
#' require(metagenomeSeq)
#' require(prokounter)
#' data(mouseData)
#' taxaBySamplesCountMatrix <- MRcounts(mouseData, norm=FALSE)
#' des <- model.matrix( ~diet,  pData(mouseData) ) #the rownames of the design matrix must contain corresponding sample names.
#' genus <- as.character(fData(mouseData)$genus)
#' pkt <- prokounter::getProkounterTrends( mat = taxaBySamplesCountMatrix,genus=genus)
#' dr.gw <- getDRforGenera( pkobj=pkt, des = des )
#' dr.gw$tab
#' @export
#' @seealso \code{\link{getDRforTaxaCollection}}, \code{\link{getDRforSampleGroups}}
getDRforGenera <- function(pkobj, des, resp.type=c("lm","chao1","ace"), glm.fam="poisson",
                           depth.inc =FALSE,
                           ... ){

  gs <- unique( as.character(pkobj$df$g) ); names(gs) <- gs

  fits <- lapply( gs, function(gi){
    so <- subset( pkobj$df, g==gi )
    tryCatch({
      getDRforSampleGroups( pkobj = list( "df"=so ), des = des, glm.fam = glm.fam, resp.type = resp.type[1],
                            depth.inc = depth.inc,
                            fg.inc = FALSE,
                            fg.sigs.only = FALSE, ...  ) #the fg.inc / fg.sigs.only do not have any effect on the genus-specific inferences as they only add global scaling terms for each genera.
    },error=function(e){
      print(e)
      NA
    })
  } )
  names(fits) <- gs

  #summarizing fits
  .summarizeMultSGFits( sgfits=fits, des=des )

}


#' Taxa collection-specific differential richness inference for 16S surveys.
#'
#' Obtains differential richness for collections of taxa and a given design matrix.
#' Each collection is viewed as a set of one or multiple genera. The function \code{\link{getDRforGenera}} is a special case of this function.
#'
#' @param pkobj the output object from prokounter::getProkounterTrends function.
#' @param des the design matrix.
#' @param genus2collection a data frame or a character matrix with the first column containing the genus
#'                         identifier and the second column the collection it has been asssigned to.
#'                         The genus identifier must be the same as those found in the original data for which the
#'                         prokounter::getProkounterTrends trends were established.
#' @param resp.type Three possible options. By default, \code{"lm"} when the association test is to be done with respect to observed number of taxa (sampled richness).
#'                  \code{"chao1"} when the association test is to be done with an asymptotic richness estimator Chao1.
#'                  \code{"ace"} when the association test is to be done with an asymptotic richness estimator ACE.
#'                  In 16S surveys, we have found these options to be equivalent because sampled richness is heavily correlated with asymptotic richness estimates.
#' @param glm.fam the glm family. Default "poisson" when estimation is performed by \code{glm}.
#'                "NB" for negative binomial, when eestimation is performed with \code{MASS::glm.nb}.
#' @param depth.inc TRUE if logged sample depth (a sample's total count) needs be incorporated as a predictor. Default FALSE as we find prokounter's trend \code{ft} to capture sampling effort well.
#' @param fg.inc TRUE by default. TRUE if prokounter genus specific effects (the fg term in the prokounter trend model) need to be incorporated in deriving the net technical effect.
#'               This is available only when the prokounter trend was derived with \code{fit.proc="ss"} option from the function \code{\link{getProkounterTrends}}.
#' @param fg.sigs.only TRUE if only significant genus effects (fg term in the prokounter trend model) needs to be incorporated.
#'                     The significance of the effects is calculated using the Benjamini-Hoschberg procedure at a default false discovery rate threshold of .1.
#'                     This FDR threshold can be controlled when running \code{\link{getProkounterTrends}} through the option \code{fg.fdr}.
#' @param ... other parameters to pass to \code{glm} or \code{glm.nb}
#' @return A list where each entry holds the association statistics (effect size, standard error, test statistic, p-value) corresponding to each column in the design matrix.
#' \itemize{
#'   \item{ fits - holds the glm fits for each collection  }
#'   \item{ tab - a list where each entry holds the association statistics (effect size, standard error, test statistic, p-value)
#'         corresponding to each column in the design matrix, row-wise for each collection.}
#' }
#'
#' @examples
#' #Obtain counts matrix and sample group information
#' require(metagenomeSeq)
#' require(prokounter)
#' data(mouseData)
#' taxaBySamplesCountMatrix <- MRcounts(mouseData, norm=FALSE)
#' des <- model.matrix( ~diet,  pData(mouseData) ) #the rownames of the design matrix must contain corresponding sample names.
#' genus <- as.character(fData(mouseData)$genus)
#' genus2collection <- fData(mouseData)[,c("genus", "phylum")]
#' pkt <- prokounter::getProkounterTrends( mat = taxaBySamplesCountMatrix,genus=genus)
#' dr.cw <- getDRforTaxaCollection( pkobj=pkt, des = des, genus2collection=genus2collection )
#' dr.cw$tab
#'
#' @export
#' @seealso \code{\link{getDRforGenera}}, \code{\link{getDRforSampleGroups}}
getDRforTaxaCollection <- function(pkobj, des, genus2collection = NULL, resp.type=c("lm","chao1","ace"), glm.fam="poisson",
                                   depth.inc =FALSE,
                                   fg.inc=TRUE,
                                   fg.sigs.only=TRUE, ... ){

  if(is.null(genus2collection)){
    collection <- pkobj$df$g
  } else {
    #first column in genus2collection correspond to genus; second
    collection <- genus2collection[match( pkobj$df$g, genus2collection[,1]),2]
  }
  gs <- unique( as.character(collection) ); names(gs) <- gs
  fits <- lapply( gs, function(gi){
    so <- pkobj$df[which(collection==gi ),]
    tryCatch({
      getDRforSampleGroups( list("df"=so), des = des, glm.fam = glm.fam, resp.type = resp.type[1],
                            depth.inc = depth.inc,
                            fg.inc = fg.inc,
                            fg.sigs.only = fg.sigs.only, ...  )
    },error=function(e){
      print(e)
      NA
    })
  } )

  names(fits) <- gs

  #summarizing fits
  .summarizeMultSGFits( sgfits=fits, des=des )
}


#' Bootstrap t confidence intervals for sample-wide differential richness inference form Prokounter.
#'
#' Developmental. The procedure aims to allow gauging reproducibility of inferences over the fitted prokounter trends (\code{ft} in manuscript) in \code{\link{getDRforSampleGroups}}.
#' We strongly recommend experimenting with the number of bootstrap samples (parameter B), until the returned confidence intervals stabilize.
#' We have found this to vary with data set size and complexities.
#' The code can be adapted for genus-level richness inference.
#' The procedure needs an ExpressionSet object (ref. example). This allows us to handle genus level calculations in an easier fashion in bootstrap samples.
#'
#' @param esetobj the ExpressionSet object
#' @param des design matrix for sample-wide differential richness inference
#' @param genusLabel \code{fData(esetobj)[[genusLabel]]} should return the genus categories. More generally, this can be any category within which we can assume most features are false.
#' @param ncore number of cores or processes to use for bootstrap
#' @param B Number of bootstrap samples. We strongly recommend experimenting with this parameter for your dataset.
#' @param alpha type I error rate. 1-alpha bootstrap t confidence intervals are constructed.
#' @param .... other parameters to pass to \code{\link{getDRforSampleGroups}}.
#'
#' @return A list where each entry holds the results for regression coefficients corresponding to each column in the design matrix.
#' \itemize{
#'   \item{ fit - holds the base result from \code{\link{getDRforSampleGroups}}   }
#'   \item{ bootstrapt.ci - the bootstrap confidence interval at the prescribed confidence level. }
#' }
#' @examples
#'  #Obtain counts matrix and sample group information
#' library(metagenomeSeq)
#' library(prokounter)
#' data(mouseData)
#' taxaBySamplesCountMatrix <- MRcounts(mouseData, norm=FALSE)
#' des <- model.matrix( ~diet,  pData(mouseData) ) #the rownames of the design matrix must contain corresponding sample names.
#' mouseEsetObj <- ExpressionSet( assayData = taxaBySamplesCountMatrix, phenoData = phenoData(mouseData), featureData = featureData(mouseData) )
#' genusLabel <- "genus" #fData(mouseEsetObj)[[genusLabel]] yields the genus annotations.
#' dr.sw.boot <- getDRforSampleGroups.bootWholeData( esetobj =  mouseEsetObj, des=des, genusLabel="genus", B = 5000 )
#' dr.sw.boot

#' @export
#' @seealso \code{\link{getDRforGenera}}, \code{\link{getDRforSampleGroups}}
getDRforSampleGroups.bootWholeData <- function(esetobj,  des, genusLabel="Genus", ncore=4, B=5000, alpha=.05, ...){
  require(doParallel)
  require(MASS)

  cl <- makeCluster( ncore )
  registerDoParallel(cl)

  bfn <- function( esetobj, des,... ){
    require(Biobase)

    tryCatch({
      si <- sample( seq(ncol( esetobj )), size = ncol(esetobj), replace = TRUE  )
      pkobj <- prokounter::getProkounterTrends( mat = exprs(esetobj[,si]),
                                                genus = as.character(fData(esetobj[,si])[[genusLabel]]),
                                                fit.proc = "ss", plt=FALSE, ...
      )
      desb <- des[si,,drop=FALSE]
      getDRforSampleGroups( pkobj = pkobj, des=desb, ... )
    },error=function(e){
      print(e)
      NA
    })
  }

  pkobj <- prokounter::getProkounterTrends( mat = exprs(esetobj),
                                            genus = as.character(fData(esetobj)[[genusLabel]]),
                                            fit.proc = "ss", plt=FALSE
  )
  giv.tab <-  getDRforSampleGroups( pkobj = pkobj, des=des, ... )$tab

  require(doParallel)
  res <- foreach( b = 1:B ) %dopar%{
    bfn( esetobj=esetobj, des=des, ...  )
  }

  stopCluster( cl )

  tokeep <- sapply( res, function(x) !is.na(x) )
  res <-res[tokeep]
  tt <- rownames( giv.tab )
  names( tt ) <- tt

  lapply( tt, function(nme){
    ts <- unlist( sapply( res, function(x){
      tryCatch({
        (x$tab[nme,1]-giv.tab[nme,1])/x$tab[nme,2]
      },error=function(e){
        NA
      })
    } ) )

    qs <- quantile( ts, c(alpha/2, 1-alpha/2), na.rm=TRUE )
    list( "fit"=giv.tab[nme,,drop=FALSE],
          "bootstrapt.ci"=c( "lower"=giv.tab[nme,1] - qs[2]*giv.tab[nme,2], "upper"=giv.tab[nme,1] - qs[1]*giv.tab[nme,2] )
          )
  }  )

}

